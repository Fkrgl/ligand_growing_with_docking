import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from construct_ligand import get_linker_atom_and_neighbor, check_molecule
from scipy.spatial.distance import cdist

def isRingAromatic(mol,bondRing):
    '''
    checks each bond of a ring for aromaticity
    :param mol: molecule that contains the ring
    :param bondRing: bond indices of ring
    '''
    for id in bondRing:
        if not mol.GetBondWithIdx(id).GetIsAromatic():
             return False
    return True

def get_aromatic_rings(mol):
    '''
    retruns the indices of all atoms of aromatic ring systems with at least one hydrogen in the molecule
    '''
    aromatic_ring_atom_idx = set()
    rings = mol.GetRingInfo()
    # check all bonds for aromaticity
    for ring in rings.BondRings():
        print(ring)
        if isRingAromatic(mol, ring):
            # get atom indices of aromatic bonds
            for bond_idx in ring:
                bond = mol.GetBondWithIdx(bond_idx)
                aromatic_ring_atom_idx.add(bond.GetBeginAtomIdx())
                aromatic_ring_atom_idx.add(bond.GetEndAtomIdx())
    return list(aromatic_ring_atom_idx)

def align_to_basefragment(mol, base_fragment, aromatic_ring_idx):
    '''
    function aligns molecule to the aromatic ring atoms of the base fragment an returns the aligned molecule
    :param mol: grown molecule
    :param base_fragment:
    :param aromatic_ring_idx:
    :return:
    '''
    # remove hydrogens
    mol = rdkit.Chem.rdmolops.RemoveHs(mol)
    base_fragment = rdkit.Chem.rdmolops.RemoveHs(base_fragment)
    mol_aligned = Chem.Mol(mol)
    # get Affine Transformation Matrix M
    alignment = rdkit.Chem.rdMolAlign.GetAlignmentTransform(mol, base_fragment,
                                                            atomMap=list(zip(aromatic_ring_idx, aromatic_ring_idx)))
    M = alignment[1]
    # built matrix from mol coords
    n = len(mol.GetAtoms())
    mol_pos = mol.GetConformer().GetPositions()
    mol_pos_ones = np.hstack((mol_pos, np.ones((n,1))))
    # transform mol onto base fragment
    transformed = M.dot(mol_pos_ones.T)
    # adapt atom coords of mol
    for i in range(n):
        mol_aligned.GetConformer().SetAtomPosition(i, transformed[:-1,i])
    return mol_aligned

def calc_group_position(mol, idx, length):
    '''
    given an atom index the function calculates the position of an atom (part of the decoration group) that
    is hypotheically substitued for the hydrogen at the atom with a specified bond length
    '''
    # add hs with coords
    mol = rdkit.Chem.rdmolops.AddHs(mol, addCoords=True)
    pos_atom = mol.GetConformer().GetAtomPosition(idx)
    # substitute one of the H atoms
    for neighbor in mol.GetAtomWithIdx(idx).GetNeighbors():
        if neighbor.GetSymbol() == 'H':
            neighbor_idx = neighbor.GetIdx()
            pos_h = mol.GetConformer().GetAtomPosition(neighbor_idx)
            # get dircetion vector of the H atom
            h_vector = pos_atom.DirectionVector(pos_h)
            # get position of new group
            group_pos = pos_atom + h_vector*length
            mol.GetConformer().SetAtomPosition(neighbor_idx, group_pos)
            return group_pos
    return


def has_spacial_neighbors(atom_pos, protein_coords, min_d=2):
    '''
    function returns True if any protein atom has a distance smaller than d to the specified atom of the ligand
    '''
    atom_pos = np.asarray(atom_pos).reshape(1, -1)
    distances = cdist(protein_coords, atom_pos, 'euclidean')
    min_dist = min(distances)[0]
    print(f'min = {min_dist:.1f}')
    if min_dist < min_d:
        return True
    else:
        return False

def add_functional_group(mol, pos, substituent_smiles, linker_atom_symbol, bond_type=Chem.rdchem.BondType.SINGLE):
    '''
    function merges a functional group to a molecle. The bond is built betwenn the nieghbor of the Au atom in
    the functional group and the specified atom in the molecule
    '''
    mol = rdkit.Chem.rdmolops.AddHs(mol, addCoords=True)
    pos_atom = mol.GetConformer().GetAtomPosition(pos)
    # substitute one of the H atoms
    for neighbor in mol.GetAtomWithIdx(pos).GetNeighbors():
        if neighbor.GetSymbol() == 'H':
            neighbor_idx = neighbor.GetIdx()
            mol.Compute2DCoords()
            substituent = Chem.MolFromSmiles(substituent_smiles)
            substituent.Compute2DCoords()
            combo = Chem.CombineMols(mol, substituent)
            linker_atom, linker_atom_neighbor = get_linker_atom_and_neighbor(combo, linker_atom_symbol)
            # make combo editable
            edcombo = Chem.EditableMol(combo)
            # add bond between linker and neighbor of au atom
            edcombo.AddBond(pos, linker_atom_neighbor.GetIdx(), order=bond_type)
            # remove AU and its bond to the linker
            mol = edcombo.GetMol()
            mol = Chem.RWMol(mol)
            mol.RemoveBond(linker_atom.GetIdx(), linker_atom_neighbor.GetIdx())
            mol.RemoveAtom(linker_atom.GetIdx())
            # remove hydrogen and its bond to the fragment
            if mol.GetAtomWithIdx(neighbor.GetIdx()).GetSymbol() != 'H':
                for n in mol.GetAtomWithIdx(pos).GetNeighbors():
                    if n.GetSymbol() == 'H':
                        neighbor = n
            mol.RemoveBond(pos , neighbor.GetIdx())
            mol.RemoveAtom(neighbor.GetIdx())
            mol.Compute2DCoords()
            # check if mol is valid
            mol = check_molecule(mol)
            if mol:
                return mol
            else:
                return None

def hasHydogenNeighbor(mol, idx):
    '''
    function checks if an atom has a hydrogen as neighbor
    '''
    for neighbor in mol.GetAtomWithIdx(idx).GetNeighbors():
        if neighbor.GetSymbol == 'H':
            return True
    return False

def decorate_ligand(mol, aromatic_atom_idx, protein_coords, bond_length, atoms_to_functionals):
    '''
    function adds functional groups to all aromatic ring atoms if they do not clash with any protein atoms
    '''
    mol_with_functionals = []
    group_atoms = ['C', 'O', 'N']
    for idx in aromatic_atom_idx:
        # check if atom has hydrogen as neighbor
        neighbors = [neighbor.GetSymbol() for neighbor in mol.GetAtomWithIdx(idx).GetNeighbors()]
        if 'H' in neighbors:
            for group in group_atoms:
                group_pos = calc_group_position(mol, idx, bond_length['C' + group])
                # check if coords of group would clash
                if has_spacial_neighbors(group_pos, protein_coords):
                    continue
                # no clashes: add functional groups to mol
                else:
                    # go through each functional group with this linker atom
                    for functional in atoms_to_functionals[group]:
                        # give multiple tries
                        while True:
                            decorated_mol = add_functional_group(mol, idx, functional, '[Au]')
                            # check if something was returned
                            if decorated_mol:
                                mol_with_functionals.append(decorated_mol)
    return mol_with_functionals


def main():
    base_fragment_path = '/home/florian/Desktop/Uni/Semester_IV/Frontiers_in_applied_drug_design/grown_molecules/root.sdf'
    grown_mol_path = '/home/florian/Desktop/Uni/Semester_IV/Frontiers_in_applied_drug_design/grown_molecules/1_0_0_0.sdf'
    base_fragment = Chem.MolFromMolFile(base_fragment_path)
    grown_mol = Chem.MolFromMolFile(grown_mol_path)
    align_to_basefragment(grown_mol, base_fragment)

if __name__ == '__main__':
    main()
