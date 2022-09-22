import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from construct_ligand import get_base_fragment_indices

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
    retruns the indices of all atoms of aromatic ring systems in the molecule
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

def main():
    base_fragment_path = '/home/florian/Desktop/Uni/Semester_IV/Frontiers_in_applied_drug_design/grown_molecules/root.sdf'
    grown_mol_path = '/home/florian/Desktop/Uni/Semester_IV/Frontiers_in_applied_drug_design/grown_molecules/1_0_0_0.sdf'
    base_fragment = Chem.MolFromMolFile(base_fragment_path)
    grown_mol = Chem.MolFromMolFile(grown_mol_path)
    align_to_basefragment(grown_mol, base_fragment)

if __name__ == '__main__':
    main()
