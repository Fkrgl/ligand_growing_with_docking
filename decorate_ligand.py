import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from construct_ligand import get_base_fragment_indices


def align_to_basefragment(mol, base_fragment):
    '''
    function alignes molecule to the base fragment an returns the aligned molecule
    '''
    mol_aligned = Chem.Mol(mol)
    # get substructure match index
    sub_idx = get_base_fragment_indices(mol, base_fragment)
    # get Affine Transformation Matrix
    alignment = rdkit.Chem.rdMolAlign.GetAlignmentTransform(mol, base_fragment, atomMap=list(zip(sub_idx, sub_idx)))
    M = alignment[1]
    # built matrix from mol coords
    n = len(mol.GetAtoms())
    mol_pos = mol.GetConformer().GetPositions()
    mol_pos_ones = np.hstack((mol_pos, np.ones((n,1))))
    # transform mol onto base fragment
    transformed = M.dot(mol_pos_ones.T)
    # adapt atom coords of mol
    for i in range(n):
        pos = mol_aligned.GetConformer().GetAtomPosition(i)
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
