from rdkit import Chem
from rdkit.Chem import AllChem
import random
from collections import Counter
import numpy as np
import sys
# +==================================================+ my packages +===================================================+
import construct_ligand
import run_plants

PLANTS = '/home/florian/Desktop/Uni/Semester_IV/Frontiers_in_applied_drug_design/PLANTS/'
RANKING_FILE = PLANTS + 'output/ranking.csv'

def main():
    # clean before new start
    run_plants.renew_ligand_file()
    scored_poses = []
    # read in mol
    mol = run_plants.get_mol_from_mol2file(PLANTS + 'ligand.mol2')
    # translate into raw mol to edit it
    mol = Chem.RWMol(mol)
    print(type(mol))
    # redock mol to get initial score
    run_plants.run_plants()
    pose_file, score = run_plants.get_best_scoring_pose(RANKING_FILE)
    scored_poses.append((score, mol))

    # construct ligand and run PLANTS
    for i in range(2):
        # set name of mol to avoid accumulation of filenname
        mol.SetProp("_Name", "1pmv_s_p")
        print(f'ITERATION {i}')
        atom_idx = construct_ligand.choose_next_position(mol)
        print(f'The next atom will be added at atom with index {atom_idx}')
        mol = construct_ligand.add_group(mol, 'C', atom_idx)
        print(atom_idx, type(atom_idx))
        print([n.GetSymbol() for n in mol.GetAtomWithIdx(int(atom_idx)).GetNeighbors()])
        # calculate valence states of all atoms
        mol.UpdatePropertyCache()
        mol = Chem.AddHs(mol)
        print([n.GetSymbol() for n in list(mol.GetAtoms())[atom_idx].GetNeighbors()])
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        # write sdf and convert in mol2
        construct_ligand.write_mol_to_sdf(mol, PLANTS + 'ligand.sdf')
        run_plants.convert_mol_files(PLANTS + 'ligand.sdf', 'sdf', 'mol2')
        # run plants on modified ligand
        run_plants.run_plants()
        pose_file, score = run_plants.get_best_scoring_pose(RANKING_FILE)
        # get best pose as rdkit molecule
        print(f'pose_file: {pose_file}')
        mol = run_plants.get_mol_from_mol2file(PLANTS + f'output/{pose_file}')
        scored_poses.append((score, mol))
        # result file becomes new ligand
        # here we need a criteria how to further modify                      ++++++++++++++++ next step ++++++++++++++++
        print(f'scored_poses: \n{scored_poses}')

    # compare first and last molecule in list
    w = Chem.SDWriter(PLANTS + 'mol_org.sdf')
    w.write(scored_poses[0][1])
    w.close()

    w = Chem.SDWriter(PLANTS + 'mol_last.sdf')
    w.write(scored_poses[-1][1])
    w.close()

if __name__ == '__main__':
    main()