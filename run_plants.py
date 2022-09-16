import subprocess
import os
import re
import pandas as pd
import shutil
from openbabel import openbabel
from rdkit import Chem
from rdkit.Chem import AllChem
import construct_ligand

PLANTS = '/home/florian/Desktop/Uni/Semester_IV/Frontiers_in_applied_drug_design/PLANTS/'
RANKING_FILE = PLANTS + 'output/ranking.csv'

def run_plants():
    '''run plants in a shell as a subprocess'''

    PLANTS = '/home/florian/Desktop/Uni/Semester_IV/Frontiers_in_applied_drug_design/PLANTS/'
    config = PLANTS + 'plantsconfig'
    os.chdir(PLANTS)
    # remove output dir if it already exists
    if os.path.exists(PLANTS + 'output'):
        shutil.rmtree(PLANTS + 'output')
    args = ['./PLANTS', '--mode', 'screen', config]
    process = subprocess.run(args)


def get_best_scoring_pose(ranking_file):
    '''Function retrieves scores of the best model'''

    ranking = pd.read_csv(ranking_file)
    pose_file = ranking.iloc[0]['LIGAND_ENTRY'] + '.mol2'
    score_norm_heavatoms = ranking.iloc[0]['SCORE_NORM_HEVATOMS']
    return pose_file, score_norm_heavatoms

def clear_output_dir():
    '''removes output dir with all subfiles'''

    if os.path.exists(PLANTS + 'output'):
        shutil.rmtree(PLANTS + 'output')

def write_mol_to_sdf(mol, path):
    # add hydrogens and optimize molecule
    mol.UpdatePropertyCache()
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    # write to file
    w = Chem.SDWriter(path)
    w.write(mol)
    w.close()

def convert_mol_files(path, in_format, out_format):
    '''converts the mol file format into another'''

    out_path = path.replace(in_format, out_format)
    print(f'in_path: {path}')
    print(f'out path: {out_path}')
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(in_format, out_format)
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, path)
    obConversion.WriteFile(mol, out_path)
    return out_path

def get_mol_from_mol2file(mol2file_path):
    sdf_path = convert_mol_files(mol2file_path, 'mol2', 'sdf')
    mol = Chem.MolFromMolFile(sdf_path)
    return mol


def renew_ligand_file():
    '''function removes all the old files from the previous run and restores the original ligand file'''

    if os.path.exists(PLANTS + 'ligand.sdf'):
        os.remove(PLANTS + 'ligand.sdf')
    if os.path.exists(PLANTS + 'ligand.mol2'):
        os.remove(PLANTS + 'ligand.mol2')
    shutil.copyfile(PLANTS + 'ligand_original.mol2', PLANTS + 'ligand.mol2')

def dock_molecule(mol):
    '''
    function optimizes the molecule, sets up the ligand.mol2 file required and runs PLANTS.
    :param mol: ligand as mol object
    :return: best pose with score
    '''
    PLANTS = '/home/florian/Desktop/Uni/Semester_IV/Frontiers_in_applied_drug_design/PLANTS/'
    RANKING_FILE = PLANTS + 'output/ranking.csv'

    mol.UpdatePropertyCache()
    mol = Chem.AddHs(mol)
    #print([n.GetSymbol() for n in list(mol.GetAtoms())[atom_idx].GetNeighbors()])
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    # write sdf and convert in mol2
    write_mol_to_sdf(mol, PLANTS + 'ligand.sdf')
    convert_mol_files(PLANTS + 'ligand.sdf', 'sdf', 'mol2')
    # run plants on modified ligand
    run_plants()
    pose_file, score = get_best_scoring_pose(RANKING_FILE)
    # get best pose as rdkit molecule
    print(f'pose_file: {pose_file}')
    mol = get_mol_from_mol2file(PLANTS + f'output/{pose_file}')
    return mol, score


def main():
    # run_plants()
    # ranking_file = '/home/florian/Desktop/Uni/Semester_IV/Frontiers_in_applied_drug_design/PLANTS/output/ranking.csv'
    # print(get_best_score(ranking_file))
    path = '/home/florian/Desktop/Uni/Semester_IV/Frontiers_in_applied_drug_design/PLANTS/ligand.sdf'
    convert_mol_files(path, 'sdf', 'mol2')

if __name__ == '__main__':
    main()