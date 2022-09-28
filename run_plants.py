import subprocess
import os
import re
import pandas as pd
import shutil
from openbabel import openbabel
from rdkit import Chem
from rdkit.Chem import AllChem
import construct_ligand


def run_plants(PLANTS):
    '''run plants in a shell as a subprocess'''

    #PLANTS = '/home/florian/Desktop/Uni/Semester_IV/Frontiers_in_applied_drug_design/PLANTS/'
    config = PLANTS + 'plantsconfig'
    os.chdir(PLANTS)
    # remove output dir if it already exists
    if os.path.exists(PLANTS + 'output'):
        shutil.rmtree(PLANTS + 'output')
    args = ['./PLANTS', '--mode', 'screen', config]
    process = subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


def get_best_scoring_poses(ranking_file, PLANTS, cutoff=0.95):
    '''
    Function retrieves pose of best model and all models with a 95% of the best models score
    '''

    #PLANTS = '/home/florian/Desktop/Uni/Semester_IV/Frontiers_in_applied_drug_design/PLANTS/output/'
    plants_output = PLANTS + 'output/'
    ranking = pd.read_csv(ranking_file)
    best_score = ranking.iloc[0]['SCORE_NORM_HEVATOMS']
    filtered_poses = ranking[ranking['SCORE_NORM_HEVATOMS'] <= best_score*cutoff][['SCORE_NORM_HEVATOMS',
                                                                                  'LIGAND_ENTRY']]
    scores = filtered_poses['SCORE_NORM_HEVATOMS'].values
    pose_files = filtered_poses['LIGAND_ENTRY'].values
    poses = [get_mol_from_mol2file(plants_output + file + '.mol2') for file in pose_files]
    return poses, scores

def clear_output_dir(PLANTS):
    '''removes output dir with all subfiles'''

    if os.path.exists(PLANTS + 'output'):
        shutil.rmtree(PLANTS + 'output')

def write_mol_to_sdf(mol, path):
    w = Chem.SDWriter(path)
    w.write(mol)
    w.close()

def convert_mol_files(path, in_format, out_format):
    '''converts the mol file format into another'''

    out_path = path.replace(in_format, out_format)
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


def renew_ligand_file(PLANTS):
    '''function removes all the old files from the previous run and restores the original ligand file'''

    if os.path.exists(PLANTS + 'ligand.sdf'):
        os.remove(PLANTS + 'ligand.sdf')
    if os.path.exists(PLANTS + 'ligand.mol2'):
        os.remove(PLANTS + 'ligand.mol2')
    shutil.copyfile(PLANTS + 'ligand_original.mol2', PLANTS + 'ligand.mol2')

def dock_molecule(mol, PLANTS):
    '''
    function optimizes the molecule, sets up the ligand.mol2 file required and runs PLANTS.
    :param mol: ligand as mol object
    :return: best pose with score
    '''

    RANKING_FILE = PLANTS + 'output/ranking.csv'

    mol.UpdatePropertyCache()
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    # write sdf and convert in mol2
    write_mol_to_sdf(mol, PLANTS + 'ligand.sdf')
    convert_mol_files(PLANTS + 'ligand.sdf', 'sdf', 'mol2')
    # run plants on modified ligand
    run_plants(PLANTS)
    # get best pose as rdkit molecule with scores
    poses, score = get_best_scoring_poses(RANKING_FILE, PLANTS)
    return poses, score


def main():
    # ranking_file = '/home/florian/Desktop/Uni/Semester_IV/Frontiers_in_applied_drug_design/PLANTS/output/ranking.csv'
    # print(get_best_score(ranking_file))
    start_mol = Chem.MolFromSmiles('C1=CC=C2C(=C1)C=CN2')
    mol, score = dock_molecule(start_mol)
    print(mol.GetConformer().GetPositions())

if __name__ == '__main__':
    main()