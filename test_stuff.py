from rdkit import Chem
from rdkit.Chem import AllChem
import os
from openbabel import openbabel
import re
import pandas as pd
from anytree import AnyNode
from multiprocessing import Pool
import multiprocessing
import glob
import argparse


def write_plantsconfig(id, lig_file):
    PLANTS = "/home/florian/Desktop/Uni/Semester_IV/Frontiers_in_applied_drug_design/PLANTS/"
    plantsconfig = PLANTS + 'plantsconfig'
    OUT = PLANTS + f'parallel_output/plantsconfig_{id}'
    # create output dir and input ligand
    tmp_out = PLANTS + f'parallel_output/{id}'
    f = open(plantsconfig, 'r')
    lines = f.readlines()
    out_dir = f'output_dir\t\t\tparallel_output/plantsconfig_{id}'
    lines[3] = out_dir
    print(lines[1])
    lines[1] = f'ligand_file\t\t\t{lig_file}'
    f.close()
    f = open(PLANTS + f"plantsconfig_{id}", 'w')
    f.writelines(lines)
    f.close()
    #os.remove(PLANTS + f"plantsconfig_{id}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', action='store_true')
    parser.add_argument('--no-index', dest='feature', action='store_false')
    parser.set_defaults(feature=True)
    args = parser.parse_args()
    print(args.feature)