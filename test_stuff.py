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


for path in glob.glob('/home/florian/Desktop/Uni/Semester_IV/Frontiers_in_applied_drug_design/grown_molecules/*.sdf'):
    print(path)