from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem.Draw
import random
from collections import Counter
import numpy as np
from openbabel import openbabel
import pandas as pd


def add_group(mol, atom_type, atom_idx):
    '''
    add a new atom to the molecule and link it to an existing atom via a single bond. For this, a hydrogen with its bond
    is removed and a new bond is added
    :param mol: molecule
    :param atom_type: atom type of the added atom
    :param atom_idx: index of atom that is linked to new atom
    :return: added molecule
    '''
    if type(mol) != Chem.rdchem.RWMol:
        mol = Chem.RWMol(mol)
    mol = Chem.AddHs(mol)
    mol = remove_bonded_hydrogen(mol, atom_idx)
    index_added_atom = mol.AddAtom(Chem.Atom(atom_type))
    mol.AddBond(int(atom_idx), int(index_added_atom), Chem.BondType.SINGLE)
    return mol

def add_linker(atom_idx, linker, mol):
    """
    adds a linker group on a specified position in the molecule
    :param atom_idx: index of atom that is bound to the linker group
    :param linker: linker group
    :param mol: molecule that is modified
    """
    if type(mol) != Chem.rdchem.RWMol:
        mol = Chem.RWMol(mol)


def remove_bonded_hydrogen(mol, atom_idx):
    '''given an atom, the function deletes an existing bonded hydrogen and the bond '''
    mol = Chem.RWMol(mol)
    for neighbor in list(mol.GetAtoms())[atom_idx].GetNeighbors():
        if neighbor.GetSymbol() == 'H':
            mol.RemoveBond(neighbor.GetIdx(), int(atom_idx))
            mol.RemoveAtom(neighbor.GetIdx())
            return mol

def choose_next_position(mol):
    '''
    choose the position in the atom on which the next fragment is linked
    :param mol:
    :return:
    '''
    # choose atom randomly and consider only carbons
    atoms = list(mol.GetAtoms())
    available_carbons = np.zeros(len(list(mol.GetAtoms())))
    atomic_symbols = np.array([atom.GetSymbol() for atom in atoms])
    available_carbons[atomic_symbols == 'C'] = 1

    for atom in atoms:
        print(atom.GetTotalNumHs())

    # check if a further atom can be added to each carbon
    # terminate if all possibilities are tried
    while len(np.where(available_carbons == 1)[0]) != 0:
        # generate list af all available carbons
        carbons = np.where(available_carbons == 1)[0]
        idx_next_carbon = random.randint(0, len(carbons)-1)
        carbon_idx = carbons[idx_next_carbon]
        # check for valenz
        n_neighbors =  len(list(atoms[carbon_idx].GetNeighbors()))
        n_hydrogens = atoms[carbon_idx].GetTotalNumHs()
        # less than four atoms and at least one hydrogen (substitute the hydrogen bby the next atom)
        if n_neighbors < 4 and n_hydrogens > 0:
            # further atom can be added even if carbon has one double bond
            return carbon_idx
        else:
            # remove carbon from list
            available_carbons[carbon_idx] = 0
            print('else block')
    print('No carbon found to further extend the chain')
    return -1

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

def load_libraries(path_frag, path_link):
    """
    loads fragment and linker library
    :param path_frag: path to fragment lib
    :param path_link: path to linker lib
    :return: list of fragment and list of linker molecules (rdkit Mol)
    """
    fragments = []
    linkers = []
    frag_lib = pd.read_csv(path_frag, sep='\t', header=None)
    linker_lib = pd.read_csv(path_link, sep='\t', header=None)
    # read in fragments
    for i in range(len(frag_lib) - 1):
        try:
            mol = Chem.MolFromSmiles(frag_lib.iloc[i].values[0])
            fragments.append(mol)
        except:
            print(f'SMILES {frag_lib.iloc[i].values[0]} cannot be translates into a molecule. Please check the SMILES')

    # read in linkers
    for i in range(len(linker_lib) - 1):
        try:
            mol = Chem.MolFromSmiles(linker_lib.iloc[i].values[0])
            linkers.append(mol)
        except:
            print(f'SMILES {linker_lib.iloc[i].values[0]} cannot be translates into a molecule. Please check the SMILES')

    return fragments, linkers

def grow_ligand_at_random(smiles, n_iterations):
    """
    grows the ligand at ra random position by one carbon each round and saves the resulting ligand in the file
    modified_compound1.sdf
    :param n_iterations: number of growing iterations
    :param smiles: base fragment as smiles
    """

    mol = Chem.MolFromSmiles(smiles)
    # translate into raw mol to edit it
    mol = Chem.RWMol(mol)
    print(type(mol))
    for i in range(n_iterations):
        print(f'ITERATION {i}')
        atom_idx = choose_next_position(mol)
        print(f'The next atom will be added at atom with index {atom_idx}')
        mol = add_group(mol, 'C', atom_idx)
        print(atom_idx, type(atom_idx))
        print([n.GetSymbol() for n in mol.GetAtomWithIdx(int(atom_idx)).GetNeighbors()])
        # calculate valence states of all atoms
        mol.UpdatePropertyCache()
        mol = Chem.AddHs(mol)
        # mol = Chem.rdmolops.AddHs(mol,addCoords=1)
        print([n.GetSymbol() for n in list(mol.GetAtoms())[atom_idx].GetNeighbors()])
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        mol = Chem.RemoveHs(mol)

    write_mol_to_sdf(mol, './results/modified_compound.sdf')

def show_indexed_mol(mol):
    '''
    molecule is drawn to help the user to enter the right index of an atom
    '''
    for at in mol.GetAtoms():
        at.SetProp('atomNote',f'{at.GetIdx()}')
    rdkit.Chem.Draw.ShowMol(mol)

def main():
    smiles = 'C1=CC=C2C(=C1)C=CN2'
    mol = Chem.MolFromSmiles(smiles)
    # grow_ligand_at_random(smiles, 10)
    # fragments, linkers = load_libraries('data/fragment_library.txt', 'data/linker_library.txt')
    show_indexed_mol(mol)

if __name__ == '__main__':
    main()