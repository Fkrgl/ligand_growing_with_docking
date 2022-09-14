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

def check_molecule(mol):
    '''
    molecule is santizised to check if it is a valid molecule
    :param mol: mol to check
    '''
    status = Chem.rdmolops.SanitizeMol(mol, catchErrors=True)
    if status == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE:
        return mol
    else:
        print(status)
        return None

def get_linker_atom_and_neighbor(mol, linker_atom_smiles):
    '''
    searches the linker atom (Au or Hg) and its neighbor. The neighbor will be later linked to an atom of the
    growing molecule
    :param mol: growing molecule
    :param linker_atom_smiles: Au or Hg linger atom
    :return: linker atom and its neighbor
    '''
    linger_idx = mol.GetSubstructMatch(Chem.MolFromSmiles(linker_atom_smiles))[0]
    print(f'match for {linker_atom_smiles}: {linger_idx}')
    linker_atom = mol.GetAtomWithIdx(linger_idx)
    linker_atom_neighbor = linker_atom.GetNeighbors()[0]
    return linker_atom, linker_atom_neighbor

def get_linker_atom_index(combo, fragment):
    '''
    Given combined but not yet bonded molecules, search for possible atoms of the fragment that can from a bond to the
    growing molecule
    :param combo: combined molecules in one object
    :param fragment: fragment (substituent) that is added to the linker
    :return: atom index of fragment atom that can from a bond to the linker group
    '''
    # because we add the fragment, it will always have the atoms with the highest index
    # this function has to be adapted for the case that more than one substructure match is found
    possible_linker_idx = combo.GetSubstructMatch(fragment)
    print(possible_linker_idx)
    print(len(possible_linker_idx))
    # select the first atom that is a carbon and has at least one hydrogen
    for linker_idx in possible_linker_idx:
        linker_atom = combo.GetAtomWithIdx(linker_idx)
        print(linker_atom)
        if linker_atom.GetSymbol() == 'C' and linker_atom.GetTotalNumHs() > 0:
            return linker_idx
    return None

def add_fragment(mol, fragment, mode, bond_type=Chem.rdchem.BondType.SINGLE):
    '''
    function adds a fragment (linker or substituent) to the growing molecule
    :param mol: growing mol
    :param fragment: mol
    :param mode: 'linker' or 'fragment'. Determines if we have a fragment linker or a linker fragment
    :param bond_type: Chem.rdchem.BondType
    :return: bonded fragment to molecule
    '''

    # select the linker atom to fuse molecules
    linker_atom_symbol = '[Au]'
    if mode == 'fragment':
        linker_atom_symbol = '[Hg]'
        # if mode == linker, atom_index is for the molecule the linker is added to, if mode == fragment
        # the atom index belongs to an fragment atom
        # a function to select an atom index is still needed
    # combine both molecules into one
    combo = Chem.CombineMols(mol,fragment)
    atom_idx = get_linker_atom_index(combo, fragment)
    print(f'linker atom index: {atom_idx}')
    # show_indexed_mol(combo)
    # get linker atom and its neighbor
    linker_atom, linker_atom_neighbor = get_linker_atom_and_neighbor(combo, linker_atom_symbol)
    print(f'linker atom index: {linker_atom.GetIdx()}\nneighbor index: {linker_atom_neighbor.GetIdx()}')
    # make combo editable
    edcombo = Chem.EditableMol(combo)
    # add bond between linker and neighbor of au atom
    edcombo.AddBond(atom_idx, linker_atom_neighbor.GetIdx(), order=bond_type)
    # remove AU and its bond to the linker
    mol = edcombo.GetMol()
    mol = Chem.RWMol(mol)
    mol.RemoveBond(linker_atom.GetIdx(), linker_atom_neighbor.GetIdx())
    mol.RemoveAtom(linker_atom.GetIdx())
    # check if mol is valid
    mol = check_molecule(mol)
    if mol:
        return mol
    else:
        return None

def main():
    smiles = 'C1=CC=C2C(=C1)C=CN2'
    mol = Chem.MolFromSmiles(smiles)
    # grow_ligand_at_random(smiles, 10)
    # fragments, linkers = load_libraries('data/fragment_library.txt', 'data/linker_library.txt')
    show_indexed_mol(mol)

if __name__ == '__main__':
    main()