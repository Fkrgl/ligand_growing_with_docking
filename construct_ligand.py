from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem.Draw
import random
from collections import Counter
import numpy as np
from openbabel import openbabel
import pandas as pd
from anytree import AnyNode

def label_base_fragment(mol):
    '''
    atoms of the base fragment are labels by a property to later decern the base fragment from similar fragments
    mol: base fragment
    '''
    for atom in mol.GetAtoms():
        atom.SetProp('base_fragment', 'True')
    return mol

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
    elements = ['C', 'N']
    possible_linker_idx = list(combo.GetSubstructMatches(fragment))
    # more that one substructure match detected
    if len(possible_linker_idx) > 1:
        # select match containing highest atom indices (this belongs to the yet unconnceted fragment)
        possible_linker_idx.sort(key=lambda x: x[0], reverse=True)
        possible_linker_idx = possible_linker_idx[0]
    # only one match
    else:
        possible_linker_idx = possible_linker_idx[0]
    # select the first atom that is a carbon and has at least one hydrogen
    for linker_idx in possible_linker_idx:
        linker_atom = combo.GetAtomWithIdx(linker_idx)
        if linker_atom.GetSymbol() in elements and linker_atom.GetTotalNumHs() > 0:
            return linker_idx
    return None

def add_fragment(mol, fragment, mode, atom_idx=None, bond_type=Chem.rdchem.BondType.SINGLE):
    '''
    function adds a fragment (linker or substituent) to the growing molecule
    :param mol: growing mol
    :param fragment: mol
    :param mode: 'linker' or 'fragment'. Determines if we have a fragment linker or a linker fragment
    :param bond_type: Chem.rdchem.BondType
    :return: bonded fragment to molecule
    '''

    # combine both molecules into one
    combo = Chem.CombineMols(mol, fragment)
    # select the linker atom to fuse molecules
    linker_atom_symbol = '[Au]'
    if mode == 'fragment':
        # for fragments, also select a atom in the fragemnt for linkage
        linker_atom_symbol = '[Hg]'
        atom_idx = get_linker_atom_index(combo, fragment)
        if atom_idx == None:
            show_indexed_mol(combo)
            show_indexed_mol(fragment)
            return
    # for linker, check if a atom index of the mol is specified
    elif mode == 'linker':
        if atom_idx is None:
            print('no atom index of the molecule is specified!')
            return
    else:
        print('mode does not exist. Please select \'linker\' or \'fragment\' as mode.')
        return
    # get linker atom and its neighbor
    linker_atom, linker_atom_neighbor = get_linker_atom_and_neighbor(combo, linker_atom_symbol)
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

def add_all_linker_fragment_combinations(mol, grow_seed, linkers, fragments, parent):
    '''
    function produces all possible base_fragment-linker-fragment combinations.
    :param mol: mol to continue growing on
    :param grow_seed: index of atom in mol at which to continue growing
    :param linkers: list of all linkers
    :param fragments: list of all fragments
    :param parent: node object of mol
    :return: all possible combination as molecules
    '''
    grown_mols = []
    nodes = []
    print(f'number of linkers: {len(linkers)}')
    print(f'number of fragments: {len(fragments)}')
    for i, linker in enumerate(linkers):
        mol_linker = add_fragment(mol, linker, 'linker', atom_idx=grow_seed, bond_type=Chem.rdchem.BondType.SINGLE)
        # check if molecule linker bond worked
        if mol_linker:
            for j, fragment in enumerate(fragments):
                mol_linker_fragment = add_fragment(mol_linker, fragment, 'fragment',
                                                   bond_type=Chem.rdchem.BondType.SINGLE)
                print(mol_linker_fragment)
                if mol_linker_fragment:
                    grown_mols.append(mol_linker_fragment)
                    new_node = AnyNode(mol=mol_linker_fragment, parent=parent)
                    nodes.append(new_node)
        else:
            continue
    parent.children = nodes
    return grown_mols

def get_possible_grow_seeds(mol, base_fragment):
    '''
    grow seeds are atoms at which we continue the growing. Function searches the base fragment atom indices
    and exculdes them from the candidate atom index list for growing. Further, growing is only continued on specific
    types (possible_atoms) and with at least one hydrogen
    :param mol: grown molecule
    :param base_fragment: base fragment that is contained in mol
    :return: indices of possible grow seeds
    '''
    possible_atoms = ['C', 'N']
    mol_atom_idx = set([atom.GetIdx() for atom in mol.GetAtoms()])
    match_idx = mol.GetSubstructMatches(base_fragment)
    for candidates in match_idx:
        sub_structure_atom = mol.GetAtomWithIdx(candidates[0])
        # check if atom is part of the base fragment
        if sub_structure_atom.HasProp('base_fragment'):
            mol_atom_idx.difference_update(candidates)
            # only take carbons with at least one hydrogen as final grow seeds
            final_grow_seeds = []
            for grow_seed_idx in mol_atom_idx:
                grow_seed = mol.GetAtomWithIdx(grow_seed_idx)
                if grow_seed.GetSymbol() in possible_atoms and grow_seed.GetTotalNumHs() > 0:
                    final_grow_seeds.append(grow_seed.GetIdx())
            return final_grow_seeds

    return None

def get_next_grow_seed(possible_grow_seeds):
    '''
    function chooses the next grow seed under all candidates. At the moment this is done at random.
    '''
    return possible_grow_seeds[random.randint(0, len(possible_grow_seeds)-1)]


def grow_molecule(n_grow_iter, base_fragment, initial_grow_seed, linkers, fragments):
    '''
    function performs n rounds of ligand growing. In each round, each linker/fragment combination is added to each of
    the current leafs.
    :param n_grow_iter: number of grow iterations
    :param base_fragment: start molecule
    :param initial_grow_seed: atom index of base fragment
    :param linkers: list with all linkers
    :param fragments: list with all fragments
    :return: grown molecules (leafs of the mol tree)
    '''

    # current molecules to continue growing (leafs of the tree)
    leafs = [base_fragment]
    k = 0
    for i in range(n_grow_iter):
        print(f'in iteration {i} we have {len(leafs)} leafs')
        grown_mols = []
        current_leafs = []
        for leaf in leafs:
            # select grow seed, grow mol on seed by one linker, fragment combo and save results in tree
            if i == 0:
                leaf_node = AnyNode(id='root',mol=leaf)
                grow_seed = initial_grow_seed
            else:
                leaf_node = AnyNode(id=f'{k}',mol=leaf)
                possible_grow_seeds = get_possible_grow_seeds(leaf, base_fragment)
                grow_seed = get_next_grow_seed(possible_grow_seeds)
            # grow all combinations for the current leaf
            grown_mols += add_all_linker_fragment_combinations(leaf, grow_seed, linkers, fragments, leaf_node)
            k += 1
        leafs = grown_mols
    print(f'end leafs {len(leafs)}')
    return leafs

def main():
    smiles = 'C1=CC=C2C(=C1)C=CN2'
    #smiles = 'c1cc(CCCO)ccc1'
    mol = Chem.MolFromSmiles(smiles)
    print(mol)
    mol = label_base_fragment(mol)
    # grow_ligand_at_random(smiles, 10)
    fragments, linkers = load_libraries('data/fragment_library.txt', 'data/linker_library.txt')
    grows = grow_molecule(2, mol, 1, linkers[:2], fragments[:2])
    print(grows)
    print(len(grows))
    # rdkit.Chem.Draw.MolsToGridImage(grows[-450:])
    # print(rdkit.Chem.Draw.MolsToGridImage(grows[-450:]))

if __name__ == '__main__':
    main()