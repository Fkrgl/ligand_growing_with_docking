from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem.Draw
import random, os, sys
import numpy as np
import pandas as pd
from anytree import AnyNode
import run_plants
from scipy.spatial.distance import cdist
import argparse
from tqdm import tqdm
from multiprocessing import Pool
from time import time
import pickle
import glob
import bisect
import shutil


PLANTS = ''
OUT_DIR = ''
BASE_FRAGMENT_NODE = None
TOP_SCORE = 0

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


def remove_bonded_hydrogen(mol, atom_idx):
    '''given an atom, the function deletes an existing bonded hydrogen and the bond '''
    mol = Chem.RWMol(mol)
    for neighbor in list(mol.GetAtoms())[atom_idx].GetNeighbors():
        if neighbor.GetSymbol() == 'H':
            mol.RemoveBond(neighbor.GetIdx(), int(atom_idx))
            mol.RemoveAtom(neighbor.GetIdx())
            return mol

def compute_3D_coordinates(mol):
    '''
    function calculates 3D coordinates for a molecule (this produces a conformer)
    '''
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def load_libraries(path_frag, path_link, has_index):
    """
    loads fragment and linker library
    :param path_frag: path to fragment lib
    :param path_link: path to linker lib
    :return: list of fragment and list of linker molecules (rdkit Mol)
    """
    fragments = []
    linkers = []
    # read in fragments
    if has_index:
        frag_lib = pd.read_csv(path_frag, sep=';', header=None)
        frag_lib[1] = frag_lib[1].fillna('999')
        index_col = frag_lib.apply(lambda x: list(map(int, list(x[1].split(',')))), axis=1)
        frag_lib[1] = index_col
        for i in range(len(frag_lib) - 1):
            try:
                mol, idx_array = frag_lib.iloc[i]
                mol = Chem.MolFromSmiles(mol)
                fragments.append([mol, idx_array])
            except:
                print(
                    f'SMILES {frag_lib.iloc[i].values[0]} cannot be translates into a molecule. Please check the SMILES')
    # if no fragment indices are supplied:
    else:
        frag_lib = pd.read_csv(path_frag, sep='\t', header=None)
        for i in range(len(frag_lib) - 1):
            try:
                mol = Chem.MolFromSmiles(frag_lib.iloc[i].values[0])
                fragments.append(mol)
            except:
                print(
                    f'SMILES {frag_lib.iloc[i].values[0]} cannot be translates into a molecule. Please check the SMILES')

    # read in linkers
    linker_lib = pd.read_csv(path_link, sep='\t', header=None)
    for i in range(len(linker_lib) - 1):
        try:
            mol = Chem.MolFromSmiles(linker_lib.iloc[i].values[0])
            linkers.append(mol)
        except:
            print(
                f'SMILES {linker_lib.iloc[i].values[0]} cannot be translates into a molecule. Please check the SMILES')

    return fragments, linkers

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
    linker_atom = mol.GetAtomWithIdx(linger_idx)
    linker_atom_neighbor = linker_atom.GetNeighbors()[0]
    return linker_atom, linker_atom_neighbor


def get_linker_atom_index(combo, fragment):
    '''
    Given combined but not yet bonded molecules, search for possible atoms of the fragment that can from a bond to the
    growing molecule
    :param combo: combined molecules in one object
    :param fragment: fragment (substituent) that is added to the linker
    :return: atom index of fragment atoms that can from a bond to the linker group
    '''

    elements = ['C', 'N']
    possible_linker_idx = list(combo.GetSubstructMatches(fragment[0]))
    # more than one substructure match detected
    if len(possible_linker_idx) > 1:
        # select match containing highest atom indices (this belongs to the yet unconnceted fragment)
        possible_linker_idx.sort(key=lambda x: x[0], reverse=True)
        possible_linker_idx = possible_linker_idx[0]
    # only one match
    else:
        possible_linker_idx = possible_linker_idx[0]

    # if no index list is given, consider all possible linker atoms of fragment
    if len(fragment) == 1:
        # select all carbon and nitrogens that have at least one hydrogen
        all_possible_linkers = []
        for linker_idx in possible_linker_idx:
            linker_atom = combo.GetAtomWithIdx(linker_idx)
            if linker_atom.GetSymbol() in elements and linker_atom.GetTotalNumHs() > 0:
                all_possible_linkers.append(linker_idx)
        return all_possible_linkers
    # if a linker list is given, return the corresponding indeces in the combo mol
    else:
        valid_positions = fragment[1]
        lowest_linker = sorted(possible_linker_idx)[0]
        valid_positions.sort()
        corresponding_index = [index + lowest_linker for index in valid_positions]
        # check if the highest corresponding index surpasses the combo atom index
        if corresponding_index[-1] > len(combo.GetAtoms()) - 1:
            return []
        else:
            return corresponding_index


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
    if mode == 'fragment':
        linker_atom_symbol = '[Hg]'
    elif mode == 'linker':
        linker_atom_symbol = '[Au]'
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

def add_all_linker_fragment_combinations(mol_node, grow_seed, linkers, fragments, node_id_parent, aromatic_idx_base,
                                         base_fragment_node, protein_coords, workers, index_list=None):
    '''
    function produces all possible base_fragment-linker-fragment combinations.
    :param mol_node: node that contains molecule to grow on
    :param grow_seed: index of atom in mol at which to continue growing
    :param linkers: list of all linkers
    :param fragments: list of all fragments
    :return: all possible combination as molecules and nodes
    '''
    bond_length = {'CO': 1.43,
                   'CC': 1.54,
                   'CN': 1.43}
    atoms_to_functionals = {'C': ['[Au]C(=O)O', '[Au]C'],
                            'O': ['[Au]O'],
                            'N': ['[Au]N']}
    nodes = []
    mol = mol_node.mol
    with Pool(workers) as pool:
        # go through all linkers
        for i, linker in enumerate(linkers):
            mol_linker = add_fragment(mol, linker, 'linker', atom_idx=grow_seed, bond_type=Chem.rdchem.BondType.SINGLE)
            # check if molecule linker bond worked
            if mol_linker:
                # go through all fragments
                for j, fragment in enumerate(fragments):
                    combo = Chem.CombineMols(mol_linker, fragment)
                    # search for linkers if no index list is specified
                    if index_list is None:
                        fragment_atom_idx = get_linker_atom_index(combo, [fragment])
                    else:
                        # only if there are specified indices (all intergers in list that are not 999)
                        if index_list[j][0] != 999:
                            fragment_atom_idx = get_linker_atom_index(combo, [fragment, index_list[j]])
                        else:
                            fragment_atom_idx = get_linker_atom_index(combo, [fragment])
                    # go through all fragment positions
                    for p, atom_idx in enumerate(fragment_atom_idx):
                        mol_linker_fragment = add_fragment(mol_linker, fragment, 'fragment', atom_idx=atom_idx,
                                                           bond_type=Chem.rdchem.BondType.SINGLE)
                        # check if fragment-linker bond worked
                        if mol_linker_fragment:
                            ### decorate ###########################
                            # 0)
                            # generate conformers for base fragment and molecule to be able to align them
                            # 1)
                            # align mol to base fragment -> get aligned mol
                            # 2)
                            # get all aromatic atoms that have at least one hydrogen
                            # 3)
                            # generate each decoration and save them as a mol

                            # get conformers
                            base_fragment_conformer = base_fragment_node.plants_pose
                            mol_linker_fragment_conformer = compute_3D_coordinates(mol_linker_fragment)
                            # align mol to base fragment
                            mol_linker_fragment_aligned = align_to_basefragment(mol_linker_fragment_conformer, base_fragment_conformer,
                                                                                aromatic_idx_base)
                            # generate all decorations
                            mol_linker_fragment_decorated = decorate_ligand(mol_linker_fragment_aligned, aromatic_idx_base,
                                                                            protein_coords, bond_length,
                                                                            atoms_to_functionals)
                            # add undecorated mol to list
                            mol_linker_fragment_decorated.append(mol_linker_fragment)
                            # check if each decorated mol is valid
                            for q, decorated_mol in enumerate(mol_linker_fragment_decorated):
                                if mol_linker_fragment_decorated:
                                    # node_id = [n_node]_[linker]_[fragment]_[fragment_atom_id]_[n_decoration]_[pose]
                                    node_id = f'{node_id_parent}_{i}_{j}_{p}_{q}_0'
                                    new_node = AnyNode(id=node_id, mol=decorated_mol, parent=mol_node, plants_pose=None,
                                                       score=None)
                                    pool.apply_async(write_node, (new_node, PLANTS + 'all_combinations/' +
                                                                  f'{node_id_parent}_{i}_{j}_{p}_{q}_0.pkl'))
            else:
                continue
        pool.close()
        pool.join()
        mol_node.children = nodes



def get_possible_grow_seeds(mol, base_fragment):
    '''
    grow seeds are atoms at which we continue the growing. Function searches the base fragment atom indices
    and exculdes them from the candidate atom index list for growing. Further, growing is only continued on specific
    types (possible_atoms) and with at least one hydrogen
    :param mol: grown molecule
    :param base_fragment: base fragment that is contained in mol
    :return: indices of possible grow seeds
    '''

    mol = Chem.RemoveHs(mol)
    base_fragment = Chem.RemoveHs(base_fragment)
    possible_atoms = ['C', 'N']
    mol_atom_idx = set([atom.GetIdx() for atom in mol.GetAtoms()])
    base_fragment_indices = get_base_fragment_indices(mol, base_fragment)
    mol_atom_idx.difference_update(base_fragment_indices)
    # only take carbons with at least one hydrogen as final grow seeds
    final_grow_seeds = []
    for grow_seed_idx in mol_atom_idx:
        grow_seed = mol.GetAtomWithIdx(grow_seed_idx)
        if grow_seed.GetSymbol() in possible_atoms and grow_seed.GetTotalNumHs() > 0:
            final_grow_seeds.append(grow_seed.GetIdx())
    return final_grow_seeds

def get_next_grow_seed(possible_grow_seeds):
    '''
    function chooses the next grow seed under all candidates. At the moment this is done at random.
    '''
    return possible_grow_seeds[random.randint(0, len(possible_grow_seeds)-1)]


def grow_molecule(n_grow_iter, initial_grow_seed, linkers, fragments, aromatic_atom_idx, protein_coords,
                  workers, keep_top=100):
    '''
    function performs n rounds of ligand growing. In each round, each linker/fragment combination is added to each
    possible atom in the current leafs (excluding the base fragment).
    :param mol_tree: molecular tree that saves all grown molecules
    :param n_grow_iter: number of grow iterations
    :param initial_grow_seed: atom index of base fragment
    :param linkers: list with all linkers
    :param fragments: list with all fragments
    '''
    global BASE_FRAGMENT_NODE
    high_scoring_nodes = []
    frag = fragments
    index_list = None
    if type(fragments[0]) == list:
        frag = []
        index_list = []
        for i in range(len(fragments)):
            frag.append(fragments[i][0])
            index_list.append(fragments[i][1])
    # begin grow with base_fragment
    k = 0
    base_fragment_node = BASE_FRAGMENT_NODE
    # write basefragment to file
    base_node_path = PLANTS + 'current_nodes/root.pkl'
    write_node(base_fragment_node, base_node_path)
    # choose best pose for base fragment
    choose_best_initial_pose(base_fragment_node)
    base_fragment = base_fragment_node.plants_pose
    base_fragment_node.mol = label_base_fragment(base_fragment_node.mol)
    for i in range(n_grow_iter):
        print(f'ITERATION {i+1}/{n_grow_iter}\n')
        print('construct all combinations ...')
        # grow each leaf molecule by all combinations
        for leaf_path in glob.glob(PLANTS + 'current_nodes/*'):
            leaf = load_node(leaf_path)
            # select grow seed, grow mol on seed by one linker_fragment combo and save results in tree
            if i == 0:
                possible_grow_seeds = [initial_grow_seed]
            else:
                possible_grow_seeds = get_possible_grow_seeds(leaf.mol, base_fragment_node.mol)
                # grow on each possible grow seed on the current leaf
            if len(possible_grow_seeds) == 0:
                continue
            for grow_seed in possible_grow_seeds:
                # grow all combinations for the current leaf
                add_all_linker_fragment_combinations(leaf, grow_seed, linkers, frag, k, aromatic_atom_idx,
                                                     base_fragment_node, protein_coords, workers, index_list)
                k += 1
            # remove the leaf from the current_leaf folder (it is now the parent of new child leafs)
            os.remove(leaf_path)
        # dock all grown molecules from this iteration and add additional poses as nodes
        print(f'Docking of iteration {i+1} is running ...')
        # parallelized docking
        pathname = PLANTS + 'all_combinations/*.pkl'
        n_passed_filter = 0
        with Pool(processes=workers) as pool:
            for passed_poses in tqdm(pool.imap_unordered(func=dock_leafs_parallel, iterable=glob.glob(pathname), chunksize=1),
                             total=len(glob.glob(pathname))):
                n_passed_filter += passed_poses
        if n_passed_filter == 0:
            # later: if no improvement in one round is made, try to continue with another node that has a good score.
            print('\nNo score improvement could be achieved.\n')
            return

        high_scoring_nodes = update_top_nodes(high_scoring_nodes)
        # #reset current leafs dir
        # if os.path.exists(PLANTS + 'current_nodes/'):
        #     shutil.rmtree(PLANTS + 'current_nodes/')
        # os.mkdir(PLANTS + 'current_nodes/')
        #
        # # get top 25 nodes
        # for leaf in top_leafs:
        #     leaf = leaf[1]
        #     write_node(leaf, PLANTS + f'current_nodes/{leaf.id}.pkl')
        #write best poses of this round in dir
        write_best_poses_to_file(high_scoring_nodes, i+1)

    return high_scoring_nodes


def update_top_nodes(high_scoring_nodes, keep_top=100):
    global TOP_SCORE
    pathname = PLANTS + 'current_nodes/*'
    for file_path in glob.glob(pathname):
        # add filtered nodes for next iteration
        with open(file_path, 'rb') as f:
            node = pickle.load(f)
            try:
                if len(high_scoring_nodes) < keep_top:
                    bisect.insort(high_scoring_nodes, (node.score, node))
                else:
                    # new node has lower score
                    try:
                        if high_scoring_nodes[-1][0] > node.score:
                            bisect.insort(high_scoring_nodes, (node.score, node))
                            del high_scoring_nodes[-1]
                    except:
                        print(f'insort failed: {type(high_scoring_nodes[-1][0])}{high_scoring_nodes[-1][0]}, {type(node.score)}{node.score}')
            except:
                print(f'not supported between instances of AnyNode and AnyNode: {type(len(high_scoring_nodes))},'
                      f' {type(keep_top)}' )
    TOP_SCORE = high_scoring_nodes[0][0]
    return high_scoring_nodes

def filter_top_leafs():
    pathname = PLANTS + 'current_nodes/*'
    top_leafs = []
    for file_path in glob.glob(pathname):
        with open(file_path, 'rb') as f:
            node = pickle.load(f)
def load_node(node_path):
    with open(node_path, 'rb') as f:
        return pickle.load(f)

def write_node(node, node_path):
    with open(node_path, 'wb') as f:
        pickle.dump(node, f)

def reset_all_folders():
    'function resets each folder to default'
    all_combinations = PLANTS + 'all_combinations'
    current_nodes = PLANTS + 'current_nodes'
    parallel_output = PLANTS + 'parallel_output'
    tmp = PLANTS + 'tmp'

    if os.path.exists(all_combinations):
        shutil.rmtree(all_combinations)
    os.mkdir(all_combinations)
    if os.path.exists(current_nodes):
        shutil.rmtree(current_nodes)
    os.mkdir(current_nodes)
    if os.path.exists(parallel_output):
        shutil.rmtree(parallel_output)
    os.mkdir(parallel_output)
    if os.path.exists(tmp):
        shutil.rmtree(tmp)
    os.mkdir(tmp)
    # remove all result folders
    for dir in glob.glob(PLANTS + 'best_poses_iter_*'):
        shutil.rmtree(dir)

def choose_best_initial_pose(base_fragment_node, cutoff=1.5):
    '''
    docks the crystal structure and finds under the initial docking poses the pose with the lowest RMSD
    '''
    global TOP_SCORE
    crystal_structure = base_fragment_node.mol
    initial_poses = None
    scores = None
    # three attempts to reach RMSD threshold
    for i in range(3):
        initial_poses, scores = run_plants.dock_molecule(crystal_structure, PLANTS)
        if min(scores) <= cutoff:
            rmsd = [calc_RMSD(crystal_structure, initial_pose) for initial_pose in initial_poses]
            idx = np.argmin(rmsd)
            # set pose and score for root node
            base_fragment_node.mol = initial_poses[idx]
            base_fragment_node.plants_pose = initial_poses[idx]
            base_fragment_node.score = scores[idx]
            TOP_SCORE = scores[idx]
            return
    # abort if threshold is not reached
    sys.exit(f'Crystal structure and docking pose of crystal structure deviate to much (RMSD > {cutoff})')


def dock_leafs_parallel(node_path):
    '''
    function performs docking for one leaf node. It is designed to work in parallel
    :param leaf_nodes: node that contains one of the current grown molecules
    '''
    passed_poses = 0
    # load node object
    leaf_node = load_node(node_path)
    poses, scores = run_plants.dock_molecule_parallel(leaf_node, PLANTS)
    # delete parent node and insert all poses as nodes
    parent = leaf_node.parent
    identifier = leaf_node.id
    mol = leaf_node.mol
    for i in range(len(poses)):
        node_id = identifier.split('_')
        node_id[-1] = str(i)
        node_id = '_'.join(node_id)
        pose_node = AnyNode(id=node_id, mol=mol, parent=parent, plants_pose=poses[i],
                            score=scores[i])
        # check for filters
        if pass_filter(pose_node):
            passed_poses += 1
            # write to file
            with open(PLANTS + f'current_nodes/{node_id}.pkl', 'wb') as pckl_file:
                pickle.dump(pose_node, pckl_file)
    # remove old node
    os.remove(node_path)
    return passed_poses

def write_best_poses_to_file(high_scoring_nodes, iter):
    '''
    writes the docking poses of the highest scoring grown molecules in the molecular tree into a file and outputs their
    molecular graph as pdf
    '''
    if len(high_scoring_nodes) == 0:
        return
    else:
        print('write best poses ...\n')
        path = PLANTS + f'best_poses_iter_{iter}/'
        ranking_file = path + 'ranking.txt'
        # check if dir exists
        if os.path.exists(path):
            shutil.rmtree(path)
        os.mkdir(path)
        # save poses to sdf
        best_poses = high_scoring_nodes
        mols = []
        legend = []
        with open(ranking_file, 'w') as f:
            f.write(f'{0}\t{BASE_FRAGMENT_NODE.id}\t{BASE_FRAGMENT_NODE.score:.4f}\n')
            for i, p in enumerate(best_poses):
                node = p[1]
                pose = node.plants_pose
                filename = str(node.id) + '.sdf'
                run_plants.write_mol_to_sdf(pose, path + filename)
                f.write(f'{i+1}\t{node.id}\t{node.score:.4f}\n')
                node.mol.Compute2DCoords()
                mols.append(Chem.RemoveHs(node.mol))
                legend.append(f'id: {node.id}\nscore: {node.score:.4f}')
        # write mol graph
        pdf_path = path + 'mols.pdf'
        write_grown_molecules_to_pdf(pdf_path, mols, legend=legend)

def get_base_fragment_indices(mol, base_fragment):
    '''
    searches for the base fragment in the mol and returns its atom indeces
    '''
    match_idx = list(mol.GetSubstructMatches(base_fragment))
    # fragment with lowest atom indeces is assumed to be base fragment
    match_idx.sort(key=lambda x: x[0])
    match_idx = match_idx[0]
    return match_idx

def calc_RMSD(base_fragment, grown_mol):
    '''
    calculates RMSD between the atoms of two molecules specified by index
    :param base_fragment
    :param grown_mol
    :param substructure_idx: indeces of atoms that belong to the base fragment in the grown molecule
    :return:
    '''
    substructure_idx = get_base_fragment_indices(grown_mol, base_fragment)
    base_fragment = Chem.RemoveHs(base_fragment)
    grown_mol = Chem.RemoveHs(grown_mol)
    n = len(base_fragment.GetAtoms())
    rmsd = 0
    for i in substructure_idx:
        p_base = base_fragment.GetConformer().GetAtomPosition(i)
        p_grow = grown_mol.GetConformer().GetAtomPosition(i)
        rmsd += p_base.Distance(p_grow)**2
    rmsd = np.sqrt(1/n * rmsd)
    return rmsd

def pass_filter(leaf_node, cut_off=5):
    '''
    function returns only true if:
    1) RMSD with the base fragment below a certain threshold
    2) a equal or better score than the base fragment
    :return: filtered leafs
    '''

    base_fragment = BASE_FRAGMENT_NODE.plants_pose
    base_fragment_score = BASE_FRAGMENT_NODE.score
    if leaf_node.score <= base_fragment_score:
        rmsd = calc_RMSD(base_fragment, leaf_node.plants_pose)
        if rmsd <= cut_off:
            return True
    return False


def read_protein_coords(protein_mol2_path):
    '''
    function reads in the protein and returns its atomic coordinates, excluding hydrogens.
    '''
    out_path = run_plants.convert_mol_files(protein_mol2_path, 'mol2', 'sdf')
    protein = Chem.MolFromMolFile(out_path, sanitize=False)
    protein = Chem.RemoveHs(protein, sanitize=False)
    protein_coords = protein.GetConformer().GetPositions()
    return protein_coords

def write_grown_molecules_to_pdf(file_path, mols, molsPerRow=4, legend=None):
    img = rdkit.Chem.Draw.MolsToGridImage(mols, molsPerRow=molsPerRow, legends=legend, returnPNG=False)
    img.save(file_path)

# ============================================= decorate aromatic rings  ============================================= #

def isRingAromatic(mol, bondRing):
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
    retruns the indices of all atoms of aromatic ring systems with at least one hydrogen in the molecule
    '''
    aromatic_ring_atom_idx = set()
    rings = mol.GetRingInfo()
    # check all bonds for aromaticity
    for ring in rings.BondRings():
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
    :param mol: grown molecule with 3D conformation
    :param base_fragment: base fragment with 3D conformation
    :param aromatic_ring_idx: indices to decorate on
    :return:
    '''
    # remove hydrogens
    mol = rdkit.Chem.RemoveHs(mol)
    base_fragment = rdkit.Chem.RemoveHs(base_fragment)
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

def calc_group_position(mol, idx, length):
    '''
    given an atom index the function calculates the position of an atom (part of the decoration group) that
    is hypotheically substitued for the hydrogen at the atom with a specified bond length
    '''
    # add hs with coords
    mol = rdkit.Chem.rdmolops.AddHs(mol, addCoords=True)
    pos_atom = mol.GetConformer().GetAtomPosition(idx)
    # substitute one of the H atoms
    for neighbor in mol.GetAtomWithIdx(idx).GetNeighbors():
        if neighbor.GetSymbol() == 'H':
            neighbor_idx = neighbor.GetIdx()
            pos_h = mol.GetConformer().GetAtomPosition(neighbor_idx)
            # get dircetion vector of the H atom
            h_vector = pos_atom.DirectionVector(pos_h)
            # get position of new group
            group_pos = pos_atom + h_vector*length
            mol.GetConformer().SetAtomPosition(neighbor_idx, group_pos)
            return group_pos
    return


def has_spacial_neighbors(atom_pos, protein_coords, min_d=2):
    '''
    function returns True if any protein atom has a distance smaller than d to the specified atom of the ligand
    '''
    atom_pos = np.asarray(atom_pos).reshape(1, -1)
    distances = cdist(protein_coords, atom_pos, 'euclidean')
    min_dist = min(distances)[0]
    if min_dist < min_d:
        return True
    else:
        return False

def add_functional_group(mol, pos, substituent_smiles, linker_atom_symbol, bond_type=Chem.rdchem.BondType.SINGLE):
    '''
    function merges a functional group to a molecle. The bond is built betwenn the nieghbor of the Au atom in
    the functional group and the specified atom in the molecule
    '''
    mol = rdkit.Chem.rdmolops.AddHs(mol, addCoords=True)
    pos_atom = mol.GetConformer().GetAtomPosition(pos)
    # substitute one of the H atoms
    for neighbor in mol.GetAtomWithIdx(pos).GetNeighbors():
        if neighbor.GetSymbol() == 'H':
            neigbor_atom = neighbor
            neighbor_idx = neighbor.GetIdx()
            mol.Compute2DCoords()
            substituent = Chem.MolFromSmiles(substituent_smiles)
            substituent.Compute2DCoords()
            combo = Chem.CombineMols(mol, substituent)
            linker_atom, linker_atom_neighbor = get_linker_atom_and_neighbor(combo, linker_atom_symbol)
            # make combo editable
            edcombo = Chem.EditableMol(combo)
            # add bond between linker and neighbor of au atom
            edcombo.AddBond(pos, linker_atom_neighbor.GetIdx(), order=bond_type)
            # remove AU and its bond to the linker
            mol = edcombo.GetMol()
            mol = Chem.RWMol(mol)
            mol.RemoveBond(linker_atom.GetIdx(), linker_atom_neighbor.GetIdx())
            mol.RemoveAtom(linker_atom.GetIdx())
            # search neighbor again, index is confused somehow (the neighbor from the for loop is no longer the right
            # atom)
            for n in mol.GetAtomWithIdx(pos).GetNeighbors():
                if n.GetSymbol() == 'H':
                    neighbor = n
            # remove hydrogen and its bond to the fragment
            mol.RemoveBond(pos, neighbor.GetIdx())
            mol.RemoveAtom(neighbor.GetIdx())
            mol = compute_3D_coordinates(mol)
            mol.Compute2DCoords()
            # check if mol is valid
            mol = check_molecule(mol)
            if mol:
                return mol
            else:
                return None

def hasHydogenNeighbor(mol, idx):
    '''
    function checks if an atom has a hydrogen as neighbor
    '''
    for neighbor in mol.GetAtomWithIdx(idx).GetNeighbors():
        if neighbor.GetSymbol == 'H':
            return True
    return False

def decorate_ligand(mol, aromatic_atom_idx, protein_coords, bond_length, atoms_to_functionals):
    '''
    function adds functional groups to all aromatic ring atoms if they do not clash with any protein atoms
    '''
    mol_with_functionals = []
    mol = Chem.AddHs(mol)

    group_atoms = ['C', 'O', 'N']
    for idx in aromatic_atom_idx:
        # check if atom has hydrogen as neighbor
        neighbors = [neighbor.GetSymbol() for neighbor in mol.GetAtomWithIdx(idx).GetNeighbors()]
        if 'H' in neighbors:
            for group in group_atoms:
                group_pos = calc_group_position(mol, idx, bond_length['C' + group])
                # check if coords of group would clash
                if has_spacial_neighbors(group_pos, protein_coords):
                    continue
                # no clashes: add functional groups to mol
                else:
                    # go through each functional group with this linker atom
                    for functional in atoms_to_functionals[group]:
                        # give max_attempts attempts to add functional group
                        attempt = 0
                        max_attempts = 10
                        while attempt < max_attempts:
                            decorated_mol = add_functional_group(mol, idx, functional, '[Au]')
                            # check if something was returned
                            if decorated_mol:
                                mol_with_functionals.append(decorated_mol)
                                break
                            attempt += 1
    return mol_with_functionals

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('l', metavar='-ligand', type=str, help='path to mol2 file for ligand')
    parser.add_argument('p', metavar='-protein', type=str, help='path to mol2 for protein')
    parser.add_argument('frag', metavar='-fragment_lib', type=str, help='path to fragment library')
    parser.add_argument('link', metavar='-linker_lib', type=str, help='path to linker library')
    parser.add_argument('i', metavar='-iter', type=int, help='number of growing iterations')
    parser.add_argument('c', metavar='-cores', type=int, help='number of cores to use')
    parser.add_argument('P', metavar='-PLANTS', type=str, help='Path to PLANTS directory')
    parser.add_argument('o', metavar='-out_dir', type=str, help='Path to output dir where result directory '
                                                                '\'grown_molecules\' and ranking.txt are located')
    parser.add_argument('--index', action='store_true', help='if specified, the program uses pre-defined indices to '
                                                             'find linker atoms of the fragment. This requires a frag-'
                                                             'ment library file with index array in the second column')
    parser.add_argument('--no-index', dest='feature', action='store_false')
    parser.set_defaults(feature=True)
    args = parser.parse_args()
    return args.l, args.p, args.frag, args.link, args.i, args.c, args.P, args.o, args.index
# ===================================================== MAIN  ======================================================== #

def main():
    # read in ligand, protein and the PLANTS path
    global PLANTS, OUT_DIR, BASE_FRAGMENT_NODE
    ligand_mol2_path, protein_mol2_path, frag_path, linker_path, iter, cores, plants, out_dir, has_index = get_arguments()
    PLANTS = plants
    OUT_DIR = out_dir
    mol = run_plants.get_mol_from_mol2file(ligand_mol2_path)
    show_indexed_mol(mol)
    protein_coords = read_protein_coords(protein_mol2_path)
    mol = label_base_fragment(mol)
    aromatic_atom_idx = get_aromatic_rings(mol)
    # read in fragment and linker libraries
    fragments, linkers = load_libraries(frag_path, linker_path, has_index)
    # initiate root node
    root = AnyNode(id='root', mol=mol, parent=None, plants_pose=None, score=None)
    BASE_FRAGMENT_NODE = root
    start_time = time()
    reset_all_folders()
    # start growing
    grow_molecule(iter, 0, [linkers[0]], [fragments[4]], aromatic_atom_idx, protein_coords, workers=cores)
    print(f'Runtime: {time()-start_time:.2f}')

if __name__ == '__main__':
    main()