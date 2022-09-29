def __init__(self, root):
    self.root = root
    self.nodes = []
    self.leafs = []
    self.nodes.append(root)
    self.leafs.append(root)


def insert_node(self, node):
    if type(node) != list:
        node = [node]
    for n in node:
        self.nodes.append((n.score, n))
    print(f'list has now length {len(self.nodes)}')

def get_best_poses(self, n_best):
    best_poses = sorted(self.nodes, key=lambda tup: tup[0])
    return best_poses[:n_best]

def write_best_poses_to_file(mol_tree):
    '''
    writes the docking poses of the highest scoring grown molecules in the molecular tree into a file
    '''
    print('write best poses ...\n')
    path = OUT_DIR + 'grown_molecules/'
    ranking_file = OUT_DIR + 'ranking.txt'
    # check if dir is empty
    if len(os.listdir(path)) > 0:
        for f in os.listdir(path):
            os.remove(os.path.join(path, f))
    # save poses to sdf
    best_poses = mol_tree.get_best_poses(n_best=100)
    with open(ranking_file, 'w') as f:
        for i, p in enumerate(best_poses):
            node = p[1]
            pose = node.plants_pose
            filename = str(node.id) + '.sdf'
            run_plants.write_mol_to_sdf(pose, path + filename)
            f.write(f'{i}\t{node.id}\t{node.score:.4f}\n')
