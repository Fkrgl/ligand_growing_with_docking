import bisect
"""
This class implemets a tree data structure that
1) saves all molecules as nodes
2) keeps track of the current leafs of the tree
"""

class Mol_Tree:
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

    def insert_leafs(self, leafs):
        if type(leafs) != list:
            leafs = [leafs]
        self.leafs += leafs

    def clear_leafs(self):
        self.leafs = []

    def get_root(self):
        return self.root

    def get_nodes(self):
        return self.nodes

    def get_leafs(self):
        return self.leafs

    def get_best_poses(self, n_best):
        best_poses = sorted(self.nodes, key=lambda tup: tup[0])
        return best_poses[:n_best]


