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
        self.best_poses = []
        self.nodes.append(root)
        self.leafs.append(root)

    def insert_node(self, node):
        if type(node) != list:
            node = [node]
        self.nodes += node
        # insert nodes in best_poses
        for n in node:
            if len(self.best_poses) < 100:
                print(self.best_poses)
                print((n.score, n))
                bisect.insort(self.best_poses, (n.score, n))
            else:
                # only sort in if score is better than score of last element
                if n.score < self.best_poses[-1][0]:
                    bisect.insort(self.best_poses, (n.score, n))
                    del self.best_poses[-1]

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

    def get_best_poses(self):
        return self.best_poses


