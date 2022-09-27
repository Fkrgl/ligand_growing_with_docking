from rdkit import Chem
from rdkit.Chem import AllChem
import os
from openbabel import openbabel
import re
import pandas as pd
from anytree import AnyNode

def alter_node(node):
    node.id = '5'

node = AnyNode(id='0')
print(node)
alter_node(node)
print(node)



