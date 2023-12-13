import os, sys
import string
import pickle
import argparse
import warnings

import numpy as np
import networkx as nx
from binarytree import tree, Node, build
from itertools import chain

import Dev_tree as dt 

sys.path.append('../../')
from python import helper as hp

sys.path.append('../')
from circuit_design import Verification as ver
from circuit_design import Register as re


parser = argparse.ArgumentParser()
parser.add_argument('path_graph_shape', type=str, help='Path to the graph of the shape')
parser.add_argument('path_shape', type=str, help='Path to the shape, in txt format')
parser.add_argument('type_mesh', type=str, help='Type of mesh')
parser.add_argument('save_path', type=str, help='Path to save the output figure and the register matrix')
parser.add_argument('verbose', type=bool, help='to print the verification process step by step or not')
args = parser.parse_args()

r = re.Register()

shape = hp.load_pickle(args.path_shape)

# the way michael's shape are saved makes duplicate
# to show connection from one shape to another, thus
# we only take the even element.
# If you do it differently, update the line below
# accordingly.
unique_cells = [x for i,x in enumerate(shape) if i%2==0]

# load the shape with the networkx library
graph = nx.read_gpickle(args.path_graph_shape)

# get the number of blocks
id_block = hp.column(shape,0)
n_blocks = np.amax(id_block)


# get the number of cells in fct
# of the building block
if args.type_mesh == 'tetrahedron':
    n_per_block = 4
elif args.type_mesh == 'rhombus':
    n_per_block = 8
    
n_cells = n_blocks * n_per_block

# get the number of division needed
# to have the right number of cells
power2_list = [2**n for n in list(range(10))]

def get_n_division(n_cells, power2_list):
    temp = np.array(power2_list)
    val = temp[temp > n_cells].min()
    n_div = np.where(temp==val)
    
    return n_div[0][0]

n_div = get_n_division(n_cells, power2_list)


# get the size of the binary tree
def get_size_tree(n_cells, power2_list):
    temp = np.array(power2_list)
    val = temp[temp > n_cells].min()
    
    return np.sum([x for x in power2_list if power2_list.index(x) <= power2_list.index(val)])

tree_size = get_size_tree(n_cells, power2_list)

# create the tree
tree = build(list(range(tree_size)))

# define a cell class to hold
# all the information we'll
# deal with in the dev tree algo
class Cell():
    def __init__(self, id_tree, id_block, coordinate, in_connectors, ex_connectors):
        self.id_tree = id_tree
        self.id_blocks = id_block
        self.in_connectors = in_connectors
        self.ex_connectors = ex_connectors
        self.coordinate = coordinate
        self.register = []
        
# fill in all the necessary info
# into our blocks
unique_id = [x for x in graph.nodes()]
blocks = {key: [] for key in unique_id}
blocks_ex_connectors = {key: [] for key in unique_id}


for x in shape:
    coord = x[2:5]
    ex_connector = x[-1]
    if coord not in blocks[x[0]]:
        blocks[x[0]].append(coord)
        blocks_ex_connectors[x[0]].append(ex_connector)
     
    
padding = [None]
for key, value in blocks_ex_connectors.items():
    if len(value) != 4:
        value.extend((4-len(value))*padding)

padding = [None]
for key, value in blocks.items():
    if len(value) != 4:
        value.extend((4-len(value))*padding)
        
               
all_cells = []

n_iter = 0
for key, block in blocks.items():
    for i,cell in enumerate(block):
        id_tree = tree.leaves[key+(n_iter*3)+(i)]
        ex_connectors = blocks_ex_connectors[key][i]
        new_cell = Cell(id_tree.value, key, cell, [key], [ex_connectors])
        all_cells.append(new_cell)
    n_iter+=1

min_tree_id = all_cells[0].id_tree

for i in range(min_tree_id):
    new_cell = Cell(i, None, None, None, None)
    all_cells.append(new_cell)


if n_div<=2: 
    expression_vector = ['A', 'B']
    expression_direction = ['U', 'D']
elif n_div<=4:
    expression_vector = ['A', 'B', 'C', 'D']
    expression_direction = ['U', 'D', 'D', 'U']
elif n_div<=8:
    expression_vector = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    expression_direction = ['U', 'D', 'D', 'U', 'U', 'D', 'D', 'U']
else:
    raise ValueError('expression vector and direction are not implemented for a needed number of division > 8' )

warnings.warn('Make sure that the expression of the register you use is in alphabetical order (A, then B, etc)')


unique_register = [
                 'c_1_BP_site_R',
                 'G_gene_D',
                 'a_1_BP_site_L',
                 'H_gene_U',
                 '0_terminator_U',
                 '0_terminator_D',
                 'A_gene_D',
                 'a_1_BP_site_L',
                 'B_gene_U',
                 'c_1_BP_site_L',
                 'd_1_BP_site_R',
                 '0_promoter_U',
                 'd_1_BP_site_L',
                 'c_2_BP_site_R',
                 'C_gene_D',
                 'a_2_BP_site_L',
                 'D_gene_U',
                 '0_terminator_U',
                 '0_terminator_D', 
                 'E_gene_D',
                 'a_2_BP_site_L',
                 'F_gene_U',
                 'c_2_BP_site_L'
                ]

idx_first_recomb = 3
superTree = dt.SuperTree(all_cells, 
                         tree, 
                         unique_register, 
                         expression_vector, 
                         expression_direction, 
                         n_div, 
                         idx_first_recomb)

    
superTree.run(args.verbose)

# saving
np.save(args.save_path, superTree.expression_from_tree)



    