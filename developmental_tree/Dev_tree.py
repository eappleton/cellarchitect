import os, sys
import numpy as np
import networkx as nx
import string
import pickle

from itertools import chain
from binarytree import tree, Node, build

sys.path.append('../../python/')
import helper as hp

class SuperTree():
    def __init__(self, all_cells, tree, register, expression_vector, 
                 expression_direction, n_div, idx_first_recomb):
        self.all_cells = all_cells
        self.tree = tree
        self.register = register
        self.expression_vector = expression_vector
        self.expression_direction = expression_direction
        self.n_div = n_div
        self.numb = False
        
        self.expression_from_tree = [None]*len(expression_vector)
        self.id_recombinase = self.get_id_from_register(self.register, 'lower')
        self.id_gene = self.get_id_from_register(self.register, 'upper')
        
        self.i_node = 0
        
    def get_id_from_register(self, unique_register, case):
        if case == 'lower':
            letters = string.ascii_lowercase
        elif case == 'upper':
            letters = string.ascii_uppercase
        else:
            raise ValueError('case must be either lower or upper, in string format')
    
        id_register_ = [x[0] for x in unique_register if x[0] in letters]
        id_register = sorted(list(set(id_register_)))[-1]
        idx = letters.index(id_register)
        
        return letters[idx:]
        
    def get_first_children_id(self, node):
        left_child = list(node.left)[0]
        right_child = list(node.right)[0]
        
        return [left_child, right_child]
        
    def create_cassette(self, cell, list_children, level_id):
        
        direction = self.expression_direction[level_id]
        recomb_id = self.id_recombinase[0]
            
        # case where the node represents a leaf
        # not that here, every elements must have
        # the same orientation as the cassette in
        # the register will be expressed for a given
        # orientation at the count corresponding to the 
        # tree leafs expression
        if list_children is None:
            for i,val in enumerate(cell.ex_connectors):
                if val is not None:
                    cell.register.append('{}_ex_gene_{}'.format(val, direction))
                    #cell.register.append('{}_gene_{}'.format(val, direction))
            for i,val in enumerate(cell.in_connectors):
                cell.register.append('{}_in_gene_{}'.format(val, direction))
                #cell.register.append('{}_gene_{}'.format(val, direction))
                                
        else:
            children_cells = [x for x in self.all_cells if x.id_tree in list_children]
            connectors_1 = children_cells[0].register
            connectors_2 = children_cells[1].register
            
            if connectors_1 == connectors_2:
                # no numb needed
                cell.register = connectors_1
            else:
                # numb needed
                self.numb = True
                register = self.create_numb_register(cell, connectors_1, connectors_2, recomb_id)
                cell.register = register
    
    def create_numb_register(self, cell, connectors_1, connectors_2, recomb_id):

        placeholder = [recomb_id  + '_' + str(self.i_node) + '_BP_site_R',
                       'spot1', 
                       '0_terminator_U',
                       '0_terminator_D',
                       'spot2',
                       recomb_id  + '_' + str(self.i_node) + '_BP_site_L']
        self.i_node+=1
           
        idx1 = placeholder.index('spot1') 
        hp.extend_at_index(placeholder, connectors_1, index=idx1)
        idx2 = placeholder.index('spot2')
        hp.extend_at_index(placeholder, connectors_2, index=idx2)
        
        register = [x for x in placeholder if x not in ['spot1', 'spot2']]
        
        return register
        
    def add_direction(self, list_, direction):
        
        new_list = []
        for x in list_:
            if 'BP' in x or '_D' in x or '_U' in x:
                new_list.append(x)
            else:
                new_list.append(x + direction)
        
        return new_list

    def update_expression_from_tree(self, level_id):
        """
        Put the cassette needed to be expressed at each
        count into expression_from_tree
        """
        
        direction = self.expression_direction[level_id]
        
        if self.numb:
            recomb_id = self.id_recombinase[0]
            cassette = [recomb_id + '_numb_' + direction]
        else:
            cassette = ['empty']
                    
        self.expression_from_tree[level_id] = cassette
        self.numb = False
            
    def clean_expression_from_tree(self):
        
        # remove none value
        temp = []
        for x in self.expression_from_tree:
            if x is not None:
                temp.append(x)
                
        self.expression_from_tree = temp
        
    def finalize_expression_from_tree(self):
        
        id_cell_zero = -1
        for i,cell in enumerate(self.all_cells):
            if cell.id_tree == 0:
                id_cell_zero = i
                break
                
        temp = self.all_cells[id_cell_zero].register
        # replace int id by upper letter for gene
        for i,val in enumerate(temp):
            if 'gene' in val:
                id_ = self.id_gene[int(val[0])]
                temp[i] = id_ + val[1:]
        
        self.expression_from_tree[self.n_div] = temp
        
    def update_register(self):
        pass
        
    def run(self, verbose=False):
        """
        Main method of the class
        Run this function to create the register
        """
                
        i_level = self.n_div        
        for level in reversed(self.tree.levels):
            if verbose:
                print('************')
                print('************')
                print('level index',i_level)
            for node in level:
                if verbose:
                    print('*****')
                    print('node',node.value)
                                    
                temp = [x for x in self.all_cells if x.id_tree==node.value]
                
                cell = temp[0]
                if node.left is None:
                    self.create_cassette(cell, None, i_level)
                else:
                    first_childrens = self.get_first_children_id(node)
                    self.create_cassette(cell, first_childrens, i_level)
                if verbose:
                    print('cell register after loop', cell.register)
                
            
            self.update_expression_from_tree(i_level)
            i_level-=1
            self.i_node = 1
            self.id_recombinase = self.id_recombinase[1:]
                        
        # cleaning of the expression vector
        self.clean_expression_from_tree()
        # add the last register from the root cell
        self.finalize_expression_from_tree()
        
    
    