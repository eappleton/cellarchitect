import os, sys
import argparse
import csv
import pickle
import numpy as np
import networkx as nx

import Connection as cn

sys.path.append('../../')
from python import helper as hp    

parser = argparse.ArgumentParser()
parser.add_argument('path_to_shape', type=str, help='Path to the shape design to mesh')
parser.add_argument('algo', type=str, help='Type of algorithm to use to optimize the connections')
parser.add_argument('algo_param', type=int, help='Parameter of the algorithm')
parser.add_argument('save_path', type=str, help='Path to save the output figure and the register matrix')
parser.add_argument('name', type=str, help='Name of the output file')
parser.add_argument('verbose', type=bool, help='to print the verification process step by step or not')
args = parser.parse_args()


graph = cn.convert_mesh_to_graph(args.path_to_shape)

if args.algo == 'MST':
    graph_optimize = nx.minimum_spanning_tree(graph)
elif args.algo == 'K_from_none':
    graph_optimize, _ = cn.k_connected_from_none(graph, 0, args.algo_param)
elif args.algo == 'K_from_mst':
    graph_mst = nx.minimum_spanning_tree(graph)
    graph_optimize, _ = cn.k_connected_from_mst(graph, graph_mst, args.algo_param)
elif args.algo == 'freely_moving_part':
    graph_mst = nx.minimum_spanning_tree(graph)
    graph_optimize = cn.moving_part(graph, graph_mst, args.algo_param)
else:
    raise ValueError('Algo must be MST, K_from_none, K_from_mst or freely_moving_part. Found %s' % (args.algo))

 
if args.verbose:
    cn.print_edges_and_vertices(graph_optimize)


if args.save_path[-1] != '/':
    args.save_path+='/'
if not os.path.exists(args.save_path):
    os.makedirs(args.save_path)
    
cn.write_output_file(graph_optimize, args.save_path, args.name)
    



