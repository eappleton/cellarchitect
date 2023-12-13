import os, sys
import argparse
import csv
import pickle
import numpy as np

import Register as re

sys.path.append('../../')
from python import helper as hp

parser = argparse.ArgumentParser()
parser.add_argument('path_circuit', type=str, help='Path to the register circuit to analyse')
parser.add_argument('path_program', type=str, help='Path to the counter regulatory program')
parser.add_argument('save_path', type=str, help='Path to save the output figure and the register matrix')
parser.add_argument('name', type=str, help='Name of the output file')
parser.add_argument('verbose', type=bool, help='to print the verification process step by step or not')
args = parser.parse_args()


with open(args.path_program, 'rb') as pp:
    regulation_program = pickle.load(pp)

with open(args.path_circuit, 'r') as pc:
    reader = csv.reader(pc)
    register = [x[0] for x in list(reader)]

r = re.Register()
    
unique_genes = [x[0] for x in register if 'gene' in x]
unique_genes.sort()

history, expression_matrix = r.read_sequence(register, regulation_program, unique_genes, verbose=args.verbose)

if args.save_path[-1] != '/':
    args.save_path+='/'
if not os.path.exists(args.save_path):
    os.makedirs(args.save_path)
    
hp.history_to_figure(history, args.name, regulation_program, dir_=args.save_path, ex_matrix=expression_matrix)
np.save(args.save_path + 'expression_matrix', expression_matrix)


