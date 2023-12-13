import os,sys
import argparse
import csv
import pickle
import numpy as np

import Verification as ver

sys.path.append('../../')
from python import helper as hp

parser = argparse.ArgumentParser()
parser.add_argument('path_circuit', type=str, help='Path to the counter circuit to analyse')
parser.add_argument('save_path', type=str, help='Path to save the output figure and the register matrix')
parser.add_argument('verbose', type=bool, help='to print the verification process step by step or not')
args = parser.parse_args()

with open(args.path_circuit, 'r') as pc:
    reader = csv.reader(pc)
    counter = [x[0] for x in list(reader)]

v = ver.Verification()

count, history, _, changes_history = v.read_sequence(counter, verbose=args.verbose)
    
if args.save_path[-1] != '/':
    args.save_path+='/'
if not os.path.exists(args.save_path):
    os.makedirs(args.save_path)

hp.history_to_figure(history, str(count), changes_history, dir_=args.save_path, verification_highlight=True)

with open(args.save_path + 'changes.txt', 'wb') as fp:
    pickle.dump(changes_history, fp)
 
if args.verbose:
    print('Counter design reach count: ', count)



