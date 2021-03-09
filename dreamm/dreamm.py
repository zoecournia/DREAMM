#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 11:10:55 2021

@author: alexis
"""

import sys
import argparse
import os
from lib.extract_features import featurizer
from lib.predict import predict


def process_command_line(argv):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description="Predict protein-membrane interfaces with DREAMM")
    
    # INPUTS
    input_args = parser.add_argument_group("Input Settings")
    input_args.add_argument("-i", "--input" , required=True, type=str,
        help="The 4 letter PDB code or the path to the PDB file")
    input_args.add_argument("-c", "--chain" , default=None, type=str, nargs="+", action='append',
        help="The chain or chains. Default is None")
    input_args.add_argument("-d", "--database" , required=True, type=str,
        help="The path to the Uniclust30_2018_8 database.")
    input_args.add_argument("-p", "--processes" , default=os.cpu_count() - 1, type=int,
        help="The number of CPUs. Default is All minus 1")
    
    args = parser.parse_args(sys.argv[1:])
    
    return args



def main(argv=None):
    args = process_command_line(argv)

    featurizer(args.input, args.chain, args.database, args.processes)
    print ('Finished feature extraction')
    
    if args.chain:
        pdb = os.path.splitext(os.path.basename(args.input))[0]
        file = pdb + "_chainA_fixed.pdb"
    else:
        file = os.path.splitext(os.path.basename(args.input))[0]
    
    predict(file)
    
    return 0 


if __name__ == "__main__":
    sys.exit(main(sys.argv))
