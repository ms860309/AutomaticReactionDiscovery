#!/usr/bin/env python
# -*- coding: utf-8 -*-


if __name__ == '__main__':
    import argparse
    import sys
    import os
    from os import path
    sys.path.append(path.join(path.dirname( path.dirname( path.abspath(__file__))),'ard'))
    sys.path.append(path.join(path.dirname( path.dirname( path.abspath(__file__))),'database'))
    from main import ARD, readInput, readXYZ, extract_bonds
    
    import logging
    from rdkit import Chem
    from rdkit import RDLogger
    import pybel
    #disable log
    RDLogger.DisableLog('rdApp.*')
    rootlogger=logging.getLogger()
    rootlogger.setLevel(logging.CRITICAL)
    ob_log_handler = pybel.ob.OBMessageHandler()
    ob_log_handler.SetOutputLevel(0)
    # Set up parser for reading the input filename from the command line
    parser = argparse.ArgumentParser(description='Automatic Reaction Discovery')
    parser.add_argument('file', type=str, metavar='infile', help='An input file describing the job options')
    parser.add_argument('reactant', type=str, metavar='infile', help='An reactant xyz input file')
    parser.add_argument('bonds', type=str, metavar='infile', help='Manual specify bonds')
    parser.add_argument('-generations',default=1, type=int, help='The network generation index',required=False)
    args = parser.parse_args()
    # Read input file
    input_file = os.path.abspath(args.file)
    ard_path = os.path.dirname(os.path.abspath(args.file))
    reactant_file = os.path.abspath(args.reactant)

    kwargs = readInput(input_file)
    
    if kwargs['manual_bonds'] == '1':
        bonds = extract_bonds(args.bonds)
        kwargs['bonds'] = bonds
    else:
        kwargs['bonds'] = []
    OBMol = readXYZ(reactant_file)

    kwargs['reac_smi']= OBMol
    # Set output directory
    output_dir = path.abspath(path.dirname(input_file))
    kwargs['output_dir'] = output_dir
    kwargs['generations'] = args.generations
    kwargs['ard_path'] = ard_path

    # Execute job
    ard = ARD(**kwargs)
    ard.executeXYZ(**kwargs)
