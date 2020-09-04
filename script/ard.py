#!/usr/bin/env python
# -*- coding: utf-8 -*-


if __name__ == '__main__':

    # standard library imports
    import sys
    import os
    from os import path
    sys.path.append(path.join(path.dirname( path.dirname( path.abspath(__file__))),'code/ard'))
    sys.path.append(path.join(path.dirname( path.dirname( path.abspath(__file__))),'code/mol_graph'))
    sys.path.append(path.join(path.dirname( path.dirname( path.abspath(__file__))),'database'))

    #third party
    import argparse
    import logging
    from rdkit import Chem
    from rdkit import RDLogger
    from openbabel import pybel

    # local application imports
    from main import ARD, readInput, readXYZ, extract_bonds, extract_constraint_index

    #disable log
    RDLogger.DisableLog('rdApp.*')
    rootlogger=logging.getLogger()
    rootlogger.setLevel(logging.CRITICAL)
    pybel.ob.obErrorLog.SetOutputLevel(0)

    # Set up parser for reading the input filename from the command line
    parser = argparse.ArgumentParser(description='Automatic Reaction Discovery')
    parser.add_argument('file', type=str, metavar='infile', help='An input file describing the job options')
    parser.add_argument('reactant', type=str, metavar='infile', help='An reactant xyz input file')
    parser.add_argument('-bonds', type=str, help='Manual specify bonds', required=False)
    parser.add_argument('-constraint', type=str, help='Manual specify constraint atom index (start from 0)', required=False)
    parser.add_argument('-generations',default=1, type=int, help='The network generation index', required=False)
    args = parser.parse_args()

    # Read input file
    input_file = os.path.abspath(args.file)
    ard_path = os.path.dirname(os.path.abspath(args.file))
    reactant_file = os.path.abspath(args.reactant)
    kwargs = readInput(input_file)
    
    # Constraint
    if kwargs['constraint'] == '1':
        index = extract_constraint_index(args.constraint)
        kwargs['constraint_index'] = index
    else:
        kwargs['constraint_index'] = None

    # Manual set up bonds
    if kwargs['manual_bonds'] == '1':
        bonds = extract_bonds(args.bonds)
        kwargs['bonds'] = bonds
    else:
        kwargs['bonds'] = None
    OBMol, reactant_graph = readXYZ(reactant_file, kwargs['bonds'])

    kwargs['reactant'] = OBMol

    if kwargs['graph'] == '1':
        kwargs['graph'] = reactant_graph
    else:
        kwargs['graph'] = None
        
    # Set output directory
    output_dir = path.abspath(path.dirname(input_file))
    kwargs['output_dir'] = output_dir
    kwargs['generations'] = args.generations
    kwargs['ard_path'] = ard_path
    kwargs['reactant_path'] = reactant_file
    
    # Execute job
    ard = ARD(**kwargs)
    ard.executeXYZ(**kwargs)
