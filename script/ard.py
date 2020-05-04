#!/usr/bin/env python
# -*- coding: utf-8 -*-


if __name__ == '__main__':
    import argparse
    import sys
    import os
    from os import path
    sys.path.append(path.join(path.dirname( path.dirname( path.abspath(__file__))),'ard'))
    sys.path.append(path.join(path.dirname( path.dirname( path.abspath(__file__))),'database'))
    from main import ARD, readInput, readXYZ
    
    import logging
    from rdkit import Chem
    from rdkit import RDLogger
    
    #disable log
    RDLogger.DisableLog('rdApp.*')
    rootlogger=logging.getLogger()
    rootlogger.setLevel(logging.CRITICAL)
    
    
    # Set up parser for reading the input filename from the command line
    parser = argparse.ArgumentParser(description='Automatic Reaction Discovery')
    parser.add_argument('file', type=str, metavar='infile', help='An input file describing the job options')
    args = parser.parse_args()

    # Read input file
    input_file = os.path.abspath(args.file)
    kwargs = readInput(input_file)
    OBMol = readXYZ(path.join(os.getcwd(),'reactant.xyz'))
    kwargs['reac_smi']= OBMol
    # Set output directory
    output_dir = path.abspath(path.dirname(input_file))
    kwargs['output_dir'] = output_dir

    # Execute job
    ard = ARD(**kwargs)
    ard.executeXYZ(**kwargs)
