#!/usr/bin/env python
# -*- coding: utf-8 -*-


if __name__ == '__main__':
    import argparse
    import sys
    import os
    from os import path
    sys.path.append(path.join(path.dirname( path.dirname( path.abspath(__file__))),'ard'))
    from main import ARD, readInput, readXYZ, add_bond

    # Set up parser for reading the input filename from the command line
    parser = argparse.ArgumentParser(description='Automatic Reaction Discovery')
    parser.add_argument('file', type=str, metavar='infile', help='An input file describing the job options')
    args = parser.parse_args()

    # Read input file
    input_file = os.path.abspath(args.file)
    kwargs = readInput(input_file)
    if kwargs['xyz'] == '1':
        OBMol = readXYZ(path.join(os.getcwd(),'reactant.xyz'))
        add_bonds = add_bond(path.join(os.getcwd(),'reactant.xyz'))
        if add_bonds[1] != '':
            kwargs['add_bonds'] = add_bonds[1]
        kwargs['reac_smi']= OBMol
        # Set output directory
        output_dir = path.abspath(path.dirname(input_file))
        kwargs['output_dir'] = output_dir

        # Execute job
        ard = ARD(**kwargs)
        ard.executeXYZ(**kwargs)
    else:
        # Set output directory
        output_dir = path.abspath(path.dirname(input_file))
        kwargs['output_dir'] = output_dir

        # Execute job
        ard = ARD(**kwargs)
        ard.execute(**kwargs)
