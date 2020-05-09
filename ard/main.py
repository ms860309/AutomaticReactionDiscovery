#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Contains the :class:`ARD` for running an automatic reaction discovery. This
includes filtering reactions, generating 3D geometries, and running transition
state searches.
"""

from __future__ import print_function

import logging
import os
import time

import pybel
from rmgpy import settings
from rmgpy.data.thermo import ThermoDatabase
from rmgpy.molecule import Molecule

import constants
import gen3D
import util
from quantum import QuantumError
from node import Node
from pgen import Generate
from network import Network

from imaginary import Imaginary
from filter_rule import _filter
import openbabel as ob
import multiprocessing as mp
###############################################################################

class ARD(object):
    """
    Automatic reaction discovery class. Filters reactions based on estimated
    thermo of reactant and products, generates force field 3D geometries, and
    runs transition state searches.
    The attributes are:

    =============== ======================== ==================================
    Attribute       Type                     Description
    =============== ======================== ==================================
    `reac_smi`      ``str``                  A valid SMILES string describing the reactant structure
    `nbreak`        ``int``                  The maximum number of bonds that may be broken
    `nform`         ``int``                  The maximum number of bonds that may be formed
    `dh_cutoff`     ``float``                Heat of reaction cutoff (kcal/mol) for reactions that are too endothermic
    `theory_low`    ``str``                  Low level of theory for pre-optimizations
    `forcefield`    ``str``                  The force field for 3D geometry generation
    `distance`      ``float``                The initial distance between molecules
    `Qclass`        ``class``                A class representing the quantum software
    `output_dir`    ``str``                  The path to the output directory
    `logger`        :class:`logging.Logger`  The main logger
    =============== ======================== ==================================

    """

    def __init__(self, reac_smi, imaginarybond=0, nbreak=3, nform=3, dh_cutoff=20.0, theory_low=None,
                 forcefield='mmff94', distance=3.5, output_dir='', **kwargs):
        self.reac_smi = reac_smi
        self.imaginarybond = int(imaginarybond)
        self.nbreak = int(nbreak)
        self.nform = int(nform)
        self.dh_cutoff = float(dh_cutoff)
        self.theory_low = theory_low
        self.forcefield = forcefield
        self.distance = float(distance)
        qprog = kwargs.get('qprog', 'gau')
        self.Qclass = util.assignQclass(qprog)
        self.output_dir = output_dir
        log_level = logging.INFO
        self.logger = util.initializeLog(log_level, os.path.join(self.output_dir, 'ARD.log'), logname='main')
        self.reactant_list = []

    def executeXYZ(self, **kwargs):
        
        reac_mol = self.reac_smi
        reac_mol.gen3D(forcefield=self.forcefield)
        network = Network(reac_mol, forcefield = self.forcefield, **kwargs)
        network.genNetwork(reac_mol)
    
###############################################################################

def readInput(input_file):
    """
    Read input parameters from a file. It is assumed that the input file
    contains key-value pairs in the form "key value" on separate lines. If a
    keyword containing the string 'geometry' is encountered, the corresponding
    geometries are read in the form (example for methane dissociation):
        geometry (
        0 1
        C                 -0.03144385    0.03144654    0.00041162
        H                  0.32521058   -0.97736346    0.00041162
        H                  0.32522899    0.53584473    0.87406313
        H                  0.32522899    0.53584473   -0.87323988
        H                 -1.10144385    0.03145972    0.00041162
        ****
        C                 -0.36061854   -0.43406458    0.80670792
        H                  0.14377652   -1.32573293    0.49781771
        H                  0.14379613    0.27926689    1.42446520
        H                  0.56523315    0.87525286   -1.46111753
        H                 -1.36941886   -0.25571437    0.49781777
        )
    If '#' is found in a line, the rest of the line will be ignored.

    A dictionary containing all input parameters and their values is returned.
    """
    # Allowed keywords
    keys = ('reac_smi', 'xyz', 'imaginarybond', 'nbreak', 'nform', 'dh_cutoff', 'forcefield', 'name',
            'nsteps', 'nnode', 'lsf', 'tol', 'gtol', 'nlstnodes', 'dh_cutoff_method',
            'qprog', 'theory')

    # Read all data from file
    with open(input_file, 'r') as f:
        input_data = f.read().splitlines()
    # Create dictionary
    input_dict = {}

    # Read geometry block
    read = False
    geometry = []
    sep_loc = -1
    for line in input_data:
        if line.strip().startswith(')'):
            break
        if read and not line.strip().startswith('#') and line != '':
            geometry.append(line)
            if line.strip().startswith('*'):
                sep_loc = len(geometry) - 1
        elif 'geometry' in line:
            read = True

    if geometry:
        if sep_loc == -1:
            raise Exception('Incorrect geometry specification')

        # Extract multiplicity, atoms, and geometries
        multiplicity = geometry[0].split()[1]
        reactant = geometry[1:sep_loc]
        reac_atoms = [line.split()[0] for line in reactant]
        reac_geo = [[float(coord) for coord in line.split()[1:4]] for line in reactant]
        product = geometry[sep_loc + 1:]
        prod_atoms = [line.split()[0] for line in product]
        prod_geo = [[float(coord) for coord in line.split()[1:4]] for line in product]

        # Create nodes
        reac_node = Node(reac_geo, reac_atoms, multiplicity)
        prod_node = Node(prod_geo, prod_atoms, multiplicity)

        # Add to dictionary
        input_dict['reactant'] = reac_node
        input_dict['product'] = prod_node


    # Extract remaining keywords and values
    for line in input_data:
        if line != '' and not line.strip().startswith('#'):
            key = line.split()[0].lower()
            if key not in keys:
                continue
            if line.split()[1] == '=':
                input_dict[key] = line.split()[2]
            else:
                input_dict[key] = line.split()[1]
    # Check if valid method was specified and default to FSM
    try:
        jobtype = input_dict['jobtype'].lower()
    except KeyError:
        input_dict['jobtype'] = 'fsm'
    except AttributeError:
        raise Exception('Invalid jobtype')
    else:
        if jobtype != 'gsm' and jobtype != 'fsm':
            raise Exception('Invalid jobtype: {}'.format(jobtype))
    return input_dict


def readXYZ(xyz):
    mol = next(pybel.readfile('xyz', xyz))
    mol = gen3D.Molecule(mol.OBMol)
    return mol
