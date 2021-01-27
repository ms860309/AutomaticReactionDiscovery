# standard library imports
import sys
import os
from os import path
sys.path.append(path.join(path.dirname(path.dirname( path.abspath(__file__))),'code/ard'))
sys.path.append(path.join(path.dirname(path.dirname( path.abspath(__file__))),'code/mol_graph'))
sys.path.append(path.join(path.dirname(path.dirname( path.abspath(__file__))),'database'))
import shutil

# local application imports
import gen3D
import pgen
from gen3D import Molecule
from subprocess import Popen, PIPE

from openbabel import openbabel as ob
from openbabel import pybel as pb

import numpy as np

def readXYZ(xyz, bonds = None):
    # extract molecule information from xyz
    mol = next(pb.readfile('xyz', xyz))
    # Manually give bond information 
    # (Because in metal system the bond information detect by openbabel usually have some problem)
    m = Molecule(pb.ob.OBMol())
    obmol = m.OBMol
    obmol.BeginModify()
    for atom in mol:
        coords = [coord for coord in atom.coords]
        atomno = atom.atomicnum
        obatom = ob.OBAtom()
        obatom.thisown = 0
        obatom.SetAtomicNum(atomno)
        obatom.SetVector(*coords)
        obmol.AddAtom(obatom)
        del obatom

    bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondOrder())
                    for bond in pb.ob.OBMolBondIter(mol.OBMol)]
    bonds.extend([(13,14,1), (13,15,1), (13,16,1), (13,17,1)])

    for bond in bonds:
        obmol.AddBond(bond[0], bond[1], bond[2])

    #obmol.PerceiveBondOrders()
    #obmol.SetTotalCharge(int(mol.charge))
    #obmol.Center()
    obmol.EndModify()
    mol_obj = gen3D.Molecule(obmol)
    return mol_obj

def extract_bonds(bonds):
    with open(bonds, 'r') as f:
        lines = f.read()
    lines = eval(lines)
    return lines

def extract_constraint_index(constraint):
    with open(constraint, 'r') as f:
        lines = f.read()
    lines = eval(lines)
    return lines

def extract_fixed_atom_index(fixed_atom):
    with open(fixed_atom, 'r') as f:
        lines = f.read()
    lines = eval(lines)
    return lines

def mopac_get_H298(reactant, product, constraint, fixed_atom, charge = 0, multiplicity = 'SINGLET'):
    """
    Create a directory folder called "tmp" for mopac calculation
    Create a input file called "input.mop" for mopac calculation
    """

    tmpdir = os.path.join(os.path.dirname(os.getcwd()), 'tmp')
    reactant_path = os.path.join(tmpdir, 'reactant.mop')
    product_path = os.path.join(tmpdir, 'product.mop')

    reac_geo, prod_geo = gen_geometry(reactant, product, constraint, fixed_atom)
    
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
    os.mkdir(tmpdir)
    with open(reactant_path, 'w') as f:
        f.write("NOSYM CHARGE={} {} {}\n\n".format(charge, multiplicity, 'PM7'))
        f.write("\n{}".format(reac_geo))

    runMopac(tmpdir, 'reactant.mop')
    reactant = getHeatofFormation(tmpdir, 'reactant.out')

    with open(product_path, 'w') as f:
        f.write("NOSYM CHARGE={} {} {}\n\n".format(charge, multiplicity, 'PM7'))
        f.write("\n{}".format(prod_geo))
    runMopac(tmpdir, 'product.mop')
    product = getHeatofFormation(tmpdir, 'product.out')
    print(float(reactant))
    print(float(product))
    print('delta H is {}'.format(float(product) - float(reactant)))

def check_bond_length(product, add_bonds):
    """
    Use reactant coordinate to check if the add bonds's bond length is too long.
    Return a 'list of distance'.
    """
    coords = [atom.coords for atom in product]
    atoms = tuple(atom.atomicnum for atom in product)
    coords = [np.array(coords).reshape(len(atoms), 3)]

    dist = []
    for bond in add_bonds:
        coord_vect_1 = coords[0][bond[0]]
        coord_vect_2 = coords[0][bond[1]]
        diff = coord_vect_1 - coord_vect_2
        dist.append(np.linalg.norm(diff))
        
    if dist == []:
        dist = [0]
    return float(max(dist))

def getHeatofFormation(tmpdir, target = 'reactant.out'):
    """
    if Error return False, which HF may be 0.0
    """
    input_path = os.path.join(tmpdir, target)
    with open(input_path, 'r') as f:
        lines = f.readlines()
    for idx, line in enumerate(lines):
        if line.strip().startswith('FINAL HEAT OF FORMATION'):
            break
    string = lines[idx].split()
    if string[0] == 'FINAL':
        HeatofFormation = string[5]
    else:
        HeatofFormation = False
    return HeatofFormation

def runMopac(tmpdir, target = 'reactant.mop'):
    input_path = os.path.join(tmpdir, target)
    p = Popen(['mopac', input_path])
    p.wait()

def gen_geometry(reactant_mol, product_mol, constraint, fixed_atom):
    print(constraint)
    reactant_mol.gen3D(constraint, forcefield='uff', method = 'ConjugateGradients', make3D=False)
    product_mol.gen3D(constraint, forcefield='uff', method = 'ConjugateGradients', make3D=False)
    print(reactant_mol.toNode())
    print('----')
    print(product_mol.toNode())
    print('----')
    arrange3D = gen3D.Arrange3D(reactant_mol, product_mol, constraint, fixed_atom)
    msg = arrange3D.arrangeIn3D()
    if msg != '':
        print(msg)
    
    print(reactant_mol.toNode())
    print('----')
    print(product_mol.toNode())
    print('----')
    gen3D.constraint_force_field(reactant_mol, constraint, forcefield='uff', method = 'ConjugateGradients')
    gen3D.constraint_force_field(product_mol, constraint, forcefield='uff', method = 'ConjugateGradients')

    prod_geo = str(product_mol.toNode()).splitlines()
    product_geometry = []
    for idx, i in enumerate(prod_geo):
        i_list = i.split()
        atom = i_list[0] + " "
        k = i_list[1:] + [""]
        if idx in constraint:
            l = " 0 ".join(k)
        else:
            l = " 1 ".join(k)
        out = atom + l
        product_geometry.append(out)
    product_geometry = "\n".join(product_geometry)

    reac_geo = str(reactant_mol.toNode()).splitlines()
    reactant_geometry = []
    for idx, i in enumerate(reac_geo):
        i_list = i.split()
        atom = i_list[0] + " "
        k = i_list[1:] + [""]
        if idx in constraint:
            l = " 0 ".join(k)
        else:
            l = " 1 ".join(k)
        out = atom + l
        reactant_geometry.append(out)
    reactant_geometry = "\n".join(reactant_geometry)

    return reactant_geometry, product_geometry


xyz_path = './reactant.xyz'
constraint_path = './constraint.txt'
fixed_atom_path = './fixed_atom.txt'

constraint = extract_constraint_index(constraint_path)
fixed_atom = extract_fixed_atom_index(fixed_atom_path)
reactant = readXYZ(xyz_path)
atoms = tuple(atom.atomicnum for atom in reactant)
rb = ((0, 1, 1), (0, 5, 1), (0, 6, 1), (0, 7, 1), (1, 2, 1), (1, 4, 1), (1, 8, 1), (2, 3, 1), (2, 9, 1), 
    (2, 10, 1), (3, 11, 1), (4, 12, 1), (5, 13, 1), (14, 15, 1), (14, 16, 1), (14, 17, 1), (14, 18, 1), (15, 20, 1), 
    (16, 22, 1), (17, 19, 1), (17, 23, 1), (18, 21, 1), (20, 24, 1), (20, 25, 1), (20, 26, 1), (21, 27, 1), (21, 28, 1), (21, 29, 1), 
    (22, 30, 1), (22, 31, 1), (22, 32, 1), (23, 33, 1), (23, 34, 1), (23, 35, 1))

product_bonds = ((0, 1, 1), (0, 5, 1), (0, 6, 1), (0, 7, 1), (1, 4, 1), (1, 8, 1), (2, 3, 2), (2, 9, 1), (1, 11, 1),
    (2, 10, 1), (4, 12, 1), (5, 13, 1), (14, 15, 1), (14, 16, 1), (14, 17, 1), (14, 18, 1), (15, 20, 1), (16, 19, 1), 
    (16, 22, 1), (17, 23, 1), (18, 21, 1), (20, 24, 1), (20, 25, 1), (20, 26, 1), (21, 27, 1), (21, 28, 1), (21, 29, 1), 
    (22, 30, 1), (22, 31, 1), (22, 32, 1), (23, 33, 1), (23, 34, 1), (23, 35, 1))

product = gen3D.makeMolFromAtomsAndBonds(atoms, product_bonds, spin=reactant.spin)
product.setCoordsFromMol(reactant)

mopac_get_H298(reactant, product, constraint, fixed_atom)