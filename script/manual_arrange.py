# standard library imports
import sys
import os
from os import path
sys.path.append(path.join(path.dirname(path.dirname( path.abspath(__file__))),'code/ard'))
sys.path.append(path.join(path.dirname(path.dirname( path.abspath(__file__))),'code/mol_graph'))
sys.path.append(path.join(path.dirname(path.dirname( path.abspath(__file__))),'database'))

# local application imports
import gen3D
import pgen
from gen3D import Molecule

from openbabel import openbabel as ob
from openbabel import pybel as pb

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
    for bond in bonds:
        obmol.AddBond(bond[0] +1, bond[1] +1, bond[2])
    obmol.SetTotalCharge(int(mol.charge))
    obmol.Center()
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

def gen_geometry(reactant_mol, product_mol, constraint):

    reactant_mol.gen3D(constraint, forcefield='uff', method = 'SteepestDescent', make3D=False)
    product_mol.gen3D(constraint, forcefield='uff', method = 'SteepestDescent', make3D=False)

    arrange3D = gen3D.Arrange3D(reactant_mol, product_mol, constraint)
    msg = arrange3D.arrangeIn3D()
    if msg != '':
        print(msg)

    reactant_mol.gen3D(constraint, forcefield='uff', method = 'SteepestDescent', make3D=False)
    product_mol.gen3D(constraint, forcefield='uff', method = 'SteepestDescent', make3D=False)

    reactant = reactant_mol.toNode()
    product = product_mol.toNode()
    print(reactant_mol.toNode())
    print(product_mol.toNode())

xyz_path = './reactant.xyz'
constraint_path = './constraint.txt'
bonds = ((0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (1, 5, 1), (1, 6, 1), (1, 7, 1), (5, 8, 1), (5, 9, 1), (5, 10, 1), (8, 11, 1), (8, 12, 1), (8, 13, 1), (14, 15, 1), (14, 16, 1), (14, 17, 1), (14, 18, 1), (15, 20, 1), (16, 22, 1), (17, 21, 1), (18, 19, 1), (18, 23, 1), (20, 24, 1), (20, 25, 1), (20, 26, 1), (21, 27, 1), (21, 28, 1), (21, 30, 1), (22, 29, 1), (22, 31, 1), (22, 34, 1), (23, 32, 1), (23, 33, 1), (23, 35, 1))
constraint = extract_constraint_index(constraint_path)
reactant = readXYZ(xyz_path, bonds=bonds)
atoms = tuple(atom.atomicnum for atom in reactant)
product_bonds = ((0, 1, 2), (0, 3, 1), (0, 4, 1), (1, 6, 1), (1, 7, 1), (2, 17, 1), (5, 8, 1), (5, 9, 1), (5, 10, 1), (5, 19, 1), (8, 11, 1), (8, 12, 1), (8, 13, 1), (14, 15, 1), (14, 16, 1), (14, 17, 1), (14, 18, 1), (15, 20, 1), (16, 22, 1), (17, 21, 1), (18, 23, 1), (20, 24, 1), (20, 25, 1), (20, 26, 1), (21, 27, 1), (21, 28, 1), (21, 30, 1), (22, 29, 1), (22, 31, 1), (22, 34, 1), (23, 32, 1), (23, 33, 1), (23, 35, 1))
product = gen3D.makeMolFromAtomsAndBonds(atoms, product_bonds, spin=reactant.spin)
product.setCoordsFromMol(reactant)

gen_geometry(reactant, product, constraint)