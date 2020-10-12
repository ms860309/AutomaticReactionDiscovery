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
        obmol.AddBond(bond[0], bond[1], bond[2])
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
    reactant_mol.separateMol()
    if len(reactant_mol.mols) > 1:
        reactant_mol.mergeMols()
    product_mol.separateMol()
    if len(product_mol.mols) > 1:
        product_mol.mergeMols()

    gen3D.constraint_force_field(reactant_mol.OBMol, constraint)
    gen3D.constraint_force_field(product_mol.OBMol, constraint)
    print(reactant_mol.toNode())
    raise
    arrange3D = gen3D.Arrange3D(reactant_mol, product_mol, constraint)
    msg = arrange3D.arrangeIn3D()
    if msg != '':
        print(msg)

    gen3D.constraint_force_field(reactant_mol.OBMol, constraint)
    gen3D.constraint_force_field(product_mol.OBMol, constraint)

    reactant = reactant_mol.toNode()
    product = product_mol.toNode()


bonds_path = './bonds.txt'
xyz_path = './reactant.xyz'
constraint_path = './constraint.txt'
bonds = extract_bonds(bonds_path)
constraint = extract_constraint_index(constraint_path)
reactant = readXYZ(xyz_path, bonds=bonds)
atoms = tuple(atom.atomicnum for atom in reactant)
product_bonds = ((0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (1, 2, 1), (1, 3, 1), (1, 5, 1), (2, 4, 1), (2, 5, 1), (2, 20, 1), (4, 5, 1), (5, 19, 1), (6, 13, 1), (7, 14, 1), (8, 9, 1), (8, 12, 1), (9, 10, 2), (9, 15, 1), (10, 11, 1), (10, 16, 1), (11, 12, 1), (11, 17, 1), (12, 13, 2), (13, 14, 1), (11, 18, 1))
product = gen3D.makeMolFromAtomsAndBonds(atoms, product_bonds, spin=reactant.spin)
product.setCoordsFromMol(reactant)

gen_geometry(reactant, product, constraint)