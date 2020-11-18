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

def mopac_get_H298(reactant, product, constraint, charge = 0, multiplicity = 'SINGLET'):
    """
    Create a directory folder called "tmp" for mopac calculation
    Create a input file called "input.mop" for mopac calculation
    """

    tmpdir = os.path.join(os.path.dirname(os.getcwd()), 'tmp')
    reactant_path = os.path.join(tmpdir, 'reactant.mop')
    product_path = os.path.join(tmpdir, 'product.mop')

    reac_geo, prod_geo = gen_geometry(reactant, product, constraint)
    
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

def gen_geometry(reactant_mol, product_mol, constraint):

    reactant_mol.gen3D(constraint, forcefield='uff', method = 'SteepestDescent', make3D=False)
    product_mol.gen3D(constraint, forcefield='uff', method = 'SteepestDescent', make3D=False)

    arrange3D = gen3D.Arrange3D(reactant_mol, product_mol, constraint)
    msg = arrange3D.arrangeIn3D()
    if msg != '':
        print(msg)

    reactant_mol.gen3D(constraint, forcefield='uff', method = 'SteepestDescent', make3D=False)
    product_mol.gen3D(constraint, forcefield='uff', method = 'SteepestDescent', make3D=False)

    print(reactant_mol.toNode())

    prod_geo = str(product_mol.toNode()).splitlines()
    product_geometry = []
    for i in prod_geo:
        i_list = i.split()
        atom = i_list[0] + " "
        k = i_list[1:] + [""]
        l = " 1 ".join(k)
        out = atom + l
        product_geometry.append(out)
    product_geometry = "\n".join(product_geometry)

    reac_geo = str(reactant_mol.toNode()).splitlines()
    reactant_geometry = []
    for i in reac_geo:
        i_list = i.split()
        atom = i_list[0] + " "
        k = i_list[1:] + [""]
        l = " 1 ".join(k)
        out = atom + l
        reactant_geometry.append(out)
    reactant_geometry = "\n".join(reactant_geometry)

    return reactant_geometry, product_geometry


xyz_path = './reactant.xyz'
constraint_path = './constraint.txt'
bonds = ((0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (1, 5, 1), (1, 6, 1), (1, 7, 1), (5, 8, 1), (5, 9, 1), (5, 10, 1), (8, 11, 1), (8, 12, 1), (8, 13, 1), (14, 15, 1), (14, 16, 1), (14, 17, 1), (14, 18, 1), (15, 20, 1), (16, 22, 1), (17, 21, 1), (18, 19, 1), (18, 23, 1), (20, 24, 1), (20, 25, 1), (20, 26, 1), (21, 27, 1), (21, 28, 1), (21, 30, 1), (22, 29, 1), (22, 31, 1), (22, 34, 1), (23, 32, 1), (23, 33, 1), (23, 35, 1))
constraint = extract_constraint_index(constraint_path)
reactant = readXYZ(xyz_path, bonds=bonds)
atoms = tuple(atom.atomicnum for atom in reactant)
product_bonds = ((0, 1, 2), (0, 3, 1), (0, 4, 1), (1, 6, 1), (1, 7, 1), (2, 17, 1), (5, 8, 1), (5, 9, 1), (5, 10, 1), (5, 19, 1), (8, 11, 1), (8, 12, 1), (8, 13, 1), (14, 15, 1), (14, 16, 1), (14, 17, 1), (14, 18, 1), (15, 20, 1), (16, 22, 1), (17, 21, 1), (18, 23, 1), (20, 24, 1), (20, 25, 1), (20, 26, 1), (21, 27, 1), (21, 28, 1), (21, 30, 1), (22, 29, 1), (22, 31, 1), (22, 34, 1), (23, 32, 1), (23, 33, 1), (23, 35, 1))
product = gen3D.makeMolFromAtomsAndBonds(atoms, product_bonds, spin=reactant.spin)
product.setCoordsFromMol(reactant)

mopac_get_H298(reactant, product, constraint)