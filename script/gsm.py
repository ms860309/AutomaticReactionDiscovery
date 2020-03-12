import pygsm
import os
import pybel as pb
def getXYZfilepath():
    """
    return xyz path
    """
    xyz = []
    abspath_pardir = os.path.abspath(os.pardir)
    rxn_path = os.path.join(abspath_pardir, 'reactions')
    dirs = os.listdir(rxn_path)
    for i in dirs:
        xyz_path = os.path.join(os.path.join(rxn_path, i), "reactant.xyz")
        xyz.append(xyz_path)
    return xyz

def genDLCobject():
    xyz_path_list = getXYZfilepath()
    mols = []
    for xyz_path in xyz_path_list:
        mol=pb.readfile("xyz",xyz_path).next()
        mols.append(mol)

genDLCobject()
