import os
import pybel
import rmgpy.molecule
from rmgpy.molecule.converter import from_ob_mol
from rmgpy.molecule.draw import MoleculeDrawer

def toRMGmol(OBMol):
    rmg_mol = from_ob_mol(rmgpy.molecule.molecule.Molecule(), OBMol)
    return rmg_mol

def readXYZ(path):
    mol = next(pybel.readfile('xyz', path))
    return mol.OBMol
    
def molDrawer():
    dirs = os.listdir('/mnt/d/reactions')
    for i in dirs:
        dir_path = os.path.join(os.path.join('/mnt/d/reactions', i), 'product.xyz')
        OBMol = readXYZ(dir_path)
        rmg_mol = toRMGmol(OBMol)
        _path = '/mnt/d/molecules/{}.png'.format(i)
        MoleculeDrawer().draw(rmg_mol, file_format='png', target=_path)
        
molDrawer()