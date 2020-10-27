# standard library imports
import os
import shutil
import time
import copy

#third party
from subprocess import Popen, PIPE
import difflib
from openbabel import pybel
from openbabel import openbabel as ob
import numpy as np

# local application imports
from node import Node
import gen3D
import util


class MopacError(Exception):
    """
    An exception class for errors that occur during mopac calculations.
    """
    pass

class Mopac(object):

    def __init__(self, reactant_mol, product_mol, mopac_method, forcefield, constraintff_alg, form_bonds, logger, count, num, constraint = None):
        self.reactant_mol = reactant_mol
        self.product_mol = product_mol
        self.mopac_method = mopac_method
        self.forcefield = forcefield
        self.constraintff_alg = constraintff_alg
        self.form_bonds = form_bonds
        self.logger = logger
        self.count = count
        self.num = num
        self.constraint = constraint

    def mopac_get_H298(self, reac_mol_copy, charge = 0, multiplicity = 'SINGLET'):
        """
        Create a directory folder called "tmp" for mopac calculation
        Create a input file called "input.mop" for mopac calculation
        """

        tmpdir = os.path.join(os.path.dirname(os.getcwd()), 'tmp')
        reactant_path = os.path.join(tmpdir, 'reactant.mop')
        product_path = os.path.join(tmpdir, 'product.mop')

        reac_geo, prod_geo, reactant_mol, product_mol = self.genInput(self.reactant_mol, self.product_mol, reac_mol_copy)
        
        if reac_geo == False and prod_geo == False:
            return False, False, False, False
        else:
            if os.path.exists(tmpdir):
                shutil.rmtree(tmpdir)
            os.mkdir(tmpdir)
            with open(reactant_path, 'w') as f:
                f.write("NOSYM 1SCF CHARGE={} {} {}\n\n".format(charge, multiplicity, self.mopac_method))
                f.write("\n{}".format(reac_geo))
            start_time = time.time()
            runMopac(tmpdir, 'reactant.mop')
            reactant = getHeatofFormation(tmpdir, 'reactant.out')

            with open(product_path, 'w') as f:
                f.write("CHARGE={} {} {}\n\n".format(charge, multiplicity, self.mopac_method))
                f.write("\n{}".format(prod_geo))
            runMopac(tmpdir, 'product.mop')
            product = getHeatofFormation(tmpdir, 'product.out')
            self.finalize(start_time, 'mopac')
            return float(reactant), float(product), reactant_mol, product_mol
    
    def genInput(self, reactant_mol, product_mol, reac_mol_copy, threshold = 5.0):
        start_time = time.time()

        # Initial optimization
        if self.constraint == None:
            Hatom = gen3D.readstring('smi', '[H]')
            ff = pybel.ob.OBForceField.FindForceField(self.forcefield)
            reactant_mol.gen3D(self.constraint, forcefield=self.forcefield, method = self.constraintff_alg, make3D=False)
            product_mol.gen3D(self.constraint, forcefield=self.forcefield, method = self.constraintff_alg, make3D=False)
        else:
            reactant_mol.gen3D(self.constraint, forcefield=self.forcefield, method = self.constraintff_alg, make3D=False)
            product_mol.gen3D(self.constraint, forcefield=self.forcefield, method = self.constraintff_alg, make3D=False)

        # Arrange
        try:
            arrange3D = gen3D.Arrange3D(reactant_mol, product_mol, self.constraint)
            msg = arrange3D.arrangeIn3D()
            if msg != '':
                print(msg)
        except:
            return False, False, False, False

        # Check reactant expected forming bond length must smaller than 4 angstrom after arrange. Default = 4
        dist = self.check_bond_length(reactant_mol, self.form_bonds) # return the maximum value in array
        # After arrange to prevent openbabel use the previous product coordinates if it is isomorphic
        # to the current one, even if it has different atom indices participating in the bonds.
        if self.constraint == None:
            ff.Setup(Hatom.OBMol)
            reactant_mol.gen3D(self.constraint, forcefield=self.forcefield, method = self.constraintff_alg, make3D=False)
            ff.Setup(Hatom.OBMol)
            product_mol.gen3D(self.constraint, forcefield=self.forcefield, method = self.constraintff_alg, make3D=False)
            ff.Setup(Hatom.OBMol)
        else:
            reactant_mol.gen3D(self.constraint, forcefield=self.forcefield, method = self.constraintff_alg, make3D=False)
            product_mol.gen3D(self.constraint, forcefield=self.forcefield, method = self.constraintff_alg, make3D=False)
        
        if dist >= threshold:
            self.logger.info('Here is the {} product.'.format(self.num))
            self.logger.info('Form bonds: {}\nDistance: {}'.format(self.form_bonds, dist))
            self.logger.info('Form bond distance is greater than threshold.')
            self.logger.info('Now finished {}/{}'.format(self.num, self.count))
            self.finalize(start_time, 'arrange')
            return False, False, False, False
        else:
            self.logger.info('\nHere is the {} product.'.format(self.num))
            self.logger.info('Structure:\n{}'.format(str(reactant_mol.toNode())))
            self.logger.info('Structure:\n{}\n'.format(str(product_mol.toNode())))
            self.logger.info('Form bonds: {}\nDistance: {}'.format(self.form_bonds, dist))    
            prod_geo = str(product_mol.toNode()).splitlines()
            product_geometry = []
            for i in prod_geo:
                i_list = i.split()
                atom = i_list[0] + " "
                k = i_list[1:] + [""]
                l = " 0 ".join(k)
                out = atom + l
                product_geometry.append(out)
            product_geometry = "\n".join(product_geometry)

            reac_geo = str(reactant_mol.toNode()).splitlines()
            reactant_geometry = []
            for i in reac_geo:
                i_list = i.split()
                atom = i_list[0] + " "
                k = i_list[1:] + [""]
                l = " 0 ".join(k)
                out = atom + l
                reactant_geometry.append(out)
            reactant_geometry = "\n".join(reactant_geometry)
            
            reactant = reactant_mol.toNode()
            reactant_mol.setCoordsFromMol(reac_mol_copy)
            self.finalize(start_time, 'arrange')
            return reactant_geometry, product_geometry, reactant, product_mol.toNode()
            
    def finalize(self, start_time, jobname):
        """
        Finalize the job.
        """
        self.logger.info('Total {} run time: {:.1f} s'.format(jobname, time.time() - start_time))

    @staticmethod
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