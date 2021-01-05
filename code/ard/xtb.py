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


class XTBError(Exception):
    """
    An exception class for errors that occur during mopac calculations.
    """
    pass

class XTB(object):

    def __init__(self, reactant_mol, product_mol, forcefield, constraintff_alg, form_bonds, logger, count, num, constraint = None, fixed_atom= None):
        self.reactant_mol = reactant_mol
        self.product_mol = product_mol
        self.forcefield = forcefield
        self.constraintff_alg = constraintff_alg
        self.form_bonds = form_bonds
        self.logger = logger
        self.count = count
        self.num = num
        self.constraint = constraint
        self.fixed_atom = fixed_atom

    def xtb_get_H298(self, reac_mol_copy, reactant_path):
        """
        Create a directory folder called "tmp" for mopac calculation
        Create a input file called "input.mop" for mopac calculation
        """

        tmpdir = os.path.join(reactant_path, 'tmp')
        reactant_path = os.path.join(tmpdir, 'reactant.xyz')
        product_path = os.path.join(tmpdir, 'product.xyz')

        reac_geo, prod_geo = self.genInput(self.reactant_mol, self.product_mol, reac_mol_copy)
        
        if reac_geo == False and prod_geo == False:
            return False, False
        else:
            if os.path.exists(tmpdir):
                shutil.rmtree(tmpdir)
            os.mkdir(tmpdir)
            os.chdir(tmpdir)
            with open(reactant_path, 'w') as f:
                f.write(str(len(str(reac_geo).splitlines())))
                f.write('\n\n')
                f.write('{}'.format(reac_geo))

            start_time = time.time()
            try:
                self.runXTB(tmpdir, 'reactant.xyz')
                reactant_energy = self.getE(tmpdir, 'reactant.xyz')
            except:
                return False, False

            with open(product_path, 'w') as f:
                f.write(str(len(str(prod_geo).splitlines())))
                f.write('\n\n')
                f.write('{}'.format(prod_geo))
            try:
                self.runXTB(tmpdir, 'product.xyz')
                product_energy = self.getE(tmpdir, 'product.xyz')
            except:
                return False, False
                
            self.finalize(start_time, 'XTB')
            return float(reactant_energy), float(product_energy)
    
    def genInput(self, reactant_mol, product_mol, reac_mol_copy, threshold = 10.0):
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
        try:  # Pass the more than 4 fragment situation
            arrange3D = gen3D.Arrange3D(reactant_mol, product_mol, self.constraint, self.fixed_atom)
            msg = arrange3D.arrangeIn3D()
            if msg != '':
                print(msg)
        except:
            self.logger.info('Here is the {} product.'.format(self.num))
            self.logger.info('Arrange fail')
            return False, False

        if self.constraint == None:
            ff.Setup(Hatom.OBMol)
            reactant_mol.gen3D(self.constraint, forcefield=self.forcefield, method = self.constraintff_alg, make3D=False)
            ff.Setup(Hatom.OBMol)
            product_mol.gen3D(self.constraint, forcefield=self.forcefield, method = self.constraintff_alg, make3D=False)
            ff.Setup(Hatom.OBMol)
        else:
            reactant_mol.gen3D(self.constraint, forcefield=self.forcefield, method = self.constraintff_alg, make3D=False)
            product_mol.gen3D(self.constraint, forcefield=self.forcefield, method = self.constraintff_alg, make3D=False)

        # Check reactant expected forming bond length must smaller than 4 angstrom after arrange. Default = 4
        # After arrange to prevent openbabel use the previous product coordinates if it is isomorphic
        # to the current one, even if it has different atom indices participating in the bonds.
        dist = self.check_bond_length(reactant_mol, self.form_bonds) # return the maximum value in array

        if dist >= threshold:
            self.logger.info('Here is the {} product.'.format(self.num))
            self.logger.info('Form bonds: {}\nDistance: {}'.format(self.form_bonds, dist))
            self.logger.info('Form bond distance is greater than threshold.')
            self.logger.info('Now finished {}/{}'.format(self.num, self.count))
            self.finalize(start_time, 'arrange')
            return False, False
        else:
            self.logger.info('\nHere is the {} product.'.format(self.num))
            self.logger.info('Structure:\n{}'.format(str(reactant_mol.toNode())))
            self.logger.info('Structure:\n{}\n'.format(str(product_mol.toNode())))
            self.logger.info('Form bonds: {}\nDistance: {}'.format(self.form_bonds, dist))

            product_geometry = product_mol.toNode()
            reactant_geometry = reactant_mol.toNode()

            reactant_mol.setCoordsFromMol(reac_mol_copy)
            self.finalize(start_time, 'arrange')
            return reactant_geometry, product_geometry
            
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

    def getE(self, tmpdir, target = 'reactant.xyz'):
        """
        Here the energy is Eh (hartree)
        """
        input_path = os.path.join(tmpdir, target)
        with open(input_path, 'r') as f:
            lines = f.readlines()
        HeatofFormation = lines[1].strip().split()[1]
        return HeatofFormation

    def runXTB(self, tmpdir, target = 'reactant.xyz'):
        input_path = os.path.join(tmpdir, target)
        outname = '{}.xyz'.format(target.split('.')[0])
        output_path = os.path.join(tmpdir, 'xtbopt.xyz')
        config_path = os.path.join(os.path.dirname(os.path.dirname(tmpdir)), 'config')
        constraint_path = os.path.join(config_path, 'constraint.inp')
        new_output_path = os.path.join(tmpdir, outname)
        if self.constraint == None:
            p = Popen(['xtb', input_path, '--opt'])
            p.wait()
            os.rename(output_path, new_output_path)
        else:
            p = Popen(['xtb', '--opt', '--input', constraint_path, input_path])
            p.wait()
            os.rename(output_path, new_output_path)