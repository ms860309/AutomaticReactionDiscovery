from node import Node
import shutil
import os
from subprocess import Popen, PIPE
import gen3D
import pybel
import psutil
import time
import logging
import util
class MopacError(Exception):
    """
    An exception class for errors that occur during mopac calculations.
    """
    pass

class mopac(object):

    def __init__(self, reac_mol, forcefield, reactant_bonds, product_bonds):
        self.reac_mol = reac_mol
        self.forcefield = forcefield
        self.reactant_bonds = reactant_bonds
        self.product_bonds = product_bonds
        log_level = logging.INFO
        process = psutil.Process(os.getpid())
        self.logger = util.initializeLog(log_level, os.path.join(os.getcwd(), 'ARD.log'), logname='main')
        self.logger.info('\nARD initiated on ' + time.asctime() + '\n')
        self.logger.info('memory usage: {}'.format(process.memory_percent()))

    def mopac_get_H298(self, InputFile, charge = 0, multiplicity = 'SINGLET', method = 'PM6'):
        """
        Create a directory folder called "tmp" for mopac calculation
        Create a input file called "input.mop" for mopac calculation
        """

        tmpdir = os.path.join(os.path.dirname(os.getcwd()), 'tmp')
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.mkdir(tmpdir)
        input_path = os.path.join(tmpdir, 'input.mop')

        geometry = self.genInput(InputFile)
        with open(input_path, 'w') as f:
            f.write(" AUX LARGE CHARGE={} {} {} GEO-OK".format(charge, multiplicity, method))
            f.write("\n{}".format(geometry))
        start_time = time.time()
        self.runMopac(tmpdir)
        result = self.getHeatofFormation(tmpdir)
        self.finalize(start_time, 'mopac')
        """
        info = psutil.virtual_memory()
        print("cpu numbers : {}".format(psutil.cpu_count()))
        print("total memory : {}".format(info.total))
        print("memory used : {}".format(psutil.Process(os.getpid()).memory_info().rss))
        print("memory used percent : {}".format(info.percent))
        """
        return float(result)
    
    def genInput(self, InputFile):
        start_time = time.time()

        reac_mol = self.reac_mol
        Hatom = gen3D.readstring('smi', '[H]')
        ff = pybel.ob.OBForceField.FindForceField(self.forcefield)
        reac_mol_copy = reac_mol.copy()
        reac_mol_copy, InputFile_copy= reac_mol.copy(), InputFile.copy()
        reac_mol.gen3D(forcefield=self.forcefield, make3D=False)
        InputFile.gen3D(forcefield=self.forcefield, make3D=False)
        
        try:
            arrange3D = gen3D.Arrange3D(reac_mol, InputFile, self.reactant_bonds, self.product_bonds)
            msg = arrange3D.arrangeIn3D()
            if msg != '':
                print(msg)
        except:
            reac_mol, InputFile = reac_mol_copy, InputFile_copy

        ff.Setup(Hatom.OBMol)  # Ensures that new coordinates are generated for next molecule (see above)
        reac_mol.gen3D(make3D=False)
        ff.Setup(Hatom.OBMol)
        InputFile.gen3D(make3D=False)
        ff.Setup(Hatom.OBMol)

        geometry = InputFile.toNode()
        geometry = str(geometry)
        self.logger.info('\nStructure:\n{}\n'.format(str(geometry)))
        geometry = geometry.splitlines()
        output = []
        for i in geometry:
            i_list = i.split()
            atom = i_list[0] + " "
            k = i_list[1:] + [""]
            l = " 1 ".join(k)
            out = atom + l
            output.append(out)
        output = "\n".join(output)
        reac_mol.setCoordsFromMol(reac_mol_copy)

        self.finalize(start_time, 'arrange')
        
        return output

    def finalize(self, start_time, jobname):
        """
        Finalize the job.
        """
        self.logger.info('Total {} run time: {:.1f} s'.format(jobname, time.time() - start_time))

    def getHeatofFormation(self, tmpdir):
        """
        if Error return False, which HF may be 0.0
        """
        input_path = os.path.join(tmpdir, "input.out")
        with open(input_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                try:
                    if line.strip().startswith('FINAL HEAT OF FORMATION'):
                        HeatofFormation = line.split()[5]
                except:
                    HeatofFormation = False
                    raise MopacError('Mopac calculation error')
        return HeatofFormation

    def runMopac(self, tmpdir):
        input_path = os.path.join(tmpdir, "input.mop")
        p = Popen(['mopac', input_path])
        p.wait()
