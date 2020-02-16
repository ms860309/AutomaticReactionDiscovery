from node import Node
import shutil
import os
from subprocess import Popen, PIPE
import gen3D
import pybel
import psutil


class mopac(object):

    def __init__(self, reac_mol, forcefield):
        self.reac_mol = reac_mol
        self.forcefield = forcefield

    def makeInputFile(self, InputFile, charge = 0, multiplicity = 'SINGLET', method = 'PM6'):
        """
        Create a directory folder called "tmp" for mopac calculation
        Create a input file called "input.mop" for mopac calculation
        """
        path = os.path.abspath("")
        tmpdir = os.path.join(path, 'tmp')
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.mkdir(tmpdir)
        input_path = os.path.join(tmpdir, 'input.mop')

        geometry = self.genInput(InputFile)
        with open(input_path, 'w') as f:
            f.write(" AUX LARGE CHARGE={} {} {}".format(charge, multiplicity, method))
            f.write("\n{}".format(geometry))

        self.runMopac(tmpdir)
        result = self.getHeatofFormation(tmpdir)
        """
        info = psutil.virtual_memory()
        print("cpu numbers : {}".format(psutil.cpu_count()))
        print("total memory : {}".format(info.total))
        print("memory used : {}".format(psutil.Process(os.getpid()).memory_info().rss))
        print("memory used percent : {}".format(info.percent))
        """
        return float(result)
    
    def genInput(self, InputFile):

        reac_mol = self.reac_mol
        Hatom = gen3D.readstring('smi', '[H]')
        ff = pybel.ob.OBForceField.FindForceField(self.forcefield)
        reac_mol_copy = reac_mol.copy()
        InputFile.gen3D(forcefield=self.forcefield, make3D=False)
        """
        arrange3D = gen3D.Arrange3D(reac_mol, InputFile)
        msg = arrange3D.arrangeIn3D()
        if msg != '':
            self.logger.info(msg)
        """
        ff.Setup(Hatom.OBMol)  # Ensures that new coordinates are generated for next molecule (see above)
        reac_mol.gen3D(make3D=False)
        ff.Setup(Hatom.OBMol)
        InputFile.gen3D(make3D=False)
        ff.Setup(Hatom.OBMol)

        geometry = InputFile.toNode()
        geometry = str(geometry)
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
        print(output)
        return output

    def getHeatofFormation(self, tmpdir):
        input_path = os.path.join(tmpdir, "input.out")
        with open(input_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.strip().startswith('FINAL HEAT OF FORMATION'):
                    HeatofFormation = line.split()[5]
        return HeatofFormation

    def runMopac(self, tmpdir):
        input_path = os.path.join(tmpdir, "input.mop")
        p = Popen(['mopac', input_path])
        p.wait()
