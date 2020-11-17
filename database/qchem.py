import os
import numpy as np


class QChemError(Exception):
    pass


class QChem(object):
    """
    Extract the information from QChem output.
    The methods of this class have only been validated for Q-Chem 5.1 DFT
    calculations.
    """

    def __init__(self, outputfile=None):
        self.logfile = outputfile

        if outputfile is None:
            self.log = None
        else:
            with open(outputfile) as f:
                self.log = f.read().splitlines()
                for line in self.log:
                    if 'fatal error' in line:
                        raise QChemError(f'Q-Chem job {outputfile} had an error!')

    def get_energy(self, first=False):
        if first:
            iterable = self.log
        else:
            iterable = reversed(self.log)
        for line in iterable:
            if 'total energy' in line:  # Double hybrid methods
                return float(line.split()[-2])
            elif 'energy in the final basis set' in line:  # Other DFT methods
                return float(line.split()[-1])
        else:
            raise QChemError(f'Energy not found in {self.logfile}')

    def get_geometry(self, first=False):
        if first:
            iterable = range(len(self.log))
        else:
            iterable = reversed(range(len(self.log)))
        for i in iterable:
            line = self.log[i]
            if 'Standard Nuclear Orientation' in line:
                symbols, coords = [], []
                for line in self.log[(i+3):]:
                    if '----------' not in line:
                        data = line.split()
                        symbols.append(data[1])
                        coords.append([float(c) for c in data[2:]])
                    else:
                        return symbols, np.array(coords)
        else:
            raise QChemError(f'Geometry not found in {self.logfile}')

    def get_frequencies(self):
        freqs = []
        for line in reversed(self.log):
            if 'Frequency' in line:
                freqs.extend([float(f) for f in reversed(line.split()[1:])])
            elif 'VIBRATIONAL ANALYSIS' in line:
                freqs.reverse()
                return np.array(freqs)
        else:
            raise QChemError(f'Frequencies not found in {self.logfile}')

    def get_normal_modes(self):
        modes = []
        for i in reversed(range(len(self.log))):
            line = self.log[i]
            if 'Raman Active' in line:
                mode1, mode2, mode3 = [], [], []
                for line in self.log[(i+2):]:
                    if 'TransDip' not in line:
                        vals = line.split()[1:]
                        mode1.append([float(v) for v in vals[:3]])
                        mode2.append([float(v) for v in vals[3:6]])
                        mode3.append([float(v) for v in vals[6:]])
                    else:
                        modes.extend([np.array(mode3), np.array(mode2), np.array(mode1)])
                        break
            elif 'VIBRATIONAL ANALYSIS' in line:
                modes.reverse()
                return modes
        else:
            raise QChemError(f'Normal modes not found in {self.logfile}')

    def get_zpe(self):
        for line in reversed(self.log):
            if 'Zero point vibrational energy' in line:
                return float(line.split()[-2]) / 627.5095  # Convert to Hartree
        else:
            raise QChemError(f'ZPE not found in {self.logfile}')


logfile='/mnt/d/irc/QDYIWPSSUWGGSS-UHFFFAOYSA-N_5/TS/ts.out'
q = QChem(outputfile=logfile)
zpe = q.get_zpe()
normal_mode = q.get_normal_modes()
freq = q.get_frequencies()
geometry = q.get_geometry()
energy = q.get_energy()
print(energy)