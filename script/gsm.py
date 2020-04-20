# standard library imports
import sys
import os
from os import path
import importlib

#third party
import argparse
import textwrap
from subprocess import Popen
#local application imports
sys.path.append(os.path.join(os.path.join(path.dirname( path.dirname( path.abspath(__file__))),'pyGSM'), 'pygsm'))



def main():
    parser = argparse.ArgumentParser(
        description="Automatically run Single-Ended String Method(SSM)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''\
                Example of use:
                --------------------------------
                python -xyzfile yourfile.xyz -add_bonds add_bonds.txt
                ''')
        )
    parser.add_argument('-xyzfile', help='XYZ file containing reactant and, if DE-GSM, product.',  required=True)
    parser.add_argument('-isomers', help='driving coordinate file', required=False)
    parser.add_argument('-coordinate_type',type=str,default='DLC',help='Coordinate system (default %(default)s)',choices=['TRIC','DLC','HDLC'])
    parser.add_argument('-lot_inp_file',type=str,default='qstart', help='external file to specify calculation e.g. qstart,gstart,etc. Highly package specific.',required=True)
    args = parser.parse_args()

    cmd = 'gsm -xyzfile {} -mode SE_GSM -package QChem -isomers {} -lot_inp_file {} -coordinate_type DLC'.format(args.xyzfile, args.isomers, args.lot_inp_file)
    p = Popen([cmd], shell = True)
    p.wait()

if __name__=='__main__':
    main()