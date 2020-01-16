# standard library imports
import sys
from os import path

# third party
import numpy as np

# local application imports
sys.path.append(path.dirname( path.dirname( path.abspath(__file__))))
from utilities import *
from coordinate_systems import rotate

ELEMENT_TABLE = elements.ElementData()

class PES(object):
    """ PES object """

    @staticmethod
    def default_options():

        if hasattr(PES, '_default_options'): return PES._default_options.copy()
        opt = options.Options() 

        opt.add_option(
                key='lot',
                value=None,
                required=True,
                doc='Level of theory object')

        opt.add_option(
                key='ad_idx',
                value=0,
                required=True,
                doc='adiabatic index')

        opt.add_option(
                key='multiplicity',
                value=1,
                required=True,
                doc='multiplicity')

        opt.add_option(
                key="FORCE",
                value=None,
                required=False,
                doc='Apply a spring force between atoms in units of AU, e.g. [(1,2,0.1214)]. Negative is tensile, positive is compresive',
                )

        opt.add_option(
                key='mass',
                value=None,
                required=False,
                doc='Mass is sometimes required'
                )

        PES._default_options = opt
        return PES._default_options.copy()

    @classmethod
    def from_options(cls,**kwargs):
        """ Returns an instance of this class with default options updated from values in kwargs"""
        return cls(cls.default_options().set_values(kwargs))

    #TODO make kwargs
    @classmethod
    def create_pes_from(cls,PES,options={},copy_wavefunction=True):
        lot = type(PES.lot).copy(PES.lot,options,copy_wavefunction)
        return cls(PES.options.copy().set_values({
            "lot":lot,
            }))

    def __init__(self,
            options,
            ):
        """ Constructor """
        self.options = options

        self.lot = self.options['lot']
        self.ad_idx = self.options['ad_idx']
        self.multiplicity = self.options['multiplicity']
        self.FORCE = self.options['FORCE']
        self._dE=1000.
        #print ' PES object parameters:'
        #print ' Multiplicity:',self.multiplicity,'ad_idx:',self.ad_idx

    @property
    def dE(self):
        return self._dE

    @dE.setter
    def dE(self,value):
        self._dE = value

    @property
    def energy(self):
        if self.lot.E:
            return self.lot.search_PES_tuple(self.lot.E,self.multiplicity,self.ad_idx)[0][2]*units.KCAL_MOL_PER_AU
        else:
            return 0.

    def fill_energy_grid2d(self,xyz_grid):

        assert xyz_grid.shape[-1] == len(self.lot.geom)*3, "xyz nneds to be 3*natoms long"
        assert xyz_grid.ndim == 3, " xyzgrid needs to be a tensor with 3 dimensions"

        energies = np.zeros((xyz_grid.shape[0],xyz_grid.shape[1]))

        # E.G DO THIS OUTSIDE PES
        # get xyz in np format
        #xyz = manage_xyz.xyz_to_np(geoms[0])
        #xyz = xyz.flatten()

        # define the underlying grid coordinate basis 
        #xvec = np.zeros((len(xyz)))
        #yvec = np.zeros((len(xyz)))
        #xvec[0]+=1.
        #yvec[1]+=1.

        # create the grid from 0 to 1
        #nx, ny = (3, 2)
        #x=np.linspace(0,1,nx)
        #y=np.linspace(0,1,ny)
        #xv,yv = np.meshgrid(x,y)

        # actually create the xyz coordinates and save as a tensor
        #xresult = np.zeros((xv.shape[0],xv.shape[1],xvec.shape[0]))
        #rc=0
        #for xrow,yrow in zip(xv,yv):
        #    cc=0
        #    for xx,yy in zip(xrow,yrow):
        #        idx = (rc,cc)
        #        xresult[rc,cc,:] = xx*xvec + yy*yvec + xyz
        #        cc+=1
        #    rc+=1


        rc=0
        for mat in xyz_grid:
            cc=0
            for row in mat:
                xyz = np.reshape(row,(-1,3))
                energies[rc,cc] = self.lot.get_energy(xyz,self.multiplicity,self.ad_idx)
                cc+=1
            rc+=1
         
        return energies

    def get_energy(self,xyz):
        #if self.checked_input == False:
        #    self.check_input(geom)
        fdE=0.
        if self.FORCE is not None:
            for i in self.FORCE:
                atoms=[i[0],i[1]]
                force=i[2]
                diff = (xyz[i[0]]- xyz[i[1]])*units.ANGSTROM_TO_AU
                d = np.linalg.norm(diff)
                fdE +=  force*d*units.KCAL_MOL_PER_AU
        return self.lot.get_energy(xyz,self.multiplicity,self.ad_idx) +fdE

    #TODO this needs to be fixed
    def get_finite_difference_hessian(self,coords):
        hess = np.zeros((len(coords)*3,len(coords)*3))
        I = np.eye(hess.shape[0])
        for n,row in enumerate(I):
            print("on hessian product ",n)
            hess[n] = np.squeeze(self.get_finite_difference_hessian_product(coords,row))
        return hess

    #TODO this needs to be fixed
    def get_finite_difference_hessian_product(self,coords,direction):
        FD_STEP_LENGTH=0.001

        # format the direction
        direction = direction/np.linalg.norm(direction)
        direction = direction.reshape((len(coords),3))

        # fd step
        fdstep = direction*FD_STEP_LENGTH
        fwd_coords = coords+fdstep
        bwd_coords = coords-fdstep

        # calculate grad fwd and bwd in a.u. (Bohr/Ha)
        grad_fwd = self.get_gradient(fwd_coords)/units.ANGSTROM_TO_AU
        grad_bwd = self.get_gradient(bwd_coords)/units.ANGSTROM_TO_AU
    
        return (grad_fwd-grad_bwd)/(FD_STEP_LENGTH*2)

    @staticmethod
    def normal_modes(
            geom,       # Optimized geometry in au
            hess,       # Hessian matrix in au
            masses,     # Masses in au 
            ):
    
        """
        Params:
            geom ((natoms,4) np.ndarray) - atoms symbols and xyz coordinates
            hess ((natoms*3,natoms*3) np.ndarray) - molecule hessian
            masses ((natoms) np.ndarray) - masses
    
        Returns:
            w ((natoms*3 - 6) np.ndarray)  - normal frequencies
            Q ((natoms*3, natoms*3 - 6) np.ndarray)  - normal modes
    
        """
    
        # masses repeated 3x for each atom (unravels)
        m = np.ravel(np.outer(masses,[1.0]*3))
    
        # mass-weight hessian
        hess2 = hess / np.sqrt(np.outer(m,m))
    
        # Find normal modes (project translation/rotations before)
        B = rotate.vibrational_basis(geom, masses)
        h, U3 = np.linalg.eigh(np.dot(B.T,np.dot(hess2,B)))
        U = np.dot(B, U3)
    
        # TEST: Find normal modes (without projection translations/rotations)
        # RMP: Matches TC output for PYP - same differences before/after projection
        # h2, U2 = np.linalg.eigh(hess2)
        # h2 = h2[6:]
        # U2 = U2[:,6:]
        # for hval, hval2 in zip(h,h2):
        #     wval = np.sqrt(hval) / units['au_per_cminv']
        #     wval2 = np.sqrt(hval2) / units['au_per_cminv']
        #     print '%10.6E %10.6E %11.3E' % (wval, wval2, np.abs(wval - wval2))
    
        # Normal frequencies
        w = np.sqrt(h)
        # Imaginary frequencies
        w[h < 0.0] = -np.sqrt(-h[h < 0.0]) 
    
        # Normal modes
        Q = U / np.outer(np.sqrt(m), np.ones((U.shape[1],)))
    
        return w, Q 
    
    def get_gradient(self,xyz):
        tmp =self.lot.get_gradient(xyz,self.multiplicity,self.ad_idx)
        grad = np.reshape(tmp,(-1,1))
        if self.FORCE is not None:
            for i in self.FORCE:
                atoms=[i[0],i[1]]
                force=i[2]
                diff = (xyz[i[0]]- xyz[i[1]])*units.ANGSTROM_TO_AU
                t = (force/d/2.)  # Hartree/Ang
                savegrad = np.copy(grad)
                sign=1
                for a in [3*(i-1) for i in atoms]:
                    grad[a:a+3] += sign*t*diff.T
                    sign*=-1
        return grad

    def check_input(self,geom):
        atoms = manage_xyz.get_atoms(self.geom)
        elements = [ELEMENT_TABLE.from_symbol(atom) for atom in atoms]
        atomic_num = [ele.atomic_num for ele in elements]
        self.checked_input =True

if __name__ == '__main__':

    QCHEM=True
    PYTC=False
    if QCHEM:
        #from .qchem import QChem
        from level_of_theories.qchem import QChem
    elif PYTC:
        from level_of_theories.pytc import PyTC 
        import psiw
        import lightspeed as ls

    filepath='../../data/ethylene.xyz'
    geom=manage_xyz.read_xyz(filepath,scale=1)   
    if QCHEM:
        lot=QChem.from_options(states=[(1,0),(1,1)],charge=0,basis='6-31g(d)',functional='B3LYP',fnm=filepath)
    elif PYTC:
        ##### => Job Data <= #####
        states = [(1,0),(1,1)]
        charge=0
        nocc=7
        nactive=2
        basis='6-31gs'

        #### => PSIW Obj <= ######
        nifty.printcool("Build resources")
        resources = ls.ResourceList.build()
        nifty.printcool('{}'.format(resources))
        
        molecule = ls.Molecule.from_xyz_file(filepath)
        geom = psiw.geometry.Geometry.build(
            resources=resources,
            molecule=molecule,
            basisname=basis,
            )
        nifty.printcool('{}'.format(geom))
        
        ref = psiw.RHF.from_options(
             geometry= geom, 
             g_convergence=1.0E-6,
             fomo=True,
             fomo_method='gaussian',
             fomo_temp=0.3,
             fomo_nocc=nocc,
             fomo_nact=nactive,
             print_level=1,
            )
        ref.compute_energy()
        casci = psiw.CASCI.from_options(
            reference=ref,
            nocc=nocc,
            nact=nactive,
            nalpha=nactive/2,
            nbeta=nactive/2,
            S_inds=[0],
            S_nstates=[2],
            print_level=1,
            )
        casci.compute_energy()
        psiw = psiw.CASCI_LOT.from_options(
            casci=casci,
            rhf_guess=True,
            rhf_mom=True,
            orbital_coincidence='core',
            state_coincidence='full',
            )

        nifty.printcool("Build the pyGSM Level of Theory object (LOT)")
        lot=PyTC.from_options(states=[(1,0),(1,1)],job_data={'psiw':psiw},do_coupling=False,fnm=filepath) 

    pes = PES.from_options(lot=lot,ad_idx=0,multiplicity=1)
    geom=manage_xyz.read_xyz(filepath,scale=1)   
    coords= manage_xyz.xyz_to_np(geom)
    print(pes.get_energy(coords))
    print(pes.get_gradient(coords))

