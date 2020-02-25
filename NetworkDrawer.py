import os
import pybel
import numpy as np
from rmgpy.molecule.molecule import Atom, Bond, Molecule
from rmgpy.molecule.draw import MoleculeDrawer
from rmgpy.molecule import element
from rmgpy.molecule.element import get_element
try:
    import cairocffi as cairo
except ImportError:
    try:
        import cairo
    except ImportError:
        cairo = None

def create_new_surface(file_format, target=None, width=1024, height=768):
    """
    Create a new surface of the specified `file_format`:
        "png" for :class:`ImageSurface`
        "svg" for :class:`SVGSurface`
        "pdf" for :class:`PDFSurface`
        "ps" for :class:`PSSurface`
    The surface will be written to the `target` parameter , which can be a
    path to save the surface to, or file-like object with a `write()` method.
    You can also optionally specify the `width` and `height` of the generated
    surface if you know what it is; otherwise a default size of 1024 by 768 is
    used.
    """
    file_format = file_format.lower()
    if file_format == 'png':
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, int(width), int(height))
    elif file_format == 'svg':
        surface = cairo.SVGSurface(target, width, height)
    elif file_format == 'pdf':
        surface = cairo.PDFSurface(target, width, height)
    elif file_format == 'ps':
        surface = cairo.PSSurface(target, width, height)
    else:
        raise ValueError(
            'Invalid value "{0}" for type parameter; valid values are "png", "svg", "pdf", and "ps".'.format(type))
    return surface

class NetworkDrawer(object):

    def __init__(self, options=None):
        self.options = {
            'structures': True,
            'fontFamily': 'sans',
            'fontSizeNormal': 12,
            'Eunits': 'kJ/mol',
            'padding': 16,
            'wellWidth': 64,
            'wellSpacing': 64,
            'Eslope': 1.5,
            'TSwidth': 16,
            'E0offset': 0.0,
        }
        if options:
            self.options.update(options)

        self.network = None
        self.left = 0.0
        self.top = 0.0
        self.right = 0.0
        self.bottom = 0.0
        self.surface = None
        self.cr = None    

    def clear(self):
        self.network = None
        self.left = 0.0
        self.top = 0.0
        self.right = 0.0
        self.bottom = 0.0
        self.surface = None
        self.cr = None

    def getEnergy(self):
        """
        PES energy is SCF energy, which could be calculate from qchem.
        """
        energy = {}
        dirs = os.listdir(os.path.join(os.path.abspath(""), 'reactions'))
        for i in dirs:
            energy_path = os.path.join(os.path.join(os.path.join(os.path.abspath(""), 'reactions'), i), "energy.out")
            with open(energy_path, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if line.strip().startswith('SCF   energy in the final basis set'):
                        SCF_Energy = line.split()[-1]
                        energy[i] = SCF_Energy
        return energy

    def getEnergyRange(self):
        """
        Get energy range and return reactant energy, minimum energy and maximum energy
        """
        energy = self.getEnergy()
        reac_energy = energy['reactant']
        min_energy =energy[min(energy, key=energy.get)]
        max_energy =energy[max(energy, key=energy.get)]

        return energy, float(reac_energy), float(min_energy), float(max_energy)

    def xyz_draw(self, file_format='png'):
        dirs = os.listdir(os.path.join(os.path.abspath(""), 'reactions'))
        for i in dirs:
            xyz_path = os.path.join(os.path.join(os.path.join(os.path.abspath(""), 'reactions'), i), "draw.xyz")
            mol = next(pybel.readfile('xyz', xyz_path))
            draw_path = os.path.join(os.path.join(os.path.join(os.path.abspath(""), 'reactions'), i), "picture.png")
            if os.path.isfile(draw_path):
                os.remove(draw_path)
            #mol.draw(show=False, filename=draw_path, update=False, usecoords=False)
            rmg_mol = self.pybel_to_rmg(mol)
            molecule_drawer = MoleculeDrawer()
            surf, content, rect = molecule_drawer.draw(rmg_mol, file_format='png', target = draw_path)
            

    def pybel_to_rmg(self, pybel_mol):
        """
        Convert Pybel molecule to RMG molecule but ignore charge,
        multiplicity, and bond orders.
        """
        mol = Molecule()
        for pybel_atom in pybel_mol:
            element = get_element(pybel_atom.atomicnum)
            atom = Atom(element=element, coords=np.array(pybel_atom.coords))
            mol.vertices.append(atom)
        for obbond in pybel.ob.OBMolBondIter(pybel_mol.OBMol):
            begin_idx = obbond.GetBeginAtomIdx() - 1  # Open Babel indexes atoms starting at 1
            end_idx = obbond.GetEndAtomIdx() - 1
            bond = Bond(mol.vertices[begin_idx], mol.vertices[end_idx])
            mol.add_bond(bond)
        return mol

    def draw_mol(self, cr, x, y, target, file_format='png'):

        xyz_path = os.path.join(os.path.join(os.path.join(os.path.abspath(""), 'reactions'), target), "draw.xyz")
        mol = next(pybel.readfile('xyz', xyz_path))
        draw_path = os.path.join(os.path.join(os.path.join(os.path.abspath(""), 'reactions'), target), "picture.png")
        if os.path.isfile(draw_path):
            os.remove(draw_path)
        #mol.draw(show=False, filename=draw_path, update=False, usecoords=False)
        rmg_mol = self.pybel_to_rmg(mol)
        molecule_drawer = MoleculeDrawer()
        surf, content, rect = molecule_drawer.draw(rmg_mol, file_format='png', target = draw_path)
        molecule_drawer.render(cr, offset=(x, y))


    def draw(self, network):
        self.xyz_draw()
        #self.network = network
        #get reactant energy e0_min and e0_max
        energy, e0_reac, e0_min, e0_max = self.getEnergyRange()
        #get total ended products
        nend_prods = len(network['products'])
        #get index which have intermediate
        intermediate = []
        for i in network['products']:
            if len(i) > 1:
                intermediate.append(network['products'].index(i))
        #get max num of intermediates to calculate length
        max_num_intermediate = len(max(network['products'], key=len))
        #determine height required for drawing
        _height = abs(e0_max - e0_min)*627.5095
        _height = 700/_height
        e_reac = 768 - abs((e0_max - e0_reac)*627.5095*_height)
        e_0000 = 768 - abs((e0_max - float(energy['0000']))*627.5095*_height)
        e_0001 = 768 - abs((e0_max - float(energy['0001']))*627.5095*_height)
        e_0002 = 768 - abs((e0_max - float(energy['0002']))*627.5095*_height)
        e_0003 = 768 - abs((e0_max - float(energy['0003']))*627.5095*_height)
        e_offest = 20
        #draw the surface
        surface = create_new_surface(file_format='pdf', target=os.path.join(os.path.abspath(""), 'network.pdf'), width=1024, height=768)
        cr = cairo.Context(surface)
        # Draw molecule
        #Some global setting
        cr.select_font_face("sans")
        cr.set_font_size(self.options['fontSizeNormal'])
        #Fill the background with white
        cr.set_source_rgba(1.0, 1.0, 1.0, 1.0)
        cr.paint()
        #draw
        #reac
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.set_line_width(2.0)
        cr.move_to(40-30, e_reac - e_offest)
        cr.line_to(40+30, e_reac - e_offest)
        cr.move_to(30, e_reac - e_offest-15)
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.set_line_width(4.0)
        cr.show_text(str(abs((e0_reac - e0_reac)*627.5095)))
        self.draw_mol(cr, 30, 110, 'reactant')
        cr.stroke()
        #intermedia (0000)
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.set_line_width(2.0)
        cr.move_to(40+472-30, e_0000 - e_offest)
        cr.line_to(40+472+30, e_0000 - e_offest)
        cr.move_to(40+472-50, e_0000 - e_offest-15)
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.set_line_width(4.0)
        cr.show_text(str(abs((float(energy['0000']) - e0_reac)*627.5095)))
        self.draw_mol(cr, 40+472-10, e_0000 - e_offest, '0000')
        cr.stroke()
        #prod (0001 0002 0003)
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.set_line_width(2.0)
        cr.move_to(40+472*2-30, e_0001 - e_offest)
        cr.line_to(40+472*2+30, e_0001 - e_offest)
        cr.move_to(40+472*2-50, e_0001 - e_offest-15)
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.set_line_width(4.0)
        cr.show_text(str((float(energy['0001']) - e0_reac)*627.5095))
        self.draw_mol(cr, 40+472*2-25, e_0001 - e_offest+10, '0001')
        cr.stroke()
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.set_line_width(2.0)
        cr.move_to(40+472*2-30, e_0002 - e_offest)
        cr.line_to(40+472*2+30, e_0002 - e_offest)
        cr.move_to(40+472*2-50, e_0002 - e_offest-15)
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.set_line_width(4.0)
        cr.show_text(str((float(energy['0002']) - e0_reac)*627.5095))
        self.draw_mol(cr, 40+472*2-25, e_0002 - e_offest+5, '0002')
        cr.stroke()
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.set_line_width(2.0)
        cr.move_to(40+472*2-30, e_0003 - e_offest)
        cr.line_to(40+472*2+30, e_0003 - e_offest)
        cr.move_to(40+472*2-50, e_0003 - e_offest-15)
        cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
        cr.set_line_width(4.0)
        cr.show_text(str((float(energy['0003']) - e0_reac)*627.5095))
        self.draw_mol(cr, 40+472*2-25, e_0003 - e_offest+10, '0003')
        cr.stroke()
        ## Draw cubic Bezier spline connecting reactants and products through TS
        cr.set_source_rgba(0.0, 0.0, 0.0, 0.5)
        cr.set_line_width(1.0)
        cr.move_to(40+30, e_reac - e_offest)
        cr.line_to(40+472*2-30, e_0001 - e_offest)
        cr.stroke()
        cr.set_source_rgba(0.0, 0.0, 0.0, 0.5)
        cr.set_line_width(1.0)
        cr.move_to(40+30, e_reac - e_offest)
        cr.line_to(40+472*2-30, e_0002 - e_offest)
        cr.stroke()
        cr.set_source_rgba(0.0, 0.0, 0.0, 0.5)
        cr.set_line_width(1.0)
        cr.move_to(40+30, e_reac - e_offest)
        cr.line_to(40+472-30, e_0000 - e_offest)
        cr.stroke()
        cr.set_source_rgba(0.0, 0.0, 0.0, 0.5)
        cr.set_line_width(1.0)
        cr.move_to(40+472+30, e_0000 - e_offest)
        cr.line_to(40+472*2-30, e_0003 - e_offest)
        cr.stroke()
        ## Draw label


network = {'reactants' : ['reactant'], 'products' : [['0001', '0003'], ['0002'], ['0003']]}
NetworkDrawer().draw(network)