from rmgpy.molecule import Molecule

import gen3D
class _filter(object):

    def __init__(self, reac_mol, product_mols):
        self.reac_mol = reac_mol
        self.prod_mols = product_mols

    def nRing(self):
        """
        input rmg molecule object and return number of carbons on the rings, if no ring return 0
        BTW,
        naphthalene = 2 rings
        Biphenyl = 1 rings
        """
        small_num_ring = []
        for idx, i in enumerate(self.prod_mols):
            i = i.toRMGMolecule()

            try:
                monoring_atoms = i.getMonocyclicRings()[0]
                length = len(monoring_atoms)
                small_num_ring.append(length)
            except:
                small_num_ring.append(0)

            try:
                polyring_atoms = i.getPolycyclicRing()
            except:
                continue
        print(small_num_ring)
        return small_num_ring