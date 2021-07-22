import unittest
import libtracy as scsi
import os

t_dir = os.path.dirname(__file__)
lattice_dir = os.path.join(t_dir, 'lattices')


class _LatticeBasicFunctionality(unittest.TestCase):
    '''can it open a file and so on
    '''

    #: lattice filename: overload in base class
    lattice_filename = None

    def setUp(self):

        lat = scsi.LatticeType()
        lat.conf.trace        = False
        lat.conf.reverse_elem = not False
        lat.conf.mat_meth     = not False

        lat.conf.H_exact    = False; lat.conf.quad_fringe    = False
        lat.conf.Cavity_on  = False; lat.conf.radiation      = False
        lat.conf.emittance  = False; lat.conf.IBS            = False
        lat.conf.pathlength = False; lat.conf.bpm            = 0
        lat.conf.Cart_Bend  = False; lat.conf.dip_edge_fudge = True

        self.lat = lat

    def test0_readFile(self):
        '''Test if reading file from path works
        '''
        self.lat.Lat_Read(self.lattice_filename)

    def test0_readInit(self):
        '''Test if reading file and init works
        '''
        lat = self.lat
        lat.Lat_Read(self.lattice_filename)
        lat.Lat_Init()


class LatticeBasicFunctionality(_LatticeBasicFunctionality):
    lattice_filename = 'simple_fodo_lat.lat'


class LatticeBasicFunctionalityB2(_LatticeBasicFunctionality):
    '''

    Todo:
        Find a better file (e.g. no split elements, no 0 drift length)
    '''
    lattice_filename = os.path.join(lattice_dir,
                                    'b2_stduser_beamports_blm_tracy_corr')


del _LatticeBasicFunctionality
# no simple example file yet
del LatticeBasicFunctionality

if __name__ == '__main__':
    unittest.main()
