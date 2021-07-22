import unittest
import libtracy as scsi


class _LatticeBasicFunctionality(unittest.TestCase):
    '''can it open a file and so on
    '''

    #: lattice filename: overload in base class
    lattice_filename = None

    def test0_readFile(self):
        '''Test if reading file from path works
        '''
        lat = scsi.LatticeType()
        lat.Lat_Read(self.lattice_filename)

    def test0_readInit(self):
        '''Test if reading file and init works
        '''
        lat = scsi.LatticeType()
        lat.Lat_Read(self.lattice_filename)
        lat.Lat_Init()


class LatticeBasicFunctionality(_LatticeBasicFunctionality):
    lattice_filename = 'simple_fodo_lat.lat'


class LatticeBasicFunctionalityB2(_LatticeBasicFunctionality):
    lattice_filename = ''


del _LatticeBasicFunctionality


if __name__ == '__main__':
    unittest.main()
