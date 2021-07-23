'''Test of constants

Rather a test that the correct names were used in the export
'''
import unittest
import libtracy as scsi


pi_est = 355 / 113
c0_est = 2.99e8
mu0_est = 4 * pi_est * 1e-7
eps0_est = 1 / c0_est**2 / mu0_est

# Wikipedia article
q_est = -1.602e-19


class ConstantsTest(unittest.TestCase):
    def test_hommax(self):
        '''HOM max variable available and reasonable
        '''
        hom_max = scsi.HOMmax
        self.assertGreaterEqual(hom_max, 0)
        self.assertLess(hom_max, 100)

    def test_speed_of_light(self):
        '''Speed of light: rough cross check

        Rather checking if existing
        '''
        self.assertAlmostEqual(scsi.c0, c0_est, delta=0.01e8)

    def test_vacuum_permeability(self):
        '''vacuum permeability cross check to 4 pi 1e-7
        '''
        self.assertAlmostEqual(scsi.mu_0, mu0_est, delta=1e-4 * mu0_est)

    def test_vacuum_permettivity(self):
        '''vacuum permetivitiy: cross check deduce from c0 and m0

        c0 = 1 / sqrt(mu0 * eps0)
        '''
        self.assertAlmostEqual(scsi.eps_0, eps0_est, delta=1e-2*eps0_est)

    def test_electron_charge(self):
        '''electron charge: positron charge to be  exact
        '''
        self.assertAlmostEqual(scsi.q_e, -q_est, delta=1e-5 * q_est)

    def test_mass_electron(self):
        '''mass electron in eV
        '''
        self.assertAlmostEqual(scsi.m_e, 511e3, delta=1)

    def test_hbar(self):
        '''Reduced Plank constant in electron volts

        estimate derived from

                   1       1  e^2
        \alpha = --------- -  -
                 4 \pi e0  c  hbar
        '''

        hbar_est_inv = 1/137 * 4 * pi_est * eps0_est * c0_est / q_est**2
        hbar_est = 1/hbar_est_inv
        hbar_est_eV = hbar_est / abs(q_est)
        self.assertAlmostEqual(scsi.h_bar, hbar_est_eV,
                               delta=1e-2 * hbar_est_eV)


if __name__ == '__main__':
    unittest.main()
