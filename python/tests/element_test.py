import os, sys

#tracy_dir = os.getenv('TRACY_LIB')
#sys.path.append(tracy_dir+'/tracy/lib')

import unittest
import tracy.lib as scsi
import abc


class _ElementTestBasis(unittest.TestCase, metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def setUp(self):
        self.elem = None
        raise ValueError('Use derived class instead')

    def test00_name(self):
        '''test if element name can be read
        '''
        elem = self.elem
        elem.Name = 'name'
        self.assertEqual(elem.Name, 'name')

    def test01_reverse(self):
        '''Can element be reverted
        '''
        self.assertEqual(self.elem.Reverse, False)
        self.elem.Reverse = True
        self.assertEqual(self.elem.Reverse, True)

    def test02_repr(self):
        repr(self.elem)

    @unittest.skip
    def test10_print(self):
        '''test if element can be printed
        '''
        name = 'print_test'
        self.elem.Name = name
        self.elem.print()
        self.assertEqual(self.elem.Name, name)


class DriftTest(_ElementTestBasis):
    def setUp(self):
        self.elem = scsi.DriftType()
        self.elem.Name = 'TestDrift'


del _ElementTestBasis


if __name__ == '__main__':
    unittest.main()
