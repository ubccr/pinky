from unittest import TestCase
from nose.tools import *
from pinky.smiles import smilin
from pinky.fingerprints import ecfp

class SmilesTestCase(TestCase):
    def test_ecfp(self):
        """Test Butyramide (example used in ECFP paper)"""
        mol = smilin('CCCC(=O)N')

        # ECFP_0
        fp = ecfp(mol, radius=0)
        assert len(fp) == 5

        # ECFP_2
        fp = ecfp(mol, radius=1)
        assert len(fp) == 11

        # ECFP_4
        fp = ecfp(mol, radius=2)
        assert len(fp) == 14

        # ECFP_6
        fp = ecfp(mol, radius=3)
        assert len(fp) == 14
