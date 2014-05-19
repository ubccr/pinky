import os
from unittest import TestCase
from nose.tools import *
from pinky.smiles import smilin
from pinky.exceptions import PinkyError

class SmilesTestCase(TestCase):
    def setUp(self):
        self.path = os.path.dirname(os.path.abspath(__file__))

    def test_smilein(self):
        smiles_strings = []
        with open("{}/smiles.txt".format(self.path), 'rb') as fh:
            for line in fh:
                smiles_strings.append(line.strip())

        for smile in smiles_strings:
            mol = smilin(smile)
            out = mol.arbsmiles()
            can = mol.cansmiles()
            for bond in mol.bonds:
                print bond.symbol, bond.bondorder, bond.bondtype, bond.fixed
            for atom in mol.atoms:
                print atom, atom.sumBondOrders()
                
            print smile, out, can
