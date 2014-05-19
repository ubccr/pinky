"""
FROWNS LICENSE

Copyright (c) 2001-2003, Brian Kelley
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met: 

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer. 
    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following disclaimer
      in the documentation and/or other materials  provided with the
      distribution. 
    * Neither the name of Brian Kelley nor the names of frowns
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
# Build a simple Molecule object given the events from the Smiles
# tokenizer.

import string
import handler
import weakref
from ..mol import Atom, Bond, Molecule

# bondlookup is of the form
# textSymbol, bondsymbol, bondorder, bondtype, equiv class, stereo

STEREO_NONE = None
STEREO_UP = "UP"
STEREO_DOWN = "DOWN"

BONDLOOKUP = {'-': ('-', 1, 1, 1, STEREO_NONE),
              '=': ('=', 2, 2, 2, STEREO_NONE),
              '#': ('#', 3, 3, 3, STEREO_NONE),
              '\\': ('\\',1, 1, 1, STEREO_DOWN),
              '/': ('/' ,1, 1, 1, STEREO_UP),
              ':':(':', 1.5, 4, 4, STEREO_NONE),              
              }

def get_symbol_aromatic(text):
    if text[0] in "cnosp":
        return string.upper(text), 1
    return text, 0

def normalize_closure(text):
    if text[:1] == "%":
        return int(text[1:])
    return int(text)

implicit_bond = -123

class DummyVFGraph:
    def __init__(self):
        self.atoms = -1
    def InsertNode(self, node):
        self.atoms += 1
        return self.atoms
    def InsertEdge(self, index1, index2, bond):
        pass
    
class BuildMol(handler.TokenHandler):

    def begin(self):
        self.closures = {}

        self.atoms = []
        self.bonds = []
        self._atom = None
        self._prev_atoms = []

        # None occurs after a '.'
        # implicit_bond means implicit single bond
        self._pending_bond = None

    def end(self):
        if len(self._prev_atoms) >= 2:
            raise AssertionError("Missing ')'")
        if self._pending_bond not in [implicit_bond, None]:
            raise AssertionError("Missing an atom after the bond")
        if self.closures:
            raise AssertionError("Missing closures for %s" %
                                 (self.closures.keys(),))
        self.mol = Molecule(self.atoms, self.bonds)
    
    def add_token(self, field, pos, text):
        getattr(self, "do_" + field)(text)

    def add_atom(self, atom):
        atoms = self.atoms
        atom.index = len(atoms)
        atoms.append(atom)
        
        if self._pending_bond == implicit_bond:
            # Implicit single or aromatic bond
            self._pending_bond = Bond()

        if self._pending_bond is not None:
            bond = self._pending_bond
            prev_atom = self._prev_atoms[-1]
            bond.atoms[:] = [prev_atom, atom]
            ##self.mol.add_bond(bond, prev_atom, atom)
            bond.atoms = [prev_atom, atom]
            atom.bonds.append(bond)
            prev_atom.bonds.append(bond)
            atom.oatoms.append(prev_atom)
            prev_atom.oatoms.append(atom)
            self.bonds.append(bond)
            
        self._pending_bond = implicit_bond
        if not self._prev_atoms:
            self._prev_atoms.append(atom)
        else:
            self._prev_atoms[-1] = atom

        #self.mol.atoms.append(atom)
        
    def do_raw_atom(self, text):
        atom = Atom()
        symbol, atom.aromatic = get_symbol_aromatic(text)
        atom.set_symbol(symbol)
        self.add_atom(atom)

    def do_open_bracket(self, text):
        self._atom = Atom()
        self._atom.has_explicit_hcount = True

    def do_weight(self, text):
        self._atom.weight = int(text)

    def do_element(self, text):
        symbol, self._atom.aromatic = get_symbol_aromatic(text)
        self._atom.set_symbol(symbol)
        
    def do_chiral_count(self, text):
        #print "setting chirality", self._atom, int(text[1:])
        self._atom.chirality = int(text[1:])

    def do_chiral_named(self, text):
        self._atom.chiral_class = text[1:3]
        self._atom.chirality = int(text[3:])

    def do_chiral_symbols(self, text):
        self._atom.chiral_class = len(text)

    def do_hcount(self, text):
        if text == "H":
            self._atom.explicit_hcount = 1
        else:
            self._atom.explicit_hcount = int(text[1:])

    def do_positive_count(self, text):
        self._atom.charge = int(text[1:])

    def do_positive_symbols(self, text):
        self._atom.charge = len(text)

    def do_negative_count(self, text):
        self._atom.charge = -int(text[1:])

    def do_negative_symbols(self, text):
        self._atom.charge = -len(text)

    def do_close_bracket(self, text):
        self.add_atom(self._atom)
        self._atom = None

    def do_bond(self, text):
        assert self._pending_bond in (implicit_bond, None)
        symbol, bondorder, bondtype, equiv_class, stereo = BONDLOOKUP[text]
        # if the bond came in as aromatic (which it
        #  CAN'T!))
        if bondtype == 4:
            assert 0, "Bond's shouldn't come in as ':'"
            fixed = 0
        else:
            fixed = 1
        bond = Bond(text, bondorder, bondtype, fixed, stereo)
        bond.equiv_class = equiv_class
        self._pending_bond = bond

    def do_dot(self, text):
        assert self._pending_bond in (implicit_bond, None)
        self._pending_bond = None

    def do_closure(self, text):
        num = normalize_closure(text)
        if self.closures.has_key(num):
            prev_atom, bond = self.closures[num]
            del self.closures[num]

            assert self._pending_bond is not None, "Can't happen"
            
            if self._pending_bond is not implicit_bond and \
               bond is not implicit_bond and \
               self._pending_bond.symbol != "-":  # according to toolkit

                # need to verify they are compatible
                prev_symbol = bond.symbol
                symbol = self._pending_bond.symbol
                if (prev_symbol == symbol) or \
                   (prev_symbol == "/" and symbol == "\\") or \
                   (prev_symbol == "\\" and symbol == "/"):
                    pass
                else:
                    raise AssertionError("bond types don't match")
            elif bond is implicit_bond and self._pending_bond is not implicit_bond:
                # see if one of the bonds is not implicit and keep it
                bond = self._pending_bond
            elif bond is implicit_bond:
                # both are implicit so make a new one
                bond = Bond()

            bond._closure = 1
            atom = self._prev_atoms[-1]
            if prev_atom is atom:
                raise AssertionError("cannot close a ring with itself")
            bond.atoms[:] = [prev_atom, atom]
            prev_atom._closure = 1
            atom._closure = 1
            ##self.mol.add_bond(bond, prev_atom, atom)
            
            bond.atoms = [prev_atom, atom]
            atom.bonds.append(bond)
            prev_atom.bonds.append(bond)
            atom.oatoms.append(prev_atom)
            prev_atom.oatoms.append(atom)
            self.bonds.append(bond)

        else:
            self.closures[num] = (self._prev_atoms[-1], self._pending_bond)
        self._pending_bond = implicit_bond
            

    def do_open_branch(self, text):
        self._prev_atoms.append(self._prev_atoms[-1])
    
    def do_close_branch(self, text):
        self._prev_atoms.pop()
