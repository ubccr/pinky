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



Cycle(atoms) -> cycle object for a ring.

   cycle.atoms -> atoms around a ring
   cycle.bonds -> bonds aroung a ring

   cycle.rotate(atom) -> rotate the ring so that atom is in front
"""
class CycleError(Exception):
    pass

class Cycle:
    def __init__(self, atoms, bonds, aromatic=0):
        """(atoms)->create a cycle object
        assumes that the atoms are in traversal order around the ring.
        That is [a1, a2, a3] means that there is a bond in the cycle
        between a1 and a2 and a2 and a3 and a3 and a1"""
        self.atoms = atoms[:]
        self.bonds = bonds[:]
        self.aromatic = aromatic
        for atom in self.atoms:
            atom.rings.append(self)
            
        for bond in self.bonds:
            bond.rings.append(self)
        
    def __len__(self):
        return len(self.atoms)

    def rotate(self, atom):
        """(atom)->start the cycle at position atom, assumes
        that atom is in the cycle"""
        try:
            index = self.atoms.index(atom)
        except ValueError:
            raise CycleError("atom %s not in cycle"%(atom))

        self.atoms = self.atoms[index:] + self.atoms[:index]
        self.bonds = self.bonds[index:] + self.bonds[:index]

    def clone(self):
        return Cycle(self.atoms, self.bonds, self.aromatic)

    def set_aromatic(self):
        """set the cycle to be an aromatic ring"""
        #XXX FIX ME
        # this probably shouldn't be here
        for atom in self.atoms:
            atom.aromatic = 1
            
        for bond in self.bonds:
            bond.aromatic = 1
            bond.bondorder = 1.5
            bond.bondtype = 4
            bond.symbol = ":"
            bond.fixed = 1

        self.aromatic = 1
