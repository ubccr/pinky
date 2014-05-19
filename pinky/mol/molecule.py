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

# sheesh, how much indirection can there be?
# I'm longing for Mitch's Graph.T right about
# now :)
from ..canonicalization import disambiguate, traverse
from ..canonicalization.equiv import compute_equiv_class


# for now these molecules will be immutable
# Fix Me, don't use _atoms use atom.parent
class Molecule(object):
    def __init__(self, atoms=None, bonds=None):
        self.atoms = atoms or []
        self.bonds = bonds or []
        index = 0
        for atom in self.atoms:
            if atom.parent:
                ValueError("Atom %s already has a parent"%atom)
            atom.parent = self
            atom.index = index
            index += 1
            if atom.chiral_class:
                atom.setchival(None, atom.chiral_class)
        index = 0
        for bond in self.bonds:
            if bond.parent:
                ValueError("Atom %s already has a parent"%atom)                
            bond.parent = self
            bond.index = index
            index += 1
            
        self.fields = {}
        self.cycles = []

        self._canonical = None
        self.dirty = 1   # a molecule is dirty
                         #  when it needs to be recanonicalized
                         #  or various properties need to be updated

        self.name = ""
        self.vfgraph = None  # place holder for a vflib matching graph
                             # must be rebuilt if None
                            
    def add_atom(self, atom):
        if atom.parent:
            raise "Already have atom", atom

        # set the parent of the atom to be this object
        # set up all the relevant datastructures for
        # the atom
        atom.parent = self
        atom.index = len(self.atoms)
        
        self.atoms.append(atom)
        self.dirty = 1
        if self.vfgraph:
            index = self.vfgraph.InsertNode(atom)
            assert index == atom.index
        self.vfgraph = None

    def add_bond(self, bond, atom1, atom2):
        if atom1.parent is not self:
            raise "Atom", atom1, "is not in the molecule"
        if atom2.parent is not self:
              raise "Atom", atom2, "is not in the molecule"          
        if bond.parent:
            raise "Bond", bond, "already has a parent"
        
        bond.parent = self
        bond.index = len(self.bonds)        
        self.bonds.append(bond)
        bond.atoms = [atom1, atom2]

        atom1.bonds.append(bond)
        atom2.bonds.append(bond)
        atom1.oatoms.append(atom2)
        atom2.oatoms.append(atom1)

        self.dirty = 1
        if self.vfgraph:
            vfgraph = self.vfgraph
            vfgraph.InsertEdge(index1, index2, bond)
            vfgraph.InsertEdge(index2, index1, bond)

    def remove_bond(self, bond):
        self.bonds.remove(bond)
        bond.destroy()
        self.dirty = 1
        self.vfgraph = None
        index = 0
        for bond in self.bonds:
            bond.index = index
            index += 1

    def remove_atom(self, atom):
        if atom.parent is not self:
            raise "Atom not in molecule"

        for bond in atom.bonds[:]:
            self.remove_bond(bond)

        self.atoms.remove(atom)        
        atom.parent = None
        atom.destroy()
        self.dirty = 1
        self.vfgraph = None
        index = 0
        for atom in self.atoms:
            atom.index = index
            index += 1
            
    def cansmiles(self, isomeric=0):
        if isomeric: draw = traverse.drawIsomeric
        else: draw = traverse.draw
        
        if self.dirty:
            # XXX FIX ME
            # Move equivalence_class setter outside?
            # recompute the equivalence classes
            for atom in self.atoms:
                atom.equiv_class = compute_equiv_class(atom)
            disambiguate.FreedDisambiguate(self)
            self._canonical, self.canonical_list = traverse.draw(self)
            self.dirty = 0

        if isomeric:
            return traverse.drawIsomeric(self)[0]
        return self._canonical

    def arbsmarts(self, isomeric=0):
        # XXX FIX ME
        # compute canonicalization smarts
        # without hydrogens and charges
        canonical, canonical_list = traverse.drawSmarts(self)
        return canonical
            
    def arbsmiles(self, isomeric=0):
        if isomeric: draw = traverse.drawIsomeric
        else: draw = traverse.draw
        
        symorders = [atom.symorder for atom in self.atoms]

        # XXX FIX ME - kind of a hack!
        # set all the symorders to the index of the
        #  atom in self.atoms
        i = 0
        for atom in self.atoms:
            atom.symorder = i
            i+= 1
            
        result, self.arb_list = draw(self)

        # restore the symorders
        for atom, order in zip(self.atoms, symorders):
            atom.symorder = order

        return result    

    def pruneToAtoms(self, atoms):
        """Prune the molecule to the specified atoms
        bonds will be removed atomatically"""
        _atoms = self.atoms[:]
        for atom in _atoms:
            if atom not in atoms:
                self.remove_atom(atom)
