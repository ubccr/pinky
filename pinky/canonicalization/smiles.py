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



Traversals - a traversal (or path) through a graph

Traversals are stored as a list of tokens representing each object in a graph.
Each token object has __str__ defined so once a list has been formed
it can be turned into a string object by

 "".join(traversal_list)

Tokens are defined by the representation desired.  See tokens.py for
a list of SmilesTokens.

One thing to remember is that closures are considered atom properties for
the purposes of traversals.  That is the AtomToken is reponsible for storing
information about closures.  This information is built up during the
depth first search of the graph.
"""

import tokens

class IDGenerator:
    def __init__(self):
        self.index = 0

    def next(self):
        self.index += 1
        return self.index
    
class SmilesTraversal:
    # override the following to create different
    # tokens for traversals
    AtomToken = tokens.Atom
    BondToken = tokens.Bond
    BranchToken = tokens.Branch
    BranchEndToken = tokens.BranchEnd

    def __init__(self, parent=None):
        if parent is None:
            atomsDone = {}
            idGenerator = IDGenerator()
            closureIdGenerator = IDGenerator()
            closures = {}
        else:
            atomsDone = parent.atomsDone
            idGenerator = parent.idGenerator
            closureIdGenerator = parent.closureIdGenerator
            closures = parent.closures

        self.atoms = []
        self.bonds = []
        self.data = []
        self.atomsDone = atomsDone
        self.idGenerator = idGenerator
        self.closureIdGenerator = closureIdGenerator
        self.closures = closures

    def addAtom(self, atom):
        #print "atom", atom,
        atomToken = self.AtomToken(atom,
                                   self.closures,
                                   self.closureIdGenerator)
        self.atomsDone[atom] = atomToken
        self.data.append(atomToken)
        self.atoms.append(atom)
        #print self.atomsDone

    def addBond(self, bond):
        #print "bond", bond
        self.data.append(self.BondToken(bond))
        # XXX HACK -> give the bond a traversal order
        self.bonds.append(bond)

    def addClosure(self, atom1, atom2, bond):
        #print "closure", atom1, atom2, bond
        atomsDone = self.atomsDone
        closures1 = atomsDone[atom1].closures
        closures2 = atomsDone[atom2].closures
        assert closures1 is not closures2
        id = self.idGenerator.next()
        
        bondToken = self.BondToken(bond)
        # only pass the bond info to the first atom
        # the other gets a None (or don't print bond info)
        closures1.append((id, bondToken))
        closures2.append((id, None))
        self.bonds.append(bond)

    def addBranch(self):
        self.data.append(self.BranchToken())

    def addBranchEnd(self):
        self.data.append(self.BranchEndToken())

    def append(self, traverse):
        """(traverse)->append the traverse to the current traverse"""
        self.data.extend(traverse.data)
        self.atoms.extend(traverse.atoms)
        self.bonds.extend(traverse.bonds)

    def __str__(self):
        # XXX HACK -> find the bond traversal order
        # and make a note of it.  This is pretty ugly
        index = 1
        for x in self.data:
            if x.__class__ == tokens.Bond:
                x.bond._traverseOrder = index
                index += 1
        return "".join(map(str,self.data))


class SmartsTraversal(SmilesTraversal):
    AtomToken = tokens.SmartsAtom
    BondToken = tokens.SmartsBond

class IsomericSmilesTraversal(SmilesTraversal):
    AtomToken = tokens.IsomericAtom
