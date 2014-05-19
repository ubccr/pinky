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


Recursively traverse a molecule building up a canonical
representation.

Each atom of a Molecule or Graph must have a attribute 'symorder' which
is a unique number.  This number guarantees only one traversal for
the graph.

Additionally each bond must have an attribute equiv_class which is
a unique value for each different type of bond.  This guarantees
proper canonicalization of bonds as well as atoms.

canonical_string = Traverse.draw(molecule)

Advanced usage:
canonical_string = Traverse.draw(molecule, TraversalType)

TraversalType controls how the traversal is represented.
SmilesTraversal is the default TraversalType but this can be
subclassed for different representations.  For example it is
easy to create a subclass to generate Tripos Line formats.
"""

from .smiles import SmilesTraversal, SmartsTraversal,  IsomericSmilesTraversal

def _traverse(atom, traverse, prevAtom,
              visitedAtoms, visitedBonds,
              atoms, bonds, Traversal, bondIndex=0):
    visitedAtoms[atom] = 1
    traverse.addAtom(atom)
    atoms.append(atom)

    bondsToTraverse = []
    traversals = []
    
    for bond in atom.bonds:
        oatom = bond.xatom(atom)
        if prevAtom is not None and oatom == prevAtom:
            # we are traversing back the way we came!
            # so don't...
            pass
        elif visitedAtoms.has_key(oatom):
            # a closure!
            traverse.addClosure(atom, oatom, bond)
            bonds.append(bond)
            visitedBonds[bond] = 1
        else:
            bondsToTraverse.append((oatom.symorder,
                                    bond.equiv_class,
                                    bondIndex,
                                    oatom,
                                    bond))
        bondIndex += 1

    if not bondsToTraverse:
        # dead end, return
        return

    bondsToTraverse.sort()

    for symorder, bondEclass, index, oatom, obond in bondsToTraverse:
        if visitedAtoms.has_key(oatom):
            # somehow, we've seen this atom so skip it
            continue
        nextTraverse = Traversal(traverse)
        traversals.append(nextTraverse)
        nextTraverse.addBond(obond)
        bonds.append(obond)
        visitedBonds[obond] = 1

        _traverse(oatom, nextTraverse, atom,
                  visitedAtoms, visitedBonds,
                  atoms, bonds, Traversal, bondIndex)

    for t in traversals[:-1]:
        traverse.addBranch()
        traverse.append(t)
        traverse.addBranchEnd()
        
    for t in traversals[-1:]:
        traverse.append(t)

def _get_lowest_symorder(atoms):
    best = atoms[0]
    for atom in atoms[1:]:
        if atom.symorder < best.symorder:
            best = atom
    return best

def draw(molecule, TraversalType=SmilesTraversal):
    """(molecule)->canonical representation of a molecule
    Well, it's only canonical if the atom symorders are
    canonical, otherwise it's arbitrary.

    atoms must have a symorder attribute
    bonds must have a equiv_class attribute"""
    result = []
    atoms = allAtoms = molecule.atoms

    visitedAtoms = {}
    #
    # Traverse all components of the graph to form
    # the output string
    while atoms:
        atom = _get_lowest_symorder(atoms)
        visitedAtoms[atom] = 1

        visitedBonds = {}
        nextTraverse = TraversalType()
        atomsUsed, bondsUsed = [], []
        _traverse(atom, nextTraverse, None,
                  visitedAtoms, visitedBonds,
                  atomsUsed, bondsUsed, TraversalType)
        atoms = []
        for atom in allAtoms:
            if not visitedAtoms.has_key(atom):
                atoms.append(atom)
        assert nextTraverse.atoms == atomsUsed
        assert nextTraverse.bonds == bondsUsed, "%s %s"%(
            nextTraverse.bonds, bondsUsed)
        

        result.append((str(nextTraverse),
                       atomsUsed, bondsUsed))

    result.sort()
    fragments = []
    for r in result:
        fragments.append(r[0])

    return ".".join(fragments), result

def drawSmarts(molecule):
    return draw(molecule, TraversalType=SmartsTraversal)

def drawIsomeric(molecule):
    return draw(molecule, TraversalType=IsomericSmilesTraversal)
