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
from ..exceptions import PinkyError

# XXX FIX ME
#   this would be much better moved into a class structure
#   a lot of data structures are used between the functions
NEVER = -1
MAYBE = 0
AROMATIC = 1
AROMATIC_ATOMS = {'N':1, 'C':1, 'S':1, 'O':1}
AROMATIC_PYROLE_ATOMS = {'C':1, 'N':1, 'O':1, 'S':1}
AROMATIC_5_RING = [(1,), (2,4), (1,4), (2,4), (1,)]

# Here's the deal.  When we get molecules out of the smiles
# parser, some properties are not known.  Mainly, we would
# like to find whether
#  1) is the bond type aromatic or single for non specified
#     bonds.  Non specified bonds are like C1CCCCC1 in the
#     smiles strings.  These are determined through the bond.fixed
#     property.  If a bond is fixed then it really is a single, double
#     triple or sacrifice fly bond :)
#  2) is an atom aromatic or not?
#
# At this point what is known is this:
#  lower case atoms are specified to be aromatic.
#  These have atom.aromatic = 1
#
#  Bonds are single, double or triple.
#   a double or triple bound can never be changed.
#   a specified single bond can never be changed
#   

# Table that indicates which atoms can be considered
# "pyrole like" and behave like the pyrole nitrogen
# in c1[nH]ccc1
#symbol charge hcount bonds
PyroleTable = {
    ('C', -1, 1, 2):1,
    ('C', -1, 0, 3):1,
    ('N', 0, 1, 2):1,
    ('N', 0, 0, 3):1,
    ('O', 0, 0, 2):1,
    ('S', 0, 0, 2):1,
    }

def getPyroleLikeAtoms(cycle):
    """cycle->return a dictionary of pyrole nitrogen-like atoms in
    a cycle or a molecule  The dictionary is keyed on the atom.handle"""    
    result = {}
    # the outgoing bonds might need to be single or aromatic
    for atom in cycle.atoms:
        lookup = (atom.symbol, atom.charge, atom.hcount, len(atom.bonds))
        if PyroleTable.get(lookup, 0):
            result[atom.handle] = atom

    return result
        
def canBeAromatic(cycle, pyroleLike):
    """(cycle)-> returns AROMATIC if a ring is conjugatable and
                                  passes the simple tests for aromaticity
                 returns MAYBE if the ring in its present form
                                can be aromatic but is not currently
                         NEVER if the ring can never be aromatic"""
    cycleLength = len(cycle)
    # *******************************************************
    #  check for kekular five membered rings
    if cycleLength == 5:
        # check atom types
        for atom in cycle.atoms:
            if not atom.symbol in AROMATIC_PYROLE_ATOMS:
                return NEVER

        # do we have exactly one pyrole nitrogen like atom?
        pyroleCount = 0
        for atom in cycle.atoms:
            if atom.handle in pyroleLike:
                pyrole = atom
                pyroleCount += 1

        
        if pyroleCount < 1 or pyroleCount > 2:
            return NEVER

        # rotate the ring so that we start on the pyrole like atom
        cycle.rotate(pyrole)
        bonds = cycle.bonds[:]
        # check the bonds for a kekular structure
        for index, bond in zip(range(len(bonds)), bonds):
            if bond.bondtype not in AROMATIC_5_RING[index]:
                return MAYBE
            
        return AROMATIC
    # *****************************************************
    #  check for kekular six membered rings
    #  kekular rings must have atoms in the AROMATIC_ATOMS
    #  groups and must belong in 6 membered rings.
    #  bonds must be conjugated
    elif cycleLength == 6:
        # XXX FIX ME -> there is a lot of problems with this
        # code I think, what about bonds that are already fixed?
        for atom in cycle.atoms:
            if not atom.symbol in AROMATIC_ATOMS:
                return NEVER
            
        bonds = cycle.bonds[:]        
                
        last = None
        switch = {1:2, 2:1}
        while bonds:
            bond = bonds.pop()
            bondtype = bond.bondtype

            if bond.bondorder == 3:
                return NEVER
            
            if last is None:
                if bond.bondtype in [1,2]:
                    last = bond.bondtype
            else:
                if last == 1 and bond.bondtype not in [2,4]:
                    return MAYBE
                elif last == 2 and bond.bondtype not in [1, 4]:
                    return MAYBE

                last = switch[last]
                if bondtype != last:
                    bond.bondorder = last

        return AROMATIC
    
    else:
        # we can never be aromatic
        return NEVER

def convert(cycle, pyroleLike, usedPyroles):
    """cycle, pyroleLike, aromatic=0-> aromatize the cycle
    pyroleLike is a lookup of the pyrole like atoms in the
    cycle.
    return 1 if the cycle was aromatized
           2 if the cycle could not be aromatized"""

    bonds = cycle.bonds
    atoms = cycle.atoms
    initialBondStates = []
    initialAtomStates = []
    _usedPyroles = {}
    for bond in bonds:
        # store the initial states but assume the
        # bond is aromatic
        initialBondStates.append((bond, bond.symbol,
                                  bond.bondorder, bond.bondtype,
                                  bond.aromatic, bond.stereo))
        # XXX FIX ME
        # until we get proper conjugation, aromatic bond orders
        # are 1.5
        bond.reset(':', bond.bondorder, 4, bond.fixed, bond.stereo)
        
    aromatized = 1
    for atom in atoms:
        initialAtomStates.append((atom, atom.aromatic))
        atom.aromatic = 1

        nonhydrogens = atom.sumBondOrders() + atom.charge

        # look for the lowest valence where we don't
        #  have to change the charge of the atom to
        #  fill the valences
        for valence in atom.valences:
            neededHydrogens = int(valence - nonhydrogens)
            if neededHydrogens >= 0:
                break
        else:
            # we can't change the aromaticity and have correct
            # valence.
            #
            # there is one special case of a five membered
            #  ring and a pyrole nitrogen like atom we need
            #  to look for.
            if len(cycle) == 5 and atom.handle in pyroleLike:
                _usedPyroles[atom.handle] = 1
            else:
                # nope, the valences don't work out so
                # we can't aromatize
                aromatized = 0
                break

    # sanity check, this should be true because of the
    # canBeAromatic routine above
    assert len(_usedPyroles) <=1, "Too many used pyroles!"
    
    cycle.aromatic = aromatized
    if not aromatized:
        for bond, symbol, order, bondtype, aromatic, stereo in initialBondStates:
            bond.reset(symbol, order, bondtype, bond.fixed, stereo)

        for atom, aromatic in initialAtomStates:
            atom.aromatic = aromatic
    else:
        # we used some pyroles, we'll have to send these to
        # the valence checker later
        usedPyroles.update(_usedPyroles)

    return aromatized

def addHydrogens(molecule, usedPyroles=None):
    """(molecule) -> add implicit hydrogens to a molecule.
    If the atom has specified valences and the atom must be
    charged then a Valence Error is raised"""
    for atom in molecule.atoms:
        # if the atom has an explicit hcount, we can't set the
        # hcount
        if atom.has_explicit_hcount:
            atom.hcount = atom.explicit_hcount
            continue
        
        if atom.valences:            
            for valence in atom.valences:
                hcount = max(0, int(valence - atom.sumBondOrders() + atom.charge))
                if hcount >= 0:
                    break
            else:
                if usedPyroles and not atom.handle in usedPyroles:
                    #print atom.symbol, atom.valences, atom.hcount, atom.charge,\
                    #      atom.sumBondOrders()
                    #print [x.bondtype for x in atom.bonds]
                    #print molecule.cansmiles()
                    raise PinkyError("Valence error in atom %s"%molecule.atoms.index(atom))
                pass

            #hcount = int(hcount)
            atom.hcount = hcount
    return molecule

def fixBonds(molecule, usedPyroles):
    bondsToBeFixed = {}

    # collect the bonds that need to be fixed
    for bond in molecule.bonds:
        if bond.fixed == 0:
            for atom in bond.atoms:
                if not atom.rings:
                    bond.fixed = 1
                    break
            else:
                bondsToBeFixed[bond.handle] = bond

    if not bondsToBeFixed:
        return molecule
    
    # seperate the bonds in cycles from the bonds outside
    # XXX FIX ME, this is a slow way to do this, perhaps
    # the bond should "know" what cycles it's in?
    # cycleBonds holds bonds in aromatic cycles
    cycleBonds = {}
    for cycle in molecule.cycles:
        if not cycle.aromatic:
            # if the cycle is not aromatic, assume
            # that the bond orders are already okay
            for bond in cycle.bonds:
                bond.fixed = 1
                if bond.handle in bondsToBeFixed:
                    del bondsToBeFixed[bond.handle]
        else:
            for bond in cycle.bonds:
                # check for bonds adjacent to pyrole like
                # atoms
                if bond.handle in bondsToBeFixed and len(cycle) == 5:
                    for atom in bond.atoms:
                        if atom in cycle.atoms and atom.handle in usedPyroles:
                            bond.reset(bond.symbol, 1, 4, 1, bond.stereo)
                            del bondsToBeFixed[bond.handle]
                            break
                cycleBonds[bond.handle] = bond

    for bond in bondsToBeFixed.values():
        if not bond.handle in cycleBonds:
            if bond.bondtype == 4:
                raise PinkyError("Aromatic Bond outside ring %s"% bond)
            else:
                # fix the bond
                # precondition, bondorder and bondtype are the
                # same!
                assert bond.bondorder == bond.bondtype
                bond.reset(bond.symbol, bond.bondorder,
                           bond.bondtype, 1, bond.stereo)
##  This below seems wrong, it doesn't work!
##        else:
##            for atom in bond.atoms:
##                # XXX FIX ME, should I check for all valences?
##                assert atom.valences
##                if atom.sumBondOrders() + atom.explicit_hcount >= atom.valences[0]:
##                    assert not bond.fixed
##                    bond.bondorder = 1
##                    assert bond.bondtype == 4 or bond.bondtype == bond.bondorder, "bondtype = %s, bondorder=%s"%(
##                        bond.bondtype, bond.bondorder)

##                    print "setting valence", bond
##                    bond.fixed = 1
##                    del bondsToBeFixed[bond.handle]
##                    break
                    
    while bondsToBeFixed:
        changed = 1
        while changed:
            changed = 0

            for bond in list(bondsToBeFixed.values()):
                for atom in bond.atoms:
                    for connectedBond in atom.bonds:
                        if connectedBond is not bond and \
                               connectedBond.handle in cycleBonds and \
                               connectedBond.fixed:
                            if connectedBond.bondorder == 2:
                                bondorder = 1
                            else:
                                bondorder = 2

                            bond.bondorder = bondorder
                            bond.fixed = 1
                            changed = 1
                            del bondsToBeFixed[bond.handle]
                            break
                    # if we fixed the bond, break out of the atom
                    # loop
                    if bond.fixed: break

        if not changed and bondsToBeFixed:
            # arbitrarily fix a bond to be a single bond
            bond = list(bondsToBeFixed.values())[0]
            bond.fixed = 1
            del bondsToBeFixed[bond.handle]
                        
    for bond in molecule.bonds:
        assert bond.fixed == 1
    return molecule
        
def aromatize(molecule, usedPyroles=None):    
    """(molecule, usedPyroles=None)->aromatize a molecular graph
    usedPyroles is a dictionary that holds the pyrole like
    atoms that are used in the conversion process.
    The following valence checker may need this information"""
    
    pyroleLike = getPyroleLikeAtoms(molecule)

    if usedPyroles is None:
        usedPyroles = {}
        
    cyclesToCheck = []
    # determine which cycles came in marked as aromatic
    # and which need to be checked form the kekular form
    #
    # if a cycle came in as aromatic, convert it
    # before going on.
    for cycle in molecule.cycles:
        for atom in cycle.atoms:
            if not atom.aromatic:
                cyclesToCheck.append(cycle)
                break
        else:
            if not convert(cycle, pyroleLike, usedPyroles):
                # XXX FIX ME
                # oops, an aromatic ring came in but
                # we can't convert it.  This is an error
                # daylight would conjugate the ring
                raise PinkyError("Bad initial aromaticity")

    # keep checking rings until something happens
    while 1:
        # assume nothing happened
        needToCheckAgain = 0

        _cyclesToCheck = []
        for cycle in cyclesToCheck:
            canAromatic = canBeAromatic(cycle, pyroleLike)
            if canAromatic == NEVER:
                # the ring can NEVER EVER be aromatic, so remove it for good
                pass
            elif canAromatic and convert(cycle, pyroleLike, usedPyroles):
                needToCheckAgain = 1
            else:
                _cyclesToCheck.append(cycle)

        cyclesToCheck = _cyclesToCheck
        if not needToCheckAgain:
            break

    # fix bonds that have no bondorder if necessary
    molecule = fixBonds(molecule, pyroleLike)
    # add implicit hydrogens
    return addHydrogens(molecule, usedPyroles)
