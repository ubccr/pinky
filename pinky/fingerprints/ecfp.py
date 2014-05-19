"""
Implementation of Extended-Connectivity Fingerprints algorithm (David Rogers and Mathew Hahn 2010). 

This is based on the MorganFingerprints.cpp code from rdkit:
  https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Fingerprints/MorganFingerprints.cpp
"""
from bitarray import bitarray
import copy

def gen_hash(lst):
    return hash(tuple(lst))

def invariants(mol):
    """Generate initial atom identifiers using atomic invariants"""
    atom_ids = {}
    for a in mol.atoms:
        components = []
        components.append(a.number)
        components.append(len(a.bonds))
        components.append(a.hcount)
        components.append(a.charge)
        components.append(a.mass)
        if len(a.rings) > 0:
            components.append(1)

        atom_ids[a.index] = gen_hash(components)

    return atom_ids

def ecfp(mol, radius=2):
    """Compute the Extended-Connectivity fingerprint for a molecule.
    
    :param mol: molecule object parsed from SMILES string
    :param radius: The number of iterations to perform. Defaults to 2 which is equivilent to ECFP4.
    :rtype: dictionary representing the molecular fingprint (atom identifiers and their counts).
    """

    atom_ids = invariants(mol)

    fp = {}
    for i in atom_ids.values():
        fp[i] = fp.get(i, 0) + 1

    neighborhoods = []
    atom_neighborhoods = [ len(mol.bonds) * bitarray('0') for a in mol.atoms]
    dead_atoms = len(mol.atoms) * bitarray('0')

    for layer in xrange(1, radius+1):
        round_ids = {}
        round_atom_neighborhoods = copy.deepcopy(atom_neighborhoods)
        neighborhoods_this_round = []

        for a in mol.atoms:
            if dead_atoms[a.index]: continue

            nbsr = []
            for b in a.bonds:
                round_atom_neighborhoods[a.index][b.index] = True
                oidx = b.xatom(a).index
                round_atom_neighborhoods[a.index] |= atom_neighborhoods[oidx]
                nbsr.append((b.bondtype, atom_ids[oidx]))

            nbsr = sorted(nbsr)
            nbsr = [item for sublist in nbsr for item in sublist]
            nbsr.insert(0, atom_ids[a.index])
            nbsr.insert(0, layer)

            round_ids[a.index] = gen_hash(nbsr)
            neighborhoods_this_round.append(
                (round_atom_neighborhoods[a.index], round_ids[a.index], a.index)
            )

        for lst in neighborhoods_this_round:
            if lst[0] not in neighborhoods:
                fp[lst[1]] = fp.get(lst[1], 0) + 1
                neighborhoods.append(lst[0])
            else:
                dead_atoms[lst[2]] = True

        atom_ids = round_ids
        atom_neighborhoods = copy.deepcopy(round_atom_neighborhoods)
    return fp
