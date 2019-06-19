# pinky - molecular fingerprint library

pinky is library for generating molecular fingerprints from SMILES strings.

pinky includes a fork of the SMILES parser and mol builder from frowns
http://frowns.sourceforge.net/ originally written by Brian Kelley with
contributions from Andrew Dalke.

This should be considered alpha software. The goal is to add more molecular
fingerprint algorithms (currently only ECFP is supported). Also, would like to
battle test the SMILES parser to provide a robust pure python SMILES parsing
library.

## Usage

```python
from pinky.smiles import smilin
from pinky.fingerprints import ecfp

mol = smilin('CCCC(=O)N')
for atom in mol.atoms:
    print atom, atom.sumBondOrders()
for bond in mol.bonds:
    print bond.symbol, bond.bondorder, bond.bondtype, bond.fixed

# Compute ECFP_4 fingerprint
fp = ecfp(mol, radius=2)
```

## License

See LICENSE file.
