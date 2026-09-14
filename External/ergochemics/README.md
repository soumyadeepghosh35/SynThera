# ergochemics
Ergonomical cheminformatics library

# Modules

1. `ergochemics.draw` | Convenient drawing functions for molecules and reactions.
2. `ergochemics.mapping` | For mapping reaction rules to reactions, generating atom-mapped reactions, and extracting reaction centers.
3. `ergochemics.standardize` | Customizable molecule and reaction standardization techniques.
4. `ergochemics.similarity` | Featurization and similarity for molecules and reactions.

# Basic Usage

Full example notebooks available [here](https://github.com/stefanpate/ergochemics/tree/main/examples).

```python
from ergochemics.draw import draw_molecule, draw_reaction
from ergochemics.mapping import operator_map_reaction, get_reaction_center
from ergochemics.standardize import (
    standardize_mol,
    standardize_smiles,
    standardize_reaction,
    hash_molecule,
    hash_reaction
)
from ergochemics.similarity import (
    MorganFingerprinter,
    ReactionFingerprinter,
    MolFeaturizer,
)
from IPython.display import SVG, display
from rdkit import Chem
```

## Draw molecules

```python
glutamic_acid_smi = 'C(CC(=O)[O-])[C@@H](C(=O)O)N'
glutamic_acid_mol = Chem.MolFromSmiles(glutamic_acid_smi)
glutamic_acid_mol.GetAtomWithIdx(4).SetProp('atomNote', ' pKa=4.15')
display(
    SVG(
        draw_molecule(
            molecule=glutamic_acid_mol,
            size=(400, 250),
            highlight_atoms=(0, 5, 6, 9),
            draw_options={
                'addAtomIndices': True, # Example of a property set to `True`
                'addBondIndices': False, # Example of a property set to `False`
                'setHighlightColour': (0.2, 0.8, 0.2, 0.9), # Example of a callable that takes a tuple argument; pass a tuple
                'useBWAtomPalette': None, # Example of a callable that takes no arguments; pass `None`
                'addStereoAnnotation': True

            },
            legend="L-Glutamate (chiral tetrahedral highlighted in green)"
        )
    )
)
```

## Draw reactions

```python
decarboxylation = 'C(CC(=O)[O-])[C@@H](C(=O)O)N>>C(CC(=O)[O-])CN.O=C=O'
display(
    SVG(
        draw_reaction(
            rxn=decarboxylation,
            sub_img_size=(300, 200),
            draw_options={
                'addAtomIndices': True,
                'comicMode': True,
            },
        )
    )
)
```

## Map rules to reactions

```python
rule = "[O:1][C:2][C:3][N:4]>>[C:3][N:4].[O:1]=[C:2]"
res = operator_map_reaction(rxn=decarboxylation, operator=rule)
print(f"Did the rule map the reaction?  {res.did_map}")
print(f"Atom mapped reaction: {res.atom_mapped_smarts}")
print(f"Indices of the atoms that matched the rule's template: {res.template_aidxs}")
```

Output:
```
>>>Did the rule map the reaction?  True
>>>Atom mapped reaction: [NH2:2][CH:1]([CH2:3][CH2:4][C:5](=[O:6])[OH:7])[C:9](=[O:10])[OH:8]>>[NH2:2][CH2:1][CH2:3][CH2:4][C:5](=[O:6])[OH:7].[O:8]=[C:9]=[O:10]
>>>Indices of the atoms that matched the rule's template: (((9, 7, 1, 0),), ((0, 1), (0, 1)))
```
## Get reaction center

```python
print(f"Separate: {get_reaction_center(res.atom_mapped_smarts, mode="separate")}")
print(f"Combined: {get_reaction_center(res.atom_mapped_smarts, mode="combined")}")
```

```
>>>Separate: (((1, 7, 9),), ((1,), (0, 1)))
>>>Combined: ((1, 7, 9), (1, 7, 8))
```

## Standardize molecules

```python
std_mol = standardize_mol(glutamic_acid_mol)
std_mol_w_stereo = standardize_mol(glutamic_acid_mol, do_remove_stereo=False)
size = (400, 200)
display(SVG(draw_molecule(glutamic_acid_mol, size=size, legend="Original")))
display(SVG(draw_molecule(std_mol_w_stereo, size=size, legend="Standardized (stereo kept)")))
display(SVG(draw_molecule(std_mol, size=size, legend="Standardized (stereo removed)")))
```

## Customizable Morgan fingerprints

```python
def constant_atom_featurizer(atom: Chem.Atom) -> list[float | int]:
    return [1.0]

std_mol = Chem.MolFromSmiles(standardize_smiles(glutamic_acid_smi))
gcc = Chem.MolFromSmiles('CC(=C)CCC(C)C(=C)C')

display(
    SVG(
        draw_molecule(
            gcc,
            size=(300, 200),
            legend="Glutamate's carbonaceous cousin"
        )
    )
)

topology_featurizer = MolFeaturizer(
    atom_featurizer=constant_atom_featurizer
)

topo_mfper = MorganFingerprinter(
    radius=2,
    length=1024,
    mol_featurizer=topology_featurizer
)

gcc_topo_mfp = topo_mfper.fingerprint(gcc)
glutamate_topo_mfp = topo_mfper.fingerprint(std_mol)

assert np.allclose(gcc_topo_mfp, glutamate_topo_mfp)
```

## Reaction-center-sensitive reaction fingerprints

```python
rxn_mfper = ReactionFingerprinter(
    radius=2,
    length=1024,
    mol_featurizer=MolFeaturizer()
)

rxnfp_w_rc_loc = rxn_mfper.fingerprint(am_rxn, use_rc=True)
```