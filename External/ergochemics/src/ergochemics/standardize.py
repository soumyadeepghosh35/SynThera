from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize
from copy import deepcopy
import hashlib

def _handle_kwargs(**kwargs) -> dict:
    default_kwargs = {
        'do_canon_taut':False,
        'do_neutralize':True,
        'do_find_parent':True,
        'do_remove_stereo':True,
        "neutralization_method": "simple",
        'max_tautomers':50,
        'quiet': True,
    }
    filtered_kwargs = {}
    for k, v in kwargs.items():
        if k not in default_kwargs:
            raise ValueError(f"Unknown keyword argument: {k}")
        filtered_kwargs[k] = v

    default_kwargs.update(filtered_kwargs)
    return default_kwargs

def standardize_mol(mol: Chem.Mol, **kwargs) -> Chem.Mol:
    '''
    Standardize a molecule using RDKit's standardization tools.

    Args
    ----
    mol:rdkit.Chem.rdchem.Mol
        Molecule to standardize.
    kwargs:dict
        Keyword arguments to pass to the standardization functions.
        - do_canon_taut:bool = False
            Whether to return canonical tautomer
        - do_neutralize:bool = True
            Whether to neutralize charges
        - do_find_parent:bool = True
            Whether to find the parent molecule
        - do_remove_stereo:bool = True
            Whether to remove stereochemistry
        - neutralization_method:str = "full"
            Method to use for neutralization, "full" or "simple".
            "full" neutralization preserves atom map numbers, "simple" does not.
            "full" neutralization may lead to Kekulization issues downstream.
        - max_tautomers:int = 50
            Maximum number of tautomers to generate
        - quiet:bool = False
            Whether to suppress warnings
    
    Returns
    -------
    mol:rdkit.Chem.rdchem.Mol
        Standardized molecule.
    '''
    kwargs = _handle_kwargs(**kwargs)
    mol = deepcopy(mol) # Defensive copy

    if kwargs['quiet']:
        _ = rdBase.BlockLogs()

    if kwargs['do_remove_stereo']:
        Chem.rdmolops.RemoveStereochemistry(mol)

    # removeHs, disconnect metal atoms, normalize the molecule, reionize the molecule
    # Also checks valency, that mol is kekulizable
    mol = rdMolStandardize.Cleanup(mol)

    # if many fragments, get the "parent" (the actual mol we are interested in)
    if kwargs['do_find_parent']:
        mol = rdMolStandardize.FragmentParent(mol)

    if kwargs['do_neutralize'] and kwargs['neutralization_method'] == "full":
        mol = neutralize_charges(mol) # Remove all charges
    elif kwargs['do_neutralize'] and kwargs['neutralization_method'] == "simple":
        mol = simple_neutralize_charges(mol)

    # Enumerate tautomers and choose canonical one
    if kwargs['do_canon_taut']:
        te = rdMolStandardize.TautomerEnumerator()
        te.SetMaxTautomers(kwargs['max_tautomers'])
        mol = te.Canonicalize(mol)
    
    return mol

def standardize_smiles(smiles:str, **kwargs) -> str:
    '''
    Standardize a molecule using RDKit's standardization tools.
    Args
    ----
    smiles:str
        SMILES string to standardize.
    kwargs:dict
        Keyword arguments to pass to the standardization functions.
        - do_canon_taut:bool
            Whether to return canonical tautomer
        - do_neutralize:bool
            Whether to neutralize charges
        - do_find_parent:bool
            Whether to find the parent molecule
        - do_remove_stereo:bool
            Whether to remove stereochemistry
        - max_tautomers:int
            Maximum number of tautomers to generate
        - quiet:bool
            Whether to suppress warnings
    
    Returns
    -------
    smiles:str
        Standardized SMILES string.
    '''
    kwargs = _handle_kwargs(**kwargs)
    mol = Chem.MolFromSmiles(smiles)
    mol = standardize_mol(
        mol,
        **kwargs
    )
    
    return Chem.MolToSmiles(mol)

def standardize_reaction(rxn: str, remove_protons: bool = False, **kwargs) -> str:
    kwargs = _handle_kwargs(**kwargs)
    '''
    Standardize a reaction using RDKit's standardization tools.

    Args
    ----
    rxn:str
        SMARTS-encoded reaction, 'reactant.reactant>>product.product'
    remove_protons:bool = False
        Whether to remove protons ([H+]) from the reaction.
    kwargs:dict
        Keyword arguments to pass to the standardization functions.
        - do_canon_taut:bool
            Whether to return canonical tautomer
        - do_neutralize:bool
            Whether to neutralize charges
        - do_find_parent:bool
            Whether to find the parent molecule
        - do_remove_stereo:bool
            Whether to remove stereochemistry
        - max_tautomers:int
            Maximum number of tautomers to generate
        - quiet:bool
            Whether to suppress warnings
    
    Returns
    -------
    rxn:str
        Standardized reaction.
    '''
    rcts, pdts = [side.split('.') for side in rxn.split('>>')]

    if remove_protons:
        rcts = [elt for elt in rcts if elt != '[H+]']
        pdts = [elt for elt in pdts if elt != '[H+]']

    rcts = [standardize_smiles(r, **kwargs) for r in rcts]
    pdts = [standardize_smiles(p, **kwargs) for p in pdts]
    
    return f"{'.'.join(rcts)}>>{'.'.join(pdts)}"

def neutralize_charges(mol: Chem.Mol) -> Chem.Mol:
    """
    Neutralize atom charges in an rdkit mol, as many as possible

    Args
    ----
    mol : rdkit.Chem.rdchem.Mol
        Molecule to neutralize.

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol
        Neutralized molecule.
    """
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    for at_idx in at_matches_list:
        atom = mol.GetAtomWithIdx(at_idx)
        chg = atom.GetFormalCharge()
        hcount = atom.GetTotalNumHs()
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(hcount - chg)
        atom.UpdatePropertyCache()
    
    return mol

def simple_neutralize_charges(mol: Chem.Mol) -> Chem.Mol:
    """
    Neutralize common charge patterns in organic molecules.

    Args
    ----
    mol : rdkit.Chem.rdchem.Mol
        Molecule to neutralize.

    Returns
    -------
    mol : rdkit.Chem.rdchem.Mol
        Neutralized molecule.
    """
    patts = (
        ("[n+;H]", "n"), # Imidazoles
        ("[N+;!H0]", "N"), # Amines
        ("[$([O-]);!$([O-][#7])]", "O"), # Carboxylic acids and alcohols
        ("[S-;X1]", "S"), # Thiols
        ("[$([N-;X2]S(=O)=O)]", "N"), # Sulfonamides
        ("[$([N-;X2][C,N]=C)]", "N"), # Enamines
        ("[n-]", "[nH]"), # Tetrazoles
        ("[$([S-]=O)]", "S"), # Sulfoxides
        ("[$([N-]C=O)]", "N"), # Amides
    )

    reactions = [
        (AllChem.MolFromSmarts(x), AllChem.MolFromSmiles(y, False)) for x,y in patts
    ]

    for (reactant, product) in reactions:
        while mol.HasSubstructMatch(reactant):
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    return mol

def fast_tautomerize(smiles: str) -> list[str]:
    '''
    Applies common tautomerization patterns and returns
    any identified tautomers. Input smiles is returned
    at the first index.

    Args
    ----
    smiles:str
        SMILES string to tautomerize.
    
    Returns
    -------
    tautomers:Iterable[str]
        List of tautomers.
    '''
    transformations = [
        "[#7H1X3&a:1]:[#6&a:2]:[#7H0X2&a:3]>>[#7H0X2:1]:[#6:2]:[#7H1X3:3]",
        # TODO: add other common patterns
    ]

    tautomer_mols = []
    for trans in transformations:
        rxn = AllChem.ReactionFromSmarts(trans)
        try:
            outputs = rxn.RunReactants((Chem.MolFromSmiles(smiles),))
        except:
            print(f"Warning: rdkit sanitization failed for: {smiles}")
            outputs = rxn.RunReactants((Chem.MolFromSmiles(smiles, sanitize=False),))

        tautomer_mols.extend([o[0] for o in outputs])

    tautomer_smiles = [Chem.MolToSmiles(m) for m in tautomer_mols]
    
    return [smiles] + list(set(tautomer_smiles))

def hash_molecule(cpd: str) -> str:
    """
    Generate a hash for a molecule based on its SMILES representation.

    Args
    ----
    smiles:str
        SMILES string of the molecule.

    Returns
    -------
    hash:str
        SHA-1 hash of the SMILES string.
    """
    return hashlib.sha1(cpd.encode('utf-8')).hexdigest()

def hash_reaction(rxn: str) -> str:
    """
    Generate a hash for a reaction based on its SMILES representation.

    Args
    ----
    rxn:str
        Reaction string in the format 'reactant.reactant>>product.product'.

    Returns
    -------
    hash:str
        SHA-1 hash of the reaction string.
    """
    lhs, rhs = [side.split('.') for side in rxn.split('>>')]
    lhs = '.'.join(sorted(lhs))
    rhs = '.'.join(sorted(rhs))
    rxn = f"{lhs}>>{rhs}"
    return hashlib.sha1(rxn.encode('utf-8')).hexdigest()