from rdkit import Chem
from typing import Callable, Iterable
import numpy as np
from rdkit.Chem import rdFMCS, rdFingerprintGenerator
from copy import deepcopy
import re
from itertools import chain
from functools import partial
from ergochemics.mapping import get_reaction_center

def rcmcs_similarity(
        rxn1: str,
        rxn2: str,
        norm: str ='max',
        average: str = 'weighted',
        enforce_ring_membership: bool = False
    ) -> float:
    '''
    Calculates the reaction center max common substructure score
    between two reactions. Accepts reactions in form 'A.B>>>C.D' 
    or iterable of reactant / product SMILES. Order of reactants / products
    in all reactions, reaction centers, and reaction center patterns must match.
    
    Args
    ----
    rxn1:str | Iterable[str]
        Atom-mapped reaction 1 in form 'A.B>>>C.D'
    rxn2:str | Iterable[str]
        Atom-mapped reaction 2 in form 'A.B>>>C.D'
    norm:str
        Normalization method for rcmcs score. If 'max', score is divided by
        # of atoms in the larger reactant. If 'min', score is divided by
        # of atoms in the smaller reactant.
    average:str
        How to average the rcmcs scores of the reactants. If 'weighted',
        the average is weighted by the number of atoms in each reactant
        or product. If 'simple', the average is the arithmetic mean.

    Returns
    -------
    rcmcs_score:float
        Reaction center max common subgraph score [0, 1]
    '''
    if norm not in ['max', 'min']:
        raise ValueError("Normalization must be 'max' or 'min'")
    
    if average not in ['weighted', 'simple']:
        raise ValueError("Average must be 'weighted' or 'simple'")
    
    lrc_idx1, rrc_idx1 = get_reaction_center(rxn1, mode='combined')
    lrc_idx2, rrc_idx2 = get_reaction_center(rxn2, mode='combined')
    lhs1, rhs1 = [Chem.MolFromSmiles(side) for side in rxn1.split(">>")]
    lhs2, rhs2 = [Chem.MolFromSmiles(side) for side in rxn2.split(">>")]

    for mol in [lhs1, rhs1, lhs2, rhs2]:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)  # Clear atom map numbers to allow string compares of canonical SMILES

    lrc1, rrc1 = [Chem.MolFragmentToSmiles(mol, atomsToUse=rc_idx) for mol, rc_idx in zip([lhs1, rhs1], [lrc_idx1, rrc_idx1])]
    lrc2, rrc2 = [Chem.MolFragmentToSmiles(mol, atomsToUse=rc_idx) for mol, rc_idx in zip([lhs2, rhs2], [lrc_idx2, rrc_idx2])]

    if (lrc1 != lrc2 or rrc1 != rrc2) and (lrc1 != rrc2 or rrc1 != lrc2):
        return 0.0  # No common reaction center structure
    
    if lrc1 == rrc2:
        tmp = lhs2
        lhs2 = rhs2
        rhs2 = tmp
        tmp_idx = lrc_idx2
        lrc_idx2 = rrc_idx2
        rrc_idx2 = tmp_idx
        tmp_smiles = lrc2
        lrc2 = rrc2
        rrc2 = tmp_smiles

    lrc_idx1 = _canonically_order_atom_idxs(lhs1, tuple(lrc_idx1), lrc1)
    rrc_idx1 = _canonically_order_atom_idxs(rhs1, tuple(rrc_idx1), rrc1)
    lrc_idx2 = _canonically_order_atom_idxs(lhs2, tuple(lrc_idx2), lrc2)
    rrc_idx2 = _canonically_order_atom_idxs(rhs2, tuple(rrc_idx2), rrc2)
    
    n_atoms = 0
    rcmcs = 0
    n_mols = 0
    for mol1, mol2, rc_idx1, rc_idx2, patt in zip([lhs1, rhs1], [lhs2, rhs2], [lrc_idx1, rrc_idx1], [lrc_idx2, rrc_idx2], [lrc1, rrc1]):
        rcmcs_i = molecule_rcmcs_score([mol1, mol2], [rc_idx1, rc_idx2], patt, norm=norm, enforce_ring_membership=enforce_ring_membership)

        n1 = mol1.GetNumAtoms()
        n2 = mol2.GetNumAtoms()

        if norm == 'max':
            atoms_i = max(n1, n2)
        elif norm == 'min':
            atoms_i = min(n1, n2)

        if average == 'weighted':
            rcmcs += rcmcs_i * atoms_i
        elif average == 'simple':
            rcmcs += rcmcs_i
        
        n_mols += 1
        n_atoms += atoms_i

    if average == 'weighted':
        return rcmcs / n_atoms
    elif average == 'simple':
        return rcmcs / n_mols

def molecule_rcmcs_score(mols: list[Chem.Mol], rcs: list[tuple[int]], patt:str, norm='max', enforce_ring_membership: bool = False) -> float:
    '''
    Args
    ----
    mols: list[Chem.Mol]
        Molecules
    rcs: list[tuple[int]]
        Reaction center atom indices
    patt:str
        Reaction center substructure pattern in SMARTS
    enforce_ring_membership:bool
        Whether to enforce that ring atoms can only match ring atoms

    Returns
    -------
    rcmcs:float
        Reaction center max common substructure score [0, 1]
    '''

    res = molecule_rcmcs(mols, rcs, patt, enforce_ring_membership)

    if res is None:
        return 0.0
    elif res['canceled']:
        return 0.0
    elif norm == 'min':
        return res['num_atoms'] / min(m.GetNumAtoms() for m in mols)
    elif norm == 'max':
        return res['num_atoms'] / max(m.GetNumAtoms() for m in mols)

def molecule_rcmcs(mols: list[Chem.Mol], rcs: list[tuple[int]], patt:str, enforce_ring_membership: bool = False) -> dict | None:
    '''
    Args
    ----
    mols: list[Chem.Mol]
        Molecules
    rcs: list[tuple[int]]
        Reaction center atom indices
    patt:str
        Reaction center substructure pattern in SMARTS
    enforce_ring_membership:bool
        Whether to enforce that ring atoms can only match ring atoms

    Returns
    -------
    res | None
        FindMCS output or None if failed or failed pre-check
    '''
    if len(mols) != len(rcs):
        raise ValueError("Number of molecules and reaction centers do not match")
    
    mols = [deepcopy(m) for m in mols]

    rc_scalar = 100

    def _replace(match):
        atomic_number = int(match.group(1))
        return f"[{atomic_number * rc_scalar}#{atomic_number}"

    def _reset(match):
        atomic_number = int(match.group(1))
        
        if atomic_number % rc_scalar == 0:
            return f"[{int(atomic_number / rc_scalar)}"
        else:
            return f"[{atomic_number}"

    patt = re.sub(r'\[#(\d+)', _replace, patt) # Mark reaction center patt w/ isotope number

    # Mark reaction center vs other atoms in substrates w/ isotope number
    for mol, rc in zip(mols, rcs):
        for atom in mol.GetAtoms():
            if atom.GetIdx() in rc:
                atom.SetIsotope(atom.GetAtomicNum() * rc_scalar) # Rxn ctr atom
            else:
                atom.SetIsotope(atom.GetAtomicNum()) # Non rxn ctr atom

    cleared, patt = _mcs_precheck(mols, rcs, patt, enforce_ring_membership) # Prevents FindMCS default behavior of non-rc-mcs

    if not cleared:
        return None

    # Get the mcs that contains the reaction center pattern
    tmp = rdFMCS.FindMCS(
        mols,
        seedSmarts=patt,
        atomCompare=rdFMCS.AtomCompare.CompareIsotopes,
        bondCompare=rdFMCS.BondCompare.CompareOrderExact,
        matchChiralTag=False,
        ringMatchesRingOnly=enforce_ring_membership,
        completeRingsOnly=False,
        matchValences=True,
        timeout=10
    )

    rcmcs_patt = Chem.MolFromSmarts(tmp.smartsString)
    rcmcs_idxs = [m.GetSubstructMatch(rcmcs_patt) for m in mols]
    smarts_string = re.sub(r'\[(\d+)', _reset, tmp.smartsString) # Remove rc scaling

    res = {
        'smarts_string': smarts_string,
        'rcmcs_idxs': rcmcs_idxs,
        'num_atoms': tmp.numAtoms,
        'num_bonds': tmp.numBonds,
        'canceled': tmp.canceled
    }

    return res

def _mcs_precheck(mols: list[Chem.Mol], rcs: list[tuple[int]], patt: str, enforce_ring_membership: bool):
    '''
    Modifies single-atom patts and pre-checks ring info
    to avoid giving FindMCS a non-common-substructure which
    results in non-reaction-center-inclusive MCSes
    '''
    if patt.count('#') == 1:
        patt = _handle_single_atom_patt(mols, rcs, patt)
    
    if enforce_ring_membership:
        cleared = _check_ring_membership(mols, rcs)
    else:
        cleared = True

    return cleared, patt

def _handle_single_atom_patt(mols: list[Chem.Mol], rcs: list[tuple[int]], patt: str):
    '''
    Pre-pends wildcard atom and bond to single-atom
    patt if mols share a neighbor w/ common isotope,
    ring membership, & bond type between
    '''
    couples = [set() for _ in range(len(mols))]
    for i, (mol, rc) in enumerate(zip(mols, rcs)):
        rc_idx = rc[0]
        for neighbor in mol.GetAtomWithIdx(rc_idx).GetNeighbors():
            nidx = neighbor.GetIdx()
            nisotope = neighbor.GetIsotope()
            in_ring = neighbor.IsInRing()
            bond_type = mol.GetBondBetweenAtoms(rc_idx, nidx).GetBondType()
            couples[i].add((nisotope, in_ring, bond_type))

    if len(set.intersection(*couples)) > 0:
        patt = '*~' + patt
    
    return patt

def _check_ring_membership(mols: list[Chem.Mol], rcs: list[tuple[int]]):
    '''
    Returns false if any "aligned" atom has distinct ring membership
    '''
    alignments = zip(*rcs)
    for elt in alignments:
        ring_membership = [mols[i].GetAtomWithIdx(idx).IsInRing() for i, idx in enumerate(elt)]

        if len(set(ring_membership)) != 1:
            return False
        
    return True

'''
Atom featurizers
'''

def dai(atom: Chem.Atom) -> list[int | float]:
    '''
    Returns Daylight atomic invariants for atom
    '''
    dai = [
        atom.GetDegree(), # Heavy atoms only
        atom.GetTotalValence() - atom.GetTotalNumHs(),
        atom.GetAtomicNum(),
        atom.GetMass(),
        atom.GetFormalCharge(),
        int(atom.IsInRing()),
        int(atom.GetIsAromatic()),
    ]

    return dai

def dai_amphoteros(atom: Chem.Atom) -> list[int | float]:
    '''
    Returns Daylight atomic invariants + amphoteros oxidation state for atom
    '''
    atomic_invariants = [
        atom.GetDegree(),
        atom.GetTotalValence() - atom.GetTotalNumHs(),
        atom.GetAtomicNum(),
        atom.GetFormalCharge(),
        int(atom.IsInRing()),
        int(atom.GetIsAromatic()),
        amphoteros_ox_state(atom)
    ]

    return atomic_invariants

def rule_default(atom: Chem.Atom) -> list[int | float]:
    '''
    Returns atom features used in reaction rule extraction
    here: https://github.com/stefanpate/coarse-grain-rxns
    '''
    atomic_invariants = [
        atom.GetDegree(),
        atom.GetTotalValence(),
        atom.GetTotalNumHs(),
        atom.GetAtomicNum(),
        atom.GetFormalCharge(),
        int(atom.IsInRing()),
        int(atom.GetIsAromatic()),
        z(atom)
    ]

    return atomic_invariants

def z(atom: Chem.Atom) -> float:
    '''
    Returns number of heteroatom neighbors for carbon,
    -1 if atom is not carbon
    '''
    if atom.GetAtomicNum() != 6:
        return -1.0
    else:
        return sum(
            float(bond.GetOtherAtom(atom).GetAtomicNum() != 6)
            for bond in atom.GetBonds()
        )
  
def amphoteros_ox_state(atom: Chem.Atom) -> float:
    '''
    Returns
    -------
    : float
        -1 if atom is not carbon
        + (# pi bonds + # heteroatom neighbors) otherwise

    Notes
    -----
    https://amphoteros.com/2013/10/22/counting-oxidation-states/
    '''
    if atom.GetAtomicNum() != 6:
        return -1.0
    else:
        return sum(
            (bond.GetBondTypeAsDouble() - 1.0) + float(bond.GetOtherAtom(atom).GetAtomicNum() != 6)
            for bond in atom.GetBonds()
        )

class MolFeaturizer:
    '''
    Featurizes molecules using a custom atom featurization scheme
    and optionally reaction center information in the case that molecules
    are reactants in a reaction.
    '''
    def __init__(self, atom_featurizer: Callable[[Chem.Atom], list[int | float]] = rule_default):
        '''
        Args
        ----
        atom_featurizer:Callable[[Chem.Atom], list[int | float]]
            Function that featurizes an atom into a list of features
        '''
        self.atom_featurizer = atom_featurizer

    def featurize(self, mol: Chem.Mol, rc: Iterable[int] = []) -> np.ndarray:
        '''
        Args
        ----
        mol: Chem.Mol
            RDKit molecule object
        rc: Iterable[int] (optional)
            List of atom indices corresponding to reaction center

        Returns
        -------
        fts: np.ndarray
            Node feature matrix of shape (num_atoms, num_features)
        
        Notes
        -----
        1. Min distance to reaction center always the last feature.
        2. If atom is not connected to an rc atom, distance is set to 10_000
        '''
        fts = []
        for atom in mol.GetAtoms():
            aidx = atom.GetIdx()
            local_fts = self.atom_featurizer(atom)
            min_dist = 10_000  # default for disconnected atoms
            for rcidx in rc:
                if aidx == rcidx:
                    min_dist = 0
                    break
                else:
                    dist = len(Chem.GetShortestPath(mol, aidx, rcidx)) - 1
                    if dist < 0:
                        continue
                    elif dist < min_dist:
                        min_dist = dist

            fts.append(local_fts + [min_dist])

        fts = np.array(fts)
        return fts
    
class MorganFingerprinter:
    '''
    Generates Morgan fingerprints using custom atom featurization scheme
    and optionally reaction center information in the case that molecules
    are reactants in a reaction.
    '''
    def __init__(self, radius: int, length: int, mol_featurizer: MolFeaturizer, allocate_ao: bool = False, **kwargs):
        '''
        Args
        ----
        radius:int
            Radius of Morgan fingerprint (max # hops to get fragments)
        length:int
            Length of Morgan fingerprint (# elements in fingerprint vector)
        mol_featurizer:MolFeaturizer
            Molecule featurizer to generate atom features
        allocate_ao:bool = False
            Whether to allocate additional output (bit info map, atom counts, atom to bits)
        '''
        self._generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=length, **kwargs)
        
        if allocate_ao:
            self._additional_ouput = rdFingerprintGenerator.AdditionalOutput()
            self._additional_ouput.AllocateBitInfoMap()
            self._additional_ouput.AllocateAtomCounts()
            self._additional_ouput.AllocateAtomToBits()
        else:
            self._additional_ouput = None

        self._fingerprint = {
            'bit': partial(self._generator.GetFingerprintAsNumPy, additionalOutput=self._additional_ouput) if allocate_ao else self._generator.GetFingerprintAsNumPy,
            'count': partial(self._generator.GetCountFingerprintAsNumPy, additionalOutput=self._additional_ouput) if allocate_ao else self._generator.GetCountFingerprintAsNumPy,
        }
        self.mol_featurizer = mol_featurizer

    def fingerprint(self, mol: Chem.Mol, output_type: str = 'bit', reaction_center: Iterable[int] = [], rc_dist_ub: int = None) -> np.ndarray:
        '''
        Fingerprints molecule.

        Args
        ----
        mol:Chem.Mol
            RDKit molecule object
        output_type:str = 'bit'
            Type of fingerprint to return. 'bit' for binary fingerprint,
            'count' for count-based fingerprint.
        reaction_center:Iterable[int] = []
            List of atom indices corresponding to reaction center
        rc_dist_ub:int = None
            Upper bound (inclusive) on bondwise distance from reaction center
            for atoms to include in fingerprint. If None, all atoms are included.

        Returns
        -------
        fp:np.ndarray
            Morgan fingerprint of molecule
        '''
        if rc_dist_ub is not None and len(reaction_center) == 0:
            raise ValueError("If providing upper bound on distance from reaction center, must also provide reaction center")

        feats = self.mol_featurizer.featurize(mol, reaction_center)

        if rc_dist_ub is not None:
            root_atoms = [
                i for i, ft in enumerate(feats)
                if ft[-1] <= rc_dist_ub
            ]       
            feats = [self.hash_features(tuple(ft.tolist())) for ft in feats]
            return self._fingerprint[output_type](mol, customAtomInvariants=feats, fromAtoms=root_atoms)
        else:
            feats = [self.hash_features(tuple(ft.tolist())) for ft in feats]
            return self._fingerprint[output_type](mol, customAtomInvariants=feats)

    def hash_features(self, atom_feats: tuple):
        return hash(atom_feats) & 0xFFFFFFFF # bitwise AND w/ 0xFFFFFFFF to get 32-bit hash expected by rdkit
    
    @property
    def bit_info_map(self) -> dict:
        '''
        Returns bit info map if allocated, else empty dict.

        Returns
        -------
        bit_info_map:dict
            Maps fingerprint bits to atom indices and fragment radius.
        '''
        if self._additional_ouput:
            return self._additional_ouput.GetBitInfoMap()
        else:
            return {}
    
    @property
    def atom_counts(self) -> tuple:
        '''
        Returns the number of bits each atom contributes to the fingerprint.

        Returns
        -------
        atom_counts:tuple
            Number of bits each atom contributes to the fingerprint.
        '''
        if self._additional_ouput:
            return self._additional_ouput.GetAtomCounts()
        else:
            return tuple()
        
    @property
    def atom_to_bits(self) -> tuple:
        '''
        Returns which bits each atom contributes to the fingerprint.

        Returns
        -------
        atom_to_bits:tuple
            Which bits each atom contributes to the fingerprint.
        '''
        if self._additional_ouput:
            return self._additional_ouput.GetAtomToBits()
        else:
            return tuple()

class ReactionFingerprinter(MorganFingerprinter):
    '''
    Convenience class to fingerprint reactions.
    '''
    def __init__(self, radius: int, length: int, mol_featurizer: MolFeaturizer, allocate_ao: bool = False, **kwargs):
        '''
        Args
        ----
        radius:int
            Radius of Morgan fingerprint (max # hops to get fragments)
        length:int
            Length of Morgan fingerprint (# elements in fingerprint vector)
        mol_featurizer:MolFeaturizer
            Molecule featurizer to generate atom features
        allocate_ao:bool = False
            Whether to allocate additional output (bit info map, atom counts, atom to bits)
        '''
        super().__init__(radius, length, mol_featurizer, allocate_ao, **kwargs)

    def fingerprint(self, rxn: str, output_type: str = 'bit', use_rc: bool = False, rc_dist_ub: int = None) -> np.ndarray:
        '''
        Fingerprints reaction

        Args
        ----
        rxn:str
            Reaction in form 'A.B>>>C.D'. Must be all-atom mapped
            if using reaction center information.
        output_type:str = 'bit'
            Type of fingerprint to return. 'bit' for binary fingerprint,
            'count' for count-based fingerprint.
        use_rc:bool = False
            Whether to use min distance to reaction center as an atom feature.
        rc_dist_ub:int = None
            Upper bound on distance from reaction center for atoms to include
            in fingerprint. If None, all atoms are included.

        Returns
        -------
        fp:np.ndarray
            Morgan fingerprint of reaction (as concatenation of left- and right-hand side fingerprints)
        '''
        if use_rc:
            lrc, rrc = get_reaction_center(rxn, mode='combined')
        else:
            lrc, rrc = [], []

        lhs, rhs = [Chem.MolFromSmiles(side) for side in rxn.split(">>")]
        lhs_fp = super().fingerprint(lhs, reaction_center=lrc, output_type=output_type, rc_dist_ub=rc_dist_ub)
        rhs_fp = super().fingerprint(rhs, reaction_center=rrc, output_type=output_type, rc_dist_ub=rc_dist_ub)
        return np.concatenate([lhs_fp, rhs_fp])
    
def _canonically_order_atom_idxs(mol: Chem.Mol, aidxs: tuple[int], canonical_smiles_pattern: str) -> tuple[int]:
    '''
    Returns atom indices ordered according to
    their order in the canonical SMILES of the
    substructure defined by canonical_smiles_pattern.
    '''
    sub_mol = Chem.MolFromSmarts(canonical_smiles_pattern)
    matches = mol.GetSubstructMatches(sub_mol)
    srt_aidxs = sorted(aidxs)

    for match in matches:
        srt_match = sorted(match)
        if srt_aidxs == srt_match:
            return tuple(match)

    raise ValueError("Atom indices do not match any substructure defined by the pattern")