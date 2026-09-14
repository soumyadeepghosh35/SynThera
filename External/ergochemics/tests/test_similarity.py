from ergochemics.similarity import MorganFingerprinter, MolFeaturizer, ReactionFingerprinter, rcmcs_similarity
import numpy as np
from rdkit import Chem

sma1 = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1.O=C(O)CCC(=O)C(=O)O"
sma2 = "NC(=O)C1=CN(C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)C=CC1.O=C(O)CCC(=O)C(=O)O"
rc1 = [3, 47, 46, 45, 5, 4, 50, 49]
rc2 = [3, 4, 5, 45, 46, 47, 49, 50]

def test_molecule_fingerprinter_without_rc():
    mfper = MorganFingerprinter(
        radius=2,
        length=2048,
        mol_featurizer=MolFeaturizer(),
    )
    fp1 = mfper.fingerprint(Chem.MolFromSmiles(sma1))
    fp2 = mfper.fingerprint(Chem.MolFromSmiles(sma2))
    assert np.allclose(fp1, fp2)  # Should be True if fingerprints are equal

def test_molecule_fingerprinter_scramble_rc_order():
    mfper = MorganFingerprinter(
        radius=2,
        length=2048,
        mol_featurizer=MolFeaturizer(),
    )
    fp1 = mfper.fingerprint(Chem.MolFromSmiles(sma1), reaction_center=rc1, rc_dist_ub=1)
    fp2 = mfper.fingerprint(Chem.MolFromSmiles(sma2), reaction_center=rc2, rc_dist_ub=1)
    assert np.allclose(fp1, fp2)  # Should be True if fingerprints are equal

sma1 = "CCO.C(=O)O"
sma2 = "C(=O)O.CCO"
rc1 = [1, 2, 3, 4]
rc2 = [0, 1, 4, 5]

def test_molecule_fingerprinter_scramble_smiles_order():
    mfper = MorganFingerprinter(
        radius=2,
        length=2048,
        mol_featurizer=MolFeaturizer(),
    )
    fp1 = mfper.fingerprint(Chem.MolFromSmiles(sma1), reaction_center=rc1)
    fp2 = mfper.fingerprint(Chem.MolFromSmiles(sma2), reaction_center=rc2)
    assert np.allclose(fp1, fp2)  # Should be True if fingerprints are equal

am_rxn = '[NH2:2][CH:1]([CH2:3][OH:4])[C:5](=[O:7])[OH:6]>>[NH2:2][CH2:1][CH2:3][OH:4].[O:6]=[C:5]=[O:7]'
rxn = 'NC(CO)C(=O)O>>NCCO.O=C=O'

def test_reaction_fingerprinter_uses_rc():
    rfper = ReactionFingerprinter(
        radius=2,
        length=2048,
        mol_featurizer=MolFeaturizer(),
    )

    # Both should return fingerprints, but they should be different
    fp_with_rc = rfper.fingerprint(am_rxn, use_rc=True)
    fp_without_rc = rfper.fingerprint(rxn, use_rc=False)
    assert not np.allclose(fp_with_rc, fp_without_rc)

am_sma = '[O:10]=[C:8]([OH:9])[c:1]1[cH:2][cH:4][cH:7][cH:5][c:3]1[OH:6]>>[OH:6][c:3]1[cH:1][cH:2][cH:4][cH:7][cH:5]1.[O:9]=[C:8]=[O:10]'
am_sma_rev = '[OH:6][c:3]1[cH:1][cH:2][cH:4][cH:7][cH:5]1.[O:9]=[C:8]=[O:10]>>[O:10]=[C:8]([OH:9])[c:1]1[cH:2][cH:4][cH:7][cH:5][c:3]1[OH:6]'

def test_rcmcs_similarity_rxn_vs_rev():
    rcmcs_score = rcmcs_similarity(am_sma, am_sma_rev)
    assert rcmcs_score == 1.0

other_sma = '[NH2:4][CH:3]([CH2:1][S:13](=[O:12])(=[O:14])[OH:15])[C:5](=[O:6])[OH:7].[O:9]=[P:8]([OH:2])([OH:10])[OH:11]>>[NH2:4][CH:3]([CH2:1][O:2][P:8](=[O:9])([OH:10])[OH:11])[C:5](=[O:6])[OH:7].[O:14]=[S:13]([OH:12])[OH:15]'

def test_rcmcs_similarity_different_rcs():
    rcmcs_score_2 = rcmcs_similarity(am_sma, other_sma)
    assert rcmcs_score_2 == 0.0

am_sma2 = '[CH3:2][CH:1]([C:8](=[O:10])[OH:9])[C:3](=[O:4])[C:5](=[O:6])[OH:7]>>[CH3:2][CH2:1][C:3](=[O:4])[C:5](=[O:6])[OH:7].[O:9]=[C:8]=[O:10]'
similar_sma = '[O:3]=[CH:2][CH2:1][C:4](=[O:6])[OH:5]>>[CH3:1][CH:2]=[O:3].[O:5]=[C:4]=[O:6]'

def test_rcmcs_similarity_similar_rxns_same_rcs():
    rcmcs_score_3 = rcmcs_similarity(am_sma2, similar_sma)
    assert rcmcs_score_3 > 0.0