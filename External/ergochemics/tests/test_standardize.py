from ergochemics.standardize import standardize_smiles, standardize_reaction, hash_reaction, hash_molecule

smi = 'O=C(O)CC[c-]1[nH]cnc1=O'
rxn = 'O=C(O)CC[c-]1[nH]cnc1=O.O=C(O)CC[c-]1[nH]cnc1=O>>O=C(O)CC[c-]1[nH]cnc1=O.O=C(O)CC[c-]1[nH]cnc1=O'

def test_standardize_reaction():
    assert hash_reaction(standardize_reaction(rxn)) == hash_reaction(standardize_reaction(rxn))

def test_standardize_smiles():
    assert hash_molecule(standardize_smiles(smi)) == hash_molecule(standardize_smiles(smi))