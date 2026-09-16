from ergochemics.mapping import get_reaction_center, operator_map_reaction

decarboxylation_rule = '[#6:1]-[#6:2]-[#8:3]>>[#6:1].[#6:2]=[#8:3]'
decarboxylation_reaction = 'OCC(N)C(=O)O>>OCC(N).O=C=O'

def test_operator_map_reaction():
    res = operator_map_reaction(decarboxylation_reaction, decarboxylation_rule)
    print(res)
    assert res.did_map
    assert res.aligned_smarts == 'NC(CO)C(=O)O>>NCCO.O=C=O'
    assert res.atom_mapped_smarts == '[NH2:2][CH:1]([CH2:3][OH:4])[C:5](=[O:7])[OH:6]>>[NH2:2][CH2:1][CH2:3][OH:4].[O:6]=[C:5]=[O:7]'
    assert res.template_aidxs == (((1, 4, 6),), ((1,), (0, 1)))

tetrahedral_chiral_inversion = '[C:1][C:2][C@H:3]([C:4])[Br:5]>>[C:1][C:2][C@@H:3]([C:4])[Br:5]'
stereo_double_bond_inversion = '[C:1]/[C:2]=[C:3]/[C:4]=[O:5]>>[C:1]/[C:2]=[C:3]\\[C:4]=[O:5]'

def test_reaction_center_tetrahedral_inversion():
    rc_tetrahedral = get_reaction_center(tetrahedral_chiral_inversion, include_stereo=True)
    assert rc_tetrahedral == (((1, 2, 3, 4), ), ((1, 2, 3, 4), ))

def test_reaction_center_stereo_double_bond_inversion():
    rc_double_bond = get_reaction_center(stereo_double_bond_inversion, include_stereo=True)
    assert rc_double_bond == (((0, 1, 2, 3), ), ((0, 1, 2, 3), ))

def test_reaction_center_non_consecutive_amns():
    rxn = "[C:1][C:3]=[O:5]>>[C:1][C:3][O:5]"
    rc = get_reaction_center(rxn)
    assert rc == (((1, 2),), ((1, 2),))

def test_mapping_unbalanced_rule():
    rule = "[#6:15]-[#6:16](-[#6:17])>>[#6:15]-[#6:16](-[#6:17])(-[#8])"
    rxn = "CCCO>>CC(O)CO"
    res = operator_map_reaction(rxn, rule)
    assert res.did_map
    assert res.atom_mapped_smarts is None
    assert res.template_aidxs[0] == ((0, 1, 2),)
    assert res.template_aidxs[1] == ((),)
    assert res.aligned_smarts == rxn


def test_mapping_explicit_hs():
    rule = "[H][C:8]([C:7])([O:9][H])[H:19]>>[C:7][C:8](=[O:9])[H:19]"
    rxn = "CCCO>>CCC=O"
    res = operator_map_reaction(rxn, rule, explicit_hs=True)
    assert res.did_map
    assert res.atom_mapped_smarts == '[H][O:3][C:2]([H])([H:4])[C:1]([H:6])([H:7])[C:5]([H:8])([H:9])[H:10]>>[H:4][C:2](=[O:3])[C:1]([H:6])([H:7])[C:5]([H:8])([H:9])[H:10]'
    assert res.template_aidxs == (((5, 1, 2, 0, 4, 6),), ((0, 1, 2),))
    assert res.aligned_smarts == '[H]OC([H])([H])C([H])([H])C([H])([H])[H]>>[H]C(=O)C([H])([H])C([H])([H])[H]'