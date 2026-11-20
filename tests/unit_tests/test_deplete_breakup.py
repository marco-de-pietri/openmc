"""Tests for inelastic breakup paths in depletion chains."""

from types import SimpleNamespace

import pytest

import openmc.data
from openmc.deplete import Chain, Nuclide, ReactionRates


@pytest.fixture
def make_chain(monkeypatch):
    def decay(name):
        z, a, _ = openmc.data.zam(name)
        return SimpleNamespace(nuclide=dict(name=name, atomic_number=z,
                                            mass_number=a, stable=True))

    monkeypatch.setattr(openmc.data, 'Decay', decay)

    def build_chain(
        records, inventory=('O16', 'C12', 'Li7', 'He4', 'H3', 'Ne20'), **kwargs
    ):
        evaluations = []
        for parent, sections in records.items():
            ev = object.__new__(openmc.data.endf.Evaluation)
            z, a, m = openmc.data.zam(parent)
            ev.target = dict(atomic_number=z, mass_number=a, isomeric_state=m)
            ev.section = {}
            ev.reaction_list = []
            for mt, lr, qm in sections:
                ev.reaction_list.append((3, mt, 2, 0))
                head = ''.join(f'{0:11}' for _ in range(6))
                cont = ''.join(f'{x:11}' for x in (qm, -mt, 0, lr, 0, 0))
                ev.section[3, mt] = f'{head}\n{cont}\n'
            evaluations.append(ev)
        return Chain.from_endf(inventory, [], evaluations, progress=False, **kwargs)

    return build_chain


@pytest.fixture
def breakup_records():
    return {
        'O16': [(4, 0, 0), (51, 0, -1), (56, 22, -4), (91, 22, -4)],
        'Li7': [(53, 33, -2)],
        'C12': [(52, 23, -7)],
        'Ne20': [(22, 0, -5)],
    }


def test_breakup_round_trip(make_chain, breakup_records, tmp_path):
    """Check breakup products and chain reduction after an XML round trip."""
    paths = {'O16': ('(n,na)', 'C12', -4, {'C12': 1, 'He4': 1}),
             'Li7': ('(n,nt)', 'He4', -2, {'He4': 1, 'H3': 1}),
             'C12': ('(n,n3a)', None, -7, {'He4': 3})}
    chain = make_chain(breakup_records, breakup_reactions=True, reactions=['(n,na)', '(n,nt)', '(n,n3a)'])
    assert chain['Ne20'].reactions == [('(n,na)', 'O16', -22, 1.0)]
    assert not chain['Ne20'].breakup_reactions
    assert len(chain.reactions) == 3
    path = tmp_path / 'chain.xml'
    chain.export_to_xml(path)
    chain = Chain.from_xml(path)
    rates = ReactionRates(['1'], chain.nuclide_dict, chain.reactions)
    for parent, (name, _, _, _) in paths.items():
        rates.set('1', parent, name, 1.0)
    matrix = chain.form_matrix(rates[0]).toarray()
    for parent, (name, target, qm, products) in paths.items():
        nuc = chain[parent]
        assert nuc.reactions == [(name, target, qm, 1.0)]
        assert nuc.breakup_reactions == [name]
        assert nuc.validate()
        for iso, row in chain.nuclide_dict.items():
            expected = -1 if iso == parent else products.get(iso, 0)
            assert matrix[row, chain.nuclide_dict[parent]] == expected
        for depth in (0, 1):
            reduced = chain.reduce([parent], level=depth)
            expected_nuclides = {parent, *products} if depth else {parent}
            assert set(reduced.nuclide_dict) == expected_nuclides
            copy = reduced[parent]
            assert copy.reactions == [(name, target if depth else None, qm, 1.0)]
            assert copy.breakup_reactions == [name]
            assert copy.breakup_reactions is not nuc.breakup_reactions


@pytest.mark.parametrize('options', [{}, {'breakup_reactions': False}])
def test_breakup_requires_opt_in(make_chain, breakup_records, options):
    chain = make_chain(breakup_records, reactions=['(n,na)'], **options)
    assert not chain['O16'].reactions
    assert chain['Ne20'].reactions == [('(n,na)', 'O16', -22, 1.0)]


def test_breakup_requires_bool(make_chain):
    with pytest.raises(TypeError, match='breakup_reactions'):
        make_chain({}, breakup_reactions='false')


def test_breakup_reaction_selection(make_chain, breakup_records):
    aggregate = make_chain({'O16': [(4, 22, -4)]}, breakup_reactions=True, reactions=['(n,na)'])
    assert not aggregate.reactions
    assert not make_chain(breakup_records).reactions
    assert not make_chain(breakup_records, breakup_reactions=True, reactions=[]).reactions
    ordinary = ['(n,2n)', '(n,3n)', '(n,4n)', '(n,gamma)', '(n,p)', '(n,a)']
    for selection in (ordinary, reversed(ordinary)):
        explicit = make_chain(breakup_records, breakup_reactions=True, reactions=selection)
        assert not explicit.reactions
    selected = make_chain(breakup_records, breakup_reactions=True, reactions=iter(['(n,na)']))
    assert selected.reactions == ['(n,na)']
    assert selected['Ne20'].reactions == [('(n,na)', 'O16', -22, 1.0)]
    assert not selected['Ne20'].breakup_reactions
    assert selected['O16'].breakup_reactions == ['(n,na)']


@pytest.mark.parametrize('conflict', [(52, 22, -3), (4, 22, -4), (22, 0, -4)])
def test_conflicting_breakup(make_chain, conflict):
    with pytest.raises(ValueError, match='O16') as exc:
        make_chain({'O16': [(56, 22, -4), conflict]}, breakup_reactions=True, reactions=['(n,na)'])
    assert '56' in str(exc.value) and str(conflict[0]) in str(exc.value)


@pytest.mark.parametrize('parent,lr,name,target,qm', [
    ('C12', 23, '(n,n3a)', None, -7_275_000),
    ('O17', 24, '(n,2na)', 'C12', -11_305_000),
])
def test_mass_balance_resolves_breakup(make_chain, parent, lr, name, target, qm):
    inventory = {'O17', 'C12', 'He4'}
    with pytest.warns(UserWarning, match='mass-balance Q'):
        chain = make_chain({parent: [(91, lr, 0), (52, lr, qm)]},
                           inventory=inventory, breakup_reactions=True, reactions=[name])
    assert chain[parent].reactions == [(name, target, qm, 1.)]


def test_ambiguous_breakup_q(make_chain):
    with pytest.raises(ValueError, match='no unique match'):
        make_chain({'C12': [(52, 23, -7_274_750), (91, 23, -7_275_000)]},
                   breakup_reactions=True, reactions=['(n,n3a)'])


def test_invalid_breakup_products(make_chain):
    with pytest.raises(ValueError, match='balance for Li7'):
        make_chain({'Li7': [(51, 23, -1)]}, breakup_reactions=True, reactions=['(n,n3a)'])
    with pytest.raises(ValueError, match='Missing breakup product C12'):
        make_chain({'O16': [(56, 22, -4)]}, inventory=['O16', 'He4'],
                   breakup_reactions=True, reactions=['(n,na)'])
    with pytest.raises(ValueError, match='Missing breakup product He4'):
        make_chain({'O16': [(56, 22, -4)]}, inventory=['O16', 'C12'],
                   breakup_reactions=True, reactions=['(n,na)'])


@pytest.mark.parametrize('marked_first', [True, False])
def test_mixed_breakup_xml(tmp_path, marked_first):
    entries = ['<reaction type="(n,na)" breakup="true"/>',
               '<reaction type="(n,na)" target="C12"/>']
    if not marked_first:
        entries.reverse()
    path = tmp_path / 'mixed.xml'
    path.write_text('<depletion_chain><nuclide name="O16">' + ''.join(entries)
                    + '</nuclide></depletion_chain>')
    with pytest.raises(ValueError, match='Mixed breakup and ordinary'):
        Chain.from_xml(path)


@pytest.mark.parametrize('name, branches', [
    ('(n,gamma)', [1.0]), ('(n,na)', []), ('(n,na)', [0.5]),
    ('(n,na)', [0.5, 0.5])])
def test_invalid_breakup_marker(name, branches):
    parent = Nuclide('O16')
    parent.breakup_reactions.append(name)
    for ratio in branches:
        parent.add_reaction(name, 'C12', 0.0, ratio, breakup=True)
    with pytest.raises(ValueError, match='Invalid breakup path'):
        parent.validate()
    assert not parent.validate(strict=False, quiet=True)


@pytest.mark.parametrize('lr, name, target, products', [
    (22, '(n,na)', 'C12', {'C12': 1, 'He4': 1}),
    (23, '(n,n3a)', 'He4', {'He4': 4}),
    (24, '(n,2na)', 'C11', {'C11': 1, 'He4': 1}),
    (25, '(n,3na)', 'C10', {'C10': 1, 'He4': 1}),
    (28, '(n,np)', 'N15', {'N15': 1, 'H1': 1}),
    (29, '(n,n2a)', 'Be8', {'Be8': 1, 'He4': 2}),
    (30, '(n,2n2a)', 'Be7', {'Be7': 1, 'He4': 2}),
    (32, '(n,nd)', 'N14', {'N14': 1, 'H2': 1}),
    (33, '(n,nt)', 'N13', {'N13': 1, 'H3': 1}),
    (34, '(n,n3He)', 'C13', {'C13': 1, 'He3': 1}),
    (35, '(n,nd2a)', 'Li6', {'Li6': 1, 'H2': 1, 'He4': 2}),
    (36, '(n,nt2a)', 'Li5', {'Li5': 1, 'H3': 1, 'He4': 2}),
])
def test_breakup_inventory(make_chain, lr, name, target, products):
    chain = make_chain({'O16': [(52, lr, -7)]}, inventory=['O16', *products],
                     breakup_reactions=True, reactions=[name])
    parent = chain['O16']
    assert parent.reactions == [(name, target, -7, 1.0)]
    assert parent.breakup_reactions == [name]
    assert parent.validate()

    rates = ReactionRates(['1'], chain.nuclide_dict, chain.reactions)
    rates.set('1', 'O16', name, 1.0)
    matrix = chain.form_matrix(rates[0]).toarray()
    for iso, row in chain.nuclide_dict.items():
        expected = -1 if iso == 'O16' else products[iso]
        assert matrix[row, chain.nuclide_dict['O16']] == expected


def test_breakup_mode_selection(make_chain):
    records = {'O16': [(mt, lr, -7) for mt, lr in
                       enumerate((0, 1, 31, 39, 40, 24, 30), start=51)]}
    inventory = ['O16', 'C11', 'Be7', 'He4']
    chain = make_chain(records, inventory, breakup_reactions=True, reactions=['(n,2na)', '(n,2n2a)'])
    assert set(chain.reactions) == {'(n,2na)', '(n,2n2a)'}
    assert chain['O16'].breakup_reactions == chain.reactions
    selected = make_chain(records, inventory, breakup_reactions=True, reactions=['(n,2n2a)'])
    assert selected['O16'].reactions == [('(n,2n2a)', 'Be7', -7, 1.0)]


def test_li7_breakup_components(make_chain, tmp_path):
    qm = -2_466_700
    sections = [(51, 0, 0)] + [(mt, 33, qm) for mt in range(52, 83)]
    chain = make_chain({'Li7': sections}, inventory=['Li7', 'He4', 'H3'],
                     breakup_reactions=True, reactions=['(n,nt)'])
    path = tmp_path / 'chain.xml'
    chain.export_to_xml(path)
    chain = Chain.from_xml(path)

    assert chain.reactions == ['(n,nt)']
    assert chain['Li7'].reactions == [('(n,nt)', 'He4', qm, 1.0)]
    assert chain['Li7'].breakup_reactions == ['(n,nt)']

    rates = ReactionRates(['1'], chain.nuclide_dict, chain.reactions)
    rates.set('1', 'Li7', '(n,nt)', 1.0)
    matrix = chain.form_matrix(rates[0]).toarray()
    parent = chain.nuclide_dict['Li7']
    for name, expected in [('Li7', -1), ('He4', 1), ('H3', 1)]:
        assert matrix[chain.nuclide_dict[name], parent] == expected
