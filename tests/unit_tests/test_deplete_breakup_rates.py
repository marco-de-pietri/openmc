"""Coupled depletion rates from parent-specific inelastic levels."""

from ctypes import c_int
from pathlib import Path
from shutil import copyfile
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

import openmc
import openmc.lib
from openmc.deplete import Chain, CoupledOperator, Nuclide
from openmc.deplete import helpers


def _make_breakup_chain(parent, reaction):
    chain = Chain()
    nuc = Nuclide(parent)
    nuc.add_reaction(reaction, None, 0., 1., breakup=True)
    chain.add_nuclide(nuc)
    return chain


@pytest.mark.parametrize(
    'flag', [22, 23, 24, 25, 28, 29, 30, 32, 33, 34, 35, 36]
)
def test_resolve_breakup_rates(monkeypatch, flag):
    """Resolve supported LR flags and reject incomplete metadata."""
    reaction = openmc.data.REACTION_NAME[flag]
    op = object.__new__(CoupledOperator)
    op.chain = _make_breakup_chain('O16', reaction)
    op.chain.add_nuclide(Nuclide('Ne20'))
    op._burnable_nucs = ['O16', 'Ne20']
    metadata = {
        4: flag, 91: flag, 51: flag, 57: 31, 58: 0, 59: 39, 60: 40, 61: 1
    }
    nuclide = SimpleNamespace(reaction_breakup_flags=metadata)
    monkeypatch.setattr(openmc.lib, 'nuclides', {'O16': nuclide})
    assert op._resolve_breakup_rates() == {('O16', reaction): (51, 91)}
    metadata[52] = None
    with pytest.raises(openmc.exceptions.DataError) as exc:
        op._resolve_breakup_rates()
    assert all(s in str(exc.value) for s in ('O16', reaction, '52'))
    metadata.update({51: 0, 52: 0, 91: 0})
    with pytest.raises(
        openmc.exceptions.DataError, match='no matching'
    ) as exc:
        op._resolve_breakup_rates()
    assert reaction in str(exc.value)


@pytest.mark.parametrize('mode, direct_nuclides', [
    ('direct', None),
    ('flux', None),
    ('flux', ['O16']),
    ('flux', ['Ne20']),
])
def test_breakup_rate_tally_columns(monkeypatch, mode, direct_nuclides):
    """Select reaction rates by nuclide as the active tally list changes."""
    monkeypatch.setattr(helpers, 'Tally', SimpleNamespace)
    monkeypatch.setattr(helpers, 'MaterialFilter', list)
    monkeypatch.setattr(helpers, 'EnergyFilter', list)
    monkeypatch.setattr(openmc.lib, 'nuclides', {
        name: SimpleNamespace(collapse_rate=lambda mt, t, e, f: mt * f[0])
        for name in ('O16', 'Ne20')})
    reaction = '(n,na)'
    if mode == 'direct':
        helper = helpers.DirectReactionRateHelper(3, 2)
    else:
        helper = helpers.FluxCollapseHelper(
            3, 2, [0.0, 1.0], [reaction], direct_nuclides
        )
    mat = SimpleNamespace(nuclides=['O16', 'Ne20'], temperature=300.)
    helper.generate_tallies(
        [mat, mat], [reaction],
        composite_reaction_mts={('O16', reaction): (51, 52)},
    )
    for names in (['Ne20'], ['O16', 'Ne20'], ['Ne20', 'O16']):
        helper.nuclides = names
        scores = helper._rate_tally.scores
        tallied_nuclides = helper._rate_tally.nuclides
        shape = (2, len(tallied_nuclides), len(scores))
        means = np.arange(1.0, 1 + np.prod(shape)).reshape(shape)
        helper._rate_tally_means_cache = means.reshape(2, -1)
        helper._flux_tally_means_cache = np.array([10., 20.])
        for mat_index in range(2):
            nuc_indices = [2 if name == 'Ne20' else 0 for name in names]
            rates = helper.get_material_rates(mat_index, nuc_indices, [1])
            for name, nuc_index in zip(names, nuc_indices):
                if name in tallied_nuclides:
                    columns = ('51', '52') if name == 'O16' else (reaction,)
                    row = means[mat_index, tallied_nuclides.index(name)]
                    expected = sum(row[scores.index(s)] for s in columns)
                else:
                    mt_sum = 51 + 52 if name == 'O16' else 22
                    expected = mt_sum * (mat_index + 1) * 10.0
                assert rates[nuc_index, 1] == expected


@pytest.mark.parametrize('mode', ['direct', 'flux'])
def test_coupled_breakup_transport(run_in_tmpdir, mode):
    """Compare coupled depletion rates with separately scored partials."""
    # Copy a single isotope and give two real partials controlled LR metadata.
    library = openmc.data.DataLibrary.from_xml()
    hydrogen = library.get_by_material('H1')['path']
    path = Path('O16.h5').resolve()
    copyfile(library.get_by_material('O16')['path'], path)
    with h5py.File(path, 'r+') as data:
        for reaction in data['O16/reactions'].values():
            mt = reaction.attrs['mt']
            if 51 <= mt <= 91:
                reaction.attrs['breakup_flag'] = int(mt in (51, 52)) * 22
        elastic = data['O16/reactions/reaction_002'].attrs
        if 'breakup_flag' in elastic:
            del elastic['breakup_flag']
    library = openmc.data.DataLibrary()
    library.register_file(path)
    library.register_file(hydrogen)
    library.export_to_xml('xs.xml')

    mat = openmc.Material()
    mat.add_nuclide('H1', 1.)
    mat.depletable = True
    mat.volume = 1.
    mat.set_density('g/cm3', 1.)
    surface = openmc.Sphere(r=2., boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-surface)
    model = openmc.Model(openmc.Geometry([cell]), openmc.Materials([mat]))
    model.materials.cross_sections = Path('xs.xml').resolve()
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 3
    model.settings.particles = 100
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([14.0e6], [1.0]),
    )
    energies = [0.0, 10.0e6, 20.0e6]
    op = CoupledOperator(
        model, _make_breakup_chain('O16', '(n,na)'),
        normalization_mode='source-rate', reaction_rate_mode=mode,
        reaction_rate_opts={'energies': energies} if mode == 'flux' else {},
    )
    for mts in ((51, 52), (52,)):
        try:
            initial = op.initial_condition()
            nuc = openmc.lib.nuclides['O16']
            metadata = nuc.reaction_breakup_flags
            assert metadata[2] is None
            assert metadata[53] == 0
            assert metadata[52] == 22
            assert metadata[51] == (22 if 51 in mts else 0)
            with pytest.raises(openmc.exceptions.OutOfBoundsError):
                openmc.lib.Nuclide(-1).reaction_breakup_flags
            mts_buffer = (c_int * 1)(-99)
            lrs_buffer = (c_int * 1)(-99)
            with pytest.raises(
                openmc.exceptions.OpenMCError,
                match='Insufficient space for reaction metadata',
            ):
                openmc.lib._dll.openmc_nuclide_reaction_breakup_flags(
                    nuc._index, 1, mts_buffer, lrs_buffer)
            assert mts_buffer[0] == lrs_buffer[0] == -99
            # Activate the initially absent marked parent after tally creation.
            initial[0][op.number.index_nuc['O16']] = 1.e22
            reference = openmc.lib.Tally()
            lib_material = openmc.lib.materials[mat.id]
            filters = [openmc.lib.MaterialFilter([lib_material])]
            if mode == 'direct':
                reference.nuclides = ['O16']
                reference.multiply_density = False
                reference.scores = ['51', '52']
            else:
                filters.append(openmc.lib.EnergyFilter(energies))
                reference.scores = ['flux']
            reference.filters = filters
            result = op(initial, 1.)
            means = reference.mean.ravel()
            if mode == 'direct':
                expected = sum(means[mt - 51] for mt in mts)
            else:
                temp = lib_material.temperature
                expected = sum(
                    nuc.collapse_rate(mt, temp, energies, means) for mt in mts
                )
            assert expected > 0.
            assert result.rates[0, 0, 0] * 1.e24 == pytest.approx(expected)
        finally:
            op.finalize()
        with h5py.File(path, 'r+') as data:
            data['O16/reactions/reaction_051'].attrs['breakup_flag'] = 0
