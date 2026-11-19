"""Breakup reaction rates from parent-specific inelastic levels."""

from ctypes import c_int
from pathlib import Path
from shutil import copyfile
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

import openmc
import openmc.lib
from openmc.deplete import (
    Chain,
    CoupledOperator,
    IndependentOperator,
    MicroXS,
    Nuclide,
    PredictorIntegrator,
    Results,
    get_microxs_and_flux,
    helpers,
)


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


@pytest.fixture
def breakup_model(run_in_tmpdir):
    """Build a fixed-source model with two O16 breakup components."""
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
    fluorine = library.get_by_material('F19')['path']
    library = openmc.data.DataLibrary()
    library.register_file(path)
    library.register_file(hydrogen)
    library.register_file(fluorine)
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
    return model


@pytest.mark.parametrize('mode', ['direct', 'flux'])
def test_coupled_breakup_transport(breakup_model, mode):
    """Compare coupled depletion rates with separately scored partials."""
    model = breakup_model
    mat = model.materials[0]
    path = Path('O16.h5')
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


@pytest.mark.parametrize('initialized', [False, True], ids=['executable', 'lib'])
def test_breakup_microxs_transport(breakup_model, initialized):
    """Compare breakup cross sections with independent tallies and collapse."""
    model = breakup_model
    mat = model.materials[0]
    mat.add_nuclide('F19', 0.5)
    mat.add_nuclide('O16', 0.5)
    chain = _make_breakup_chain('O16', '(n,na)')
    names = ['F19', 'O16']
    reactions = ['(n,gamma)', '(n,na)']
    energies = [0., 9.e6, 15.e6, 20.e6]
    mesh = openmc.RegularMesh()
    mesh.dimension = (2, 1, 1)
    mesh.lower_left = (-1., -1., -1.)
    mesh.upper_right = (2., 1., 1.)
    filters = [openmc.MeshFilter(mesh), openmc.EnergyFilter(energies)]
    ref = openmc.Tally(
        filters=filters, scores=['51', '52', *reactions], nuclides=names
    )
    ref.multiply_density = False
    ref_flux = openmc.Tally(filters=filters, scores=['flux'])
    model.tallies = [ref, ref_flux]
    original = model.tallies
    kwargs = {
        'nuclides': names,
        'reactions': reactions,
        'energies': energies,
        'chain_file': chain,
        'include_model_tallies': True,
    }
    if initialized:
        model.init_lib()
    try:
        fluxes, direct = get_microxs_and_flux(
            model, mesh, path_statepoint='direct.h5', **kwargs)
    finally:
        if initialized:
            model.finalize_lib()
    assert model.tallies is original
    assert list(original) == [ref, ref_flux]
    with openmc.StatePoint('direct.h5') as sp:
        rates = sp.tallies[ref.id].get_reshaped_data()
        flux = sp.tallies[ref_flux.id].get_reshaped_data()[:, :, 0, 0]
    np.testing.assert_allclose(fluxes, flux)
    assert np.all(flux[:, -1] == 0.)
    assert np.all(rates[:, 1, 1, :2] > 0.)
    for i, table in enumerate(direct):
        assert table.nuclides == names
        assert table.reactions == reactions
        for n, name in enumerate(names):
            for j, reaction in enumerate(reactions):
                rate = rates[i, :, n, 2 + j]
                if (name, reaction) == ('O16', '(n,na)'):
                    rate = rates[i, :, n, :2].sum(axis=-1)
                expected = np.divide(rate, flux[i], out=np.zeros_like(rate),
                                     where=flux[i] != 0.)
                np.testing.assert_allclose(table[name, reaction], expected)

    if initialized:
        return
    _, collapsed = get_microxs_and_flux(
        model, mesh, reaction_rate_mode='flux', **kwargs
    )
    with openmc.lib.TemporarySession(model=model):
        for i, table in enumerate(collapsed):
            spectrum = flux[i] / flux[i].sum()
            for name in names:
                nuc = openmc.lib.nuclides[name]
                for reaction in reactions:
                    mts = (openmc.data.REACTION_MT[reaction],)
                    if (name, reaction) == ('O16', '(n,na)'):
                        mts = (51, 52)
                    expected = sum(
                        nuc.collapse_rate(mt, 293.6, energies, spectrum) for mt in mts
                    )
                    assert table[name, reaction][0] == pytest.approx(expected)

    for selected in (['O16'], names):
        opts = {'reactions': ['(n,na)']}
        if selected != names:
            opts['nuclides'] = selected
        totals, hybrid = get_microxs_and_flux(
            model, mesh, reaction_rate_mode='flux', reaction_rate_opts=opts, **kwargs)
        for i, table in enumerate(hybrid):
            assert totals[i] == pytest.approx(flux[i].sum())
            for name in names:
                for reaction in reactions:
                    expected = collapsed[i][name, reaction][0]
                    if reaction == '(n,na)' and name in selected:
                        direct_avg = np.average(
                            direct[i][name, reaction], weights=flux[i]
                        )
                        assert direct_avg != pytest.approx(expected)
                        expected = direct_avg
                    assert table[name, reaction][0] == pytest.approx(expected)


@pytest.mark.parametrize(('flag', 'diagnostic'), [
    pytest.param(None, 'missing LR.*51', id='missing-flag'),
    pytest.param(0, 'no matching', id='no-match'),
])
def test_breakup_microxs_errors(breakup_model, flag, diagnostic):
    """Check both MicroXS entry points report invalid breakup metadata."""
    model = breakup_model
    chain = _make_breakup_chain('O16', '(n,na)')
    kwargs = dict(nuclides=['O16'], chain_file=chain)
    original = model.tallies
    with h5py.File('O16.h5', 'r+') as data:
        levels = data['O16/reactions']
        levels['reaction_052'].attrs['breakup_flag'] = 0
        attrs = levels['reaction_051'].attrs
        if flag is None:
            del attrs['breakup_flag']
        else:
            attrs['breakup_flag'] = flag
    match = r'O16 \(n,na\).*' + diagnostic
    with pytest.raises(openmc.exceptions.DataError, match=match):
        get_microxs_and_flux(model, model.materials, **kwargs)
    with pytest.raises(openmc.exceptions.DataError, match=match):
        MicroXS.from_multigroup_flux([0., 20.e6], [1.], model=model, **kwargs)
    assert model.tallies is original
    table = MicroXS.from_multigroup_flux(
        [0., 20.e6], [1.], reactions=['(n,gamma)'], model=model, **kwargs)
    assert table['O16', '(n,gamma)'][0] > 0.
    with openmc.lib.TemporarySession(model=model):
        # O16 has invalid metadata and U235 exists only in the global index.
        table = MicroXS.from_multigroup_flux(
            [0., 20.e6], [1.], chain_file=chain, nuclides=['F19', 'U235'], model=model)
        assert table['U235', '(n,na)'] == 0.
        assert table['F19', '(n,na)'][0] > 0.


def test_breakup_microxs_independent_restart(breakup_model, monkeypatch):
    """Compare independent depletion and restart with the expected inventories."""
    model = breakup_model
    chain = Chain()
    parent = Nuclide('O16')
    parent.add_reaction('(n,na)', 'C12', 0., 1., breakup=True)
    for nuc in (parent, Nuclide('C12'), Nuclide('He4')):
        chain.add_nuclide(nuc)
    energies = [0., 10.e6, 20.e6]
    kwargs = dict(chain_file=chain, nuclides=['O16'], temperature=600., model=model)
    with openmc.lib.TemporarySession(model=model):
        nuc = openmc.lib.load_nuclide('O16')
        expected = sum(nuc.collapse_rate(mt, 600., energies, [0.25, 0.75])
                       for mt in (51, 52))
        table = MicroXS.from_multigroup_flux(energies, [1., 3.], **kwargs)
        zero = MicroXS.from_multigroup_flux(energies, [0., 0.], **kwargs)
    assert expected > 0.
    assert table['O16', '(n,na)'] == pytest.approx(expected)
    assert zero['O16', '(n,na)'] == 0.
    table.to_hdf5('microxs.h5')
    table = MicroXS.from_hdf5('microxs.h5')
    fuel = openmc.Material(components={'O16': 1.}, density=0.01,
                           density_units='atom/b-cm')
    fuel.depletable = True
    fuel.volume = 1.
    flux = 2.
    source = 1.e23
    monkeypatch.setattr(openmc.deplete.pool, 'USE_MULTIPROCESSING', False)

    def integrate(directory, steps, previous=None):
        operator = IndependentOperator(
            [fuel], [np.array([flux])], [table], chain,
            normalization_mode='source-rate', prev_results=previous,
        )
        operator.output_dir = directory
        integrator = PredictorIntegrator(operator, steps, source_rates=source)
        integrator.integrate(write_rates=True)
        return Results(operator.output_dir / 'depletion_results.h5')

    full = integrate('full', [1., 2.])
    first = integrate('restart', [1.])
    restart = integrate('restart', [2.], first)
    decay_constant = expected * flux * source * 1.e-24
    times = np.array([0., 1., 3.])
    remaining = 1.e22 * np.exp(-times * decay_constant)
    for result in (full, restart):
        t, rate = result.get_reaction_rate(fuel, 'O16', '(n,na)')
        np.testing.assert_allclose(t, times)
        np.testing.assert_allclose(rate, remaining * decay_constant)
        for name in ('O16', 'C12', 'He4'):
            atoms = result.get_atoms(fuel, name)[1]
            expected_atoms = remaining if name == 'O16' else 1.e22 - remaining
            np.testing.assert_allclose(atoms, expected_atoms)
