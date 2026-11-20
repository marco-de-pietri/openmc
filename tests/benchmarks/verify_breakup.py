"""Opt-in ENDF/B-VIII.1 breakup benchmark. See BREAKUP_DEPLETION_PHASE_4.md."""

import hashlib
from io import StringIO
import os
from pathlib import Path

import h5py
import numpy as np
import pytest
from scipy.integrate import trapezoid

import openmc
import openmc.lib
from openmc.data import Evaluation
from openmc.data.endf import get_cont_record, get_head_record
from openmc.deplete import (
    Chain, CoupledOperator, IndependentOperator, MicroXS, Nuclide,
    PredictorIntegrator, Results, get_microxs_and_flux,
)
from openmc.deplete.microxs import _model_breakup_mts


ENERGIES = np.array([1., 8.e6, 16.e6, 20.e6])
SPECTRUM = np.array([1., 3., 2.])
TEMPERATURE = 293.6
RTOL = 1.e-10


def check(actual, expected):
    np.testing.assert_allclose(actual, expected, rtol=RTOL, atol=0.)


def reference_levels(evaluation):
    levels = {}
    for mf, mt in evaluation.section:
        if mf == 3 and 51 <= mt <= 91:
            records = StringIO(evaluation.section[mf, mt])
            get_head_record(records)
            levels[mt] = get_cont_record(records)[3]
    return levels


def check_metadata(index, endf, parent, lr, expected_mts):
    path = Path(openmc.data.DataLibrary.from_xml(index).get_by_material(parent)['path'])
    for file in (endf, path):
        with file.open('rb') as stream:
            print(parent, file, hashlib.file_digest(stream, 'sha256').hexdigest())
    evaluation = Evaluation(endf)
    assert evaluation.info['library'] == ('ENDF/B', 8, 1)
    assert evaluation.gnds_name == parent
    levels = reference_levels(evaluation)
    mts = sorted(mt for mt, flag in levels.items() if flag == lr)
    with h5py.File(path) as f:
        observed = {int(r.attrs['mt']): r.attrs.get('breakup_flag')
                    for r in f[parent]['reactions'].values() if 51 <= r.attrs['mt'] <= 91}
    print(parent, 'MF3', levels, 'HDF5', observed, 'reference MTs', mts)
    assert observed == levels
    if expected_mts is not None:
        assert mts == expected_mts
    return mts, path


def group_cross_sections(data, mts):
    """Integrate linear-linear partials over each group's threshold overlap."""
    group_integrals = np.zeros(len(ENERGIES) - 1)
    for mt in mts:
        partial = data[mt].xs['294K']
        for i, (low, high) in enumerate(zip(ENERGIES[:-1], ENERGIES[1:])):
            left = max(low, partial.x[0])
            right = min(high, partial.x[-1])
            if left < right:
                interior = partial.x[(partial.x > left) & (partial.x < right)]
                grid = np.r_[left, interior, right]
                group_integrals[i] += trapezoid(partial(grid), grid)
    return group_integrals / np.diff(ENERGIES)


def fixed_source_model(parent, index):
    mat = openmc.Material(components={parent: 1.}, density=0.01, density_units='atom/b-cm')
    mat.depletable = True
    mat.volume = 1.
    surface = openmc.Sphere(r=(3 / (4 * np.pi))**(1 / 3), boundary_type='vacuum')
    model = openmc.Model(openmc.Geometry([openmc.Cell(region=-surface, fill=mat)]),
                         openmc.Materials([mat]))
    model.materials.cross_sections = index
    mat.temperature = TEMPERATURE
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 5
    model.settings.particles = 1000
    model.settings.seed = 53
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(), energy=openmc.stats.Discrete([14.e6, 18.e6], [.5, .5]))
    return model


def check_inventory(material, chain, table, parent, reaction, products, reference):
    flux = 2.
    source = 1.e23
    initial_atoms = 1.e22
    times = [0., 1., 3.]
    rate_constant = reference * flux * source * 1.e-24
    remaining = initial_atoms * np.exp(-rate_constant * np.array(times))

    def run(directory, steps, previous=None):
        op = IndependentOperator(
            [material], [np.array([flux])], [table], chain,
            normalization_mode='source-rate', prev_results=previous)
        op.output_dir = directory
        PredictorIntegrator(op, steps, source_rates=source).integrate(write_rates=True)
        return Results(op.output_dir / 'depletion_results.h5')

    full = run('full', [1., 2.])
    restarted = run('split', [2.], run('split', [1.]))
    inventories = {parent: remaining}
    for name, count in products.items():
        inventories[name] = count * (initial_atoms - remaining)
    for label, results in [('full', full), ('split', restarted)]:
        t, rate = results.get_reaction_rate(material, parent, reaction)
        print(parent, label, 'times', t, 'event rates', rate, 'expected', remaining * rate_constant)
        check(t, times)
        check(rate, remaining * rate_constant)
        for name, expected in inventories.items():
            t, atoms = results.get_atoms(material, name)
            print(name, label, 'atoms', atoms, 'expected', expected)
            check(t, times)
            check(atoms, expected)


@pytest.mark.parametrize('parent,filename,lr,reaction,target,products,expected_mts', [
    ('Li7', 'n-003_Li_007.endf', 33, '(n,nt)', 'He4', {'He4': 1, 'H3': 1}, list(range(52, 83))),
    ('C12', 'n-006_C_012.endf', 23, '(n,n3a)', None, {'He4': 3}, [*range(52, 63), 91]),
    ('O16', 'n-008_O_016.endf', 22, '(n,na)', 'C12', {'C12': 1, 'He4': 1}, None),
])
def test_evaluated_breakup(parent, filename, lr, reaction, target, products, expected_mts,
                           run_in_tmpdir, monkeypatch):
    index = Path(os.environ['BREAKUP_XS_XML'])
    endf = Path(os.environ['BREAKUP_ENDF_DIR']) / filename
    mts, path = check_metadata(index, endf, parent, lr, expected_mts)

    chain = Chain()
    nuc = Nuclide(parent)
    nuc.add_reaction(reaction, target, 0., 1., breakup=True)
    chain.add_nuclide(nuc)
    for name in products:
        chain.add_nuclide(Nuclide(name))

    model = fixed_source_model(parent, index)
    mat = model.materials[0]
    if not mts:
        with pytest.raises(openmc.exceptions.DataError, match=f'{parent}.*no matching'):
            _model_breakup_mts(model, chain, [parent], [reaction])
        op = CoupledOperator(model, chain, normalization_mode='source-rate')
        try:
            with pytest.raises(openmc.exceptions.DataError, match=f'{parent}.*no matching'):
                op.initial_condition()
        finally:
            op.finalize()
        return
    selected = _model_breakup_mts(model, chain, [parent], [reaction])[parent, reaction]
    print(parent, 'selected MTs', selected)
    assert list(selected) == mts

    filters = [openmc.MaterialFilter(mat), openmc.EnergyFilter(ENERGIES)]
    production = 'H3-production' if parent == 'Li7' else 'He4-production'
    ref = openmc.Tally(filters=filters, nuclides=[parent],
                       scores=[*map(str, mts), production])
    ref.multiply_density = False
    flux = openmc.Tally(filters=filters, scores=['flux'])
    model.tallies = [ref, flux]
    with openmc.StatePoint(model.run(cwd='baseline', output=False)) as sp:
        baseline = sp.tallies[ref.id].mean.copy()
    _, direct = get_microxs_and_flux(
        model, [mat], nuclides=[parent], reactions=[reaction], chain_file=chain,
        energies=ENERGIES, include_model_tallies=True, path_statepoint='marked.h5',
        run_kwargs={'output': False})
    with openmc.StatePoint('marked.h5') as sp:
        rates = sp.tallies[ref.id].mean
        fluence = sp.tallies[flux.id].mean.ravel()
    print(parent, production, rates[:, 0, -1], 'baseline', baseline[:, 0, -1])
    check(rates, baseline)
    assert rates[:, 0, -1].sum() > 0.
    partial_sum = rates[:, 0, :-1].sum(axis=1)
    assert partial_sum.sum() > 0.
    expected_direct = np.divide(partial_sum, fluence, out=np.zeros_like(fluence),
                                where=fluence > 0.)
    print(parent, 'direct XS', direct[0][parent, reaction], 'reference', expected_direct)
    check(direct[0][parent, reaction], expected_direct)

    group_xs = group_cross_sections(openmc.data.IncidentNeutron.from_hdf5(path), mts)
    for mode, rate_sum in [('direct', partial_sum.sum()), ('flux', group_xs @ fluence)]:
        op = CoupledOperator(
            model, chain, normalization_mode='source-rate', reaction_rate_mode=mode,
            reaction_rate_opts={'energies': ENERGIES} if mode == 'flux' else {})
        try:
            initial = op.initial_condition()
            assert list(op._resolve_breakup_rates()[parent, reaction]) == mts
            result = op(initial, 1.)
            check(openmc.lib.tallies[ref.id].mean, baseline[:, 0, :])
            coupled = result.rates.get(str(mat.id), parent, reaction)
            expected_rate = rate_sum * 1.e-24 / mat.volume
            print(parent, mode, 'coupled rate', coupled, 'reference', expected_rate)
            check(coupled, expected_rate)
        finally:
            op.finalize()

    reference = np.average(group_xs, weights=SPECTRUM)
    collapsed = MicroXS.from_multigroup_flux(
        ENERGIES, SPECTRUM, model=model, chain_file=chain, nuclides=[parent],
        reactions=[reaction], temperature=TEMPERATURE)
    print(parent, 'collapsed XS', collapsed[parent, reaction], 'reference', reference)
    assert reference > 0.
    check(collapsed[parent, reaction], reference)
    monkeypatch.setattr(openmc.deplete.pool, 'USE_MULTIPROCESSING', False)
    check_inventory(mat, chain, collapsed, parent, reaction, products, reference)
