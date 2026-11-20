# Phase 4 — Verify breakup rates and inventories

**Goal.** Trace a breakup event from its evaluated neutron data through the transport tally to the final atom counts. Use an independent reference at each step so a wrong MT list shared by both depletion operators cannot make a comparison pass.

## Starting point

Phases 1–3 are committed through `c3c4aad2b`. The [chain tests](tests/unit_tests/test_deplete_breakup.py) already cover supported LR flags, conflicting representations, XML round trips, chain reduction, and matrix columns. The [breakup-rate tests](tests/unit_tests/test_deplete_breakup_rates.py) use controlled O16 flags to exercise coupled transport, direct and collapsed MicroXS generation, direct overrides, and independent depletion with restart. Retain those portable tests. Phase 4 adds evaluated-data comparisons, light-product multiplicity over time, and fixed-composition transport checks.

## 1. Verify component discovery against the evaluation

Read MT and LR directly from the original ENDF/B-VIII.1 neutron MF3 records. Select MT51–91 with the relevant LR and record the exact list for each parent. Compare it with the LR attributes in the processed neutron HDF5 file and the MTs selected by depletion. Build the reference independently of the production selector. Record the evaluation and processed file identifiers or checksums alongside the result.

| Parent and channel | Component check | Inventory per event |
|---|---|---|
| Li7 `(n,nt)` | Verify MT52–82 with LR33 | −1 Li7, +1 He4, +1 H3 |
| C12 `(n,n3a)` | Verify MT52–62 and MT91 with LR23 | −1 C12, +3 He4 |
| O16 `(n,na)` | Derive the LR22 list from the evaluation. MT56/LR22 is a controlled fixture in the overall plan. | −1 O16, +1 C12, +1 He4 |

Treat a discrepancy as a result to investigate. Report the parent, MT, and LR at each stage. A missing HDF5 LR attribute and an explicit LR0 must remain distinguishable. Keep the evaluated-data benchmark separately runnable with explicit local input paths so the ordinary test suite can run with its prescribed NNDC data.

## 2. Verify transport rates at fixed composition

- **Direct scoring.** Use a fixed-source model with energies above the breakup thresholds. Tally each reference MT separately for the parent and compare their sum with the conventional breakup entry. For direct MicroXS, compare each group against its component-rate sum divided by the same domain’s flux. Include at least one positive breakup rate.
- **Flux collapse.** Use a specified temperature and a nonuniform group spectrum. Independently integrate each reference partial cross section against the same within-group flux convention used by OpenMC. Sum the partials and compare with the reported one-group cross section. Record energy boundaries, spectrum, temperature, and normalization.
- **Ordinary transport tallies.** Run the same fixed-composition model and random seed with and without the marked depletion path. Compare the individual component and particle-production tallies, including a nonzero production score. Keep tally filters, estimators, and source normalization identical. The existing mixed-parent and direct-override tests continue to cover selection between ordinary and breakup scores.

## 3. Verify atom counts and saved rates

For each of Li7 and C12, isolate the breakup path in a reduced chain, with stable products, fixed flux and source normalization, and no competing losses. Obtain the per-parent event rate constant λ from the independent component-rate reference. At elapsed time t, the expected parent inventory is N(t) = N(0) exp(−λt). Each product gains its multiplicity times N(0) − N(t).

Compare the independent depletion history with these values at each saved step. For C12, check that `Results.get_reaction_rate` returns λN(t) under `(n,n3a)` while He4 gains three atoms per event. For Li7, check both He4 and H3. Include a run split at a step boundary with `prev_results` and compare its saved times, rates, and inventories with the uninterrupted schedule. Reuse the O16 restart test pattern.

## Verification and completion

Run the focused chain, reaction-metadata, breakup-rate, MicroXS, independent-operator, and restart tests. Run the affected transport and depletion regressions with the prescribed NNDC data, a strict floating-point build, and `OMP_NUM_THREADS=2`. Pass `--mpi` for an MPI-enabled build. Run the ENDF/B-VIII.1 benchmarks separately with the matching evaluated and processed data.

Record the reference and observed MT lists, rates, and inventories, the comparison tolerances, and the test commands. Use numerical tolerances for deterministic comparisons and account for tally uncertainty when comparing independent Monte Carlo runs. Identify any benchmark still awaiting compatible data. Generating and publishing a distributable chain/library pair belongs to [Phase 5](BREAKUP_DEPLETION_PLAN.md#phase-5--generate-and-publish-a-compatible-data-pair).

## Verification results

The [opt-in benchmark](tests/benchmarks/verify_breakup.py) reads each original MF3 section and compares every MT51–91 LR with the processed file, including zeroes. It exercises both the file-based MicroXS selector and the loaded-library coupled selector. The local installed ENDF/B-VIII.1 HDF5 library lacked LR attributes. For this run, the three original evaluations were processed separately with NJOY 2016.78 and `IncidentNeutron.from_njoy(..., temperatures=[293.6], heatr=False)`. The Li7 input was the original evaluation preserved locally as `n-003_Li_007.endf.bak`.

| Parent | MF3, processed HDF5 and depletion selection | Remaining evaluated levels |
|---|---|---|
| Li7 | MT52–82 inclusive, LR33 | MT51 has LR0 |
| C12 | MT52–62 inclusive and MT91, LR23 | MT51 has LR0 |
| O16 | Empty LR22 selection. Both selectors diagnose the marked channel as having no matching levels. | MT51–57 and MT91 all have LR0, including MT56. The evaluation contains a separate MF3 MT22. |

| Parent | Evaluated file SHA-256 | Processed file SHA-256 |
|---|---|---|
| Li7 | `14d44ff4aacc5d433ed1b7dbe2a8800f90cdb337a5d739c90add5367661c83f8` | `e24c20bc9ce74c7bcfa97061bd219d52e413d580a956d1eaf3ac789891b3f382` |
| C12 | `249b431cec84be40b957f8499862ceb4380e7d6054492b6500f4d3d255a907fb` | `054736277a342db8694a6c087ad2d1ea176fadaf74ef2eb653a4f0ea2e69dedc` |
| O16 | `850acd8a1c1f103a63864437281605abc8d4d3c6b8f444041f24514ad45a71c9` | `15ddebb271378cb6821a998141e5e46d9226071403bd94f57770b5f36c23eafa` |

The transport cases used the same fixed composition, source, seed and ordinary tallies in each run. Component and particle-production means agreed group by group, including positive H3 production for Li7 and He4 production for C12. Energy boundaries were 1 eV, 8 MeV, 16 MeV and 20 MeV. The collapse reference integrated the linear partial cross sections on their native 293.6 K grids, with constant flux per eV within each group. Coupled flux mode used the measured group flux. The prescribed-spectrum MicroXS comparison used group-integrated weights of 1/6, 3/6 and 2/6. All numerical comparisons used relative tolerance `1e-10` and zero absolute tolerance. The fixed-seed comparisons reused the same neutron histories.

| Parent | Direct MicroXS by group (b) | Direct coupled rate at one source neutron/s (1/s per atom) | Flux-mode coupled rate (1/s per atom) | Prescribed-spectrum MicroXS (b) |
|---|---|---|---|---|
| Li7 | 0.09830604, 0.30245123, 0.22514720 | 1.64013425e-25 | 1.73546691e-25 | 0.2639507582 |
| C12 | 0, 0.20979067, 0.30766188 | 1.59776364e-25 | 1.32031781e-25 | 0.1649692346 |

Each entry below gives the common rounded values from the analytic reference, uninterrupted run and restart at 1 s. The inventory calculation used a volume of 1 cm³, flux of 2 per cm² per source neutron and source rate of 1e23 neutrons/s.

| Parent | Quantity | 0 s | 1 s | 3 s |
|---|---|---:|---:|---:|
| Li7 | Parent atoms | 1e22 | 9.48579049e21 | 8.53533528e21 |
| Li7 | H3 atoms and He4 atoms, each | 0 | 5.14209506e20 | 1.46466472e21 |
| Li7 | Saved (n,nt) events/s | 5.27901516e20 | 5.00756319e20 | 4.50581644e20 |
| C12 | Parent atoms | 1e22 | 9.67544513e21 | 9.05759428e21 |
| C12 | He4 atoms | 0 | 9.73664611e20 | 2.82721717e21 |
| C12 | Saved (n,n3a) events/s | 3.29938469e20 | 3.19230155e20 | 2.98844879e20 |

### Reproducing the checks

Use the Python bindings and executable from the same strict-FP build. The verification build used RelWithDebInfo, MPI and two OpenMP threads. The focused and regression command passed 96 tests and skipped nine multiprocessing variants disabled under MPI.

```bash
export OMP_NUM_THREADS=2
export OPENMC_CROSS_SECTIONS=/absolute/path/to/nndc_hdf5/cross_sections.xml
python -m pytest -q --mpi \
  tests/unit_tests/test_{deplete_breakup,data_reaction_metadata,deplete_breakup_rates,deplete_microxs,deplete_independent_operator,deplete_restart}.py \
  tests/regression_tests/{particle_production_fission,deplete_with_transport,deplete_no_transport,microxs}
```

Run the evaluated benchmark explicitly. Supply absolute paths to the directory containing the original neutron ENDF files and the index for the matching LR-preserving HDF5 files with a 293.6 K table. All three cases passed. The `-s` output includes file hashes, complete MT/LR maps and reference and observed values at each saved step.

```bash
export BREAKUP_ENDF_DIR=/absolute/path/to/endfb-viii.1-endf/neutron
export BREAKUP_XS_XML=/absolute/path/to/processed/cross_sections.xml
OPENMC_CROSS_SECTIONS="$BREAKUP_XS_XML" OMP_NUM_THREADS=2 \
  python -m pytest -sq --mpi tests/benchmarks/verify_breakup.py
```
