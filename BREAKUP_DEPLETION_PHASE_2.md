# Phase 2 — Resolve neutron inelastic breakup rates in coupled depletion

For each participating parent and reaction marked by Phase 1, return one event rate per target atom in the existing `ReactionRates` column. Obtain the contributing partial cross sections from the continuous-energy neutron library loaded for this run. The chain supplies the parent and conventional reaction name. The transport library supplies the MTs and their LR flags.

## Resolution rules

1. Use the parent’s `Nuclide.breakup_reactions` and the existing `BREAKUP_MODES` mapping. The supported LR values are 22–25, 28–30, and 32–36. LR1, LR31, LR39, and LR40 remain outside scope.
2. For each marked pair, inspect that parent’s loaded MT51–91 sections, inclusive. Collect the distinct MTs whose LR maps to the requested reaction name, in increasing MT order. An explicit LR0 is an ordinary level. MT4 is the aggregate and contributes no separate component.
3. If any candidate section lacks LR metadata, its contribution cannot be classified. Raise an error before transport, even if another component matches. Include parent, conventional reaction, and affected MTs. If all candidates are classified but none match, report the parent and requested reaction. Never silently substitute zero or the conventional MT for a marked path.
4. Unmarked pairs use the existing ordinary score and collapse logic, including fission. A conventional name can refer to different rate sources on different parents. Require LR metadata only for parents with participating marked paths.

For example, a Li7 library carrying LR33 on MT52–82 supplies 31 components for its marked `(n,nt)` column. Discovery uses the metadata in each library. A reduced chain can retain the marker while dropping its target, so resolution must also work for a null target.

## Implementation

### Expose the loaded reaction metadata

`Reaction` already stores `mt_` and `breakup_flag_`. The HDF5 reader in `src/reaction.cpp` leaves the latter at −1 when the attribute is absent. Add a read-only query alongside the nuclide functions in `src/nuclide.cpp`, declare it in `include/openmc/capi.h`, and bind it on `openmc.lib.Nuclide` in `openmc/lib/nuclide.py`. Return the loaded MT/LR pairs, preserving the distinction between missing LR and zero. Use the existing C API index checks, error codes, and ctypes error handler. Return Python-owned values so callers can safely retain the result after finalization.

### Resolve once per library initialization

In `CoupledOperator.initial_condition`, prepare the mapping after library initialization and before the rate helpers create their tallies. Scope discovery to marked parents eligible for neutron reaction rates, including eligible nuclides that can appear at a later depletion step. Load those nuclides through the existing library interface if necessary. Resolve once per parent, then pass the mapping to the coupled rate helper. Keep the mapping on the operator or helper instance and rebuild it on reinitialization. Under MPI, each rank needs the same score order even when its local materials differ.

### Direct tally mode

In `DirectReactionRateHelper.generate_tallies`, build a deterministic union of ordinary scores and required partial MT scores. Retain the conventional score wherever an ordinary pair needs it. Deduplicate scores by transport identity and record their positions in the actual tally. Keep the existing material filter and per-atom scoring with `multiply_density = False`.

In `get_material_rates`, index the tally by material, current tally nuclide, and transport score. The current flat iteration assumes one tally score per output reaction. With expanded scores, sum only the component means assigned to that parent and output column. Assign an ordinary pair its existing single score. Preserve `nuc_index` and `rx_index` when writing the output, including when the active nuclide list changes between steps.

### Flux-collapse mode and direct overrides

In `FluxCollapseHelper.get_material_rates`, use the existing material temperature, group boundaries, and tallied group flux in each call to `Nuclide.collapse_rate`. For a marked pair, sum the per-atom rates from the resolved MTs. For ordinary pairs, keep the existing single-MT calculation.

When `reaction_rate_opts` requests direct scoring for selected reactions or nuclides, expand the override tally using the same mapping as direct mode. Index its results using the actual number of transport scores. A selected pair takes its complete rate from the override tally once. Other pairs use flux collapse. Exercise both an explicit override nuclide list and the default that follows the active nuclides.

Both modes return conventional reaction columns in chain order. The existing operator handles normalization. Atom multiplicities enter when the depletion matrix consumes the event rate. Integration with standalone MicroXS generation and independent depletion belongs to Phase 3.

## Acceptance tests

| Case | Expected result |
|---|---|
| Loaded C API and Python metadata | A fixture with matching LR, LR0, and an omitted LR attribute returns three distinguishable values. Invalid nuclide indices follow the existing error convention. |
| Component discovery | Parameterize the twelve supported flags. Include MT51, MT91, unrelated levels, excluded flags, and aggregate MT4. Verify exact component lists. A legacy library remains usable for unmarked parents. |
| Incomplete data | A missing LR on a candidate fails even alongside a matching component. A fully classified parent with no matches also fails before scoring. Check the diagnostic and exercise a reduced chain with a marked null target. |
| Parent-local direct rates | Give marked O16 two LR22 partials with distinct means and unmarked Ne20 an ordinary MT22 mean. Set unrelated bins to nonzero values. Assert that each parent receives only its assigned rate. Repeat with a different nuclide ordering and more than one material. |
| Flux and overrides | Compare a marked column with the sum of individual component collapses for the same temperature and flux. Select only one parent for a direct override and give its direct and collapsed rates distinct values to detect double counting. |
| Lifecycle and compatibility | Exercise a marked nuclide becoming active after initialization and two successive library initializations with different component lists. Run existing ordinary-chain coupled-depletion tests, including MPI where available. |

Run focused Python helper and binding tests first, then a coupled-operator test for each rate mode using a small compatible neutron library with known LR metadata. Verify individual component tallies or collapses independently of the resolver’s output. ENDF/B-VIII.1 benchmark comparisons and time-dependent inventory verification belong to [Phase 4](BREAKUP_DEPLETION_PHASE_4.md).
