# Phase 3 — Generate breakup cross sections for independent depletion

**Goal.** Make a marked reaction usable through both MicroXS generation APIs. For each requested parent, reaction and energy group, the conventional reaction column holds the microscopic cross section for **one event**, summed over that parent's matching inelastic components. A user can pass this table and the corresponding flux directly to `IndependentOperator`.

## Assessment of the committed phases

| Commit | Implemented behavior | Consequence for Phase 3 |
|---|---|---|
| Phase 1, `8f0cadb48` | The chain generator records the breakup marker, target and secondaries. XML and chain reduction preserve the marker. Matrix tests cover parent loss and product multiplicity. | The chain already describes what each event produces. For C12 breakup, one event removes one C12 and produces three He4. The factor of three belongs in inventory production, never in the MicroXS value. |
| Phase 2, `aeb1d35fe` | The coupled operator discovers parent-specific partial MTs from the loaded library. Direct and flux helpers use that mapping, including direct overrides. Tests exercise both transport modes and library reinitialization. | Reuse the discovery and score-selection rules. Discovery currently lives on `CoupledOperator`, while standalone MicroXS generation constructs its own tallies and runs its own collapse loop. |

The remaining gap is in the **producer** of the cross sections. Today `get_microxs_and_flux` assigns conventional names directly to tally scores, and `MicroXS.from_multigroup_flux` collapses one conventional MT per reaction. The independent rate helper already multiplies the supplied group cross sections by the group flux and sums them.

## Implementation

1. **Share component selection.** Extract the selection logic from `CoupledOperator._resolve_breakup_rates` into a private helper that accepts the requested parent, marked reaction names and MT-to-LR mapping. Keep data loading with the caller. For each participating marked pair, inspect MT51–91, distinguish LR0 from missing metadata, and retain the Phase 2 diagnostics for incomplete data or no matching levels. Intersect the markers with the requested reactions and nuclides so an unrelated marker cannot block a subset request.
2. **Resolve against the library used for transport.** Direct MicroXS tally scores must be known before `model.run`, including when transport runs in a separate process. Read the marked parents' reaction attributes from the neutron HDF5 files selected by the model's cross-section index and pass those mappings to the shared selector. Reuse this library selection for the subsequent model-based flux collapse. In standalone `from_multigroup_flux`, read the loaded nuclide's `reaction_breakup_flags` in the collapse session. Resolve once per applicable library and reuse the mapping across domains. Preserve caller tallies on success and on errors.
3. **Direct MicroXS.** Use the Phase 2 score builder with the actual tally nuclides. For each domain and parent, read the conventional score or sum its mapped partials separately in each energy group. Divide by that domain's group flux, retaining zero for zero-flux groups. Build the returned `MicroXS` with the requested nuclide and reaction ordering. Partial MT labels belong only to the internal tally.
4. **Flux collapse and overrides.** In `from_multigroup_flux`, sum the individual component collapses at the requested temperature using the same normalized spectrum. Keep the existing one-group output and zero-total-flux behavior. In `get_microxs_and_flux`, a direct override uses the direct component sum over the full energy range divided by the integrated flux. Replace only the selected nuclide/reaction entries with these direct values, once. Preserve other entries from the flux collapse.
5. **Independent depletion and documentation.** Document the per-event convention in the MicroXS and independent-depletion user documentation, including externally supplied tables. Keep existing reaction names, units, group shapes, CSV/HDF5 formats and results keys. Use the current independent operator, matrix assembly and restart interfaces.

## Acceptance tests

| Case | Check |
|---|---|
| Direct, multiple groups and domains | Give a marked parent two distinct partial rates and an unmarked parent the same conventional reaction name. Compare every output group with independently tallied rates divided by its own domain flux. Include a zero-flux group. |
| Flux collapse | For a nonuniform supplied spectrum, compare the one-group MicroXS with the sum of separate component collapses. Check an ordinary pair and an all-zero spectrum. Exercise both standalone and model-based generation, including a model library different from the global default. |
| Direct overrides | Select a subset of reactions and one parent, then repeat with the default override nuclides. Use distinguishable direct and collapsed values to catch a missed override or double counting. |
| Library and selection errors | Require a diagnostic naming the parent, reaction and offending MT for missing LR on a candidate. Check the no-match case and a request excluding the marked pair. Exercise direct generation with executable transport and with an initialized library. |
| Independent consumption and restart | Round-trip a generated MicroXS table, then use it with a known flux and source normalization. Check the parent loss and product gains against the combined rate. Save rates with `integrate(write_rates=True)`, read the conventional column with `Results.get_reaction_rate`, and compare a restart using `prev_results` with an uninterrupted run on the same chain and step schedule. |

Run the focused MicroXS, breakup-rate, independent-operator and restart tests. Retain the existing ordinary-chain cases as compatibility checks. The evaluated-library benchmarks and broader inventory comparisons are covered by [Phase 4 of the overall plan](BREAKUP_DEPLETION_PLAN.md).
