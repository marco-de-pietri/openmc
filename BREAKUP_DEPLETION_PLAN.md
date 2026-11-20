# Neutron inelastic breakup in OpenMC depletion

Generate a depletion chain and a continuous-energy neutron HDF5 library from compatible evaluations. The chain describes the change in inventory. The loaded transport library identifies the partial cross sections that supply the event rate. A user selects the published pair and runs depletion through the ordinary operator APIs.

## Scope exclusions

LR1, LR31, LR39, and LR40 are outside the scope of all phases. LR1 requires
MF6 product identities and multiplicities to interpret general breakup.
LR31 (gamma emission), LR39 (internal conversion), and LR40 (internal pair
formation) describe electromagnetic de-excitation. Generating depletion paths
for these four flags and changing secondary-particle transport for these
processes are outside this work.

## Phase 1 — Generate conventional chain reactions

**Goal.** Represent each supported inventory change using an ordinary reaction name, residual target, and a Boolean breakup marker. For neutron inelastic scattering, support the twelve particle-specific flags LR22–25, LR28–30, and LR32–36. Interpret LR together with the primary MT. Preserve the marker through XML and reduction. Have the C++ chain reader skip conventional MT-to-target registration for marked paths.

**Outcome.** The generated chain describes O16 producing C12 and He4, Li7 producing He4 and H3, and C12 producing three He4 atoms. The detailed object contract, generation rules, and tests are in the [Phase 1 plan](BREAKUP_DEPLETION_PHASE_1.md).

## Phase 2 — Resolve event rates from transport data

**Goal.** At initialization, resolve each participating marked parent and reaction name to the neutron inelastic MTs (51–91) with the corresponding LR in that parent’s loaded evaluation. Expose loaded MT and LR metadata through the C API and Python bindings as needed. Cache the mapping for the lifetime of that library initialization. Distinguish an explicit LR of zero from unavailable metadata. Fail before scoring if a candidate component lacks the metadata needed for classification or if a marked channel has no matches. Name the nuclide and reaction in the diagnostic.

**Outcome.** Direct mode tallies the union of required scores and sums only the assigned component means for each parent. Flux mode collapses those components individually against the flux and sums their rates. Direct overrides in flux mode use the same mapping once per selected pair. Keep operator normalization, material and nuclide indexing, and conventional rate-array ordering. A marker on one parent has no effect on another parent with the same reaction name. See the [Phase 2 implementation and acceptance tests](BREAKUP_DEPLETION_PHASE_2.md).

## Phase 3 — Apply rates and support independent depletion

**Goal.** Feed one combined event rate into the ordinary matrix column for each parent and channel. Reuse existing target and secondary accounting. At rate r, Li7 contributes −r Li7, +r He4 and +r H3. C12 contributes −r C12 and +3r He4.

**Outcome.** Use the Phase 2 resolver in direct and flux-based MicroXS generation, including multigroup flux collapse and direct overrides. Store the combined cross section under the conventional name. Independent depletion with supplied MicroXS consumes those combined values. Document this convention for external cross sections and verify result lookup and restart with the same chain.

See the [Phase 3 implementation and acceptance tests](BREAKUP_DEPLETION_PHASE_3.md).

## Phase 4 — Verify transport and inventory behavior

**Goal.** Check discovery, rates, and inventories independently. Use portable MF3 fixtures for ordinary levels, supported LR flags, missing metadata, and conflicting same-parent representations. Validate the component sets against the original ENDF/B-VIII.1 evaluations. For Li7 check MT52–82 and for C12 check MT52–62 plus MT91. Exercise O16 MT56/LR22 as a fixture; use the evaluated component list when testing a real library.

**Outcome.** Compare direct rates to separately tallied component sums and collapsed rates to component-wise flux integrals. Test two parents sharing a reaction name with different rate sources. Verify analytic constant-rate inventories, independent depletion, and restart. At fixed composition, verify that ordinary component and particle-production tallies retain their values. Run affected regression tests with the project’s prescribed data, strict floating-point build, and MPI options as applicable.

See the [Phase 4 verification plan and acceptance checks](BREAKUP_DEPLETION_PHASE_4.md).

## Phase 5 — Generate and publish a compatible data pair

**Goal.** Use the supported OpenMC chain generator in openmc-dev/data and regenerate neutron HDF5 files through an ENDF-aware processing path that preserves LR. Record evaluation releases and processing versions. Inspect LR attributes in the distributed files and validate full and reduced chains against those files.

**Outcome.** Publish tested chain/library combinations with the minimum reader version. Include an end-to-end depletion check using the distributed artifacts. General MF6 energy-dependent yields and other breakup modes require separate physics and rate-accounting work.

## References

- [ENDF-6 formats manual, Tables 9 and 10](https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf#page=38)
- [OpenMC PR 4146, preservation of evaluated reaction metadata](https://github.com/openmc-dev/openmc/pull/4146)
- [ENDF/B-VIII.1 chain generator](https://github.com/openmc-dev/data/blob/master/depletion/generate_endf81_chain.py)
