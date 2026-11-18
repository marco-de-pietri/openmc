# Phase 1 — Generate depletion chain entries for neutron inelastic breakup

**Goal.** Generate one conventional depletion path per representable breakup channel. The path describes the atoms produced per event. Its marker tells a rate consumer to discover the relevant inelastic components from the paired neutron HDF5 library. Phase 1 delivers the generated Python objects, XML representation, and safe loading through both chain readers.

## Scope exclusions

LR1, LR31, LR39, and LR40 are outside the scope of all phases. LR1 requires
MF6 product identities and multiplicities to interpret general breakup.
LR31 (gamma emission), LR39 (internal conversion), and LR40 (internal pair
formation) describe electromagnetic de-excitation. Generating depletion paths
for these four flags and changing secondary-particle transport for these
processes are outside this work.

## Supported neutron inelastic modes

Apply this mapping to MF3 sections with MT51–91. The offsets are relative to the parent. Outgoing neutron counts include the primary scattered neutron. Inventory secondaries are in addition to any residual target.

| LR | OpenMC reaction | Outgoing neutrons | Inventory secondaries | ΔA | ΔZ |
|---:|---|---:|---|---:|---:|
| 22 | `(n,na)` | 1 | He4 | −4 | −2 |
| 23 | `(n,n3a)` | 1 | 3 He4 | −12 | −6 |
| 24 | `(n,2na)` | 2 | He4 | −5 | −2 |
| 25 | `(n,3na)` | 3 | He4 | −6 | −2 |
| 28 | `(n,np)` | 1 | H1 | −1 | −1 |
| 29 | `(n,n2a)` | 1 | 2 He4 | −8 | −4 |
| 30 | `(n,2n2a)` | 2 | 2 He4 | −9 | −4 |
| 32 | `(n,nd)` | 1 | H2 | −2 | −1 |
| 33 | `(n,nt)` | 1 | H3 | −3 | −1 |
| 34 | `(n,n3He)` | 1 | He3 | −3 | −2 |
| 35 | `(n,nd2a)` | 1 | H2 + 2 He4 | −10 | −5 |
| 36 | `(n,nt2a)` | 1 | H3 + 2 He4 | −11 | −5 |

## Generated objects and XML

Keep the four fields of `ReactionTuple(type, target, Q, branching_ratio)`. Add a parent-local `Nuclide.breakup_reactions` list of marked reaction names. For the following paths the branch ratio is one and Q is the evaluated common MF3 mass-difference Q in eV. The table shows the breakup contribution to each nuclide’s objects.

| Parent and evaluated mode | Tuple type | Tuple target | Entry in breakup_reactions | Atoms produced per event |
| --- | --- | --- | --- | --- |
| O16, neutron inelastic with LR22 | (n,na) | C12 | (n,na) | One C12 target and one He4 secondary |
| Li7, neutron inelastic with LR33 | (n,nt) | He4 | (n,nt) | One He4 target and one H3 secondary |
| C12, neutron inelastic with LR23 | (n,n3a) | None | (n,n3a) | Three He4 secondaries |

The corresponding reaction elements below omit evaluated Q values for illustration. Each belongs under its own parent nuclide.

```xml
<reaction type="(n,na)" target="C12" breakup="true"/>
<reaction type="(n,nt)" target="He4" breakup="true"/>
<reaction type="(n,n3a)" breakup="true"/>
```

The existing global reaction table supplies secondary species and multiplicities. The chain-wide reaction list continues to contain distinct conventional names. The marker belongs to a parent because another nuclide can use the same name for a conventional transport MT. Source MT lists and their LR values reside in the paired transport library.

## Implementation by function

1. **Read and interpret MF3 in Chain.from_endf().** Extend the local parent/MT record from QI alone to QM, QI, and LR. Reuse the MF3 CONT reader. Keep QI handling for ordinary paths. In a small neutron-reaction interpretation helper, recognize the twelve LR values above on MT51–91. Use the same LR-to-name mapping for generation and marker validation. Compute the residual by adding the existing openmc.data.DADZ offsets to the parent mass and charge. These offsets already account for all outgoing neutrons and charged fragments. Use REACTIONS[name].secondaries for inventory production, preserving repeated entries. Convert a physical residual to its nuclide name. Zero remaining mass and charge gives a null target. Reject inconsistent balances and diagnose missing inventory products in the full generated chain.
2. **Select and group paths in Chain.from_endf().** Use the requested conventional reaction names for both ordinary and breakup paths. The existing default reaction list and any explicit list follow the same rule. An empty list selects no breakup paths. Group only exclusive components with the same parent, type, target, products, and QM. Exclude aggregate MT4. Reject incompatible Q values, competing ordinary and marked paths of the same parent/type, or an ambiguous aggregate/component representation. Report the parent and conflicting MTs. A single rate slot and marker cannot distinguish these cases.
3. **Construct and persist the marker.** Initialize the set in Nuclide.__init__(). Give Nuclide.add_reaction() an optional positional or keyword breakup flag defaulting to false. Register the tuple and marker together. In Nuclide.to_xml_element() write breakup="true" on marked entries; in Nuclide.from_xml() treat an absent marker as false. Chain.add_nuclide() continues to register names and invalidate its cache. Exercise Chain.export_to_xml() and Chain.from_xml() as public entry points.
4. **Validate and reduce.** In Nuclide.validate() check that every marked name denotes exactly one supported, unit-branch path on that parent. Reject mixed marked and ordinary entries during insertion or XML loading before path-level information is lost. Chain._follow() can use the ordinary target and known secondaries. In Chain.reduce() carry the marker onto the copied path, including at a depth cutoff where the target is dropped to preserve parent destruction. Construct an independent marker list for each copied nuclide.
5. **Respect the marker in the C++ reader.** In ChainNuclide::ChainNuclide() in src/chain.cpp inspect the marker before mapping the name to an MT for the transport product table. Skip marked entries in this table, including those with a target. For O16, mapping the name alone would incorrectly register C12 under MT22. Runtime component discovery for rate calculations belongs to Phase 2.

## Acceptance tests

Generate each example through Chain.from_endf() from portable MF3 fixtures including unrelated inelastic levels. Assert the tuples, markers, Q values, path counts, product multiplicities, default selection, and explicit selection. Cover mixed parents sharing a reaction name, complete breakup, differing QM, competing representations, invalid balances, and legacy XML. Export and reload the chain, then reduce it at depths zero and one and check the retained paths, markers, and reachable products. Load the target-bearing O16 example through C++ and verify that it adds no MT22 product mapping while an ordinary unmarked reaction still registers normally.


Parameterize portable MF3 tests over all twelve LR values with independently specified residuals and inventory multiplicities. For example, synthetic O16 fixtures for LR24, LR25, and LR30 must yield C11 + He4, C10 + He4, and Be7 + 2 He4, respectively. Check the depletion matrix column at unit rate for one parent loss and the expected atom counts. Include a mixed fixture with multiple supported modes on one parent, an ordinary inelastic level, and excluded flags. Verify that explicit selection retains only the requested breakup modes.
