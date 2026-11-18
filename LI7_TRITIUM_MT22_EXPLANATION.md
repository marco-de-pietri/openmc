#pjkWhy Li-7 tritium production can be represented as MT22

## What the ENDF/B-VIII.1 file actually says

The Li-7 ENDF file `n-003_Li_007.endf` does not contain a separate MF=3,
MT=22 cross-section table. Instead, it divides the tritium-producing reaction
among the discrete inelastic channels MT=52 through MT=82.

The general description in MF=1, MT=451 states:

```text
(3) (N,NT) X/S (SUM OF MT52 THRU MT82)
```

The more detailed description states:

```text
MT=52-82 (N,NPRIME)ALPHA-T CROSS SECTION.  EXCEPT FOR MT56,
        EXCITATION ENERGY BINS WERE CONSTRUCTED SUCH THAT CALCU-
        LATED NEUTRON EMISSION SPECTRA AGREE WITH THE MEASURE-
        MENTS BY LI80 AT 5.96, 9.83, AND 14.1 MEV INCIDENT ENER-
        GY. INTERPOLATION OF CROSS SECTIONS BETWEEN ENERGIES
        BASED ON LINEAR INTERPOLATION OF COMPUTED CM EMISSION
        SPECTRA, SUITABLY CONTRACTED OR EXPANDED TO CONSERVE
        ENERGY.  ABOVE 14.1 MEV, THE EXP EMISSION SHAPE WAS USED
        TO INFER THE X/S, WITH THE SHAPES AGAIN SUITABLY EX-
        PANDED TO CONSERVE ENERGY.
        AT ALL INCIDENT ENERGIES, THE SUM OF MT=52-82 WAS
        CONSTRAINED TO EQUAL THE (N,NT) CROSS SECTION FROM THE
        VARIANCE-COVARIANCE ANALYSIS.
```

Therefore, the statement made explicitly by this evaluation is

$$
\sigma_{(n,nt)}(E) = \sum_{i=52}^{82}\sigma_i(E).
$$

Here, $E$ is the incident-neutron energy. The equality applies at every
incident energy, rather than only after averaging over an energy group or a
particular neutron spectrum.

ENDF-6 encodes the individual MT numbers in the section records and lists the
sections individually by MF and MT in the directory at the end of MF=1,
MT=451. A parser can identify the tritium-producing breakup bins by selecting
the lithium-7 inelastic MF=3 sections with `LR=33`. In this evaluation, that
selection gives MT=52 through MT=82. The descriptive text in MF=1, MT=451
explicitly confirms that their cross sections sum to the total (n,nt) cross
section.

MT=51 must not be included in this sum. MT=51 describes ordinary inelastic
scattering to the first bound excited state of Li-7, followed by gamma
emission. It does not describe the alpha-plus-tritium breakup. The tritium
production sum begins at MT=52, exactly as specified by the evaluation.

## What the reaction means physically

The combined MT=52--82 channels describe the final reaction

$$
{}^7\mathrm{Li} + n \rightarrow n + \alpha + {}^3\mathrm{H}.
$$

In simple terms:

1. A neutron strikes Li-7.
2. A neutron is present in the outgoing products.
3. The remaining Li-7 system breaks into an alpha particle, He-4, and a
   triton, H-3.
4. Each MT number from 52 through 82 represents a different discrete
   excitation-energy bin used by the evaluation.
5. Adding all those bins gives the total tritium-production cross section for
   this reaction mechanism.

The many MT numbers do not represent many different final radioactive
products of interest for depletion. They are a detailed way of describing how
the outgoing neutron energy and the breakup process are distributed. For an
inventory calculation concerned with the amount of H-3 produced, their sum is
the relevant quantity.

## What data are stored under MT=52 through MT=82

These MT sections contain more than 31 unrelated cross-section curves. They
form a set of **pseudo-levels**, or excitation-energy bins, used to describe
the continuum breakup of Li-7 into an alpha particle and a triton while a
neutron leaves the collision.

Ordinarily, MT=51, 52, 53, and so on describe inelastic neutron scattering to
particular levels of the target nucleus. In this evaluation, MT=52--82 are
used as bins for the `(n,n')alpha-t` continuum. Each MT identifies a different
amount of excitation energy and therefore a different outgoing-neutron energy
range. The sum over all bins is the total tritium-production cross section.

The [ENDF-6 format](https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf)
organizes records as **material (MAT) → file (MF) → section (MT) → records**.
Here, the material is the lithium-7 target. Within that material, each logical
ENDF file groups a different kind of data, and each section in a file has an MT
number. The same MT can appear under several MFs.

In this evaluation, each of the three files below contains a section for every
MT from 52 through 82. For example, MF=3/MT=52 gives the cross section for
the MT=52 bin, while MF=4/MT=52 gives the neutron angular distribution for
that same bin.

- **MF=3:** the reaction cross section as a function of incident-neutron
  energy;
- **MF=4:** the angular distribution of the emitted neutron;
- **MF=33:** covariance information associated with the cross section.

There are no MF=5 or MF=6 sections for MT=52--82 in this file. Consequently,
there is no independent continuous outgoing-energy distribution stored in
MF=5 or correlated energy-angle distribution stored in MF=6. The pseudo-level
energy assigned to each MT, together with the reaction kinematics and MF=4
angular data, represents the outgoing neutron.

### MF=3: cross-section data

Every MF=3 section contains a tabulated pair of arrays:

$$
\{E_j,\sigma_i(E_j)\},
$$

where $i$ is the MT number, $E_j$ is an incident-neutron energy in eV, and
$\sigma_i$ is the partial cross section in barns. The interpolation flag is
linear-linear, so the cross section is linearly interpolated between adjacent
tabulated points.

The ENDF fields distinguish two energy quantities:

- `QM = -2.4667 MeV` is common to MT=52--82. It is the overall mass-difference
  Q value used for the alpha-plus-tritium breakup channel.
- `QI` is different for every MT. Its negative value defines the effective
  pseudo-level or excitation-bin energy shown in the table below.

Each of these MF=3 sections carries the `LR=33` residual-breakup flag.
The MT identifies the inelastic neutron-scattering bin, and `LR=33` indicates
that the excited lithium-7 nucleus then emits a triton, leaving helium-4.

The following table inventories the numerical MF=3 content. “Level/bin” is
`-QI`; “threshold” is the first tabulated incident energy, where the partial
cross section is zero; “points” is the number of tabulated energy/cross-section
pairs. The maximum is included only as a convenient description of the curve,
not as another ENDF field.

| MT | Level/bin (MeV) | Threshold (MeV) | MF=3 points | Maximum cross section (b) |
|---:|---:|---:|---:|---:|
| 52 | 2.75 | 3.1454 | 95 | 0.0209333 |
| 53 | 3.25 | 3.7172 | 88 | 0.0527220 |
| 54 | 3.75 | 4.2891 | 82 | 0.0634225 |
| 55 | 4.25 | 4.8610 | 75 | 0.0749295 |
| 56 | 4.63 | 5.2956 | 73 | 0.1752000 |
| 57 | 4.75 | 5.4329 | 68 | 0.0595125 |
| 58 | 5.25 | 6.0048 | 61 | 0.0449164 |
| 59 | 5.75 | 6.5767 | 58 | 0.0384606 |
| 60 | 6.25 | 7.1485 | 55 | 0.0336449 |
| 61 | 6.75 | 7.7204 | 52 | 0.0295095 |
| 62 | 7.25 | 8.2923 | 48 | 0.0263222 |
| 63 | 7.75 | 8.8642 | 45 | 0.0240611 |
| 64 | 8.25 | 9.4361 | 42 | 0.0225851 |
| 65 | 8.75 | 10.008 | 38 | 0.0217271 |
| 66 | 9.25 | 10.580 | 36 | 0.0207316 |
| 67 | 9.75 | 11.152 | 34 | 0.0203447 |
| 68 | 10.25 | 11.724 | 32 | 0.0198962 |
| 69 | 10.75 | 12.295 | 30 | 0.0197100 |
| 70 | 11.25 | 12.867 | 28 | 0.0193212 |
| 71 | 11.75 | 13.439 | 26 | 0.0180677 |
| 72 | 12.25 | 14.011 | 23 | 0.0166212 |
| 73 | 12.75 | 14.583 | 21 | 0.0151414 |
| 74 | 13.25 | 15.155 | 19 | 0.0136963 |
| 75 | 13.75 | 15.727 | 17 | 0.0123440 |
| 76 | 14.25 | 16.299 | 15 | 0.0111268 |
| 77 | 14.75 | 16.871 | 13 | 0.0100950 |
| 78 | 15.25 | 17.442 | 11 | 0.00917586 |
| 79 | 15.75 | 18.014 | 8 | 0.00835206 |
| 80 | 16.25 | 18.586 | 6 | 0.00761483 |
| 81 | 16.75 | 19.158 | 4 | 0.00597473 |
| 82 | 17.25 | 19.730 | 2 | 0.00191626 |

Most bins are spaced by 0.5 MeV. MT=56, at 4.63 MeV, is the exception because
it represents the known 4.63 MeV Li-7 level and was evaluated separately.
The file documentation explicitly says that MT=56 used a separate
variance-covariance analysis. MT=52--55 and MT=57--82 provide the surrounding
constructed continuum bins.

The threshold is higher than the listed level/bin energy because incident
energy must also provide recoil energy. As the bin energy rises, the threshold
rises. The highest bins contain fewer tabulated points because the evaluation
ends at 20 MeV; for example, MT=82 begins at 19.73 MeV and therefore contains
only its zero-at-threshold point and its value at 20 MeV.

### Why the particle-specific LR codes are legacy

ENDF-6 retains the particle-specific `LR` codes for compatibility with
ENDF/B-V. These codes use a fixed lookup table. For example, `LR=33`
specifies that the residual nucleus emits a triton after the primary
particle emission. The `LR` field remains part of ENDF-6.
[ENDF-6 Formats Manual, Section 0.4 and Appendix B.2](https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf)

For inelastic neutron scattering, adding the primary scattered neutron to
the particles specified by LR gives the reaction designated by the same MT
number. [ENDF-6 Formats Manual, Section 0.4](https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf)

| Code | Additional particles specified by LR | Reaction designated by the same MT |
|---:|---|---|
| 22 | α | (n,nα) |
| 23 | 3α | (n,n3α) |
| 24 | n + α | (n,2nα) |
| 25 | 2n + α | (n,3nα) |
| 28 | p | (n,np) |
| 29 | 2α | (n,n2α) |
| 30 | n + 2α | (n,2n2α) |
| 32 | d | (n,nd) |
| 33 | t | (n,nt) |
| 34 | ³He | (n,n³He) |
| 35 | d + 2α | (n,nd2α) |
| 36 | t + 2α | (n,nt2α) |

Appendix B.2 also defines residual de-excitation flags outside this
LR-to-MT correspondence. `LR=31` means gamma-only decay, `LR=39` means
internal conversion, and `LR=40` means electron–positron pair formation
during de-excitation. With `MT=51, LR=40`, for example, MT=51 identifies
inelastic scattering to the first excited state and LR=40 specifies how
that state decays. [ENDF-6 Formats Manual, Section 0.4 and Appendix B.2](https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf)

ENDF-6 also provides a general representation for breakup reactions.
`LR=1` indicates that MF=6 identifies the products and their multiplicities
explicitly. MF=6 can also describe their energy and angular distributions.
[ENDF-6 Formats Manual, Chapter 6](https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf)

In this lithium-7 evaluation, a parser can use `LR=33` in MF=3 to identify
the triton-producing breakup bins, MT=52 through MT=82.

### MF=4: outgoing-neutron angular data

MF=4 describes the direction of the emitted neutron relative to the incident
neutron. For MT=52--71, the file stores Legendre coefficients at selected
incident energies. These coefficients describe angular anisotropy; a zero set
corresponds to isotropic emission, while nonzero coefficients bias emission
toward or away from the incident-neutron direction.

The evaluation comments say:

```text
MT=52-71  LEGENDRE COEFFICIENTS OBTAINED FROM ANALYSIS OF LI80
          EXPERIMENT DESCRIBED IN MF3 ABOVE (EXCEPT MT56). COEF-
          FICIENTS FOR MT=56 OBTAINED BY DRAWING SMOOTH CURVES
          THROUGH FITTED VALUES FROM EXPERIMENTS LISTED ABOVE
          UNDER MF3/MT56, ESPECIALLY HO68 AND HO79.
MT=72,82  ISOTROPY ASSUMED.
```

Here `MT=72,82` means the inclusive range MT=72 through MT=82. Thus:

- MT=52--55 and MT=57--71 use angular shapes derived mainly from the Li80
  neutron-emission measurements;
- MT=56 uses separately fitted angular data associated with the 4.63 MeV
  level;
- MT=72--82 assume isotropic neutron emission.

For example, MT=52 contains angular distributions at its threshold and at
5.96, 9.83, 14.1, and 20 MeV. The threshold coefficients are zero. The
nonzero Legendre coefficients at the other energies describe the measured or
inferred angular shape. Similar data are supplied for the other anisotropic
bins. MT=56 has a denser set of incident-energy points because it was treated
separately.

### MF=33: covariance data

MF=33 describes uncertainty and energy-to-energy correlation in the partial
cross sections. The file contains an MF=33 section for every MT from 52
through 82, but the evaluator describes them as **lump indicators only**.
They are grouped into covariance lumps rather than supplied as 31 fully
independent experimental covariance evaluations:

- MT=52--55 are combined in lump MT=853;
- MT=56 has its own lump MT=854;
- MT=57--61 are combined in lump MT=855;
- MT=62--66 are combined in lump MT=856;
- MT=67--71 are combined in lump MT=857;
- MT=72--76 are combined in lump MT=858;
- MT=77--82 are combined in lump MT=859.

These lumped partial-reaction covariances were adjusted so their combination
is consistent with the covariance of the total `(n,nt)` reaction. This matters
for uncertainty propagation, but it does not alter the central cross-section
identity used here:

$$
\sigma_{(n,nt)}(E)=\sum_{i=52}^{82}\sigma_i(E).
$$

The evaluator adds an important qualification: the MT=852--859 covariance
matrices sum only approximately to the MT=4 covariance matrix. The file
recommends using these lumps only to calculate neutron-emission-spectrum
uncertainties. They should not be treated as a precise covariance matrix for
the aggregate MT22 cross section without additional covariance processing.

### What is retained when only the MT22 sum is needed

For `MicroXS.from_multigroup_flux`, the necessary quantity is the incident-
energy-dependent MF=3 sum. The operation is performed pointwise before flux
collapse:

$$
\sigma_{22}(E)=\sigma_{52}(E)+\sigma_{53}(E)+\cdots+\sigma_{82}(E).
$$

The MF=4 angular distributions are unnecessary for this use because no
particles are being transported and `MicroXS.from_multigroup_flux` only
collapses reaction cross sections. The MF=33 covariance data are also not used
by that method. Consequently, creating an aggregate MT22 cross-section entry
for this purpose retains the required tritium-production rate while omitting
transport-only angular detail and uncertainty information from the aggregate
record.

## Why MT22 is used for the combined cross section

In ENDF terminology, MT=22 denotes the `(n,na)` reaction. The letters after
the comma describe the emitted light particles: an outgoing neutron (`n`) and
an alpha particle (`a`). Applying mass-number and charge conservation to Li-7
leaves H-3 as the residual nucleus:

$$
\begin{aligned}
A &: 7 + 1 = 1 + 4 + 3,\\
Z &: 3 + 0 = 0 + 2 + 1.
\end{aligned}
$$

Thus, for Li-7, MT22 corresponds to

$$
{}^7\mathrm{Li}(n,n\alpha){}^3\mathrm{H}.
$$

The effective MT22 cross section can therefore be constructed as

$$
\boxed{\sigma_{22}(E) \equiv
       \sum_{i=52}^{82}\sigma_i(E)}.
$$

The symbol $\equiv$ is intentional. It means that MT22 is being defined as
an aggregate representation for the intended OpenMC calculation. The original
ENDF/B-VIII.1 file does not itself contain an MT22 record or explicitly print
the equation `sigma_22 = sum(MT52 ... MT82)`.

The justification has two parts:

- The ENDF evaluation explicitly says that MT=52--82 sum to its total
  `(n,nt)` cross section.
- The physical final state of that total reaction is neutron + alpha + H-3,
  which is represented by MT22 `(n,na)` in ENDF reaction notation.

Therefore, assigning the sum to MT22 is a conversion of the evaluated partial
channels into a single tritium-production channel. It is not a quotation of an
existing MT22 value from `n-003_Li_007.endf`.

## Relationship to `Li7g.asc`

The EAF file `Li7g.asc` includes a separate entry:

```text
3   22   N,NA        JENDL-3.2/A; ADD   [12]
```

This shows that the EAF library chose MT22 as the representation of the
Li-7 `(n,na)` channel. It supports using the same reaction identifier for a
single tritium-production channel in OpenMC.

It does not mean that the numerical MT22 values in `Li7g.asc` are necessarily
identical to the sum of MT=52--82 in `n-003_Li_007.endf`. The two files come
from different evaluated-data sources and processing histories. If the goal is
to retain the ENDF/B-VIII.1 cross section, the numerical data should be formed
by summing MT=52--82 from `n-003_Li_007.endf`; `Li7g.asc` supplies an example
of the desired MT22 organization.

## Meaning for OpenMC `MicroXS.from_multigroup_flux`

`openmc.deplete.MicroXS.from_multigroup_flux` collapses the cross section for
the requested reaction MT. It does not know that the Li-7 MT=52--82 channels
should be combined into an H-3 production rate. If `(n,na)` is requested while
the OpenMC HDF5 file has no MT22 reaction, the resulting microscopic cross
section is zero.

Adding an aggregate MT22 reaction to the processed Li-7 HDF5 file gives
OpenMC the required mapping:

$$
\text{OpenMC `(n,na)`}\longrightarrow \text{MT22}
\longrightarrow \sum_{i=52}^{82}\sigma_i(E).
$$

The added reaction should be marked as redundant in the HDF5 data. This tells
OpenMC that MT22 is a sum of reactions already present in the file and prevents
it from being added again to the total transport cross section. The aggregate
remains available to `collapse_rate` and consequently to
`MicroXS.from_multigroup_flux`.

For a depletion or inventory calculation, the depletion-chain entry should
then map Li-7 `(n,na)` to H-3. In that setup, OpenMC uses the collapsed MT22
cross section as the Li-7 destruction and H-3 production rate.

## Precise conclusion

The defensible statement is:

> The ENDF/B-VIII.1 Li-7 evaluation explicitly defines its total `(n,nt)`
> cross section as the sum of MT=52 through MT=82. Because these channels have
> the physical final state neutron + alpha + H-3, that sum can be represented
> as an aggregate MT22 `(n,na)` cross section for OpenMC inventory calculations.

The statement that should be avoided is:

> The ENDF/B-VIII.1 file explicitly provides MT22 as the sum of MT=52 through
> MT=82.

That second statement is too strong: the file supplies the partial channels
and their stated `(n,nt)` sum, while the MT22 record is added during the
conversion for OpenMC.
