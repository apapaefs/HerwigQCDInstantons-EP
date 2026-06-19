# Herwig QCD Instantons

This repository contains a Herwig plugin for generating QCD-instanton-like
multi-parton final states,

```text
g g -> quark-antiquark pairs + extra gluons
```

together with a many-body phase-space generator and a small Rivet analysis for
particle-level shape checks. The default run cards are set up for the
Khoze-Krauss-Schott (KKS) model of arXiv:1911.09726.

The implementation is deliberately phenomenological. It does not calculate a
full instanton matrix element event by event. Instead, it enumerates allowed
partonic final states, generates flat high-multiplicity phase space, and assigns
weights from either a simple multiplicity model or a tabulated KKS-inspired
partonic cross section.

## Components

- `MEInstanton.{h,cc}`: Herwig `BlobME` plugin. Defines the process list,
  matrix-element weight, quark/gluon multiplicity choices, factorization-scale
  option, and colour connections.
- `Phasespace/MamboPhasespace.{h,cc}`: many-body phase-space generator used by
  `MEInstanton`.
- `LHC-Instanton*.in`: example Herwig input cards.
- `Rivet/QCD_INSTANTON_KKS.cc`: particle-level Rivet analysis inspired by the
  observables in KKS Section 4.

## Installing `MamboPhasespace` into a Herwig source tree

Let `HERWIG_SRC` be the top-level Herwig source directory, for example the
directory that contains `configure` and `MatrixElement/`.

1. Copy the phase-space source files:

```sh
cp Phasespace/MamboPhasespace.cc "$HERWIG_SRC/MatrixElement/Matchbox/Phasespace/"
cp Phasespace/MamboPhasespace.h  "$HERWIG_SRC/MatrixElement/Matchbox/Phasespace/"
```

2. Edit `$HERWIG_SRC/MatrixElement/Matchbox/Phasespace/Makefile.am`:

Add `MamboPhasespace.h` to `ALL_H_FILES` and `MamboPhasespace.cc` to
`ALL_CC_FILES`.

For an already-configured build tree, either rerun the autotools/configure step,
or make the same additions in the generated `Makefile.in` and `Makefile`.

3. Rebuild and install Herwig:

```sh
cd "$HERWIG_SRC"
make -j"$(sysctl -n hw.ncpu 2>/dev/null || nproc)"
make install
```

4. Check that Herwig can instantiate the new class:

```sh
cat >/tmp/mambo-phasespace-smoke.in <<'EOF'
cd /Herwig/MatrixElements/Matchbox/Phasespace
create Herwig::MamboPhasespace MamboPS
EOF

Herwig read /tmp/mambo-phasespace-smoke.in
```

The final command should exit without an error.

## Building `MEInstanton` as a Herwig contrib plugin

After `MamboPhasespace` has been added to Herwig and Herwig has been rebuilt,
copy this repository into Herwig's `Contrib` directory:

```sh
cp -R /path/to/HerwigQCDInstantons "$HERWIG_SRC/Contrib/HerwigQCDInstantons"
cd "$HERWIG_SRC/Contrib"
bash make_makefiles.sh
```

Build and install the plugin:

```sh
cd "$HERWIG_SRC/Contrib/HerwigQCDInstantons"
make
make install
```

This installs `Instantons.so` into the Herwig plugin directory.

Check that Herwig can load both pieces:

```sh
cat >/tmp/instanton-plugin-smoke.in <<'EOF'
cd /Herwig/MatrixElements/Matchbox/Phasespace
create Herwig::MamboPhasespace MamboPS
cd /Herwig/MatrixElements
create Herwig::MEInstanton MEInstanton Instantons.so
set MEInstanton:Phasespace /Herwig/MatrixElements/Matchbox/Phasespace/MamboPS
EOF

Herwig read /tmp/instanton-plugin-smoke.in
```

## Minimal Herwig setup

The instanton plugin needs a phase-space object and must be inserted as the
matrix element for the `SubProcess` object:

```text
cd /Herwig/MatrixElements/Matchbox/Phasespace
create Herwig::MamboPhasespace MamboPS
set MamboPS:CouplingData PhasespaceCouplings

cd /Herwig/MatrixElements
create Herwig::MEInstanton MEInstanton Instantons.so
insert SubProcess:MatrixElements[0] MEInstanton
set MEInstanton:Phasespace /Herwig/MatrixElements/Matchbox/Phasespace/MamboPS
```

The example cards then configure PDFs, sampler settings, mass cuts, the shower,
and optional Rivet output.

## Process Content

`MEInstanton` registers gluon-gluon initiated final states:

```text
g g -> q qbar pairs + N extra gluons
```

The quark pairs are one pair per active flavour, ordered as
`d dbar`, `u ubar`, `s sbar`, `c cbar`, and optionally `b bbar`.

The number of extra gluons is controlled by Herwig's inherited
`NAdditional` interface. The code copies this value into `ngluon_max` at
initialization and registers processes with 0 through `NAdditional` extra
gluons.

For example:

```text
set MEInstanton:NAdditional 10
```

allows final states with 0, 1, ..., 10 extra gluons. You can raise or lower this
without changing the C++ code, but the phase-space dimension grows as
`3*nFinal - 4`, so larger values can make integration much slower.

## Matrix-Element Models

The main model switch is:

```text
set MEInstanton:MEModeling PureMultiplicity
set MEInstanton:MEModeling KKS
```

### `PureMultiplicity`

This is a toy model. The base matrix element is flat and only the relative
weights of different extra-gluon multiplicities are changed. The available
multiplicity parametrisations are described below.

### `KKS`

This is the default in the supplied run cards. It uses hardcoded interpolation
tables as functions of `sqrt(shat)` for:

- `1/rho`, used for the optional factorization scale;
- `alpha_s(1/rho)`, currently diagnostic only;
- the mean number of gluons;
- `sigmahat`, the tabulated partonic instanton cross section.

In KKS mode, the event weight is built from:

1. a PDF reweighting from the default scale `shat` to the selected
   factorization scale;
2. a Poisson weight for the actual extra-gluon multiplicity, using the
   interpolated mean number of gluons;
3. the interpolated `sigmahat`, converted from pb to `GeV^-2`;
4. a division by the generated phase-space Jacobian;
5. an overall `2*shat` factor used by the implementation.

The interpolation table covers approximately:

```text
10.7 GeV <= sqrt(shat) <= 2895.5 GeV
```

The code returns zero above the upper edge. The example cards set an explicit
`MHatMax` of `2895.5*GeV`.

## Gluon Multiplicity Parameters

These settings affect `PureMultiplicity` mode:

```text
set MEInstanton:MultiplicityParametrisation Poisson
set MEInstanton:MultiplicityParametrisation Gaussian
set MEInstanton:MultiplicityParametrisation Flat
set MEInstanton:MultiplicityParametrisation UserDefined
```

- `Poisson`: multiplies the flat weight by
  `PoissonMean^ngluon * exp(-PoissonMean) / ngluon!`.
- `Gaussian`: multiplies the flat weight by
  `exp(-(ngluon - GaussianParamA)^2/GaussianParamB)
  / sqrt(pi*GaussianParamB)`.
- `Flat`: leaves all allowed gluon multiplicities equally weighted by the
  matrix element.
- `UserDefined`: currently has no additional function implemented, so it is
  equivalent to `Flat` unless code is added.

The numerical parameters are:

```text
set MEInstanton:PoissonMean 3.0
set MEInstanton:GaussianParamA 5.0
set MEInstanton:GaussianParamB 200.0
```

In `KKS` mode, these toy multiplicity parameters are ignored. KKS mode always
uses the interpolated KKS mean gluon multiplicity.

Important: the allowed gluon range is truncated by `NAdditional`, and the
Poisson distribution is not renormalized after truncation. If `NAdditional` is
too small, the high-multiplicity tail is removed rather than redistributed.

## Quark-Pair Options

The quark-pair mode is controlled by:

```text
set MEInstanton:QuarkPairs Fixed
set MEInstanton:QuarkPairs Variable
```

### `Fixed`

Uses `NQuarkPair` quark-antiquark pairs:

```text
set MEInstanton:NQuarkPair 4
```

The practical range is 1 to 5, corresponding to the active flavours
`d,u,s,c,b`. The supplied cards use 4 as the base value.

### `Variable`

Registers both 4-pair and 5-pair final states. In KKS mode the matrix-element
weight is multiplied by `0.5` for this option, so that the two quark-pair
choices are averaged rather than simply summed.

The example cards use:

```text
set MEInstanton:NQuarkPair 4
set MEInstanton:QuarkPairs Variable
```

## Factorization-Scale Options

This switch is relevant for KKS mode:

```text
set MEInstanton:FactorizationScale InvRho
set MEInstanton:FactorizationScale sHat
```

- `InvRho`: uses `(1/rho)^2` from the KKS interpolation table.
- `sHat`: uses the partonic invariant mass squared.

The implementation first evaluates the gluon PDFs at the default scale `shat`,
then reweights by the ratio of PDFs at the selected factorization scale to PDFs
at `shat`.

The renormalization-like scale returned by `MEInstanton::scale()` is still
`sHat()`. The `FactorizationScale` switch affects the PDF reweighting in KKS
mode.

## Colour-Connection Options

The colour model is selected with:

```text
set MEInstanton:ColourConnections Simple
set MEInstanton:ColourConnections Random
set MEInstanton:ColourConnections Random2
set MEInstanton:ColourConnections Random3
```

- `Simple`: connects the incoming gluons as a singlet, connects most quark
  pairs directly, and pairs final-state gluons in a simple deterministic
  pattern. If the number of extra gluons is odd, the last quark pair is
  connected through one gluon.
- `Random`: randomizes final-state colour connections while keeping the
  incoming gluons as a singlet.
- `Random2`: randomizes colour connections with one incoming gluon colour line
  included in the map, allowing more initial-final colour structure.
- `Random3`: randomizes colour connections with both incoming gluon colour
  lines included. This is the most permissive option and is used in the example
  KKS cards.

These are phenomenological colour assignments for showering and hadronization,
not exact instanton colour amplitudes.

## Typical KKS Settings

A compact KKS setup looks like:

```text
set MEInstanton:NQuarkPair 4
set MEInstanton:QuarkPairs Variable
set MEInstanton:NAdditional 10

set MEInstanton:MEModeling KKS
set MEInstanton:FactorizationScale InvRho
set MEInstanton:ColourConnections Random3

create ThePEG::Cuts /Herwig/Cuts/ResetCuts
set /Herwig/EventHandlers/EventHandler:Cuts /Herwig/Cuts/ResetCuts
set /Herwig/Cuts/ResetCuts:MHatMin 30*GeV
set /Herwig/Cuts/ResetCuts:MHatMax 2895.5*GeV
```

Use a larger `NAdditional` if you want more of the high-gluon tail at large
`sqrt(shat)`, and use a lower value for faster pilot studies.

## Typical Toy-Multiplicity Settings

For a quick non-KKS multiplicity study:

```text
set MEInstanton:NQuarkPair 4
set MEInstanton:QuarkPairs Fixed
set MEInstanton:NAdditional 8

set MEInstanton:MEModeling PureMultiplicity
set MEInstanton:MultiplicityParametrisation Poisson
set MEInstanton:PoissonMean 4.0
set MEInstanton:ColourConnections Simple
```

This generates a flat matrix element with the selected relative distribution of
extra gluons.

## PDFs, Cuts, and Sampler Settings

The example cards use:

```text
set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 13000.0
set /Herwig/Partons/LHAPDF:PDFName NNPDF31_nnlo_as_0118
set /Herwig/Partons/LHAPDF:RangeException Freeze
```

and the flat bin sampler with `AlmostUnweighted Yes`. The `MHatMin` and
`MHatMax` cuts are important in KKS mode because the tabulated model is only
defined over a finite `sqrt(shat)` range.

The supplied Table 3-style cards use:

- `LHC-Instanton-Rivet-Table3-Low.in`: `MHatMin 30*GeV`, 10,000 events.
- `LHC-Instanton-Rivet-Table3-High.in`: `MHatMin 500*GeV`, 1,000 events.

## Rivet Analysis

The local Rivet plugin can be built with:

```sh
cd Rivet
make
```

The direct Herwig+Rivet cards load it with:

```text
read snippets/Rivet.in
insert /Herwig/Analysis/Rivet:Paths 0 Rivet
insert /Herwig/Analysis/Rivet:Analyses 0 QCD_INSTANTON_KKS
set /Herwig/Analysis/Rivet:Filename LHC-Instanton-Rivet.yoda
```

The analysis fills particle-level observables:

- charged tracks with `pT > 0.5 GeV` and `|eta| < 2.5`;
- anti-kt `R = 0.4` jets with `pT > 20 GeV` and `|eta| < 4.5`;
- track and jet reconstructed mass proxies;
- multiplicity, scalar `ST`, average pairwise `Delta phi`, and sphericity in
  low-mass track and high-mass jet windows.

It is intended for fast particle-level shape comparisons, not as a detector or
Delphes reproduction of the KKS plots.

## Caveats and Checks

- `alpha_s(1/rho)` is tabulated and printed by the interpolation test, but the
  current KKS weight does not use it directly. The interpolated `sigmahat` table
  should be treated as already containing the KKS rate information.
- `NAdditional` truncates the gluon multiplicity. Poisson multiplicity weights
  are not renormalized over the truncated range.
- Increasing `NAdditional` quickly increases integration dimensionality and can
  reduce sampler efficiency.
- The colour options are modelling choices. Compare colour options if
  hadronization-sensitive observables are central to the study.
- The Rivet analysis is normalized to unit area and is meant for shape checks.
  Validate rates and shapes separately before using the output for phenomenology.
