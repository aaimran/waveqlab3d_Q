# Independent projected-cell responses: implementation and validation

`anelastic-cQ-cg-t` and `anelastic-fQ-cg-t` implement a common projected-cell
backend with independent readers, configuration flags, and state allocations.
Both retain eight additive relaxation mechanisms once per complete disjoint
2×2×2 fine-node group. They support traditional, upwind, and upwind_drp spatial
order 6 within the restricted envelope below.

## Numerical and implementation design

The arithmetic mean of the eight physical symmetric strain rates drives every
cell mechanism. Its summed stress memory is returned unchanged to every member
node: restriction divides by eight; stress return does not. The additive/full
coefficient family and corrected unrelaxed moduli are identical in cell and full
regions. No harmonic nodal-CG coefficients or strengths multiplied by eight are
used. Material fitting is checked independently on at least 1025 frequencies,
including reference frequency and in-band transition endpoints.

For homogeneous Cartesian interior groups this restriction/prolongation is a
volume-weighted adjoint pair. Nonnegative bulk/shear increments, positive times,
and positive relaxed bulk/shear stiffness are required. The six-component
stress/strain tensor product uses the proper shear factor two; independent
energy checks use Mandel coordinates. The interior passivity construction is
not a new proof of the existing full boundary/SAT/PML discretization.

Each cell has one owner, selected by the rank holding its first global node.
Other ranks send their strain slots and receive the six-component summed stress.
Routes include direct edge/corner peers; there is no per-stage global reduction
of cell values and no replication of all eight memories at participants.
Startup route checks compare global anchors/slots. Every stage checks for
missing or duplicate derivative writes. Cell sums and mechanism sums use fixed
slot/mechanism order, including all-local cells.

The common block-rate call is bracketed by stage begin/finalize. Full nodes update
locally with physical/PML-corrected strain; cell nodes record strain. Finalization
exchanges slot values/current memory stresses and adds cell derivatives before
boundary/interface SAT enforcement. Memory residuals are scaled by A, then
memories are updated by B*dt exactly once, using the existing RK convention.

## Runtime envelope

- N=8; other counts are rejected in this initial release.
- One or two uniform, axis-aligned Cartesian material/metric blocks.
- Locked interfaces; no topography, dynamic faults, or heterogeneous blocks.
- Qs0/Qp0>=400 with Qp0>=Qs0 and a feasible passive joint bulk/shear fit.
- sqrt(3)<=Vp/Vs<=2; minimum/maximum grid-spacing ratio>=0.5.
- CFL<=0.25 and at least 24.1 fine-grid points per shortest fitted S wavelength
  over the fit band, using the largest physical grid spacing and the fitted
  relaxed-shear speed as a conservative lower bound for phase speed.
- Full memories at boundaries, interfaces, PML and incomplete/excluded groups.
- Automatic guards: 12 nodes for traditional/upwind 6, 16 for upwind_drp 6,
  plus PML thickness on affected faces.
- `max_fit_error=0.01` by default, hard maximum 0.02. The supplied smooth fQ
  examples explicitly request 0.02 because this eight-mechanism spectrum does
  not pass 0.01 for those settings. Failed fits remain errors.

The examples use fmax=5 to satisfy the stricter projected spatial-resolution
floor on their small grids. Both sharp/smooth target definitions are retained;
gamma=0 has the constant-target coefficient identity. The supplied short-pulse input examples are integration/decomposition smoke tests.
For physical attenuation studies, the fit band must cover the important source
frequencies and the mesh must satisfy the resolution floor for that band.
This support envelope is
a conservative tested starting point, not a proof for every continuous target,
transition width, source spectrum, or material ratio in those intervals.

## Storage and resource accounting

For double precision and N=8, memory values plus RK residuals occupy
`768*(n_owned_cells+n_full_nodes)` bytes globally. Relative to 768 bytes per
full-memory node, projected nodes use 96 bytes each for these pools: an eightfold
raw interior attenuation-memory reduction.

The implementation deliberately keeps compact eight-slot cell-strain workspace
for deterministic summation and missing/duplicate-write detection. All-local
cells therefore use 384 strain-workspace bytes and 48 stress-feedback bytes per
cell, before maps/flags and any split-cell communication. Pool plus this real
workspace is 150 bytes per projected node, giving about 5.12-fold reduction
against full attenuation pools alone. The raw eightfold claim does not include
this workspace, full buffers, maps, solver fields, geometry, or MPI startup/requests.

Initialization prints pool bytes, real strain/feedback/communication workspace,
and index/route/coverage array bytes separately. It allocates state directly;
there is no deep-copy initialization of the large attenuation pools. Startup
peer construction has temporary route/buffer copies, and each stage has small
MPI request workspace. These are not included in the printed persistent-array
categories. Full simulation peak RSS and wall-time speedup still require a
representative large-grid benchmark; the small supplied examples have a large
full-buffer fraction and do not demonstrate eightfold total savings.

## Completed scientific checks

The initial exploratory 1008-case audit is retained in
`docs/cg-t-projection-audit.json`. The registered expanded gate compares the
120-DOF projected periodic-cell operator with the fitted additive full reference,
using actual production stencils and checking optimized coefficient arrays.
Its 9000 supported cases cover all three sixth-order operators, four material
cases, ten directions, four spacing-aspect choices, fit-band/transition samples,
and 24.1/32/64 nominal S points per wavelength. All eigenbranches are included in
stability/RK checks; physical branch selection requires two S and one P branch.

Measured expanded extrema:

- Projected/full temporal modal-Q difference: 1.942458% (limit 3%).
- Phase difference: 0.011535% (limit 0.5%).
- Directional/polarization Q spread: 1.514379% (limit 2%).

Material Q and temporal modal Q are distinct finite-loss measurements; they are
reported separately rather than equated. Origin-shift checks cover all eight
origins for representative constant and frequency-dependent cases. The target
material fit is independently checked by the production additive adapter.

The centered elastic stencil still has checkerboard/null modes. Projection fixes
access to the complete relaxation spectrum for physical branches but does not
remove elastic checkerboards; modes with zero projected strain can remain
undamped. No growing branches were found in the tested sweeps. This is not a
claim that arbitrary grid-scale sources will reproduce full-memory damping.

The independent nine-field packet test uses a full/projected/full slab and evolves
all eight cell mechanisms, with P and S packets, constant/frequency targets, and
all three stencils. All six gates passed the added-reflection limit 0.001 of
incident amplitude. The largest measured ratio was 0.000144875 (0.0144875%).
These normal-incidence cases do not establish all oblique transition waveforms.

## Production and MPI checks

The tensor/RK unit uses nonzero initial memories and prior RK residuals, spatially
varying prescribed strain, and an independent 3×3 tensor reference over three
stages. It checks full nodes, projected owner-cell forcing, stress redistribution,
RK scaling/update, legacy-allocation isolation, and repeated cleanup.

MPI variants use 1/2/4/8 ranks with global slot ordering and direct edge/corner
routes. A one-cell variant on 1/2/3/4/8 ranks covers zero-cell owners, entirely
full-only ranks, and a single cell spanning all eight ranks. An explicit check
requires seven peers on the eight-rank owner so this is an actual corner-routing
test rather than a nominal rank-count smoke test.

All six response/operator examples produced identical nonzero receiver samples
on 1/2/3/4 ranks in the completed checks. The receiver gate compares samples,
not merely rounded final maxima. Further smoke tests compare final states across
one/two ranks and one/eight ranks for one-block examples.

Longer PML/interface regressions add PML on outer x faces, a free y face, and a
second block with different density and P/S speeds. They run long enough for
waves to reach the boundary/PML buffers and compare nonzero final field/memory
states on one/two ranks. These are integration/stability/decomposition checks,
not a broad quantified oblique-PML or transmission-accuracy reference suite.

Negative regressions cover unsupported order, low Q, topography, CFL, fault
coupling, buffers excluding every projected cell, insufficient resolution,
incompatible P/S targets, unsupported mechanism count/coefficient policy, and
the hard fit-error cap. No full-only fallback occurs.

## Reproduction

```sh
cmake -S src -B build -DCMAKE_BUILD_TYPE=Release \
  -DWQL3D_SCIENTIFIC_PYTHON=/path/to/python-with-numpy
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure -R 'cgt_'
ctest --test-dir build --output-on-failure -E 'cgt_'
cmake --install build --prefix .
```

A fresh Debug directory with GNU enables bounds checking, backtraces, NaN
initialization and invalid/zero/overflow traps. Run the coefficient, material,
corner/zero-cell tests and at least one one-/eight-rank solver case in that build.
MPI tests use the configured launcher and flags. Scientific gates require NumPy;
legacy binary-comparison tests using a separate default interpreter are not
registered if that interpreter lacks NumPy.

Broader lower-Q support, more mechanism counts, curvilinear/heterogeneous
projection, quantified oblique boundary/PML accuracy, full convergence/source
alias studies, and large-grid resource benchmarks remain follow-up work. The
initial implementation does not claim these extensions are validated.


## Build/test completion record

The clean Release build completed successfully. The 57 previously registered
full-memory/nodal-CG regressions all passed after integration. The new projected
response tests passed coefficient/tensor/RK, 1–4-rank receiver comparisons for
all six response/operator combinations, six P/S transition gates, six longer
PML/interface tests, eight-rank corner/zero-owner routes, and invalid-input checks.
Final one-/eight-rank solver examples explicitly require a 2×2×2 rank topology
and seven cell peers, so they exercise actual solver corner communication.

The fresh Debug build enabled bounds checks, NaN initialization and floating-point
traps. The original direct singleton MPI kernel launch hit an Apple Metal-driver
instruction during MPI startup, before the cell kernel ran; debugger inspection
located the failure outside solver code. The MPI-initializing unit now uses the
configured launcher even for one rank. All corrected Debug coefficient/kernel
and zero-owner/corner MPI tests passed, without disabling numerical traps.


Final completion: all 95 registered Release tests passed across the validation
runs (57 prior full/nodal regressions and 38 projected-response tests). A clean
Release rebuild was followed by 13 focused kernel/input/corner-solver tests and
the updated 9000-case gate; all passed. The final Debug run passed all 11 selected
coefficient/kernel/MPI/one-eight-rank traditional solver checks. The spectral
report also includes six positive-storage/nonpositive-dissipation matrix checks,
48 origin shifts, and 108 low-Q research stability cases; lower Q is not enabled.
The installed executables are `bin/waveqlab3d` and `bin/pre_wql3d`.

Reports are under `build/test_runs/cgt_bloch_scientific_gate/`,
`build/test_runs/cgt_transition_*/`, and the six
`build/test_runs/{cQ,fQ}_cgt_*_histories/` directories. The Release build log is
`build/cgt-clean-release.log`. Debug build diagnostics were captured in
`/tmp/cgt-debug-build.log`; regenerate the Debug build to reproduce them.
