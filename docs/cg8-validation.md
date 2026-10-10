# Independent CG8 implementation and validation

The executable now provides `anelastic-cQ8-cg` and `anelastic-fQ8-cg` with
independent readers/configuration/state. Both use eight mechanisms in complete
period-two Cartesian cells, one active mechanism per interior node, and
separately fitted eight-mechanism full buffers.

## Supported operator envelope

The numerical gates do not justify transferring the staggered-grid paper's
low-Q accuracy claims to this collocated solver. The initial runtime therefore
enforces:

- Upwind order 4 or 6, or upwind_drp order 6, with CFL <= 0.25.
- One or two axis-aligned uniform Cartesian blocks, without topography.
- Uniform loaded material within each block; blocks may differ across locked interfaces.
- Qs0/Qp0 >= 400, with Qp0 >= Qs0 and passive, feasible bulk/shear fits.
- sqrt(3) <= Vp/Vs <= 2 and grid-spacing aspect ratios in [0.5,1].
- At least 16.1 points per shortest S wavelength over the requested fitting band,
  using fitted phase velocities and the largest physical grid spacing.
- Full buffers at all physical faces, block interfaces, and throughout PML plus a guard strip.
- No dynamic faults, curved meshes, heterogeneous blocks, or other FD stencils.

The coefficient library admits lower Q for research tests, but the executable
rejects it outside this envelope. This distinction is deliberate: scalar fits
alone passed at Q=20/50 while discrete directional attenuation failed the gates.
Other bands/transition settings must still pass both the coarse and full-buffer
fit, including convergence and a maximum relative material-Q error <=5%.
The support envelope describes a conservative initial implementation, not a
proof over every possible continuously varying input parameter.

## Constitutive and storage design

The common backend fits nonnegative bulk/shear strengths jointly against P/S
quality-factor targets, using bounded damped Gauss-Newton with a projected
stationarity test and central finite differences where feasible. Full-buffer
strengths are fitted independently rather than obtained by dividing finite-loss
coarse strengths by eight. Both use an exact reference phase-speed correction.

Interior storage is `eta_cg(6,n_cg)` and its RK accumulator. Buffer storage is
`eta_full(6,8,n_full)` and its accumulator. No legacy attenuation arrays are
allocated. Pool construction writes directly into the selected state instance,
so initialization does not deep-copy the large memory pools.

Memory-variable bytes (excluding maps, fields, coefficients, and geometry) are
`96*(n_cg+8*n_full)`, versus `768*(n_cg+n_full)` for the full eight-mechanism model.
The supplied 41^3 block examples contain 8000 CG nodes and 60921 buffer nodes,
using 47555328 memory-variable bytes per block. These small examples have a large
buffer fraction; they demonstrate correctness rather than maximum memory savings.

Whole-cell eligibility and MODULO parity are calculated from global block
indices. An MPI split can cut a cell without changing its assignment. Incomplete
physical-edge cells use full buffers. Geometry uniformity is checked from loaded
metrics, not merely from the input mesh label. The timestep is reduced after
material normalization to account for corrected unrelaxed wave speeds.

## Scientific tests

`tests/cg8/scientific_gate.py` reads the actual production interior derivative
coefficients, builds the 3-D periodic-cell velocity/stress/memory operator, and
compares coarse physical branches against the separately fitted full reference.
It includes parasitic branches, RK amplification, multiple frequencies,
polarizations, directions, spacing aspect ratios, and all eight origin shifts.
It explicitly checks that the centered stencil and low-Q examples expose failures.

Measured quantities are modal decay/phase quantities compared with the full
reference. They are not silently identified with the material definition
`Q=Re(M)/Im(M)` at finite loss. The coefficient tests separately verify material-Q
fit and reconstructed reference phase speeds.

The original sweep contained 9072 periodic-cell cases, including 2016
supported upwind-4 cases across seven directions and four spacing-aspect choices.
The sixth-order extension expands this to 15120 cases, including 6048 supported
upwind-4/upwind-6/upwind-DRP-6 cases. Directional spread is assessed separately
for each operator, and all eight origin shifts are checked for each supported
operator. The harness also compares the pointwise and optimized sixth-order
interior coefficient arrays and checks the optimized negative-adjoint pair.
The original upwind-4 supported extrema were (expanded extrema are recorded below):

- Maximum coarse/full modal-Q difference: 0.920128%.
- Maximum phase difference: 0.003081%.
- Maximum directional/polarization spread: 1.336845%.

The reproducible JSON is
`build/test_runs/cg8_bloch_scientific_gate/cg8-bloch-report.json`.

`tests/cg8/transition_wave_test.py` independently advances all nine fields and
six memory components on a 3-D periodic grid with a full/CG/full slab and an
incident shear packet. The coarse pattern is resolved in both transverse
directions; this is not a scalar effective-modulus calculation. Measured
extra reflected amplitudes relative to the full reference were:

- Constant target: 0.0001124 of incident amplitude (0.01124%).
- Frequency-dependent target: 0.0001335 (0.01335%).

Both are below the proposed 0.1% switch-reflection limit. Current values are
recorded in `build/test_runs/cg8_transition_wave_gate/cg8-transition-wave-report.json`.
This establishes the tested interior-switch cases; it does not claim that every
boundary condition, source spectrum, or oblique PML incidence is characterized.

`tests/cg8/mpi_histories.py` compares actual nonzero receiver histories on
1, 2, 3, and 4 ranks with an origin that cuts cells at odd partition boundaries.
Both response examples produced identical samples in the completed clean-build
checks. Rounded final maxima are used only for additional smoke regressions.

`cg8_material_unit` checks normal/shear tensor forcing, nonzero initial memory,
RK scaling/update, compact ownership, cleanup, complete-cell parity, and a split
that cuts a cell. `cg8_coefficients_unit` checks constant/frequency fits,
normalization, readers, compatibility failures, and convergence failures.

PML regressions exercise integration and buffer allocation using the existing
full-response PML strain correction. Broad oblique-incidence PML accuracy,
free-surface reflection suites, dynamic faults, and material-discontinuity
reference-waveform suites remain follow-up validation work. Dynamic faults and
heterogeneous blocks are not enabled. The earlier implementation plan remains
an expansion/validation roadmap beyond this initial support envelope.

## Legacy migration

`anelastic-fQ8` is now full-layout only. Its runtime coarse selector, parity
branches, harmonic initialization, and eight-slot coarse-state flag are removed.
`coarse_grain=2` and `coefficient_method='withers-2015'` return guidance to the
new `anelastic-fQ8-cg` response. An explicit old `coarse_grain=0` is accepted as
an ignored deprecated parser field for this migration release. It cannot activate
coarse behavior. Other full response models retain their existing behavior.

Pure published-table harmonic calculations remain mathematical regression data;
they are not an alternative old runtime mode.

## Reproducing the checks

The scientific gates require Python with NumPy. MPI receiver tests use the same
configured launcher/flags as the solver tests. No SciPy dependency is required.

```sh
cmake -S src -B build -DCMAKE_BUILD_TYPE=Release \
  -DWQL3D_SCIENTIFIC_PYTHON=/path/to/python-with-numpy
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
cmake --install build --prefix .
```

For a bounds-checked build, configure a separate directory with
`-DCMAKE_BUILD_TYPE=Debug` and run at least the CG8 material/coefficient tests and
one-/two-block MPI smoke tests. A clean build is required when upgrading the
Fortran derived-type layout; retain no old module/object files.

Runnable examples:

- `inputfile/test_anelastic_cQ8_cg_dynamic.in`
- `inputfile/test_anelastic_fQ8_cg_dynamic.in`
- `inputfile/test_anelastic_cQ8_cg_nb1.in`

## Completed build checks

The clean Release build used the existing `.zshrc` rebuild function. All 47
registered Release tests passed, including both scientific gates, receiver
histories, compact tensor/RK tests, unsupported-input rejection, migration, and
legacy full-response regressions. Legacy Python binary comparisons requiring
NumPy on the default legacy interpreter were not registered; the separate
scientific interpreter with NumPy ran the CG gates.

The Debug build enabled bounds checks and floating-point traps. The CG material
and coefficient units passed. A pre-existing absent-optional-argument access in
MPI shutdown was exposed and corrected with an explicit PRESENT guard. The
Debug one-/two-rank one-block CG runtime test then passed with bounds checks
and floating-point traps enabled.

## Sixth-order extension

Both independent CG responses accept `fd_type='upwind', order=6` and
`fd_type='upwind_drp', order=6`. No change to the attenuation equations or
mechanism assignment is required: the production derivative paths already call
the common compact constitutive update. The preflight whitelist is extended only
for these two operators. All previous material, geometry, resolution, timestep,
and interface restrictions still apply.

The automatic guard uses boundary closure width + stencil reach + two nodes,
rounded upward to an even width: 12 for upwind 6 and 16 for upwind DRP 6.
PML thickness is added on the affected faces. The 41^3 examples contain 4096
coarse nodes for upwind 6 and 512 for upwind DRP 6, with respectively 50178816
and 52587264 bytes of memory variables per block. Wider guards reduce savings
on these small grids.

Runnable sixth-order examples are the four
`inputfile/test_anelastic_{cQ,fQ}8_cg_{upwind,upwind_drp}6.in` files.
The source's separate moment-tensor order remains 4 in each example; it is not
changed when selecting the spatial operator.

Traditional order 6 remains supported by the existing full-memory responses,
but is explicitly rejected for nodal CG. For propagation along x, the centered
y/z derivative matrices vanish on a period-two pattern at zero transverse
wavenumber. The four transverse parity planes decouple and each samples only
its two x-assigned relaxation mechanisms. Changing stencil order or relaxing
mode-selection thresholds does not restore the complete fitted spectrum.
The rejection regression covers traditional order 6 explicitly.

A possible extension would store eight mechanisms per coarse cell, restrict
fine-grid strain to that cell, and return the summed memory stress through an
energy-consistent prolongation. That is a different spatial approximation with
new MPI communication and accuracy/stability requirements. It is not implemented
or claimed as validated by this extension. No silent full-memory fallback occurs.

Sixth-order Release runtime validation completed all 11 focused tests: four
1–4-rank receiver-history comparisons, four one-/two-rank PML regressions,
unsupported-input rejection, and the two compact-material/coefficient unit tests.
All four receiver-history comparisons had zero sample differences. The material
unit was rebuilt and rerun after adding the sixth-order guard assertions.

The independent slab-transition results (extra reflection / incident amplitude)
were 0.000114217 and 0.000132225 for upwind 6, and 0.000113910 and 0.000132332
for upwind DRP 6, for constant and frequency-dependent targets respectively.
All are below 0.001. The Release executables are installed in `bin/`.

The expanded 15120-case Bloch sweep passed all gates for 6048 supported cases.
Overall maximum modal-Q difference is 0.920373%, phase difference 0.003082%,
and directional/polarization spread 1.341926%. All three independent transition
gates passed. Thus all 15 focused Release tests passed for this extension.
