# README #

WaveQLab3D is a code for 3D seismic wave propagation and earthquake rupture dynamics. It solves the elastic wave equation in curvilinear coordinates (i.e., complex geometries) with a possibly nonplanar frictional fault interface. The current version supports off-fault viscoplasticity, spatially variable elastic properties, and several friction laws (including rate-and-state and slip-weakening). The code is under development and is available under the MIT license. Authors include Kenneth Duru, Sam Bydlon, Eric Dunham, and Kyle Withers with parallelization by Hari Radhakrishnan.

Build with CMake 3.18 or newer, a Fortran compiler, and an MPI installation
compatible with that compiler. The configuration uses compiler identity rather
than the MPI wrapper name and supports GNU, Intel/IntelLLVM, NVHPC/PGI, Cray,
and Flang toolchains, subject to the precision checks performed at configure time.
Only macOS with GNU Fortran has been verified locally; Linux and Windows need
validation with their native compiler and MPI installation. This covers build
configuration; some runtime output routines still invoke Unix directory commands.

```sh
cmake -S src -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel 4
ctest --test-dir build -C Release --output-on-failure
cmake --install build --config Release --prefix .
```

Installation places the executables in `bin`, including `.exe` on Windows.
For a multi-configuration generator, `--config Release` chooses the build
configuration. To choose an MPI compiler wrapper explicitly, pass
`-DCMAKE_Fortran_COMPILER=mpifort` on the first configuration. To omit test
executables and registration, pass `-DBUILD_TESTING=OFF`.

MPI tests use the launcher and process-count/pre/post flags discovered by CMake.
A cluster that requires an `mpirun` wrapper can set
`-DWQL3D_MPI_RUNNER=/path/to/mpirun`; launcher options can be supplied through
`MPIEXEC_PREFLAGS` and `MPIEXEC_POSTFLAGS` as CMake lists. Legacy binary comparison
tests additionally require Python 3 with NumPy; select that interpreter using
`-DPython3_EXECUTABLE=/path/to/python3`.

The legacy source requires both default REAL and DOUBLE PRECISION to occupy
eight bytes. Configuration checks this without running a program, including
when cross compiling. Other compilers can supply their precision options through
`-DWQL3D_REAL8_FLAG="..."`. User-provided CMake optimization and debug flags are
preserved.

Supported attenuation response options currently include `anelastic`, `anelastic-Q`, `anelastic-Q8`, `anelastic-cQ8-b2`, `anelastic-cQ`, `anelastic-cQ8-cg`, `anelastic-fQ8-cg`, `anelastic-fQ`, `anelastic-fQ8`, `anelastic-Qf`, `constant-Q-4M`, `constant-Q-8M`, `frequency-Q-4M`, and `frequency-Q-8M`.

For the fixed eight-mechanism constant-Q response, prefer explicit P- and
S-wave quality factors in anelastic_Q8_list:

    &problem_list
      response = 'anelastic-Q8'
    /

    &anelastic_Q8_list
      Qs0  = 50.0
      Qp0  = 50.0
      fref = 1.0
    /

Qs0 and Qp0 must be supplied together and must be positive. The stored eight-mechanism weights are normalized
spectral-shape coefficients; the RHS scales them once by the local inverse Q.

For a two-block model with a different constant Q pair in each block, use the
additive `anelastic-cQ8-b2` response. It requires `nblocks=2`; array element 1
applies to block 1 and element 2 applies to block 2:

    &problem_list
      response = 'anelastic-cQ8-b2'
      nblocks = 2
    /

    &anelastic_cQ8_b2_list
      Qs0 = 40.0, 100.0
      Qp0 = 80.0, 180.0
      fref = 1.0
      fmin = 0.05
      fmax = 20.0
      weight_method = 'fixed-q50'
    /

The reference frequency, approximation band, and eight-mechanism weight method
are shared by the blocks. Existing `anelastic-Q8` inputs and behavior are unchanged.

The independent fitted constant-Q response is selected with `anelastic-cQ`.
It requires two blocks, accepts 4 through 8 relaxation mechanisms, and keeps
its configuration and runtime memory separate from the fixed Q8 responses:

    &problem_list
      response = 'anelastic-cQ'
      nblocks = 2
    /

    &anelastic_cQ_list
      Qs0 = 40.0, 100.0
      Qp0 = 80.0, 180.0
      fref = 1.0
      fmin = 0.05
      fmax = 20.0
      n_mechanisms = 6
      coefficient_policy = 'nnls-block-ps'
      nnls_samples = 256
      nnls_objective = 'relative-q'
      nnls_tolerance = 1.0e-10
      max_fit_error = 0.10
    /

Supported coefficient policies are `nnls-shared`, `nnls-block`,
`nnls-block-ps`, and `fixed-q50`. The fixed table requires eight mechanisms;
the NNLS policies support 4, 5, 6, 7, or 8. `max_fit_error` is the maximum
allowed relative error over the requested frequency band and should be chosen
to match the accuracy required by the simulation.
The response supports `fd_type='traditional'`, `fd_type='upwind'`, and
`fd_type='upwind_drp'` with the order combinations accepted by preflight.

The independent frequency-dependent response is selected with `anelastic-fQ`.
It requires exactly two blocks and keeps its flags, fitted strengths, and memory
arrays separate from all existing attenuation responses. All 4 through 8
mechanisms are applied at every grid point; coarse graining and published Withers
strength tables are not supported by this response.

```fortran
&problem_list
  response = 'anelastic-fQ'
  nblocks = 2
/
&anelastic_fQ_list
  Qs0 = 40.0, 120.0
  Qp0 = 69.3, 155.9
  gamma = 0.6
  f_transition = 1.0
  fref = 1.0
  fmin = 0.05
  fmax = 20.0
  n_mechanisms = 8
  coefficient_policy = 'nnls-block-ps'
  relaxation_policy = 'band'
  nnls_objective = 'relative-q'
  nnls_samples = 256
  nnls_tolerance = 1.0e-10
  nnls_max_iterations = 200000
  max_fit_error = 0.10
/
```

Each block has its own P/S Q plateau, with both values at least 15. The shared
law is `Q(f)=Q0` below `f_transition`, and
`Q(f)=Q0*(f/f_transition)**gamma` above it; `gamma` must be in `[0,0.9]`.
The default `transition_policy='sharp'` preserves this law. With
`transition_policy='smooth'`, Q remains constant up to the lower transition
boundary, then joins the power law at the upper boundary.

`fmin` and `fmax` set the fitting band, while `fref` sets the storage-modulus
normalization frequency and must lie within that band. The transition may lie
outside the band. P and S strengths are fitted independently in each block and
used once in the memory equations.

`relaxation_policy='band'` places mechanisms logarithmically including both band
edges. `relaxation_policy='fq8-table'` is available only with eight mechanisms
and uses the existing fQ8 relaxation times scaled by `f_transition`, while still
fitting conventional nonnegative full-layout strengths. To compare with the
existing conventional full-layout fQ8 fit, use this policy with
`fmin=0.1*f_transition` and `fmax=10*f_transition`.

Optional transition controls are:

```fortran
  transition_policy = 'smooth'       ! default: 'sharp'
  transition_lower_ratio = 0.8       ! default: 0.8
  transition_upper_ratio = 1.2       ! default: 1.2
```

The boundaries are the ratios multiplied by `f_transition`. The smooth target
uses a cubic Hermite bridge in log-frequency and log-Q, matching the plateau's
zero slope and the power law's `gamma` slope. Q and its first derivative are
continuous at both joins. `Qs0` and `Qp0` remain plateau values; with smoothing,
Q at `f_transition` is generally slightly above that plateau.

Ratios must be finite with `0 < lower < 1 < upper`. For a monotonic cubic bridge,
smooth mode additionally requires `lower >= upper**(-2)`. Very asymmetric
intervals violating this condition are rejected because they would introduce
a dip below the plateau. Boundaries may lie outside the fitting band, but must
remain finite and positive. The transition settings are shared by both blocks
and P/S targets. Relaxation times and reference-modulus normalization are not
changed by smoothing.

This is inspired by the paper's 0.8-1.2 transition interval, not an exact
implementation of its stated intermediate power law. The smooth example is
`inputfile/test_anelastic_fQ_smooth.in`.

The fit must converge within `nnls_max_iterations` and pass the relative Q error
bound on a dense grid, explicitly checking the transition frequency and, for
smooth mode, both joins that fall within the fitting band. Unsupported policies,
negative strengths, and non-positive relaxed moduli fail initialization.
Fewer mechanisms may need a narrower frequency band or a larger explicit error
bound: the example with four mechanisms has approximately 25% maximum Q error,
while eight mechanisms give approximately 7.2%. `gamma=0` gives the constant-Q
limit; comparisons against cQ allow for the two fitters' numerical tolerances.
The relaxation timestep bound is `2*min(tau)`.

The response supports the existing upwind, traditional, and upwind DRP paths,
including the interior PML correction and MPI decomposition. A complete runnable
example is `inputfile/test_anelastic_fQ_dynamic.in`.

The independent eight-mechanism coarse responses are `anelastic-cQ8-cg` and
`anelastic-fQ8-cg`, with `&anelastic_cQ8_cg_list` and `&anelastic_fQ8_cg_list`.
Both use compact one-mechanism interior storage and independently fitted full
buffers. The supported operators are upwind order 4 or 6 and upwind_drp order 6,
with CFL <=0.25, uniform
axis-aligned Cartesian material blocks, locked interfaces, Qs0/Qp0 >=400 with
Qp0 >=Qs0, sqrt(3)<=Vp/Vs<=2, spacing aspect ratios >=0.5, and at least 16.1 points
per shortest S wavelength in the fit band. Unsupported inputs fail explicitly.

Traditional order 6 remains available in the full-memory responses. It is rejected
for these nodal CG responses: the centered derivative leaves transverse period-two
planes uncoupled for axis-aligned waves and fails the eight-mechanism wave-response
gate. Supporting centered CG requires a separately validated memory-coupling scheme.
Automatic full-buffer guards are 10 nodes for upwind 4, 12 for upwind 6, and 16
for upwind_drp 6, plus PML thickness on the relevant faces. Order-6 examples are
`inputfile/test_anelastic_cQ8_cg_upwind6.in`,
`inputfile/test_anelastic_fQ8_cg_upwind6.in`, and their `upwind_drp6` counterparts.

The new names fix N=8 and period two; there is no `coarse_grain` option. Supply
one Q entry per active block. Common controls include `fmin/fmax/fref`,
`relaxation_policy='band'` or `'withers-times'`, `fit_samples=256`,
`fit_tolerance=1e-9`, `fit_max_iterations=500`, `max_fit_error=0.05` (maximum 0.05),
`boundary_policy='full-buffer'`, `buffer_layers=-1` for automatic sizing, and
`pattern_origin=1,1,1`. `coefficient_policy` is `constrained-effective-fit`.
Frequency CG additionally accepts the existing fQ gamma/transition controls,
including sharp/smooth targets. Some target/band combinations cannot satisfy the
error bound and are rejected rather than silently adjusted.

`anelastic-fQ8` and its deprecated alias are now always full-layout conventional
NNLS. Old coarse inputs and raw Withers coefficient-method requests must migrate
to the new response/configuration. Explicit `coarse_grain=0` is a deprecated
ignored compatibility field. Published coefficients are not copied unchanged
into the new effective-response fit.

See [the CG8 validation report](docs/cg8-validation.md) for measured scientific
checks, limitations, memory accounting, and commands. Runnable examples are
`inputfile/test_anelastic_cQ8_cg_dynamic.in` and
`inputfile/test_anelastic_fQ8_cg_dynamic.in`.

Station output columns default to `t vx vy vz`. Their order can be changed in
the `output_list` namelist; for example:

    &output_list
      output_seismograms = T,
      station_output_order = 't vz vx vy'
    /

`station_output_order` is case-insensitive, accepts spaces or commas between
names, and must contain each of `t`, `vx`, `vy`, and `vz` exactly once.

A leading integer station number can optionally be included on every station
list row and used in the output filename:

    &output_list
      output_seismograms = T,
      station_number_in_list = T,
      station_number_in_filename = T
    /

    !---begin:station_list---
    1   0.693d0   0.000d0   0.000d0
    2   5.543d0   0.000d0   0.000d0
    3  10.392d0   0.000d0   0.000d0
    !---end:station_list---

This produces names such as `fname_station-1.dat`. Both options default to
false. `station_number_in_filename = T` requires
`station_number_in_list = T`.

Station files can be written directly in `station_file_directory` instead of
its `block1` and `block2` subdirectories, and stations can be restricted to one
block:

    &output_list
      output_seismograms = T,
      station_use_block_subdirectories = F,
      common_stations_blocks = 'both'
    /

`common_stations_blocks` accepts `block1`, `block2`, or `both` and defaults to `both`.
`station_use_block_subdirectories` defaults to true. When `both` is selected,
station files that use physical coordinates or station numbers in their names
receive a `_block1` or `_block2` suffix, preventing common-plane stations from
overwriting one another.

Optional commented headers and station metadata can be written at the start of
each station `.dat` file:

    &output_list
      station_add_header = T,
      station_add_metadata = T
    /

For `station_output_order = 't vz vx vy'`, the beginning of a numbered station
file is:

    # station_number: 1
    # x y z:  6.9300000000000000E-001  0.0000000000000000E+000  0.0000000000000000E+000
    # grid_i j k: 58 1 51
    # grid_x y z:  7.0000000000000000E-001  0.0000000000000000E+000  0.0000000000000000E+000
    # mapping_distance:  7.0000000000000000E-003
    # t vz vx vy

Requested and mapped physical coordinates, grid indices, and mapping distance
are included when metadata is enabled. The station number line is included when
`station_number_in_list = T`. Both preamble options default to false.

The station-to-grid mapping printed during startup is controlled separately
from the boxed station summary:

    &output_list
      output_station_mapping = F
    /

It defaults to true. When enabled, every matched station is printed on one
clearly delimited line ending in a semicolon, for example:

    station 1: distance= 7.000000E-003, indices=(58 1 51), grid_xyz=(...), requested_xyz=(...);

`output_station_info` only controls the boxed configuration summary;
`output_station_mapping` controls these individual station mapping lines.


The independent projected-cell responses are `anelastic-cQ-cg-t` and
`anelastic-fQ-cg-t`, configured with `&anelastic_cQ_cg_t_list` and
`&anelastic_fQ_cg_t_list`. They support traditional, upwind, and upwind_drp
spatial order 6. Each complete 2×2×2 node group shares all eight additive
relaxation mechanisms, driven by the arithmetic mean of physical strain rate;
the summed memory stress is returned to every node. This differs from the
one-mechanism-per-node layout of the existing `anelastic-cQ8-cg` and `anelastic-fQ8-cg` responses.

Initial support requires N=8, one or two uniform axis-aligned Cartesian material
blocks, locked interfaces, CFL<=0.25, Qs0/Qp0>=400 with Qp0>=Qs0,
sqrt(3)<=Vp/Vs<=2, spacing aspect ratios>=0.5, and at least 24.1 fine-grid
points per shortest fitted S wavelength over the band. Full-memory buffers
protect boundaries, interfaces, incomplete groups, and PML. Automatic guards
are 12 nodes for traditional/upwind 6 and 16 for upwind_drp 6, plus PML thickness.
Unsupported inputs fail explicitly; there is no full-only fallback.

Common controls are Qs0/Qp0 per active block, fmin/fmax/fref,
`n_mechanisms=8`, `coefficient_policy='additive-passive-fit'`,
`relaxation_policy='band'` or `'withers-times'`, `fit_samples=256`,
`fit_tolerance=1e-9`, `fit_max_iterations=500`, `max_fit_error=0.01`
(hard maximum 0.02), `cell_origin=1,1,1`, `boundary_policy='full-buffer'`,
and `buffer_layers=-1` for automatic sizing. Frequency-dependent response
controls retain gamma and the sharp/smooth transition options. The fitting
error is independently checked on a dense frequency grid and at transition
edges. Eight mechanisms cannot fit every target to the default 1%; the supplied
fQ examples explicitly request 2%, and failed fits remain errors.

The raw interior attenuation memories and RK residuals use one-eighth of full
storage. Cell-strain workspace, maps, MPI routes, full buffers, and other solver
arrays reduce the total memory saving. Initialization prints these categories
separately. The elastic stencil remains sixth order; the projected attenuation
approximation must not be described as a sixth-order full scheme.

Examples are `inputfile/test_anelastic_{cQ,fQ}_cg_t_{traditional,upwind,upwind_drp}6.in`
and `inputfile/test_anelastic_{cQ,fQ}_cg_t_nb1.in`. The implementation/validation
record is `docs/cg-t-validation.md`; the development plan is
`cg-t-implementation-plan.md`.
