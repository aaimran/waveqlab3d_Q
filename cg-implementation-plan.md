# Independent eight-mechanism coarse-grained responses: implementation plan

Prepared 2026-10-09 for the collocated WaveQLab3D solver.

## 1. Scope and verification status

Implement two independent public responses:

- `anelastic-cQ8-cg`: constant-Q targets, eight relaxation mechanisms distributed over a period-two Cartesian pattern.
- `anelastic-fQ8-cg`: frequency-dependent Q targets, the same spatial layout, and the existing fQ sharp/smooth transition choices.

Remove executable coarse-graining paths from every other response. Existing full-layout responses retain their full constitutive models; the old fQ8 coarse input contract receives an explicit migration diagnostic.

This document is a plan. No solver implementation or input behavior has been changed while preparing it. Code facts have been checked against the current source, and representative baseline tests are recorded below. The proposed collocated coarse method is NOT already scientifically validated. Release requires the spectral, discrete-wave, interface, and memory gates specified here.

Related analysis: [cg-analysis.md](cg-analysis.md).

### Evidence checked against the current source

| Finding | Evidence | Consequence |
|---|---|---|
| Runtime coarse selection exists only in fQ8 | `anelastic_fq8_model.f90`, `coarse_grained_Qf8` in `datatypes.f90` | Remove this feature from fQ8 after introducing the new responses |
| fQ8 coarse storage still has eight slots per node | `material.f90`, `init_anelastic_Qf8_properties` | Do not reuse these arrays for the new implementation |
| There are 32 legacy coarse branches in the traditional/DRP source and two in the upwind helper source | `JU_xJU_yJU_z6.f90`, `RHS_Interior.f90` | Both derivative families require a migration audit |
| Existing independent fQ has 32 traditional/DRP call sites and 16 upwind dispatch call sites | Same two sources | Use this integration pattern; avoid copying constitutive equations into every stencil |
| Derivative reach and boundary closure width already have a central description | `decomposition_safety.f90`, `get_stencil_requirements` | Derive initial buffer sizes from this description rather than hardcoded per-operator assumptions |
| Subdomain ranges retain global block indices | `mpi3dcomm.f90`, `decompose1d` and `C%mq/mr/ms` | Mechanism assignment can be independent of MPI ownership |
| fQ target and error evaluation share the transition configuration | `anelastic_fq_model.f90` | Extract/reuse pure target mathematics without changing sharp/smooth behavior |
| Time stepping uses low-storage Runge-Kutta | `fields.f90`, `time_step.f90` | Both compact memory values and derivative accumulators must participate in every stage |
| Existing MPI exchange primarily exchanges velocity/stress fields | `domain.f90`, `block.f90` | The nodal coarse formulation should not introduce neighbor-memory communication unless an actual kernel requires it |

Additional mathematical cross-checks using the published table entries:

- At gamma=0 and Q0=20, Table 2 gives individual coarse strengths between 0.20645 and 0.4525, but their sum is 2.15965. A sum-less-than-one test incorrectly rejects this scalar coarse configuration; the local one-mechanism condition is the relevant test.
- At gamma=0.6 and Q0=50, the first interpolated Table 2 strength is -0.000688. Published polynomial interpolation therefore cannot be presumed to preserve nonnegativity at every input value.
- The current fQ8 conventional fit produces full-layout strengths, whereas coarse nodes require independently calibrated coarse strengths. The presently accepted conventional/coarse combination cannot be migrated by copying its coefficients unchanged.

Baseline verification completed on 2026-10-09: all six representative tests passed in 6.24 seconds:

- `fq_coefficients_unit`
- `fq_transition_unit`
- `fq8_effective_response_unit`
- `cq_coefficients_unit`
- `fq_dynamic_decomposition`
- `fq8_coarse_grain_2_dynamic_decomposition`

CTest currently registers 35 tests. Only the six audit-baseline tests above were rerun while preparing this plan. They establish a starting regression baseline, not correctness of the new collocated formulation.

### Scientific sources and their limits

- [Withers, Olsen, and Day (2015), DOI 10.1785/0120150020](https://doi.org/10.1785/0120150020): eight-mechanism frequency-dependent attenuation, harmonic coarse response, reference phase-velocity correction, published tables.
- [Day (1998)](https://steveday.sdsu.edu/PUBLISHED/Day_1998_coarse-grained.pdf): redistribution of relaxation mechanisms and its subwavelength limitation. Its demonstrated staggered-grid result is not a validation of this collocated discretization.
- [Kristek and Moczo (2003)](https://www.nuquake.eu/Publications/Kristek_Moczo_BSSA_2003.pdf): special treatment of material discontinuities and neighboring anelastic functions. It supports the need for interface care; its particular staggered-grid construction is not copied into our solver.

The implementations, fitting algorithms, thresholds, and buffer rules proposed below are engineering choices to be verified. The papers do not guarantee these particular choices.

## 2. Public contract and supported scope

### Independent response identity

Use separate readers, parameter types, configuration entries, and material-state instances for cQ8-cg and fQ8-cg. Neither response is an alias of cQ, fQ, Q8, or fQ8. Initializing a new CG response must not enable any legacy attenuation flag or allocate any legacy memory array.

Shared pure target helpers, coefficient mathematics, layout construction, and constitutive kernels are allowed. Independence means separate public configuration and runtime ownership, not unnecessary duplication of identical algorithms.

Proposed modules:

- `anelastic_cq8_cg_model.f90`: cQ8-cg configuration and constant target wrapper.
- `anelastic_fq8_cg_model.f90`: fQ8-cg configuration and frequency-dependent target wrapper.
- `anelastic_cg8_model.f90`: shared constrained coefficient fitting, response evaluation, normalization, and fit reports.
- `anelastic_cg8_layout.f90`: global parity, whole-cell eligibility, boundary/PML masks, and compact indexing.
- `anelastic_cg8_material.f90`: initialization, common strain kernel, buffer kernel, RK-state lifecycle, and diagnostics.

### Input schema

Use `&anelastic_cQ8_cg_list` and `&anelastic_fQ8_cg_list`.

Common fields:

| Field | Proposed contract |
|---|---|
| `Qs0`, `Qp0` | Arrays of length `nblocks`; required, finite, initially >=15 |
| `fref` | Positive frequency inside the fitting band |
| `fmin`, `fmax` | Finite frequencies with `0 < fmin < fmax` |
| `relaxation_policy` | `band` initially; optional `withers-times` uses the published eight times, not published strengths |
| `coefficient_policy` | Initially only `constrained-effective-fit` |
| `fit_samples` | Initially 256 logarithmic samples; validate sufficient size |
| `fit_tolerance`, `fit_max_iterations` | Positive convergence controls with explicit non-convergence failure |
| `max_fit_error` | Initially 0.05 relative material-response Q error; fail if not achieved |
| `boundary_policy` | Initially only `full-buffer` |
| `buffer_layers` | `-1` for automatic stencil-derived minimum; explicit values must meet the minimum |
| `pattern_origin` | Three integers defining global parity; useful for pattern-shift validation |

Additional fQ8-cg fields:

- `gamma`, `f_transition`.
- `transition_policy='sharp'` by default, or `smooth`.
- `transition_lower_ratio=0.8`, `transition_upper_ratio=1.2`, with the existing finite/monotonicity checks.

Exactly eight mechanisms and period two are fixed by the response names. Do not expose `n_mechanisms`, `coarse_grain`, or alternative cell sizes in these namelists.

Support one or two blocks, matching the solver's existing block limit. Allocate/read the required number of Q entries; validate only the active entries. This permits inexpensive one-block scientific tests without imposing the full cQ/fQ two-block limitation on new responses.

No unchecked table-exact strengths, silent clipping, or single-frequency-only refit will be production coefficient policies. Published tables remain available as test data and feasible initial guesses.

### Initial geometry support

Production support begins with uniform Cartesian blocks whose material properties are constant within each block. Allow differing materials between blocks, with full buffers on both sides of the interface. Detect material uniformity from the loaded/interpolated fields, not merely from the input's material-source label.

Reject curved/topographic grids and unsupported heterogeneous blocks with actionable diagnostics. Extending CG to those grids requires physical-volume/SBP-weighted sampling and additional calibration. Do not silently treat unequal-volume nodes as equal-volume samples.

Each FD family/order is enabled only after its discrete-wave validation gate passes. An unsupported combination must fail preflight, rather than borrowing another stencil's measured wavelength limit.

## 3. Constitutive model, coefficient calibration, and physical validity

### Targets

cQ8-cg uses constant P/S Q targets. fQ8-cg uses the same sharp/smooth pure target functions as full fQ. With gamma=0, identical bands, times, input velocities, and fitter settings must give the same cQ8-cg and fQ8-cg coefficients and normalized moduli.

Use identical relaxation times for both bulk and shear channels. Default band placement can match the full response's endpoints; alternative published times are scaled by f_transition. Report actual times and fit band explicitly.

### Fit bulk/shear channels while retaining a P/S input interface

A physically passive isotropic material is naturally parameterized by positive bulk modulus K and shear modulus mu. Each coarse node has one relaxation time and two channel strengths, wK and wMu.

For a homogeneous scalar channel X, the initial coarse-cell model is:

\[
R_{X,k}(\omega)=1-\frac{w_{X,k}}{1+i\omega\tau_k},
\qquad
H_X(\omega)=\left[\frac18\sum_{k=1}^{8}R_{X,k}(\omega)^{-1}\right]^{-1}.
\]

Use this harmonic model as an initial homogenization assumption, not an exact elastic/collocated theorem. Keep the scalar-response diagnostics separate from measured directional elastic responses throughout fitting and validation. Start from effective shear and bulk channels, then form a candidate effective P modulus from K_eff + 4*mu_eff/3. The final P/S response must be checked against the actual discrete elastic operator in multiple directions and polarizations.

Jointly fit the channel strengths to the supplied P/S Q targets. This prevents a superficially valid independent P/S fit from implying a negative bulk relaxation. At weak loss, the diagnostic relation is:

\[
Q_K^{-1}\approx
\frac{(K+4\mu/3)Q_P^{-1}-(4\mu/3)Q_S^{-1}}{K}.
\]

It is a useful compatibility screen, not a substitute for finite-loss checks. Physically incompatible input targets must produce an explicit error; do not adjust the requested Q values silently.

### Numerical fitter

Implement a deterministic bounded nonlinear least-squares fit of effective Q. A proposed implementation is damped Gauss-Newton with bounded/projected steps and a feasible line search, with regularization for ill-conditioned normal equations. Verify the optimizer on synthetic known-coefficient recovery cases and finite-difference derivative checks before using it for production fitting.

Constraints for every mechanism/channel:

- finite positive tau;
- `0 <= wK,wMu < 1-epsilon_margin`;
- positive unrelaxed and relaxed bulk/shear moduli;
- finite normalized responses and positive dissipation throughout the validation band.

Do not constrain the SUM of eight coarse strengths below one. Each node contains only one coarse mechanism per channel. Conversely, the separate full-buffer coefficient set must satisfy its additive relaxed-modulus conditions.

Use `8*lambda` only as a weak-loss seed when converting full strengths. It is not a general finite-loss coefficient construction. If the seed violates local constraints, construct another feasible seed and fit; do not present clipping as the final calibrated result.

The nonlinear fitter must use the actual effective response objective. Do not relabel the existing linear full-layout NNLS as a coarse fit.

### Reference phase velocity

For a scalar effective response R at fref, use:

\[
M_u=\rho c_{ref}^{2}\left[\operatorname{Re}(R^{-1/2})\right]^2.
\]

Derive the coupled bulk/shear normalization so that both prescribed P and S phase velocities are recovered. Shear normalization can be solved directly; positive bulk normalization and the P target generally need a coupled solve/refit. Require convergence and positivity, then verify reconstructed reference velocities independently.

Do not use the storage-only correction `rho*c_ref**2/Re(R)` as an exact phase-velocity normalization. This plan changes normalization only in the new CG responses and their internally owned full buffers; correcting legacy full responses is a separate task.

### Independent full-buffer model

Fit a separate additive eight-mechanism bulk/shear response to the SAME Q targets and reference velocities. It uses its own strengths and normalization, owned by the new response.

Do not obtain full-buffer coefficients by blindly dividing finite-loss coarse strengths by eight. Coarse and additive effective responses are different. Both sets must meet their own fit/phase-velocity criteria before hybrid simulations begin.

## 4. Compact storage and constitutive integration

Use a common state type with separate instances, for example `cq8_cg` and `fq8_cg`, only one allocated for the selected response.

Proposed per-block storage:

- immutable coarse and full-buffer coefficient sets;
- deterministic node-to-state lookup with regime and compact index;
- coarse memory `eta_cg(6,n_cg)` and derivative accumulator of the same size;
- buffer memory `eta_full(6,8,n_full)` and matching derivative accumulator;
- coefficient-group identifiers rather than redundant dense per-node eight-element strength arrays;
- owned-node storage wherever all verified callers operate on owned nodes.

Audit array lower bounds carefully: global block indices must not be confused with local compact positions. Use MODULO for parity, including negative indices/origins.

The common strain kernel accepts physical strain-rate components already corrected for the selected derivative/PML path. It then:

1. Looks up the node's regime and active mechanism.
2. Subtracts the appropriate memory stress-rate contribution once.
3. Updates one coarse mechanism or all eight buffer mechanisms in bulk/shear form.
4. Accumulates memory derivatives using the existing RK convention.

Avoid duplicate constitutive equations in the large stencil source. Keep ordinary, traditional/DRP, and PML entry points thin and verify that each applies attenuation exactly once.

RK scaling, updates, initialization, cleanup, and shutdown finite checks must operate on both compact pools. Never allocate legacy etaQ8, etaQf8, etacQ, or etafQ arrays for a new CG response.

### Memory accounting

For six components, eight-byte reals, and separate value/derivative arrays, memory-variable bytes are:

\[
B_{hybrid}=96(n_{cg}+8n_{full}),\qquad B_{full}=768(n_{cg}+n_{full}).
\]

The ratio approaches eight only for pure compact CG interior memory. Include lookup maps, coefficient tables, ghost storage, buffers, and other solver arrays separately in total-memory reports.

With full-buffer fraction b, the memory-variable saving factor is `8/(1+7*b)`. For b=0.2, it is approximately 3.33, not eight. Report actual counts and bytes at initialization and verify allocations in tests.

## 5. Global layout, boundaries, interfaces, and PML

### Pattern

Use global block indices and a documented origin:

\[
k=1+\operatorname{modulo}(I-o_I,2)
+2\operatorname{modulo}(J-o_J,2)
+4\operatorname{modulo}(K-o_K,2).
\]

One complete eligible supercell contains all eight mechanisms exactly once. MPI ownership boundaries may split a supercell without changing its assignment or eligibility.

### Eligibility mask and full buffers

Construct the mask in global block coordinates before rank-local compact allocation. It must not depend on MPI partition shape.

A supercell is coarse only if ALL its nodes are eligible; otherwise the entire supercell uses the full-buffer model. Incomplete cells at physical domain edges also use full storage. This preserves complete interior sampling and supports odd global grid dimensions.

Initial automatic exclusion thickness: boundary closure width plus halo reach plus one period-two guard cell, rounded outward to a complete-cell boundary. Derive closure/reach from `get_stencil_requirements`. Treat this as an initial conservative choice whose sufficiency is tested, not a proven universal buffer thickness.

Exclude:

- physical faces and their SBP/SAT closure regions;
- both sides of computational block interfaces initially, even for matching materials;
- both sides of faults, maintaining separate histories;
- every PML node and a non-PML guard strip near PML entrance;
- material-discontinuity neighborhoods if heterogeneous support is added later;
- partial/unsupported cells.

Reject a requested CG configuration if no meaningful eligible interior remains, explaining the geometry/buffer requirement. Do not silently run an entirely full-layout simulation under a CG name. An internal full-only mode may exist in the scientific test harness, not as a public substitute for an unsupported CG configuration.

### Physical boundaries and interfaces

Retain existing traction, free-surface, characteristic, fault, and coupling conditions. Apply them to total stress; do not erase memory histories at a free surface.

Use independently calibrated full coefficients and normalized moduli in buffers. The effective target must agree with the CG interior, but the discrete full/CG switch can still scatter waves. Measure this effect and widen/revise the buffer if necessary.

Do not blend coarse and full strengths arbitrarily across partial cells. Any later blending law needs its own constitutive/passivity and reflection validation.

For continuous duplicated block-face nodes, record a consistent logical origin mapping, although both sides initially remain in full buffers. For a fault, identical physical coordinates do not imply shared memory state.

PML must receive the already-established derivative correction before the full-buffer memory update. Keep attenuation active in PML; do not treat PML damping as a replacement for physical Q.

Curved meshes and general material variation remain rejected until a weighted layout and interface-specific formulation have passed separate gates.

## 6. Mandatory collocated-grid scientific gate

Before claiming production support, analyze the actual discrete operator for the candidate periodic CG medium. The scalar harmonic material fit alone is insufficient.

### Periodic-cell/Bloch harness

Construct the semi-discrete velocity-stress-memory operator using the actual interior derivative coefficients, variable constitutive factors, density, and periodic phase factors. Include the stencil's full reach, even when it exceeds the two-node cell. Do not substitute a generic second-order centered derivative.

For each supported FD family/order, sweep wavevectors and identify physical P/S branches, extra branches, growth rates, and checkerboard behavior. Include the RK amplification polynomial for the proposed timestep. All mode branches, not merely the desired propagating branch, must satisfy the chosen stability condition.

If the existing centered/upwind/DRP family fails this gate, keep that combination unsupported. Do not mask a failure with undocumented filtering or dissipation. A filtered or projected-memory alternative would be a new numerical method requiring separate planning.

### Plane-wave harness

Use homogeneous periodic or sufficiently large Cartesian tests, initially without boundaries/PML. Measure temporal decay and phase speed for modes; for receiver-based spatial attenuation, fit the complex wavenumber consistently with the chosen Q definition. Modulus Q and an apparent spatial-decay Q are not interchangeable at finite loss.

Test:

- cQ and fQ, including sharp/smooth transitions and gamma=0;
- low Q near 15/20, moderate Q, and high Q;
- compatible unequal P/S targets;
- axes, face diagonals, body diagonals, and both S polarizations;
- all eight parity-origin shifts;
- wavelength resolution sweep, e.g. 4, 6, 8, 12, 16, 24, 32 grid spacings;
- different timestep fractions and long propagation distances;
- reference frequencies away from the transition as well as at it.

Measure the elastic solver's own phase/dissipation errors separately. Upwind numerical dissipation must not be misidentified as physical attenuation. Report both raw measured behavior and the elastic baseline; an elastic correction alone is not proof of a valid CG model.

Determine a measured minimum usable wavelength for every enabled stencil. Do not use the paper's four-grid-spacing discussion as our acceptance threshold.

### Proposed acceptance thresholds

These are initial engineering targets, not paper guarantees:

- Material-response fit: maximum relative Q error <= configured bound, initially 5%, on a dense grid including fmin/fmax and fQ transition joins.
- Scalar/joint normalization reconstruction: relative reference-velocity residual <= 1e-8 in the material model.
- Discrete-wave measured Q: within 5% of the target over the declared usable band, with independently controlled measurement uncertainty.
- CG-induced phase-speed error relative to a scientifically normalized full reference: <=0.5% over that band.
- Directional/polarization and pattern-origin differences: <=2% in Q and <=0.5% in phase speed.
- No positive-growth physical or parasitic branch beyond numerical eigenvalue tolerance; all RK amplification factors within their allowed bound.

If these targets cannot be met for a configuration, narrow its supported band, recalibrate, or reject it. Do not silently relax thresholds to make tests pass.

## 7. Source integration map

| File/module | Planned changes | Required verification |
|---|---|---|
| New model modules | Independent input schemas, shared effective fitter and target helpers | Parser, constrained recovery, gamma-zero equivalence, finite checks |
| `simulation_config.f90` | Separate cQ8-cg/fQ8-cg config and presence flags | Broadcast round-trip of every field |
| `input_preflight.f90` | Register names, validate blocks/material/geometry/operator, migration diagnostics | Invalid combinations fail with explicit codes |
| `datatypes.f90` | New owned state type/instances; delete old fQ8 coarse flag | State exclusivity and compact allocation counts |
| `block.f90` | Initialize the selected independent material after grid/material setup | No legacy allocations/flags activated |
| `domain.f90` | Pass configs, timestep constraints, summaries, shutdown checks and cleanup | One/two-block behavior, MPI and lifecycle checks |
| `fields.f90` | RK accumulator scaling and memory updates for both pools | Stage-by-stage comparison against a small reference ODE |
| `RHS_Interior.f90` | Add thin CG dispatch at all upwind paths; remove old parity branches | Ordinary/PML point parity and constitutive tests |
| `JU_xJU_yJU_z6.f90` | Add shared CG strain calls in traditional/DRP paths; remove old branches | Every supported path applies one update |
| `material.f90` / `anelastic_fq8_model.f90` | Make legacy fQ8 full-only; remove harmonic coarse runtime initialization | Existing full fQ8 waveforms unchanged |
| `src/CMakeLists.txt`, tests, inputs, README | Register modules/tests, migrate examples, document supported bands and limits | Clean Release build, required tests, install |

Use one independent constitutive reference test to compare bulk/shear updates against direct tensor equations for normal and shear strains, including nonzero initial memory. Do not rely only on source-text call counts.

## 8. Remove coarse graining from all other responses

Perform this migration atomically with introduction of the new responses.

1. Remove `coarse_grained_Qf8`, old one-mechanism selection branches, conditional coarse modulus initialization, coarse parameter broadcasts, and legacy coarse summaries.
2. Make `anelastic-fQ8` and its deprecated alias full-layout only. Keep its conventional NNLS behavior and existing full-layout reference correction unchanged in this migration.
3. Reject `coarse_grain=2` in the old fQ8 contract with a diagnostic directing users to `anelastic-fQ8-cg` and the new namelist. Do not silently rename the response or reuse old coefficients.
4. In the release introducing the new CG responses, accept explicit `coarse_grain=0` only as a deprecated, ignored compatibility field. It must not select any runtime layout; no CG flag or branch remains. Update repository full-layout examples to omit it. Remove this parser-only shim in the following release, with a targeted migration error for the obsolete field. Coarse values other than zero are rejected immediately.
5. Reject old `coefficient_method='withers-2015'` under full-only fQ8 with an explanation that the published coarse weights require the new response and new calibration. Do not divide its weights by eight silently.
6. Keep cQ, fQ, Q4/Q8, constant-Q-4M/8M, and other frequency-Q variants full-layout. The fQ `relaxation_policy='fq8-table'` remains valid because it selects times only, not coarse weights/layout.
7. Migrate the old coarse example/test to the new fQ8-cg response; retain explicit migration-rejection tests for the previous contract.
8. Retain pure harmonic-response helpers only where used by the new CG model or scientific tests. Move them out of a legacy runtime response module if necessary.
9. Search all active sources, input examples, tests, documentation, and alias paths to confirm that executable CG behavior is reachable only through the two new response names. Ignore historical `.bak` copies for build coverage, but label/remove stale documentation references.

## 9. Verification matrix beyond the scientific gate

### Configuration and migration

- Valid one/two-block cQ8-cg/fQ8-cg inputs.
- Missing/wrong Q array extents; nonfinite parameters; invalid times/bands; infeasible P/S targets.
- Unsupported mesh/material/operator modes and domains lacking CG interior.
- Smooth-transition options and all MPI-broadcast fields.
- Old coarse response/options rejected with actionable guidance; old full inputs unchanged.

### Layout and allocation

- All eight mechanisms occur once in every eligible supercell.
- Negative/shifted origins use correct MODULO semantics.
- Odd global dimensions and partial physical-edge cells enter full buffers.
- MPI partitions with odd starts/extents, and partitions cutting supercells, give identical global assignments/masks.
- Distinct physical interfaces and material blocks do not share states.
- Memory-pool counts and bytes match the formula; no hidden eight-slot interior arrays.
- Destruction/reinitialization is safe; no legacy attenuation flag or array is active.

### Boundary/interface/PML behavior

- Free-surface traction and reflected P/S amplitudes/phases versus a full reference.
- Characteristic boundaries and full/CG entrance reflection tests.
- Continuous block join: changing block location/parity does not alter the modeled solution beyond tolerance.
- Material contrast and Q contrast tested separately and together; reflection/transmission against a trusted full/reference solution.
- Locked and dynamic fault tests once their integration is enabled.
- PML incidence along axes and at oblique angles; physical attenuation remains active.
- Buffer thickness sensitivity and source/receiver positions near the full/CG switch.

Initial target for extra full/CG-switch reflection: amplitude <=1e-3 of incident amplitude over the declared band, evaluated with controlled baseline/measurement error. Report when this is too strict or unachievable; do not advertise the interface as reflection-free.

### MPI and performance

- Compare actual receiver/field samples and memory summaries on 1, 2, 3, and 4 ranks with compatible decompositions; rounded final maxima alone are insufficient.
- Use explicit tolerances scaled to signal norms; investigate decomposition-dependent differences.
- Time and memory benchmarks on matched full/CG problems after correctness gates pass. Report total memory, attenuation-state memory, update cost, wall time, and buffer fraction separately.
- Verify the relaxation limit and elastic CFL using corrected unrelaxed moduli. A smaller memory footprint does not permit an assumed larger timestep.

## 10. Implementation phases and stop conditions

1. **Record baseline and migration inventory.** Freeze representative full and legacy coarse outputs, allocation counts, and available compiler/test dependencies.
2. **Build scientific constitutive and discrete-wave harnesses.** Implement constrained harmonic fitting, passive bulk/shear normalization, and an actual-stencil periodic-cell analysis before promoting the runtime mode.
3. **Introduce independent configuration and compact layout/state.** Verify masks, ownership, whole cells, and memory counts without legacy allocation reuse.
4. **Integrate cQ8-cg first.** Use homogeneous Cartesian interiors plus full buffers. Pass scalar/joint fit, tensor ODE, RK, discrete-wave, and MPI gates.
5. **Add fQ8-cg through the shared kernel.** Pass gamma-zero equivalence, sharp/smooth targets, band/transition checks, and frequency-dependent plane-wave tests.
6. **Validate buffers, block interfaces, and PML.** Enable each feature/operator only after its gate passes; keep unsupported combinations explicitly rejected.
7. **Remove legacy CG runtime paths and migrate repository inputs/tests.** Perform the public migration only once the replacement's supported configurations are concrete and validated.
8. **Run clean Release and appropriate debug/bounds builds, existing full-response regressions, install checks, and performance measurements.** Publish the measured support matrix and known limits.

A failure of the collocated Bloch/plane-wave gate is a stop condition for that numerical combination, not a reason to transfer staggered-grid claims. A failure of joint P/S passivity or fit accuracy is a rejected input, not a reason to hide negative strengths. Excessive full/CG-switch scattering requires revised buffers/calibration or a different method.

The plan's endpoint is not merely two recognized response strings. It is two independent, compact, scientifically characterized responses; explicit removal of CG from other response paths; and a documented support/migration contract with reproducible tests.


## Implementation status (2026-10-09)

The independent responses, compact pools, separate constrained coarse/full fits,
phase normalization, global whole-cell layout, RK/PML integration, migration,
and test harnesses have been implemented. See [docs/cg8-validation.md](docs/cg8-validation.md)
for the measured support envelope and reproducible checks.

The collocated gates rejected the original broad low-Q/stencil proposal. Initial
runtime support is consequently restricted to upwind order 4, high-Q Cartesian
locked-interface configurations with the documented velocity/aspect/resolution
limits. Curved/heterogeneous grids, dynamic faults, centered/other stencils, and
low-Q use are explicitly rejected. Broader free-surface/material-interface and
oblique-PML reference-waveform characterization, larger-scale performance
benchmarks, and extension of the measured parameter envelope remain roadmap
items; they are not claimed as completed or supported by the papers alone.


### Sixth-order operator extension

The nodal CG8 support envelope now includes upwind order 6 and upwind DRP order 6,
in addition to upwind order 4, subject to the existing material/geometry limits.
Automatic full-buffer guards are 12 and 16 nodes respectively, plus PML thickness.
Traditional centered order 6 fails the period-two wave-response gate and remains
rejected for CG; it is available in the full-memory responses. A projected-cell
memory scheme would require a separate implementation and validation effort.
See `docs/cg8-validation.md` for the operator checks and reproducible regressions.
