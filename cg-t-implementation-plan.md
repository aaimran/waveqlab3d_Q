# Plan: independent anelastic-cQ-cg-t and anelastic-fQ-cg-t

Status: initial production implementation completed; validation results and limits are recorded in `docs/cg-t-validation.md`. This document retains the original roadmap, including broader follow-up validation.
Repository: `/Users/aimran/Documents/Code Development/waveqlab3d_Q`.

## 1. Objective and first release

Implement two independent projected-cell attenuation responses:

- `anelastic-cQ-cg-t`: constant-Q target.
- `anelastic-fQ-cg-t`: frequency-dependent target, including the existing sharp and smooth transition functions.

The suffix `-t` identifies the traditional-compatible projected-cell backend. It does not select the old nodal CG implementation. The first release targets exactly `traditional`, `upwind`, and `upwind_drp`, all with spatial `order=6`. Each operator is enabled only after its own acceptance gates pass.

Use complete, disjoint 2×2×2 fine-node groups, with all eight mechanisms stored once per group. Velocity and stress remain at every fine node. Keep full eight-mechanism memories in boundary/interface/PML buffers and incomplete groups. Do not blend nodal CG and projected CG within a response.

Start with N=8. The requested names deliberately omit N, but this does not authorize unverified mechanism counts. Expose `n_mechanisms=8` and reject other counts in the initial implementation; structure the backend around N-sized arrays so a later coefficient-library extension can support more counts. The eightfold grouping reduction is independent of N, and N mechanisms are not assigned to individual parity positions.

Preserve all current full-memory responses and both existing nodal CG8 responses. Give the new responses distinct configuration, allocations, diagnostics, examples, and tests. Share pure coefficient/target mathematics through explicit interfaces where appropriate, without activating legacy state flags.

## 2. What was checked before making this plan

### Repository cross-checks

The current `cg8_state` in `src/anelastic_cg8_types.f90` stores `eta_cg(6,n_cg)` and its RK accumulator, with one mechanism chosen by nodal parity. That layout cannot represent the proposed cell evolution without a new ownership model.

`src/anelastic_cg8_material.f90` currently fits two coefficient sets: harmonic coarse strengths and additive full-buffer strengths. Its coarse coefficients must NOT be reused for projected cells. The new method uses the additive/full constitutive spectrum in both regions.

`src/elastic.f90` calls `RHS_Center` followed by `RHS_Near_Boundaries`; `src/block.f90:set_rates_block` wraps that call. The projected update needs a finalization step after BOTH derivative regions have contributed strain.

`src/time_step.f90` stages are: exchange fine fields; scale residuals by A; calculate block rates; exchange interface fields; enforce boundary/interface SAT conditions; update fields by B*dt. The new memories must follow exactly the same residual convention and stage timing.

`src/mpi3dcomm.f90:exchange_all_neighbors` exchanges six faces using datatypes whose transverse ranges are owned nodes. It does not establish a general edge/corner exchange for cell reduction. Reusing it unmodified for a cell spanning up to eight ranks would omit contributions.

`src/decomposition_safety.f90` gives (boundary closure, stencil reach): traditional 6=(6,3), upwind 6=(6,4), upwind DRP 6=(8,5). The existing conservative guard formula rounds closure+reach+2 upward to an even count, giving 12, 12, and 16 nodes respectively.

### Mathematical cross-checks

The weighted restriction/prolongation identities below were checked numerically on a nonuniform positive-weight cell: adjoint defect 2.22e-16, RP defect zero, and positive relaxed-energy quadratic form for a representative 80% strength sum. These are algebra checks, not a proof that all solver boundaries satisfy an energy estimate.

### Preliminary actual-stencil calculation

An independent 120-DOF periodic-cell prototype was run using coefficients extracted from production sixth-order stencils and production additive/full material fits. It contains 72 velocity/stress DOFs plus 48 cell-memory DOFs. Results are in `docs/cg-t-projection-audit.json`; the reproducible research harness is `tests/cg_t/projection_audit.py`.

The 1008 cases cover three stencils, four material cases, seven directions, nominal S frequencies 0.1/1/5, and 16/24/32/64 fine-grid points per nominal S wavelength. The material cases include equal Q=400 and Qs/Qp=400/600, Vp/Vs=2 or sqrt(3), and constant or smooth frequency-dependent targets. This sweep uses isotropic spacing only; it does not cover all fit-band endpoints, shifts, interfaces, PML, or runtime code.

| Operator | Cases | Missing clean physical branches | Maximum modal-Q difference from full | Maximum phase difference from full |
|---|---:|---:|---:|---:|
| traditional 6 | 336 | 0 | 4.24136% | 0.023779% |
| upwind 6 | 336 | 0 | 3.96780% | 0.023753% |
| upwind DRP 6 | 336 | 0 | 3.96926% | 0.023757% |

At PPW>=24, maximum modal-Q differences fall to 1.95773%, 1.74231%, and 1.74258%. The maximum positive eigenvalue real part over the sweep is 4.30e-12 (roundoff-sized), and maximum RK amplification is 1+6.7e-16 at the tested timestep.

This resolves the missing-spectrum issue in this limited prototype. It does not establish production support. It also confirms that projection has a larger spatial attenuation error than the current validated nodal upwind CG in comparable tests.

## 3. Constitutive model and projection

Let E(v) be the symmetric fine-grid strain rate computed by the selected production derivative. Let C_U be the unrelaxed stiffness tensor and ΔC_k its nonnegative bulk/shear relaxation increment. Each cell carries six stress-memory components η_k for every k=1..N.

Use the equations:

```
v_dot       = rho^-1 Div(sigma)
sigma_dot   = C_U E(v) - P sum_k eta_k       [projected-cell nodes]
eta_k_dot   = (DeltaC_k R E(v) - eta_k)/tau_k
```

At full-buffer nodes use the ordinary local equations with R=P=I and one complete mechanism set per node.

Do not average already-computed stress rates, derivatives of material coefficients, SAT forces, source terms, or PML auxiliary fields. Restrict the physical symmetric strain rate obtained from the velocity gradients. Full-buffer PML nodes use the same PML-corrected strain as the existing full attenuation implementation.

### Restriction and prolongation

For fine-node quadrature weights w_i>0 in cell c:

```
W_c = sum_i w_i
(R e)_c = sum_i w_i e_i / W_c
(P q)_i = q_c
H_f P = R^T H_c,   H_c(c,c)=W_c,   R P = I
```

For initial uniform Cartesian interior cells, w_i are equal and R is the arithmetic mean. P replicates the SAME combined cell memory stress to all eight nodes. There is no extra division by eight in P and no multiplication by eight in the mechanism strengths.

The weights must correspond to the elastic discretization's energy inner product, including volume/Jacobian and SBP weights where applicable. The initial guards exclude closure nodes from projected cells; nevertheless, identify and document each operator's actual norm and derivative adjoint relation before claiming a global stability result. A generic field-sum diagnostic is not an SBP norm.

Six-component tensor contractions must include the shear factor two. Use a Mandel basis in algebra/energy tests, or the equivalent diagonal tensor-product weight diag(1,1,1,2,2,2). Runtime storage can retain physical stress components to match the solver. Explicitly test conversion; a plain six-vector Euclidean dot product is insufficient.

### Passivity conditions and derivation to complete

Require tau_k>0, ΔC_k positive semidefinite, and C_U-sum_k ΔC_k positive definite, separately for bulk and shear channels. Zero increments are allowed and must not require inversion.

Introduce conceptual relaxed-strain variables z_k with z_k_dot=(R e-z_k)/tau_k and η_k=ΔC_k(R e-z_k)/tau_k. For a uniform material projected region, the candidate storage energy is:

```
1/2 <v,rho v>_f
+ 1/2 <e,C_U e>_f
- 1/2 <R e,sum_k DeltaC_k R e>_c
+ 1/2 sum_k <R e-z_k,DeltaC_k(R e-z_k)>_c
```

Together with the norm adjoint identity, its material dissipation is
`-sum_k <R e-z_k, DeltaC_k(R e-z_k)>_c/tau_k`.
Positive relaxed stiffness and contractive restriction give nonnegative storage for homogeneous disjoint cells. Full nodes contribute their ordinary local storage. Derive the actual elastic boundary flux/SAT contributions and confirm them with the production derivative pair. This does not prove unconditional explicit-RK stability or PML stability.

The projection principle is informed by energy-norm-compatible transfer literature, particularly [Kozdon and Wilcox](https://arxiv.org/abs/1410.5746) and [Mattsson and Carpenter](https://arxiv.org/abs/0902.2791). Those papers do not validate this attenuation method; the proposed interior equations and their tests are a new adaptation.

## 4. Coefficients and reference material

Use ONE additive/full coefficient family per material block for projected cells and full buffers. This removes a constitutive/unrelaxed-modulus jump at the layout switch, although spatial projection can still reflect waves.

For initial N=8, factor a narrowly scoped additive fitting interface from `build_cg8_coefficients(...,coarse=.false.)`, or expose an adapter that calls only that branch. Do not call the harmonic/coarse fitter, copy raw staggered-grid coarse tables, divide strengths by eight, or reinterpret the old nodal effective fit as a cell fit. The physical target function can share `fq_target_q`; set gamma=0 for cQ.

The new coefficient type stores dynamic N arrays for tau and bulk/shear strengths, corrected positive C_U moduli, reconstructed reference velocities, passivity margins, maximum dense-grid target error, and fit convergence status. Retain joint P/S fitting and exact full-spectrum reference phase-speed normalization. The same normalized moduli are written at all nodes and material ghosts in a block.

Initial requested material-fit limit: 1% default, hard maximum 2%; reject targets that fail it. Validate on a denser independent frequency grid than the optimizer, with explicit samples near transition endpoints and fref. A finite-loss material-Q target is not identical to a temporal modal-Q measurement: keep material-fit and spatial/modal gates separate and report both.

Retain sharp and smooth fQ transitions without changing their definitions. Preserve the constant-Q identity at gamma=0, using identical relaxation and fitting settings. Reject incompatible P/S targets or a negative bulk dissipation; do not clip silently.

## 5. Configuration and initial safety envelope

Use independent namelists `&anelastic_cQ_cg_t_list` and `&anelastic_fQ_cg_t_list` and distinct parameter wrappers. Proposed common fields:

| Field | Initial behavior |
|---|---|
| Qs0, Qp0 | One value per active block; initially one or two blocks |
| n_mechanisms | 8, other counts rejected initially |
| fmin, fmax, fref | Positive ordered frequencies; fref inside fit band |
| relaxation_policy | `band` or `withers-times`, only where the additive fit passes |
| coefficient_policy | `additive-passive-fit`, no harmonic option |
| fit_samples, fit_tolerance, fit_max_iterations | 256, 1e-9, 500 initially; retain convergence checks |
| max_fit_error | 0.01 default, <=0.02 allowed |
| cell_origin | Global block-node origin, default 1,1,1; periodicity fixed at 2 |
| boundary_policy | `full-buffer` only |
| buffer_layers | -1 automatic, or validated even-round-up explicit width |

fQ additionally accepts gamma, f_transition, transition_policy, transition_lower_ratio, and transition_upper_ratio with the existing definitions/defaults. cQ does not expose frequency-power controls.

Initially require uniform axis-aligned Cartesian material/metrics per block, no topography, locked interfaces, Qs0/Qp0>=400, Qp0>=Qs0, sqrt(3)<=Vp/Vs<=2, spacing aspect ratio >=0.5, CFL<=0.25, and at least 24.1 fine-grid points per shortest fitted S wavelength over the band using maximum physical spacing. These are proposed conservative limits, not inherited validation guarantees. They must survive the expanded tests before runtime enablement.

Do not lower the Q restriction merely because the scalar additive fit works at lower Q. Future lower-Q/heterogeneous/curvilinear support needs separate gates. Reject dynamic faults, unsupported orders, incompatible coefficients, nonfinite inputs, excessive buffers leaving no cells, and insufficient resolution with dedicated CFG-CGT/RUN-CGT diagnostics and a nonzero exit status. No silent full-only fallback.

## 6. Whole-cell eligibility, boundaries, and interfaces

Compute the cell's first global node using safe nested MODULO parity arithmetic; identify it by (block ID, first global q/r/s). Use 64-bit IDs/counts with checked products. Every eligible cell contains exactly eight distinct owned fine nodes globally, even when distributed across ranks.

A cell is eligible only if ALL eight nodes lie beyond the full-memory guard and outside PML. Incomplete physical-edge cells and cells intersecting excluded regions are full. Cells never cross physical block/material interfaces. Different blocks can have different uniform materials and coefficient sets; use the existing full-layout SAT coupling at their interfaces.

Automatic starting guard widths: traditional6=12, upwind6=12, DRP6=16, plus PML thickness on each applicable face. Whole-cell eligibility must expand the final exclusion consistently when that thickness is odd. These widths are conservative starting choices; transition/boundary tests can require larger ones, never smaller without evidence.

Initially do not group across internal heterogeneity, geometry jumps, fault interfaces, or blocks. The existing uniformity checks must inspect loaded material and metrics collectively, not just input labels. An eventual heterogeneous extension would require material-dependent projection and a new energy argument.

## 7. Ownership, pools, and MPI

Give each eligible coarse cell ONE owner: the rank owning its lexicographically first global fine node. Ownership is based on the block Cartesian decomposition and global indices, not rank-number ordering. Other ranks owning any of its nodes are participants. A cell can span faces, edges, or corners and touch up to eight ranks.

Store on each cell owner:

- `eta_cell(6,N,n_owned_cells)` and `deta_cell`.
- An eight-slot mapping from global fine nodes to strain contributions.
- Compact accumulated/projected six-component strain and combined memory stress.

Store on each fine-node owner:

- Node-to-cell/slot mapping for its projected nodes and full-node compact indices.
- `eta_full(6,N,n_full)` and `deta_full` for its full nodes only.
- Six-component outgoing strain and received stress data for split cells, not replicated N-mechanism memories.

Precompute owner/participant routing from the gathered per-rank owned ranges within EACH block communicator. Validate exactly eight unique slots per cell and identical eligibility decisions. Interior cells require no communication.

Use a sparse owner/participant peer exchange (preposted receives plus sends, then wait) with separate deterministic tags for strain and stress phases. Message descriptors are fixed at initialization; avoid per-stage all-to-all metadata, global cell arrays, or all-reducing every cell's data. Corner ranks communicate directly; the existing face-only exchange is insufficient.

Gather eight node strain contributions in fixed global slot order before weighted summation. Do not first sum arbitrary rank-local subsets if strict decomposition reproducibility is required. The owner combines memories in fixed k order and sends one six-component stress sum to participants. Transport owner-cell IDs only in startup descriptors; stage packets use stable precomputed offsets.

All required peers participate even if a rank owns zero cells or only full nodes. Check counts, MPI errors, empty sends, ranks with no projected contributions, and unmatched tags. Complete requests before RK updates and before deallocation. An ownership change with decomposition is acceptable; the physical sums must remain the same.

## 8. RK stage integration

Add a projected stage lifecycle around the existing elastic rate calculation:

1. Existing fine-field halo exchange and A scaling occur unchanged. Scale BOTH new memory residual pools by A exactly once.
2. Begin the projected stage: sum current owner-cell η in k order and initiate its distribution; reset only strain-contribution workspace and coverage flags, never residual accumulators.
3. Run existing center and near-boundary elastic derivatives. A new strain dispatch does a direct local full-memory update at buffer nodes. At projected nodes it records the six physical strain-rate components in the correct cell/slot, without a nodal mechanism update.
4. Finalize after both derivative regions: finish stress distribution and subtract that cell sum from stress residuals at each projected owned node exactly once. Gather all eight strain slots, form R E, and add `(ΔC_k R E-η_k)/tau_k` to the owner-cell residuals exactly once per k.
5. Finish communication before returning to existing boundary/interface SAT enforcement. Sources and SAT terms remain fine-grid terms; do not project them as if they were strain.
6. Existing B*dt update advances fine fields and both new memory pools, with each cell advanced only by its owner.

Use the production convention `Deta = A*Deta + eta_dot`, `eta += B*dt*Deta`; do not insert dt twice. Do not use the updated fine field or next-stage η while calculating the current stage. Coverage flags in debug/tests must detect duplicate or missing strain writes, especially around the center/boundary partition.

No dense six-component full-volume staging array should be the permanent default. All-local cells can use fixed-slot local accumulation (or a canonical-order gather); do not rely on kernel traversal order matching slot order. For reproducibility across decompositions, use the same canonical summation everywhere; split cells need compact per-node slot data to retain deterministic sums. Use a debug/reference gather implementation for correctness comparison, then retain the sparse compact production route.

## 9. File and integration work packages

Proposed new files (paths relative to the repository root):

| File | Responsibility |
|---|---|
| src/anelastic_cg_t_types.f90 | Independent projected state, pools, route descriptors |
| src/anelastic_cg_t_model.f90 | Additive coefficient adapter, tensor/passivity mathematics |
| src/anelastic_cq_cg_t_model.f90 | cQ reader and wrapper |
| src/anelastic_fq_cg_t_model.f90 | fQ reader and wrapper |
| src/anelastic_cg_t_layout.f90 | Global cell IDs, eight-slot maps, complete-cell eligibility |
| src/anelastic_cg_t_comm.f90 | Sparse split-cell gather/distribution and cleanup |
| src/anelastic_cg_t_material.f90 | Initialization, full/record dispatch, stage finalize, RK/stats |
| tests/cg_t/* | Independent mathematical, production-kernel, MPI, and scientific tests |
| inputfile/test_anelastic_{cQ,fQ}_cg_t_*.in | All three sixth-order examples and negative inputs |
| docs/cg-t-validation.md | Measured limits, failures, reproducible commands |

Integrate in `datatypes.f90` through distinct allocatable `cq_cg_t` / `fq_cg_t` state members; `simulation_config.f90` / `input_preflight.f90` through distinct readers, broadcasts, and allowlists; `block.f90` through initialization and the projected stage bracket; `fields.f90` through RK scaling/update; `domain.f90` through relaxation timestep, normalized-speed CFL recalculation, statistics and cleanup.

Audit every traditional/upwind/DRP6 pointwise and optimized interior hook in `JU_xJU_yJU_z6.f90`, plus ordinary/PML near-boundary hooks in `RHS_Interior.f90`. Use one explicit new dispatch interface. Do not tie new behavior to `allocated(cq8_cg)` or reuse the old one-mechanism helper.

Update `src/CMakeLists.txt` for both solver and preprocessor dependency lists, new coefficient/unit executables, and isolated CTest work directories. Follow the existing per-target Fortran module directories and MPI::MPI_Fortran linkage. Scientific tests need the configured NumPy interpreter; runtime MPI tests use the generated launcher configuration. Do not claim unregistered optional tests passed.

## 10. Validation gates before enabling a response/operator

### Gate A: algebra and constitutive kernels

Verify weighted adjoint identity, RP=I, constant-strain preservation, positive storage and nonpositive dissipation in the correct tensor inner product. Include random positive volume weights, zero mechanism strengths, independently varying bulk/shear strengths, and stress memories initially nonzero. Compare cell forcing and feedback to eight identical full-memory nodes under uniform strain. Check factor-of-eight errors explicitly.

Cross-check additive coefficients against the production full branch and independent dense-frequency reconstruction. Verify reference P/S velocities and passive strength sums. Check gamma=0 identities, both transitions, zero bulk loss, incompatible Q pairs, failed convergence, and hard maximum fit error.

Verify RK residual scaling and update against an independent small ODE reference over multiple steps; include nonzero previous-stage residuals. Exercise both owner/full pools and cleanup twice.

### Gate B: full periodic-cell spectral sweep

Build an independent periodic-cell operator from actual production derivative coefficients; confirm pointwise/optimized forward/backward arrays and weighted adjoint pairing. Include all 120 eigenbranches, not only desired modes. Expand the prototype to:

- All three operators separately; axes, face/body diagonals and intermediate directions.
- Both P and S polarizations, mode continuity tracking and near-degenerate eigenspace checks. The simple preliminary coherence selector alone is insufficient for final centered-stencil validation.
- Q=20/50/100/200 research cases and Q>=400 candidate-supported cases; several Qp/Qs and Vp/Vs values within and outside proposed limits.
- All eight cell origins, isotropic and aspect ratios 0.5/0.75/1 with axis permutations.
- Fit-band endpoints, logarithmic interior frequencies, fref and densely sampled transition edges for sharp/smooth targets and both relaxation policies.
- PPW=8/12/16/24/32/48/64/128, including candidate-runtime boundary 24.1.
- Actual RK polynomial and full projected/full-reference timestep; stable material ODEs alone are not enough.

Require no positive physical growth beyond scale-aware floating-point tolerance, no growing parasitic branch, and RK amplification <=1 within roundoff tolerance. Resolve/check centered checkerboard branches; projection does not remove the elastic centered stencil's pre-existing null modes. Modes with zero projected strain can receive no attenuation from the cell memories even when they evolve elastically. This is not a growth instability, but their source excitation and late-time contribution must be quantified against full memory; an eigenvalue-growth test alone cannot establish acceptable waveforms.

Provisional accuracy ceilings at enabled inputs: projected/full modal-Q difference <=3%, phase difference <=0.5%, and directional/polarization Q spread <=2% after separating the full operator's baseline dispersion. Report coarse/full directional effects and total target error separately. Track branch correspondence so a missing or wrong branch cannot pass by omission. Every operator must satisfy gates independently.

### Gate C: production serial wave tests

First compare the production projected kernel to the independent small-cell reference. Then use P/S packets and monochromatic excitation with uniform full and projected media, plus full/projected/full slabs. Test normal/oblique incidence, all origins, several resolutions, and interfaces moved relative to cell origin.

Require measurable signals and correct attenuation/phase against the corresponding full-memory production reference; fit decay only after separating source transients. Require additional switch-reflected amplitude <=0.1% of incident amplitude in specified separated windows. Test cQ and fQ, including transition-crossing broadband signals. A shear-only axis packet is insufficient.

Check actual convergence: elastic derivatives remain sixth order, but piecewise-constant memory projection may reduce attenuation-related convergence. State observed order; do not advertise a sixth-order complete CG scheme merely from the derivative setting.

### Gate D: MPI ownership and runtime equivalence

Test 1,2,3,4,8 ranks; construct cases forcing q/r/s cuts, edge cuts, corner cuts, and odd-length partitions. Include ranks with zero owned cells. Avoid relying on the automatic decomposition picking the intended cut; make the test decomposition explicit or assert its ranges.

Compare actual nonzero receiver histories AND matched owner-cell memories/global invariants to the serial reference. Require differences <=1e-10*signal_peak+1e-14 (and report bitwise equality when obtained), independent of origin and ownership changes. Match cell data by global ID, not local pool index. Check total cells*8 equals projected node count and every node is covered once. Assert projected cells exist in every smoke test claiming to exercise CG.

Run under bounds-checking/floating-point-trap Debug builds as well as Release. Where available use MPI timeout/checking tools to detect deadlocks and unmatched messages, without adding a new dependency requirement.

### Gate E: boundaries, locked interfaces and PML

Test free/characteristic surfaces, nonzero normal and oblique P/S waves, two differing uniform material blocks, impedance contrast, PML thickness/origin parity changes, and fields reaching buffers. Compare with matched full-layout additive coefficients. Record receiver, reflected/transmitted amplitude, energy and late-time growth; final finite maxima alone are smoke checks.

Require PML/full-buffer constitutive updates to match their independent tensor reference. Validate adequate guards with refinement and increasing-buffer comparisons. Distinguish numerical CG errors from existing full-solver/PML errors. Do not claim a new full PML energy proof from the interior projection identity.

### Gate F: resources and regression

For N mechanisms and double precision, memory values plus RK residuals use:

```
full everywhere: 96*N*(n_projected_nodes+n_full_nodes) bytes
projected pools: 96*N*(n_projected_nodes/8+n_full_nodes) bytes
```

For N=8 this is 768 bytes/full node versus 96 bytes/projected node before maps/workspace. Log owner-cell/full-node counts, pool bytes, maps, communication/workspace bytes, and peak allocation/RSS separately. Do not count participant replicas as owner memories or promise eightfold total-process savings. Avoid transient deep-copy initialization.

Benchmark sufficiently large blocks with a meaningful projected interior, at fixed resolution and rank counts. Report stage time, communication time, imbalance, peak memory and receiver errors against full memory and existing nodal CG8. Speedup is not an acceptance assumption. Run existing full and nodal-CG regressions after integration to catch accidental dispatch/allocation changes.

## 11. Ordered implementation sequence and stop conditions

1. Finalize the tensor projection/passivity derivation and expanded independent spectral harness. Stop before runtime enablement if a requested operator fails; retain the failure evidence.
2. Add independent readers/configuration and additive coefficient adapter, with malformed-input/passivity tests. Do not allocate solver state yet.
3. Implement global whole-cell layout, serial compact pools, stage collection/finalization and RK lifecycle. Validate uniform-strain and serial spectral/wave cases.
4. Implement sparse MPI owner/participant routing, including edge/corner cells, deterministic slot sums and zero-cell ranks. Validate before adding PML/interface examples.
5. Integrate all three derivative paths, direct full buffers, existing PML-corrected strain, timestep recalculation and global stats. Verify center/boundary write coverage.
6. Run complete gates A–F, clean Release/Debug builds and legacy regression suite. Resolve failures rather than weakening mode filters or accuracy thresholds.
7. Enable only the operator/input envelope that passed. Add examples, README support matrix and measured validation report. Install executables after the validated build.

If PPW=24.1 is inadequate, increase the runtime PPW floor based on demonstrated bounds. If simple cell projection remains unacceptable, investigate a higher-order, energy-compatible restriction/prolongation pair with overlapping support as a separate design: it changes ownership, guards, accuracy and storage/workspace. Do not silently substitute such an operator or retune coefficients against a single propagation direction.

## 12. Completion definition

Both response names have independent readers/state, one complete additive mechanism spectrum per eligible coarse cell, correct full-buffer updates and MPI-independent stage evolution. Each enabled stencil has published numerical limits and passing gates. No legacy response behavior changes. Documentation separates raw attenuation savings from total memory and includes exact reproducing commands, registered/skipped tests, compiler/build mode and remaining limitations.

The plan is cross-checked against present source, projection identities, relevant primary literature, and the limited prototype. Production MPI, boundary/PML, broader spectral accuracy and runtime passivity checks are REQUIRED future work; they are not described as already verified.

The original [Withers, Olsen and Day paper](https://kbolsen.sdsu.edu/PUBL_dir/withers_Qf_15.pdf) provides the frequency-target/coarse-graining context. This projected collocated-cell scheme differs from its distributed staggered-grid mechanism arrangement, so its published coarse effective coefficients and accuracy claims are not transferred to these new responses.


## 13. Reproducing the preliminary verification

These commands rebuild only the existing coefficient-export test, then run the
saved research harness. They do not enable or implement the new responses.
Run from the repository root with the configured NumPy interpreter:

```sh
cmake --build build --target cg8_coefficients_test --parallel 4
mkdir -p build/test_runs/cg8_coefficients_unit
(cd build/test_runs/cg8_coefficients_unit && ../../cg8_coefficients_test)
/Users/aimran/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 \
  tests/cg_t/projection_audit.py \
  --coefficients-dir build/test_runs/cg8_coefficients_unit \
  --output docs/cg-t-projection-audit.json
```

The JSON includes every case, per-operator extrema and the weighted projection
algebra check. It is an exploratory calculation rather than a registered
production test. Default coefficient exports depend on the existing configured
build; explicit paths make an alternate build possible.
