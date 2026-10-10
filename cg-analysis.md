# Coarse graining on a collocated grid

Coarse graining can be used on a collocated grid, but the paper's coefficients and accuracy guarantees cannot simply be transferred unchanged. The key requirement is that waves sample the distributed relaxation mechanisms correctly through our discrete spatial operator.

Coarse graining leaves the velocity/stress grid unchanged. It reduces the spatial sampling of the attenuation memory variables.

## Full versus coarse representation

For a full model, each point contains the entire relaxation spectrum:

\[
M_{\mathrm{full}}(\omega)
=M_u\left[1-\sum_{k=1}^{N}
\frac{\lambda_k}{1+i\omega\tau_k}\right].
\]

For the paper's coarse model, each point contains one relaxation mechanism:

\[
M_j(\omega)
=M_u\left[1-\frac{w_{k(j)}}{1+i\omega\tau_{k(j)}}\right].
\]

The full spectrum emerges from wave propagation through neighboring points. Consequently, coarse attenuation depends on the spatial discretization and wavelength, as well as on the local constitutive equations. This dependence is central to [Day's original analysis](https://steveday.sdsu.edu/PUBLISHED/Day_1998_coarse-grained.pdf).

The paper discussed here is Withers, Olsen, and Day (2015), *Memory-Efficient Simulation of Frequency-Dependent Q*, Bulletin of the Seismological Society of America, 105(6), 3129–3142, [doi:10.1785/0120150020](https://doi.org/10.1785/0120150020). The local reference is [2015_Withers_Memory_Efficient_Frequency_Dependent_Q.pdf](</Users/aimran/Documents/ChatGPT/A_Paper/Papers/Used/2015_Withers_Memory_Efficient_Frequency_Dependent_Q.pdf>).

## Incorporating the paper's approach into our collocated solver

### 1. Assign mechanisms over a global 2 × 2 × 2 pattern

On a uniform Cartesian grid, a candidate assignment is

\[
k(I,J,K)=1+\operatorname{mod}(I,2)
+2\operatorname{mod}(J,2)+4\operatorname{mod}(K,2),
\]

where I, J, K are global node indices with a fixed origin.

Each node stores six stress-memory components for its assigned mechanism. Its relaxation time and bulk/shear strengths come from that mechanism.

The pattern must remain identical when MPI decomposition changes. Ghost points and duplicate interface nodes must receive consistent assignments.

Collocation makes this bookkeeping straightforward because all stress components share a position. It does not, by itself, establish accuracy.

### 2. Fit coefficients for the coarse representation

In the weak-loss, equal-volume limit, the starting relationship is

\[
w_k \approx N\lambda_k.
\]

Thus, eight-mechanism full-layout strengths cannot be used directly as one-mechanism-per-node strengths.

For stronger attenuation, the paper uses a harmonic effective-modulus model:

\[
M_{\mathrm{eff}}(\omega)
=\left[\frac1N\sum_{k=1}^{N}\frac1{M_k(\omega)}\right]^{-1}.
\]

That provides a starting calibration. However, this scalar harmonic expression should not be treated as an exact description of every P/S wave direction on our collocated elastic grid.

The recommendation is to fit against this effective response initially, then validate—and, if necessary, recalibrate—against the actual collocated solver.

### 3. Check the collocated operator's interaction with the pattern

This is the most important additional step.

Period-two coefficients introduce variation at grid-scale wavenumbers. Some centered collocated derivative operators have weakly controlled or stationary alternating-sign modes. Upwind operators behave differently.

Therefore, test each supported derivative family separately:

- P waves and both S polarizations;
- propagation along axes and diagonals;
- different wavelengths measured in grid points;
- several shifts of the mechanism pattern.

A small periodic-cell, or Bloch-wave, analysis would reveal whether the pattern introduces extra branches, directional attenuation, or unstable modes. Plane-wave simulations would verify the result.

These are requirements inferred for our discretization; they are not validations supplied by the Withers paper.

### 4. Use the correct normalization and validity conditions

The reference-modulus correction must use the coarse effective response, including its phase-velocity correction.

For a scalar one-mechanism local modulus, require nonnegative strengths with w_k < 1. Applying the full-layout condition sum(lambda_k) < 1 to the sum of coarse weights is inappropriate.

For elasticity, also check the resulting bulk and shear relaxation, rather than assuming independently acceptable P/S fits guarantee a passive constitutive model.

### 5. Allocate compact memory

A true coarse implementation stores one active mechanism per node:

- six memory arrays;
- six derivative arrays for our Runge–Kutta scheme;
- no eight-slot mechanism dimension.

The existing fQ8 coarse mode selects one mechanism but still allocates eight slots. That needs to change to obtain the memory benefit.

### 6. Introduce it in homogeneous interiors first

Initially, restrict coarse mode to uniform Cartesian regions with slowly varying material properties.

Free surfaces, faults, sharp material interfaces, PML, and strongly distorted meshes require additional validation. In particular, computationally equal-sized cells can have unequal physical volumes on curvilinear grids, so the equal-volume factor N is no longer automatically appropriate.

Material discontinuities are a documented difficulty: [Kristek and Moczo (2003)](https://www.nuquake.eu/Publications/Kristek_Moczo_BSSA_2003.pdf) developed a modified coarse formulation using neighboring anelastic functions and consistent material averaging.

## Another option for a collocated grid

We could instead store all eight mechanisms on a coarser memory grid, while retaining the fine velocity/stress grid:

1. Restrict fine-grid strain rates to coarse memory locations.
2. Update all mechanisms there.
3. Transfer the summed memory contribution back to the fine stress grid.

With eight fine nodes per coarse memory location, this also reduces memory-variable storage by roughly eight.

This is a different method from the paper's nodal redistribution. Its advantage is that every coarse location retains the complete temporal spectrum. Its cost is spatial filtering and transfer operations. Restriction and prolongation should be adjoint under the appropriate discrete energy weights, and interfaces need special treatment; arbitrary interpolation could spoil stability or accuracy.

## Consequences of full versus coarse graining

| Property | Full mechanisms at each node | Coarse memory representation |
|---|---|---|
| Local attenuation spectrum | Complete at every point | Recovered over neighboring points or coarse locations |
| Memory-variable storage | Proportional to mechanism count | Potentially much smaller |
| Attenuation update work | Updates every mechanism | Fewer updates; transfer costs may apply |
| Total simulation speed | Higher attenuation cost | Improvement depends on how much runtime attenuation consumes |
| Short wavelengths | Limited mainly by the underlying solver and temporal fit | Additional attenuation-sampling limit |
| Directional behavior | No artificial mechanism-pattern variation | Pattern can introduce directional errors |
| Sharp material changes | Local coefficients are straightforward | Neighboring sampling can mix different materials |
| Stability/timestep | Relaxation and elastic constraints | Savings do not automatically permit a larger timestep |
| Validation | Spectral fit plus solver checks | Requires spectral, spatial, directional, and interface checks |

For the period-two method, Day identifies artifacts around wavelengths of four grid spacings and shorter. That is not a safe operating threshold for our solver: the acceptable minimum wavelength must be measured for its particular stencils. High-order accuracy can make coarse graining the limiting factor even when elastic propagation remains accurate. [Day (1998)](https://steveday.sdsu.edu/PUBLISHED/Day_1998_coarse-grained.pdf).

## Recommended direction

For this code, keep full fQ as the reference, develop compact coarse storage in homogeneous Cartesian interiors, and establish its usable wavelength range before extending it to faults, interfaces, and curved grids. The acceptance criterion should be matching measured attenuation and phase velocity, not merely matching shutdown diagnostics across MPI ranks.

## Mechanism counts, spatial patterns, boundaries, and interfaces

The number of mechanisms N and the spatial period are separate choices. Eight mechanisms fit naturally into eight positions in a 2 × 2 × 2 pattern, but N = 5 does not require a five-node cell, nor does N = 16 necessarily require a larger spatial period.

For our collocated solver, preserve a short, approximately balanced spatial pattern and change how many mechanisms are assigned to each position. The layouts below are design candidates requiring validation, not established results for this solver.

### Sampling-weight consistency

Suppose mechanism k occurs at a fraction p_k of the positions. In the homogeneous, equal-volume, weak-loss limit, its coarse strength starts from

\[
w_k=\frac{\lambda_k}{p_k},
\]

where lambda_k is its full-layout strength.

For a cell containing M positions, with mechanism k appearing m_k times,

\[
p_k=\frac{m_k}{M},
\qquad
w_k=\frac{M}{m_k}\lambda_k.
\]

The multiplier is determined by occurrence frequency—not automatically by N. This generalizes Day's redistribution argument; finite-loss coefficients still require calibration against the effective response and our spatial operator. [Day (1998)](https://steveday.sdsu.edu/PUBLISHED/Day_1998_coarse-grained.pdf).

### Candidate layouts for 3–16 mechanisms

The following keeps a repeating 2 × 2 × 2 pattern with eight positions.

| N | Assignment within the eight positions | Maximum mechanisms per node |
|---:|---|---:|
| 3 | Repeat mechanisms with occurrence counts `3,3,2` | 1 |
| 4 | Repeat each mechanism twice | 1 |
| 5 | Occurrence counts `2,2,2,1,1` | 1 |
| 6 | Occurrence counts `2,2,1,1,1,1` | 1 |
| 7 | Occurrence counts `2,1,1,1,1,1,1` | 1 |
| 8 | Each mechanism occurs once | 1 |
| 9 | Two mechanisms at one position; one at seven | 2 |
| 10 | Two at two positions; one at six | 2 |
| 11 | Two at three positions; one at five | 2 |
| 12 | Two at four positions; one at four | 2 |
| 13 | Two at five positions; one at three | 2 |
| 14 | Two at six positions; one at two | 2 |
| 15 | Two at seven positions; one at one | 2 |
| 16 | Two mechanisms at every position | 2 |

For N <= 8, repeated occurrences use the corresponding multiplier. For example, with three mechanisms:

\[
w_1=\frac83\lambda_1,\qquad
w_2=\frac83\lambda_2,\qquad
w_3=4\lambda_3.
\]

For 9 <= N <= 16, each mechanism occurs once within the cell, so the initial multiplier is eight—even though some nodes carry two mechanisms.

Which mechanisms are repeated or paired should depend on their strengths and relaxation times. Avoid concentrating the strongest relaxation contributions at the same node, and test multiple spatial arrangements for directional bias.

Another possibility for N < 8 is to assign each mechanism once and leave the other positions without memory variables. Day explicitly discusses such sparse layouts. They offer greater savings with packed storage, but require larger local strengths and additional checks. [Day (1998)](https://steveday.sdsu.edu/PUBLISHED/Day_1998_coarse-grained.pdf).

### Why not simply use a cell with N positions?

That is possible, but geometric factorization can produce poor patterns:

- N = 4: `2×2×1` is compact but privileges one direction.
- N = 6: `3×2×1` has unequal periods.
- N = 16: `4×2×2` permits one mechanism per node, but doubles the longest period.
- Prime counts such as 5, 7, 11, or 13 give highly elongated cells if restricted to exactly N rectangular positions.

Longer periods require longer resolved wavelengths and can increase directional errors. Allowing two mechanisms per node is often preferable to enlarging the pattern.

### Main design considerations

1. **Choose N from spectral accuracy first.** Fit the requested Q law and bandwidth. Then choose a spatial distribution. A convenient spatial cell should not dictate an inadequate temporal spectrum.
2. **Check local constitutive validity after redistribution.** With several mechanisms at node j, the scalar relaxed-modulus requirement becomes

   \[
   \sum_{k\in A_j}w_{j,k}<1.
   \]

   Elasticity additionally requires valid bulk/shear relaxation. A fit that is acceptable in the full layout can become invalid after concentrated coarse weighting.
3. **Check spatial accuracy independently of spectral fit accuracy.** Test wavelength, propagation direction, P/S polarization, stencil family, and pattern origin. Correct average weights alone do not establish correct attenuation on a collocated grid.
4. **Use physical/discrete integration weights on nonuniform grids.** Node counts represent sampling fractions only on an appropriate equal-volume interior grid. Curvilinear geometry and SBP boundary weights complicate this assumption.
5. **Allocate only the active slots.** For N = 16, compact two-slot storage gives an eightfold reduction of memory-variable storage relative to full 16-slot storage. For N = 9, padded two-slot storage gives a reduction of 9/2 = 4.5; packed storage can save more. These factors do not describe total application memory or runtime.

### Physical boundaries

A periodic cell is truncated, boundary stencils change, and boundary nodes may carry different integration weights. The interior sampling argument therefore no longer applies automatically.

For an initial implementation, use full mechanisms in a boundary buffer, with coarse storage in the interior. The buffer must cover the boundary operator's reach and should include complete coarse cells; its sufficient thickness must be tested.

Do not compensate for a missing mechanism by arbitrarily multiplying the surviving boundary weights. Nor should memory variables simply be zeroed at a free surface: the traction condition applies to the total stress.

The full/coarse transition itself needs reflection and attenuation tests.

### Interfaces and faults

Distinguish three cases:

| Interface type | Treatment |
|---|---|
| MPI partition inside one material | Continue the same global pattern; it is not a physical boundary |
| Joined computational blocks representing continuous material | Keep assignments consistent at shared physical positions |
| Material discontinuity or fault | Keep attenuation histories separate on each side; initially use full-mechanism buffers |

A coarse cell should not combine different materials merely to complete its mechanism inventory. Doing so can blur attenuation contrasts and alter reflection/transmission. Specialized interface-aware formulations exist, but they require more than rearranging mechanism indices. [Kristek and Moczo (2003)](https://www.nuquake.eu/Publications/Kristek_Moczo_BSSA_2003.pdf).

For PML, initially use full mechanisms with the existing PML correction, and test the transition into that region.

### Recommended implementation sequence

For our first implementation, support N = 8 in homogeneous Cartesian interiors, with compact one-slot storage and full-mechanism boundary/interface buffers. After establishing attenuation, phase velocity, and stability over a measured wavelength range, extend to repeated assignments for N < 8 and two-slot assignments for 9–16.

The current `anelastic-fQ` supports only 4–8 mechanisms; extending its coefficient model to 3 and 9–16 would be separate work.


## Implementation follow-up (2026-10-09)

The descriptions above of the former fQ8 coarse runtime are historical. Coarse
runtime behavior has moved to independent `anelastic-cQ8-cg` and `anelastic-fQ8-cg`
responses; other responses are full-layout. The initial collocated implementation
is deliberately limited to its measured high-Q/upwind support envelope. See
[docs/cg8-validation.md](docs/cg8-validation.md) for results and limitations.


### Sixth-order operator extension

The nodal CG8 support envelope now includes upwind order 6 and upwind DRP order 6,
in addition to upwind order 4, subject to the existing material/geometry limits.
Automatic full-buffer guards are 12 and 16 nodes respectively, plus PML thickness.
Traditional centered order 6 fails the period-two wave-response gate and remains
rejected for CG; it is available in the full-memory responses. A projected-cell
memory scheme would require a separate implementation and validation effort.
See `docs/cg8-validation.md` for the operator checks and reproducible regressions.
