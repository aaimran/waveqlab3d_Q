# README #

WaveQLab3D is a code for 3D seismic wave propagation and earthquake rupture dynamics. It solves the elastic wave equation in curvilinear coordinates (i.e., complex geometries) with a possibly nonplanar frictional fault interface. The current version supports off-fault viscoplasticity, spatially variable elastic properties, and several friction laws (including rate-and-state and slip-weakening). The code is under development and is available under the MIT license. Authors include Kenneth Duru, Sam Bydlon, Eric Dunham, and Kyle Withers with parallelization by Hari Radhakrishnan.

Supported attenuation response options currently include `anelastic`, `anelastic-Q`, `anelastic-Q8`, `anelastic-Qf`, `constant-Q-4M`, `constant-Q-8M`, `frequency-Q-4M`, and `frequency-Q-8M`.

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

Qs0 and Qp0 must be supplied together and must be positive. Existing inputs
that omit them remain supported through the legacy parameter c, which defines
Qs = c*Vs and Qp = 2*Qs. The stored eight-mechanism weights are normalized
spectral-shape coefficients; the RHS scales them once by the local inverse Q.
