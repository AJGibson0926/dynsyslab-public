# Fourier coefficients under periodic time changes

`dynsyslab/math/PeriodicClock.hpp` transforms a periodic observable when its
available solution uses an auxiliary time. All algorithms live in the native
mathematics library; the existing compiled `mathematics` CLI and optional
`dynsyslab_client.mathematics.query` expose them without Python numerics.

Let \(\alpha=\omega_s s\), \(dt/ds=w(\alpha)>0\), and let \(w_0\) be the
clock mean. Physical frequency is \(\Omega=\omega_s/w_0\). With
\(\beta(0)=0\), \(\beta'=w/w_0-1\), physical phase is \(\alpha+\beta\).
The Fourier change-of-variable identity is
\[
 a_n=\frac1{2\pi}\int_0^{2\pi}g(\alpha)\frac{w(\alpha)}{w_0}
          e^{-in(\alpha+\beta(\alpha))}\,d\alpha.
\]

`reparameterizePeriodicFourier` evaluates this by compensated periodic
trapezoidal summation. It accepts finite Fourier polynomials with identical
positive auxiliary frequencies. Clock coefficients must have exact Hermitian
symmetry and a positive real mean. It checks a sufficient positivity bound
using the larger of the Fourier-l1 bound and a sample/Lipschitz bound.
These are exact-arithmetic sufficient bounds for the supplied polynomial,
evaluated in floating point; they are not interval certificates and say
nothing about an omitted clock tail. The grid must resolve the supplied
harmonics. Refine its size to assess coefficient quadrature error.

The query operation is `periodic_clock`, action `transform`, with
`clock_rate`, `channels`, `samples`, and `maximum_harmonic`. Each series uses
the usual `angular_frequency` and `[harmonic,real,imaginary]` coefficient rows.
Results include transformed channels, `phase_correction`, `physical_period`,
`mean_clock_rate`, `minimum_sampled_clock_rate`, and
`clock_polynomial_lower_bound`; `certified` is false.

Action `coefficient_expansion` instead takes `series`, `clock_rate`, `harmonic`
and `order` and computes the finite convolution sum
\[
 \sum_{k=0}^K\frac{(-in)^k}{k!}[g(w/w_0)\beta^k]_n.
\]
The returned remainder bound uses the real-exponential integral remainder
and Fourier-l1 majorants. It excludes floating-point cancellation and omitted
input coefficients. Positivity of the clock is a caller precondition for
this algebraic evaluator. Orders are capped at 128 and convolution work at
the shared Fourier-algebra limit. Large harmonic/phase corrections favor
quadrature over Taylor summation.

`special-function`, function `quartic_orbit`, additionally exposes the existing
native elliptic quartic evaluator at a batch of `elapsed` values. Supply five
ascending polynomial `coefficients`, `initial` and `velocity`; these must obey
the first integral. `parameter_period` refers only to the quartic parameter.
The thin Python adapter is `periodic.evaluate_quartic_orbit`.

These are standard change-of-variable, Fourier-algebra and Taylor algorithms;
see Katznelson, *An Introduction to Harmonic Analysis*, chapter I, and the
Sundman examples in [Kudryashov–Sinelshchikov](https://arxiv.org/abs/1409.6699).
Tests compare the transformed Kepler anomaly coefficients against independent
Fourier–Bessel formulas, check non-even clocks, constant observables, changing
mean clock scales, Taylor remainders and invalid inputs.
