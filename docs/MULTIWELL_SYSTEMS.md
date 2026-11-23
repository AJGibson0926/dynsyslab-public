# Native multiwell systems

`triple_well_duffing` and `three_leaf_hamiltonian_clover` are native two-state,
dimensionless systems. They use the standard Solver, observable catalog,
expression derivatives, input interface, archive and Viewer routes.

## Triple-well Duffing

With state `(position, velocity)=(x_1,x_2)`, the equation is

\[
 \dot x_1=x_2,\qquad
 \dot x_2=-x_1(x_1^2-a^2)(x_1^2-1)-d x_2+g u.
\]

Parameters: `inner_saddle` a=0.5, `damping` d=0, `input_gain` g=1.
Require 0<a<1 and d>=0. The conservative potential is
V(x_1)=x_1^6/6-(1+a^2)x_1^4/4+a^2*x_1^2/2.
The three minima are at 0 and +/-1; the two barriers are at +/-a.
Only d=0 and zero input conserve H=x_2^2/2+V.

Presets: `conservative_chain`, `left_libration`, `right_libration`,
`exterior_orbit`, `equal_depth_wells`. The last uses a=1/sqrt(3).
Preset descriptions apply to their default parameters; parameter overrides
can change which well a supplied initial condition occupies.

Observables include position, velocity, force (zero-input acceleration),
kinetic/potential/total energy, zero-input energy rate, input power per unit
input, squared position and its rate. Derivative links distinguish exact
all-input identities from zero-input Lie derivatives. The energy balance is
H-dot=-d*x_2^2+g*x_2*u. State labels never identify energy as conserved under
arbitrary forcing.

For the classical triple-well context see Chacón and Díaz Bejarano,
[J. Sound and Vibration 186 (1995),269–278](https://doi.org/10.1006/jsvi.1995.0448).
The exact polynomial normalization above is the API definition.

## Three-leaf Hamiltonian clover

State `(x_1,x_2)` consists of canonical phase-plane coordinates. **x_2 is not
the time derivative of x_1.** With H=(x_1^2+x_2^2)^2-x_1^3+3*x_1*x_2^2,

\[
 \dot x_1=4x_1^2x_2+4x_2^3+6x_1x_2,\qquad
 \dot x_2=3x_1^2-4x_1^3-3x_2^2-4x_1x_2^2+g u.
\]

The only configurable coefficient is `input_gain` g=1. Zero input gives the
canonical Hamiltonian flow. Its three centers lie at radius 3/4 and angles
0, 2pi/3 and 4pi/3. H=0 is the three-petal separatrix based at the degenerate
origin. These public model facts follow directly from the polynomial.

Presets: `right_leaf`, `upper_left_leaf`, `lower_left_leaf`, `exterior_orbit`,
`separatrix`. Observables include both state rates, H, both H partial
derivatives, squared phase-plane radius and its rate, real/imaginary parts of
(x_1+i*x_2)^3, zero-input H rate, and input power per unit input. Coordinate
presets expose `state`, `state_and_rate`, and `cubic_coordinates`.
The cubic coordinates are observables; the catalog does not claim they form
a globally invertible coordinate chart.

The public equation source is A. Kirichuka, *Periodic solutions for Hamiltonian
systems* (2012), section 4.6, printed p.97, equations (23)–(24), Figure 4.8:
[institutional paper](https://www.lumii.lv/uploads/sadirb./2012/LUMII_2012_A_Kirichuka.pdf).
The implementation is independent of the separate educational code pack.

## API and validation

Use `createSystem(SystemConfiguration)` with either identifier, or the normal
compiled `systems`, `observable`, `validate` and `run` CLI routes. No new
solver or archive format is introduced. Python clients use the existing
native workflow adapters.

`dynsyslab_multiwell_system_tests` checks every preset, Jacobians against
finite differences, observable derivative links, the controlled energy
balance, unforced conservation, equilibrium values, and invalid parameters.
Private spectral research is not embedded in these public model definitions.
