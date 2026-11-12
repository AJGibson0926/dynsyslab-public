# Bubble-preset provenance and limits

All values are rounded *starting points*, kept explicit so that studies can be
re-run with measured values.  They are not intended as a property database.

- Water density/viscosity near ambient conditions: NIST's [water-viscosity
  correlation](https://srd.nist.gov/JPCRD/jpcrd121.pdf) and [water transport
  tables](https://srd.nist.gov/JPCRD/jpcrd175.pdf).
- The RPE form and its Newtonian, incompressible assumptions follow Brennen's
  [Cavitation and Bubble Dynamics](https://media.library.caltech.edu/CaltechBOOK%3A1995.001/chap2.htm)
  and the microbubble review at [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC2818980/).
- Glycerol density at 20 C is documented by [PubChem](https://pubchem.ncbi.nlm.nih.gov/compound/glycerol).
- Dimethicone density and viscosity range are documented by
  [PubChem](https://pubchem.ncbi.nlm.nih.gov/compound/dimethyl%20polysiloxane).
- Blood viscosity and surface-tension ranges are summarized in this
  [biomechanics review](https://www.sciencedirect.com/science/article/pii/S0010482508001832).

The saline, glycerol/water, silicone-air interfacial-tension, CO2, and honey
entries are nominal engineering values/ranges chosen for parameter exploration.
They must be replaced by formulation- and temperature-specific measurements for
quantitative work.  Honey is non-Newtonian; the supplied entry is expressly a
stress-test surrogate, not a validated RPE model.

## Default RPE scaling convention

The `air_water_20c` default is the project reference case: \(R_0=10^{-5}\) m,
\(\kappa=1.3333\), \(p_0=10^5\) Pa, \(p_V=2330\) Pa, \(\rho_L=998\)
kg/m3, \(\sigma=0.0725\) N/m, and \(\mu_L=0.001\) Pa s. The equilibrium gas
pressure is \(p_{G0}=p_0-p_V+2\sigma/R_0\). Its frequency and characteristic
time are

\[
f_0=\frac{\sqrt{3\kappa p_{G0}-2\sigma/R_0}}
           {2\pi R_0\sqrt{\rho_L}}, \qquad T_0=\frac{1}{f_0}.
\]

The independent dimensionless time is \(\tau=t/T_0=f_0t\). The pressure scale
used consistently in the dimensionless RPE is
\(P_*=\rho_L(R_0/T_0)^2\); every pressure term is divided by \(P_*\).
The public scalar input is external acoustic pressure \(p_a\), so it enters the
pressure jump as \(-p_a\): positive input compresses the bubble. This is stored
as `external-acoustic-pressure/v2`. Archives without an explicit field retain
the earlier `legacy-positive-acceleration/v1` behavior during reconstruction.
`CL` and `CG` are not parameters of this incompressible Rayleigh--Plesset model.

## Shared ODE preset conventions

All named ODE presets are nominal, reproducible starting points rather than
claims that one parameter set is universally chaotic or physiological. Use the
common loader as `dynsyslib.load_system(name, frequency_hz=..., amplitude=...,
preset=...)`. Explicit construction keywords override preset values; the
recommended driver values stored with a preset are informational only.

- **Duffing / van der Pol / Hénon--Heiles:** standard nondimensional benchmark
  regimes: double/single wells, weak/relaxation limit cycles, and low-energy or
  near-escape Hamiltonian initial conditions.
- **Rössler:** `classic_chaos` uses \(a=b=0.2,c=5.7\), Rössler's original
  continuous-chaos selection. `common_alternative` records the frequently used
  \(a=b=0.1,c=14\) choice.
- **Lorenz:** `lorenz_1963` is \((\sigma,\rho,\beta)=(10,28,8/3)\).
  `arbabi_mezic_hdmd_2017` uses the same unforced parameters, matching the
  Arbabi--Mezić illustration. It samples \(z_1\) every 0.01 seconds and forms
  a \(10000\times500\) Hankel matrix from a random trajectory; therefore the
  preset does not assert a paper-specific initial condition. The final preset
  is a continuation starting point near the classical transition region.
- **Chua:** both presets use the canonical dimensionless double-scroll values
  \(\alpha=15.6,\beta=28,m_0=-1.143,m_1=-0.714\), with opposite lobe seeds.
- **Hindmarsh--Rose / Hodgkin--Huxley:** these provide standard slow-fast
  bursting/tonic-spiking regimes and classic squid-axon rest/DC-bias conditions.
- **FPUT / Kuramoto:** alpha and beta fixed-end 32-site chains; weak-coupling
  and synchronizing mean-field phase populations.

## References for non-RPE presets

1. E. N. Lorenz, *Deterministic Nonperiodic Flow*, J. Atmos. Sci. 20 (1963),
   doi:10.1175/1520-0469(1963)020<0130:DNF>2.0.CO;2.
2. H. Arbabi and I. Mezić, *Ergodic theory, Dynamic Mode Decomposition and
   Computation of Spectral Properties of the Koopman Operator* (2017),
   arXiv:1611.06664, Sec. 4.
2. O. E. Rössler, *An equation for continuous chaos*, Phys. Lett. A 57 (1976),
   397–398.
3. L. O. Chua, *Chua's circuit*, Scholarpedia 2(10):1488 (2007), with the
   matching dimensionless double-scroll implementation at chuacircuits.com.
4. J. L. Hindmarsh and R. M. Rose, *A model of neuronal bursting using three
   coupled first order differential equations*, Proc. R. Soc. Lond. B 221
   (1984).
5. A. L. Hodgkin and A. F. Huxley, *A quantitative description of membrane
   current and its application to conduction and excitation in nerve*, J.
   Physiol. 117 (1952), 500–544.
6. The FPUT report (Fermi, Pasta, Ulam; Tsingou's MANIAC computations) and
   Y. Kuramoto, *Chemical Oscillations, Waves, and Turbulence* (1984).
