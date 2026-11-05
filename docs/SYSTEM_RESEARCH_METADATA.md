# System research metadata and paper catalogues

Every native system descriptor owns a `research` section, including ODE, delay,
PDE and lattice systems. ODE instances additionally expose configured facts and
guides through `ODESystem::dynamicalKnowledge()`. These values contain versioned
facts, applicability conditions, source references, and a paper catalogue. An empty catalogue means no references have been registered;
it does not mean the system has no relevant literature.

Bibliographic records contain a stable ID, title, authors, publication year,
venue, DOI, URLs, topics, findings, coordinate conventions, applicability,
verification notes, and related fact IDs. The verification field distinguishes
full-text inspection from a reference found through another paper. For example,
RPE's Chang-Chen record is dated **1986**, not the 1960s or 1970s, and records that
its original full text has not yet been inspected.

## List and add references

The compiled CLI supports:

```text
dynsyslab-cli.exe research list-rpe.json --archive-root PATH
dynsyslab-cli.exe research add-paper.json --archive-root PATH
```

List request:

```json
{"schema":"dynsyslab.research-query/v1","action":"list","system":{"name":"rpe"}}
```

An add request uses the same schema, system and an explicit paper record:

```json
{
  "schema": "dynsyslab.research-query/v1",
  "action": "add",
  "system": {"name": "rpe"},
  "paper": {
    "id": "chang-chen-1986",
    "title": "Growth of a gas bubble in a viscous fluid",
    "authors": ["Hsueh-Chia Chang", "Liang-Heng Chen"],
    "year": 1986,
    "venue": "Physics of Fluids 29(11), 3530-3536",
    "doi": "10.1063/1.865830",
    "urls": ["https://doi.org/10.1063/1.865830"],
    "topics": ["viscosity", "Melnikov method", "cavitation threshold"]
  }
}
```

This example is already built in and returns `already_present`. New references
return `added`. Reusing an ID for another DOI is an error. DOI matching is
case-insensitive and accepts the `https://doi.org/` prefix. Entries without a DOI
use their stable ID. A duplicate is returned without overwriting its existing
notes; editing/removing existing references is not exposed by this first protocol.

C++ clients use `evaluateResearchCatalogueQuery` and `loadSystemKnowledge` from
`<dynsyslab/workflow/SystemKnowledge.hpp>`. The C ABI provides
`dsl_workflow_research_json` with the existing result-handle protocol. The
request does not fetch websites or download papers. URLs are descriptive data.

## Persistence and numerical meaning

Additions use the canonical immutable NPZ writer with a distinct `.dslpaper`
extension under `system_metadata/paper_catalogues/<system-id SHA-256>/` in the
selected archive root. Records include a schema, system identity, creation time,
and content digest. The system identifier, not the current parameter values,
selects the bibliography, so references survive preset/parameter changes. An
OS-managed writer lock prevents concurrent writers from creating conflicting
identities. These files do not enter Viewer's trajectory tree or affect its Mask.

Bibliographies live in the dynamical systems, not their solution archives. Solve
archives contain numerical data, configuration and analytic-comparison identity
without paper links or bibliography snapshots. Existing archives remain immutable.
The phase sidecar records numerical guide checks alongside the system cache.
Adding a paper does not invalidate a numerical cache; changing a built-in
executable guide or its metadata does change the cache key.

Paper additions are advisory metadata. They never install code, executable
expressions, equilibrium locations, or manifold curves. Native system guide
callbacks supply candidate seeds. The sidecar checks the actual zero-input field
and Jacobian, reports seed residuals in `knowledge_checks`, and checks applicable
invariants along sampled manifold curves. Published parameter-specific results
are not presumed to hold for other parameters or coordinate systems.

## RPE references currently registered

- Chang and Chen (1986), *Growth of a gas bubble in a viscous fluid*,
  [DOI](https://doi.org/10.1063/1.865830). Bibliography verified against
  [Chang's publication list](https://academicweb.nd.edu/~changlab/Publications.html);
  Melnikov use identified in the introduction and reference 4 of
  [Harkin et al.](https://arxiv.org/pdf/physics/9911072).
- Hegedus and Kullmann (2012), *Basins of attraction in a harmonically excited
  spherical bubble model*, [publisher](https://pp.bme.hu/me/article/view/1243).
  Computes saddle-periodic-orbit manifolds and attraction basins in a forced,
  viscous, isothermal regime below vapor pressure.
- Harkin, Nadim and Kaper (1999), *On acoustic cavitation of slightly subcritical
  bubbles*, [arXiv](https://arxiv.org/abs/physics/9911072).
- Prosperetti (1975), *Nonlinear oscillations of gas bubbles in liquids:
  transient solutions and the connection between subharmonic signal and
  cavitation*, [author copy](https://pages.jh.edu/aprospe1/publications/PapersPublished/GasBubbles/JASA_NonlinOscilTrans.pdf).
  The separatrices in this analysis are in averaged amplitude coordinates.
- Kudryashov and Sinelshchikov (2014), *Analytical solutions of the Rayleigh
  equation for empty and gas-filled bubble*, [arXiv](https://arxiv.org/abs/1409.6699).
- Mancas and Rosu (2016; preprint 2015), *Evolution of spherical cavitation
  bubbles: parametric and closed-form solutions*, [arXiv](https://arxiv.org/abs/1508.01157).

## Heat and delay references

The heat-equation system owns the MIT heat-equation/Fourier-series course reference and records continuum versus discrete boundary-mode conventions. Mackey–Glass owns its [1977 model paper](https://doi.org/10.1126/science.267326) and the [2015 delayed-exponential reference](https://doi.org/10.1186/s13662-015-0538-z). The native scalar constant-history formula is an independently checked reformulation; the latter paper studies a second-order pure-delay oscillator. Neither bibliography is copied into PDE or DDE solution archives.

## Lattice-system applicability

Ising and Potts system catalogues record the relevant equilibrium statistical-mechanics literature. These results concern equilibrium ensembles or thermodynamic limits; they are not deterministic trajectories of the two quantities recorded during stochastic sweeps. The catalogue keeps this distinction explicit and does not register a fabricated analytic sweep trajectory or apply the ODE phase-space sidecar to those observations.

The double-pendulum catalogue includes Enolskii–Pronine–Richter (2003) and Bocci–Mingari Scarpello (2022, arXiv v2). Native metadata records checked equilibrium candidates, energy/angular-momentum applicability and the reduced relative separatrix. See [the gravity-free reference](GRAVITY_FREE_DOUBLE_PENDULUM.md).


## Remaining ODE catalogue reviews

Hodgkin–Huxley metadata reconstructs stationary gates from a candidate voltage and searches the scalar native current balance. The source convention is converted to the native millivolt voltage and seconds clock (gating rates carry a factor of 1000). A fixed-voltage exponential is attached only to a verified invariant limit of the free-voltage equations.

Kuramoto metadata records finite-array phase-sum cancellation, equal-frequency synchrony and the exact two-oscillator Adler relation. Global phase shifts give a neutral direction, so the locking guide is not an isolated-hyperbolic-point claim. Its catalogue adds [Ha–Ko–Park–Zhang (2016), Example 4.2](https://doi.org/10.4171/EMSS/17), plus the [Gupta–Campa–Ruffo review](https://arxiv.org/abs/1403.2083) consulted through the Class Hub. The native two-oscillator reference uses a projective matrix exponential and retains winding through phase slips.

Spring-pendulum energy uses the actual bob separation and spring extension. Vertical rods are equilibrium seeds only when their spring torque vanishes. The synchronous equal-length, unstretched-spring invariant set supports an exact nonlinear pendulum reference; small-angle normal modes are not substituted for the geometric model.

The custom Marmottant and Reddy–Szeri native classes expose configured knowledge directly. Shell metadata restricts smooth phase-space work to memory-independent branches. Reddy–Szeri metadata records its piecewise gas exponent, the free axial coordinate and exact stationary-radius shape limit. The latter is independently derived from the implemented closure; [Reddy–Szeri's 2002 shape-stability paper](https://doi.org/10.1063/1.1483840) is separately catalogued with explicit bibliographic-only verification from the author's publication list.

All candidate searches are bounded by the requested domain. Sign-change scans are seed generators, not exhaustive root certificates. Sidecar field-residual and stability checks remain authoritative, and all paper records remain in system metadata rather than solution archives.
