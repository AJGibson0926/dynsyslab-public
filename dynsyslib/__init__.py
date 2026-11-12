"""Plug-and-play dynsyslab system library.

Use ``load_system`` for ODE entries; select delay and lattice entries from their
own registries because their numerical contracts differ by design.
"""
from .drivers import SinusoidalDriver, TabulatedDriver, ZeroDriver
from .interfaces import DelaySystem, LatticeSystem, LibrarySystem, ODESystem
from .lib import (SystemPreset, BUBBLE_PRESETS, BubblePreset, RayleighPlessetSystem, DuffingSystem,
                  MARMOTTANT_PRESETS, MarmottantSystem, REDDY_SZERI_PRESETS, ReddySzeriSystem, VanDerPolSystem,
                  RosslerSystem, LorenzSystem, ChuaCircuit, HenonHeilesSystem, HindmarshRoseSystem,
                  HodgkinHuxleySystem, FPUTSystem, KuramotoSystem)
from .lib import NonlinearPendulumSystem, DoublePendulumSystem, SpringCoupledPendulumsSystem
from .ode import ODE_LIBRARY, SYSTEM_LIBRARY, SYSTEM_PRESETS, load_system, resolve_archived_system
from .delay import MackeyGlassSystem
from .lattice import Ising3DModel, IsingModel, PottsModel
from .observables import (Observable, ObservableCatalog, coerce_observable,
                          coerce_coordinate_observable, observable_from_summary,
                          state_observable_catalog)
from .analytic_expressions import (AnalyticExpression, state_variables,
                                   exp, log, sin, cos, tanh, sqrt)
from .legacy_dynbatch_csv import (LegacyControlHistory, LegacyDynBatchCSV,
                                  LegacySystem, LegacyTrajectory,
                                  export_legacy_dynbatch_csv,
                                  read_legacy_dynbatch_csv)

DELAY_LIBRARY = {"mackey_glass": MackeyGlassSystem}
STATISTICAL_LIBRARY = {"ising": IsingModel, "ising_3d": Ising3DModel, "potts": PottsModel}

__all__ = ["load_system", "resolve_archived_system", "ODE_LIBRARY", "SYSTEM_LIBRARY", "SYSTEM_PRESETS", "DELAY_LIBRARY", "STATISTICAL_LIBRARY", "SystemPreset", "BUBBLE_PRESETS",
           "RayleighPlessetSystem", "DuffingSystem", "MARMOTTANT_PRESETS","MarmottantSystem",
           "REDDY_SZERI_PRESETS","ReddySzeriSystem", "VanDerPolSystem", "RosslerSystem", "LorenzSystem", "ChuaCircuit",
           "HenonHeilesSystem", "HindmarshRoseSystem", "HodgkinHuxleySystem", "FPUTSystem", "KuramotoSystem", "NonlinearPendulumSystem", "DoublePendulumSystem", "SpringCoupledPendulumsSystem",
           "SinusoidalDriver", "ZeroDriver", "TabulatedDriver", "ODESystem", "DelaySystem", "LatticeSystem", "LibrarySystem",
           "Observable", "ObservableCatalog", "state_observable_catalog",
           "coerce_observable", "coerce_coordinate_observable", "observable_from_summary",
           "AnalyticExpression", "state_variables", "exp", "log", "sin", "cos", "tanh", "sqrt"]
__all__ += ["LegacyControlHistory", "LegacyDynBatchCSV", "LegacySystem", "LegacyTrajectory",
            "export_legacy_dynbatch_csv", "read_legacy_dynbatch_csv"]
