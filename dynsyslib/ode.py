"""ODE plug-in registry and shared system-construction overhead.

Individual equations of motion live in their own modules beside ``rpe.py``.
This module deliberately owns only registry lookup and common driver attachment.
"""
from __future__ import annotations
from .drivers import SinusoidalDriver, ZeroDriver
from .lib.rpe import BUBBLE_PRESETS, BubblePreset, RayleighPlessetSystem
from .lib.duffing import DUFFING_PRESETS, DuffingSystem
from .lib.marmottant import MARMOTTANT_PRESETS, MarmottantSystem
from .lib.reddy_szeri import REDDY_SZERI_PRESETS, ReddySzeriSystem
from .lib.van_der_pol import VAN_DER_POL_PRESETS, VanDerPolSystem
from .lib.rossler import ROSSLER_PRESETS, RosslerSystem
from .lib.lorenz import LORENZ_PRESETS, LorenzSystem
from .lib.chua import CHUA_PRESETS, ChuaCircuit
from .lib.henon_heiles import HENON_HEILES_PRESETS, HenonHeilesSystem
from .lib.hindmarsh_rose import HINDMARSH_ROSE_PRESETS, HindmarshRoseSystem
from .lib.hodgkin_huxley import HODGKIN_HUXLEY_PRESETS, HodgkinHuxleySystem
from .lib.fput import FPUT_PRESETS, FPUTSystem
from .lib.kuramoto import KURAMOTO_PRESETS, KuramotoSystem
from .lib.nonlinear_pendulum import NONLINEAR_PENDULUM_PRESETS, NonlinearPendulumSystem
from .lib.double_pendulum import DOUBLE_PENDULUM_PRESETS, DoublePendulumSystem
from .lib.spring_coupled_pendulums import SPRING_COUPLED_PENDULUM_PRESETS, SpringCoupledPendulumsSystem

ODE_LIBRARY = {
    "rpe": RayleighPlessetSystem,
    "duffing": DuffingSystem,
    "marmottant": MarmottantSystem,
    "reddy_szeri": ReddySzeriSystem,
    "van_der_pol": VanDerPolSystem,
    "rossler": RosslerSystem,
    "lorenz": LorenzSystem,
    "chua": ChuaCircuit,
    "henon_heiles": HenonHeilesSystem,
    "hindmarsh_rose": HindmarshRoseSystem,
    "hodgkin_huxley": HodgkinHuxleySystem,
    "fput": FPUTSystem,
    "kuramoto": KuramotoSystem,
    "nonlinear_pendulum": NonlinearPendulumSystem,
    "double_pendulum": DoublePendulumSystem,
    "spring_coupled_pendulums": SpringCoupledPendulumsSystem,
}
SYSTEM_LIBRARY = ODE_LIBRARY
SYSTEM_PRESETS = {"rpe": BUBBLE_PRESETS, "duffing": DUFFING_PRESETS,
                  "marmottant": MARMOTTANT_PRESETS,"reddy_szeri": REDDY_SZERI_PRESETS,
                  "van_der_pol": VAN_DER_POL_PRESETS,
                  "rossler": ROSSLER_PRESETS, "lorenz": LORENZ_PRESETS, "chua": CHUA_PRESETS,
                  "henon_heiles": HENON_HEILES_PRESETS, "hindmarsh_rose": HINDMARSH_ROSE_PRESETS,
                  "hodgkin_huxley": HODGKIN_HUXLEY_PRESETS, "fput": FPUT_PRESETS, "kuramoto": KURAMOTO_PRESETS,
                  "nonlinear_pendulum": NONLINEAR_PENDULUM_PRESETS, "double_pendulum": DOUBLE_PENDULUM_PRESETS,
                  "spring_coupled_pendulums": SPRING_COUPLED_PENDULUM_PRESETS}

def load_system(name: str, *, frequency_hz: float | None=None, amplitude: float=0., driver=None, **kwargs):
    """Construct an ODE with a sinusoidal, zero, tabulated, or custom driver."""
    cls=SYSTEM_LIBRARY[name]
    if driver is None: driver=ZeroDriver() if frequency_hz is None else SinusoidalDriver(frequency_hz,amplitude)
    preset_name=kwargs.pop("preset", None)
    if isinstance(preset_name,str):
        preset_values=SYSTEM_PRESETS[name][preset_name].construction_kwargs()
        preset_values.update(kwargs)  # Explicit construction keywords take precedence.
        kwargs=preset_values
    elif preset_name is not None:
        kwargs["preset"]=preset_name
    return cls(driver=driver,**kwargs)

def resolve_archived_system(system_name: str, parameters: dict | None=None):
    """Recreate an ODE library object from dynsyslab archive provenance.

    This is deliberately a library-side operation: archives retain a portable
    name and parameter snapshot, and consumers may opt in to reconstruction.
    """
    parameters=dict(parameters or {})
    key=next((item for item,cls in ODE_LIBRARY.items() if item==system_name or getattr(cls,"name",None)==system_name),None)
    if key is None:
        # Names are generally instance fields, so compare against a default
        # construction only as a compatibility fallback for older archives.
        for item,cls in ODE_LIBRARY.items():
            try:
                if cls(driver=ZeroDriver()).name==system_name: key=item; break
            except Exception as error:
                import warnings
                warnings.warn(f"Could not construct {item!r} while resolving archived system {system_name!r}: {error}",RuntimeWarning,stacklevel=2)
    if key is None: raise KeyError(f"No dynsyslib ODE is registered for archived system '{system_name}'.")
    parameters.pop("driver",None); parameters.pop("name",None); parameters.pop("_dynsyslab_class",None); parameters.pop("_dynsyslab_module",None)
    # Archives predating the declared pressure convention were produced by
    # the positive-acceleration implementation. Preserve that behavior when
    # reconstructing their executable system/KROM lineage.
    if key=="rpe" and "pressure_input_convention" not in parameters:
        parameters["pressure_input_convention"]="legacy-positive-acceleration/v1"
    preset=parameters.pop("preset",None)
    if isinstance(preset,dict) and preset.get("key") in SYSTEM_PRESETS.get(key,{}):
        base=SYSTEM_PRESETS[key][preset["key"]].construction_kwargs(); base.update({name:value for name,value in parameters.items() if name not in preset})
        if isinstance(base.get("preset"),str): return load_system(key,driver=ZeroDriver(),**base)
        return SYSTEM_LIBRARY[key](driver=ZeroDriver(),**base)
    reconstructed_bubble_preset=False
    if isinstance(preset,dict) and key in {"marmottant","reddy_szeri"}:
        # These bubble variants have system presets of their own but retain a
        # nested physical BubblePreset in solution provenance.
        fields={name:preset[name] for name in ("key","description","temperature_c","rho_kg_m3",
            "mu_pa_s","sigma_n_m","vapor_pressure_pa","gas","kappa","caution") if name in preset}
        parameters["preset"]=BubblePreset(**fields); reconstructed_bubble_preset=True
    if preset is not None and not isinstance(preset,str) and not reconstructed_bubble_preset: parameters["preset"]=preset
    elif isinstance(preset,str): parameters["preset"]=preset
    if isinstance(parameters.get("initial_state_override"),list): parameters["initial_state_override"]=tuple(parameters["initial_state_override"])
    return load_system(key,driver=ZeroDriver(),**parameters)
