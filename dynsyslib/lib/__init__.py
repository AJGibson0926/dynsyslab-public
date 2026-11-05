"""Individual ODE system definitions for dynsyslib."""
from .presets import SystemPreset
from .rpe import BUBBLE_PRESETS, BubblePreset, RayleighPlessetSystem
from .duffing import DUFFING_PRESETS, DuffingSystem
from .marmottant import MARMOTTANT_PRESETS, MarmottantSystem
from .reddy_szeri import REDDY_SZERI_PRESETS, ReddySzeriSystem
from .van_der_pol import VAN_DER_POL_PRESETS, VanDerPolSystem
from .rossler import ROSSLER_PRESETS, RosslerSystem
from .lorenz import LORENZ_PRESETS, LorenzSystem
from .chua import CHUA_PRESETS, ChuaCircuit
from .henon_heiles import HENON_HEILES_PRESETS, HenonHeilesSystem
from .hindmarsh_rose import HINDMARSH_ROSE_PRESETS, HindmarshRoseSystem
from .hodgkin_huxley import HODGKIN_HUXLEY_PRESETS, HodgkinHuxleySystem
from .fput import FPUT_PRESETS, FPUTSystem
from .kuramoto import KURAMOTO_PRESETS, KuramotoSystem
from .nonlinear_pendulum import NONLINEAR_PENDULUM_PRESETS, NonlinearPendulumSystem
from .double_pendulum import DOUBLE_PENDULUM_PRESETS, DoublePendulumSystem
from .spring_coupled_pendulums import SPRING_COUPLED_PENDULUM_PRESETS, SpringCoupledPendulumsSystem

__all__=["SystemPreset","BUBBLE_PRESETS","BubblePreset","RayleighPlessetSystem","DUFFING_PRESETS","DuffingSystem",
         "MARMOTTANT_PRESETS","MarmottantSystem","REDDY_SZERI_PRESETS","ReddySzeriSystem","VAN_DER_POL_PRESETS","VanDerPolSystem",
         "RosslerSystem","LorenzSystem","ChuaCircuit","HenonHeilesSystem","HindmarshRoseSystem",
         "HodgkinHuxleySystem","FPUTSystem","KuramotoSystem","NonlinearPendulumSystem","DoublePendulumSystem","SpringCoupledPendulumsSystem","ROSSLER_PRESETS","LORENZ_PRESETS","CHUA_PRESETS","HENON_HEILES_PRESETS","HINDMARSH_ROSE_PRESETS","HODGKIN_HUXLEY_PRESETS","FPUT_PRESETS","KURAMOTO_PRESETS","NONLINEAR_PENDULUM_PRESETS","DOUBLE_PENDULUM_PRESETS","SPRING_COUPLED_PENDULUM_PRESETS"]
