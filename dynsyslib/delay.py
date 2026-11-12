"""Delay-system library kept separate from the ODE tangent campaign backend."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .interfaces import DelaySystem

@dataclass(frozen=True)
class MackeyGlassSystem(DelaySystem):
    """x' = beta*x(t-tau)/(1+x(t-tau)^n) - gamma*x(t), solved by method of steps."""
    beta: float=.2; gamma: float=.1; delay: float=17.; exponent: float=10.; history_value: float=1.2; name: str="mackey-glass"
    def integrate(self, duration: float, dt: float=.05, checkpoint_steps: int|None=None, on_checkpoint=None):
        if dt <= 0 or self.delay < dt: raise ValueError("Require dt > 0 and delay >= dt.")
        steps, lag = int(round(duration/dt)), int(round(self.delay/dt))
        x = np.full(steps+lag+1, self.history_value, dtype=float)
        for i in range(lag, steps+lag):
            xd = x[i-lag]
            x[i+1] = x[i] + dt*(self.beta*xd/(1+xd**self.exponent)-self.gamma*x[i])
            done=i-lag+1
            if on_checkpoint and checkpoint_steps and (done % checkpoint_steps == 0 or done == steps):
                start=max(0,done-checkpoint_steps)
                on_checkpoint(np.arange(start,done+1)*dt, x[lag+start:lag+done+1])
        t = np.arange(steps+1)*dt
        return t, x[lag:]
