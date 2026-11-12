"""Small serializable expression trees for analytic state observables.

The expression layer is intentionally narrower than a computer-algebra system:
it represents the scalar operations used by DynSysLab model observables, can
differentiate them exactly, and can be archived without executable code.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import cos as _cos, exp as _exp, log as _log, sin as _sin, tanh as _tanh
from typing import Any

import numpy as np


@dataclass(frozen=True)
class AnalyticExpression:
    """A scalar real expression of one dynamical state vector."""

    operation: str
    state_dimension: int
    arguments: tuple["AnalyticExpression", ...] = ()
    scalar: float | None = None
    state_index: int | None = None

    def __post_init__(self):
        dimension = int(self.state_dimension)
        if dimension < 1:
            raise ValueError("Analytic-expression state_dimension must be positive.")
        object.__setattr__(self, "state_dimension", dimension)
        if any(argument.state_dimension != dimension for argument in self.arguments):
            raise ValueError("Analytic-expression arguments must share one state dimension.")
        if self.operation == "constant" and (self.scalar is None or not np.isfinite(self.scalar)):
            raise ValueError("Analytic-expression constants must be finite.")
        if self.operation == "state" and (
            self.state_index is None or not 0 <= int(self.state_index) < dimension
        ):
            raise ValueError("Analytic-expression state index is out of range.")
        supported = {"constant", "state", "add", "subtract", "multiply", "divide",
                     "power", "negative", "exp", "log", "sin", "cos", "tanh"}
        if self.operation not in supported:
            raise ValueError(f"Unsupported analytic-expression operation {self.operation!r}.")

    @classmethod
    def constant(cls, value: float, state_dimension: int) -> "AnalyticExpression":
        return cls("constant", state_dimension, scalar=float(value))

    @classmethod
    def state(cls, state_dimension: int, index: int) -> "AnalyticExpression":
        return cls("state", state_dimension, state_index=int(index))

    @property
    def is_constant(self) -> bool:
        return self.operation == "constant"

    def evaluate(self, state: np.ndarray) -> float:
        x = np.asarray(state, dtype=float).reshape(-1)
        if x.size != self.state_dimension:
            raise ValueError("Analytic-expression state has the wrong dimension.")
        op = self.operation
        if op == "constant": return float(self.scalar)
        if op == "state": return float(x[int(self.state_index)])
        values = tuple(argument.evaluate(x) for argument in self.arguments)
        if op == "add": return values[0] + values[1]
        if op == "subtract": return values[0] - values[1]
        if op == "multiply": return values[0] * values[1]
        if op == "divide": return values[0] / values[1]
        if op == "power": return values[0] ** float(self.scalar)
        if op == "negative": return -values[0]
        if op == "exp": return _exp(values[0])
        if op == "log": return _log(values[0])
        if op == "sin": return _sin(values[0])
        if op == "cos": return _cos(values[0])
        if op == "tanh": return _tanh(values[0])
        raise RuntimeError("Unreachable analytic-expression operation.")

    def derivative(self, state_index: int) -> "AnalyticExpression":
        """Return the exact symbolic partial derivative for one state axis."""
        index = int(state_index)
        if not 0 <= index < self.state_dimension:
            raise IndexError("Analytic-expression derivative index is out of range.")
        zero = self._constant(0.0); one = self._constant(1.0)
        op = self.operation
        if op == "constant": return zero
        if op == "state": return one if self.state_index == index else zero
        a = self.arguments[0]; da = a.derivative(index)
        if op == "negative": return -da
        if op == "exp": return exp(a) * da
        if op == "log": return da / a
        if op == "sin": return cos(a) * da
        if op == "cos": return -sin(a) * da
        if op == "tanh": return (one - tanh(a) ** 2) * da
        if op == "power":
            exponent = float(self.scalar)
            return self._constant(exponent) * (a ** (exponent - 1.0)) * da
        b = self.arguments[1]; db = b.derivative(index)
        if op == "add": return da + db
        if op == "subtract": return da - db
        if op == "multiply": return da * b + a * db
        if op == "divide": return (da * b - a * db) / (b ** 2)
        raise RuntimeError("Unreachable analytic-expression derivative.")

    def summary(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "encoding": "analytic-state-expression/v1",
            "operation": self.operation,
            "state_dimension": self.state_dimension,
        }
        if self.arguments: payload["arguments"] = [argument.summary() for argument in self.arguments]
        if self.scalar is not None: payload["scalar"] = float(self.scalar)
        if self.state_index is not None: payload["state_index"] = int(self.state_index)
        return payload

    @classmethod
    def from_summary(cls, summary: dict[str, Any]) -> "AnalyticExpression":
        payload = dict(summary)
        if payload.get("encoding") != "analytic-state-expression/v1":
            raise ValueError("Unsupported analytic-expression encoding.")
        return cls(
            str(payload["operation"]),
            int(payload["state_dimension"]),
            tuple(cls.from_summary(item) for item in payload.get("arguments", ())),
            None if "scalar" not in payload else float(payload["scalar"]),
            None if "state_index" not in payload else int(payload["state_index"]),
        )

    def _constant(self, value: float) -> "AnalyticExpression":
        return AnalyticExpression.constant(value, self.state_dimension)

    def _coerce(self, value: float | "AnalyticExpression") -> "AnalyticExpression":
        if isinstance(value, AnalyticExpression):
            if value.state_dimension != self.state_dimension:
                raise ValueError("Cannot combine expressions with different state dimensions.")
            return value
        return self._constant(float(value))

    def _binary(self, operation: str, other, *, reverse: bool = False):
        right = self._coerce(other); left = self
        if reverse: left, right = right, left
        if operation == "add":
            if left.is_constant and left.scalar == 0: return right
            if right.is_constant and right.scalar == 0: return left
        if operation == "subtract" and right.is_constant and right.scalar == 0: return left
        if operation == "multiply":
            if (left.is_constant and left.scalar == 0) or (right.is_constant and right.scalar == 0):
                return self._constant(0)
            if left.is_constant and left.scalar == 1: return right
            if right.is_constant and right.scalar == 1: return left
        if operation == "divide" and right.is_constant and right.scalar == 1: return left
        if left.is_constant and right.is_constant:
            if operation == "add": result=left.scalar+right.scalar
            elif operation == "subtract": result=left.scalar-right.scalar
            elif operation == "multiply": result=left.scalar*right.scalar
            else: result=left.scalar/right.scalar
            return self._constant(result)
        return AnalyticExpression(operation, self.state_dimension, (left, right))

    def __add__(self, other): return self._binary("add", other)
    def __radd__(self, other): return self._binary("add", other, reverse=True)
    def __sub__(self, other): return self._binary("subtract", other)
    def __rsub__(self, other): return self._binary("subtract", other, reverse=True)
    def __mul__(self, other): return self._binary("multiply", other)
    def __rmul__(self, other): return self._binary("multiply", other, reverse=True)
    def __truediv__(self, other): return self._binary("divide", other)
    def __rtruediv__(self, other): return self._binary("divide", other, reverse=True)

    def __neg__(self):
        if self.is_constant: return self._constant(-float(self.scalar))
        if self.operation == "negative": return self.arguments[0]
        return AnalyticExpression("negative", self.state_dimension, (self,))

    def __pow__(self, exponent: float):
        value = float(exponent)
        if not np.isfinite(value): raise ValueError("Analytic-expression powers must be finite constants.")
        if value == 0: return self._constant(1)
        if value == 1: return self
        if self.is_constant: return self._constant(float(self.scalar) ** value)
        return AnalyticExpression("power", self.state_dimension, (self,), scalar=value)


def state_variables(state_dimension: int) -> tuple[AnalyticExpression, ...]:
    return tuple(AnalyticExpression.state(state_dimension, index) for index in range(int(state_dimension)))


def _unary(operation: str, value: AnalyticExpression) -> AnalyticExpression:
    if not isinstance(value, AnalyticExpression):
        raise TypeError("Analytic unary functions require an AnalyticExpression.")
    if value.is_constant:
        functions = {"exp": _exp, "log": _log, "sin": _sin, "cos": _cos, "tanh": _tanh}
        return AnalyticExpression.constant(functions[operation](float(value.scalar)), value.state_dimension)
    return AnalyticExpression(operation, value.state_dimension, (value,))


def exp(value: AnalyticExpression) -> AnalyticExpression: return _unary("exp", value)
def log(value: AnalyticExpression) -> AnalyticExpression: return _unary("log", value)
def sin(value: AnalyticExpression) -> AnalyticExpression: return _unary("sin", value)
def cos(value: AnalyticExpression) -> AnalyticExpression: return _unary("cos", value)
def tanh(value: AnalyticExpression) -> AnalyticExpression: return _unary("tanh", value)
def sqrt(value: AnalyticExpression) -> AnalyticExpression: return value ** 0.5


__all__ = ["AnalyticExpression", "state_variables", "exp", "log", "sin", "cos", "tanh", "sqrt"]
