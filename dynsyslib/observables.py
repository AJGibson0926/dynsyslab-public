"""Physical/model observables used at Koopman method boundaries.

An :class:`Observable` is a coordinate map ``y = h(x)``.  It is deliberately
different from a SINDy/EDMD feature dictionary: an observable decides which
quantities are modeled, while a feature dictionary decides how functions of
those quantities are represented.
"""
from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from typing import Any, Callable, Sequence

import numpy as np


ArrayFunction = Callable[[np.ndarray], np.ndarray]


@dataclass(frozen=True)
class Observable:
    """A named state-to-output map with an explicit Jacobian contract."""

    name: str
    output_dimension: int
    state_dimension: int
    function: ArrayFunction
    jacobian_function: ArrayFunction | None = None
    descriptor: dict[str, Any] | None = None
    validation_state: tuple[float, ...] | None = None
    value_type: str = "real"
    derivative_observables: tuple[tuple["Observable", ...], ...] | None = None

    def __post_init__(self):
        if not str(self.name).strip():
            raise ValueError("Observable name must be nonempty.")
        if int(self.output_dimension) != self.output_dimension or self.output_dimension < 1:
            raise ValueError("Observable output_dimension must be a positive integer.")
        if int(self.state_dimension) != self.state_dimension or self.state_dimension < 1:
            raise ValueError("Observable state_dimension must be a positive integer.")
        object.__setattr__(self, "output_dimension", int(self.output_dimension))
        object.__setattr__(self, "state_dimension", int(self.state_dimension))
        if self.value_type not in {"real", "complex"}:
            raise ValueError("Observable value_type must be 'real' or 'complex'.")
        object.__setattr__(self, "descriptor", dict(self.descriptor or {"kind": "callable"}))
        probe = (np.zeros(self.state_dimension, dtype=float) if self.validation_state is None
                 else np.asarray(self.validation_state, dtype=float).reshape(-1))
        if probe.size != self.state_dimension or np.any(~np.isfinite(probe)):
            raise ValueError("Observable validation_state has the wrong dimension or is non-finite.")
        dtype=complex if self.value_type=="complex" else float
        value = np.asarray(self.function(probe), dtype=dtype).reshape(-1)
        if value.size != self.output_dimension or np.any(~np.isfinite(value)):
            raise ValueError("Observable function returned an invalid probe value.")
        if self.derivative_observables is not None:
            derivatives=tuple(tuple(row) for row in self.derivative_observables)
            if len(derivatives)!=self.output_dimension or any(len(row)!=self.state_dimension for row in derivatives):
                raise ValueError("Observable derivative_observables must form an output-by-state matrix.")
            if any(not isinstance(item,Observable) or item.output_dimension!=1 or item.state_dimension!=self.state_dimension
                   for row in derivatives for item in row):
                raise ValueError("Each stashed observable partial derivative must be scalar and state aligned.")
            object.__setattr__(self,"derivative_observables",derivatives)
        if self.jacobian_function is not None:
            jacobian = np.asarray(self.jacobian_function(probe), dtype=dtype)
            if jacobian.shape != (self.output_dimension, self.state_dimension):
                raise ValueError("Observable Jacobian has the wrong shape.")
            if np.any(~np.isfinite(jacobian)):
                raise ValueError("Observable Jacobian probe contains non-finite values.")
        elif self.derivative_observables is not None:
            jacobian=np.asarray([[item.evaluate(probe).item() for item in row]
                                 for row in self.derivative_observables],dtype=dtype)
            if jacobian.shape!=(self.output_dimension,self.state_dimension) or np.any(~np.isfinite(jacobian)):
                raise ValueError("Observable derivative observables produced an invalid Jacobian probe.")

    @property
    def has_analytic_jacobian(self) -> bool:
        return self.jacobian_function is not None or self.derivative_observables is not None

    def evaluate(self, states: np.ndarray) -> np.ndarray:
        """Evaluate one state or a row-major batch of states."""
        values = np.asarray(states, dtype=float)
        one = values.ndim == 1
        if one:
            values = values.reshape(1, -1)
        if values.ndim != 2 or values.shape[1] != self.state_dimension:
            raise ValueError(
                f"Observable states must have shape (samples, {self.state_dimension})."
            )
        dtype=complex if self.value_type=="complex" else float
        result = np.vstack([np.asarray(self.function(row), dtype=dtype).reshape(-1) for row in values])
        if result.shape != (values.shape[0], self.output_dimension):
            raise ValueError("Observable function returned an inconsistent output shape.")
        if np.any(~np.isfinite(result)):
            raise ValueError("Observable evaluation produced non-finite values.")
        return result[0] if one else result

    def jacobian(
        self,
        state: np.ndarray,
        *,
        backend: str = "analytic",
        relative_step: float = 1e-6,
    ) -> np.ndarray:
        """Return ``J_h(x)`` analytically or by explicit finite difference."""
        x = np.asarray(state, dtype=float).reshape(-1)
        if x.size != self.state_dimension or np.any(~np.isfinite(x)):
            raise ValueError("Observable Jacobian state has the wrong dimension or is non-finite.")
        if backend == "analytic":
            if not self.has_analytic_jacobian:
                raise NotImplementedError(
                    "This observable has no analytic Jacobian; explicitly request "
                    "backend='finite_difference' to approximate it."
                )
            dtype=complex if self.value_type=="complex" else float
            if self.jacobian_function is not None:
                result = np.asarray(self.jacobian_function(x), dtype=dtype)
            else:
                result=np.asarray([[item.evaluate(x).item() for item in row]
                                   for row in self.derivative_observables],dtype=dtype)
        elif backend == "finite_difference":
            if relative_step <= 0 or not np.isfinite(relative_step):
                raise ValueError("relative_step must be finite and positive.")
            dtype=complex if self.value_type=="complex" else float
            result = np.empty((self.output_dimension, self.state_dimension), dtype=dtype)
            for index in range(self.state_dimension):
                step = relative_step * max(1.0, abs(x[index]))
                plus = x.copy(); minus = x.copy()
                plus[index] += step; minus[index] -= step
                result[:, index] = (self.evaluate(plus) - self.evaluate(minus)) / (2.0 * step)
        else:
            raise ValueError("backend must be 'analytic' or 'finite_difference'.")
        if result.shape != (self.output_dimension, self.state_dimension):
            raise ValueError("Observable Jacobian has the wrong shape.")
        return result

    def summary(self) -> dict[str, Any]:
        return {
            "encoding": "state-observable/v1",
            "name": self.name,
            "state_dimension": self.state_dimension,
            "output_dimension": self.output_dimension,
            "analytic_jacobian": self.has_analytic_jacobian,
            "value_type": self.value_type,
            "descriptor": dict(self.descriptor),
            "validation_state": None if self.validation_state is None else list(self.validation_state),
            "derivative_observables": None if self.derivative_observables is None else [
                [item.summary() for item in row] for row in self.derivative_observables],
        }

    @classmethod
    def from_expressions(
        cls,
        name: str,
        expressions,
        *,
        state_dimension: int | None = None,
        descriptor: dict[str,Any] | None = None,
        validation_state: tuple[float,...] | None = None,
        derivative_axis_names: Sequence[str] | None = None,
        stash_derivatives: bool = True,
    ) -> "Observable":
        """Build an archive-safe observable and exact analytic partials.

        The supplied expression trees are differentiated symbolically.  The
        resulting scalar partial derivatives are retained as observables, so a
        coordinate Jacobian remains inspectable rather than becoming a hidden
        implementation detail.
        """
        from .analytic_expressions import AnalyticExpression
        values=tuple(expressions)
        if not values or any(not isinstance(item,AnalyticExpression) for item in values):
            raise TypeError("Observable.from_expressions requires analytic expression objects.")
        dimension=int(state_dimension or values[0].state_dimension)
        if any(item.state_dimension!=dimension for item in values):
            raise ValueError("Observable expressions must share one state dimension.")
        axes=tuple(derivative_axis_names or (f"x{index}" for index in range(dimension)))
        if len(axes)!=dimension:
            raise ValueError("derivative_axis_names must match state_dimension.")
        partials=tuple(tuple(expression.derivative(index) for index in range(dimension))
                       for expression in values)
        def evaluate(x,items=values):
            return np.asarray([item.evaluate(x) for item in items],dtype=float)
        def jacobian(x,items=partials):
            return np.asarray([[item.evaluate(x) for item in row] for row in items],dtype=float)
        base_descriptor={**dict(descriptor or {}),"kind":"analytic-expression",
            "expressions":[item.summary() for item in values],"derivative_axis_names":list(axes)}
        derivative_matrix=None
        if stash_derivatives:
            rows=[]
            for output_index,row in enumerate(partials):
                prefix=name if len(values)==1 else f"{name}_{output_index}"
                rows.append(tuple(cls.from_expressions(
                    f"{prefix}_{axes[index]}_derivative",(partial,),state_dimension=dimension,
                    descriptor={"kind":"analytic-expression-derivative","source_observable":name,
                                "output_index":output_index,"state_index":index,"state_axis":axes[index]},
                    validation_state=validation_state,derivative_axis_names=axes,stash_derivatives=False)
                    for index,partial in enumerate(row)))
            derivative_matrix=tuple(rows)
        return cls(name,len(values),dimension,evaluate,jacobian,base_descriptor,
                   validation_state,"real",derivative_matrix)

    def flow_derivative(self, rhs: Callable[[float,np.ndarray],np.ndarray], *,
                        name: str | None=None, time: float=0.) -> "Observable":
        """Return the material derivative ``dh/dt = J_h(x) f(t,x)``."""
        source=self
        return Observable(name or f"{self.name}_time_derivative",self.output_dimension,
            self.state_dimension,
            lambda x:source.jacobian(x)@np.asarray(rhs(float(time),np.asarray(x,dtype=float)),dtype=float),
            descriptor={"kind":"observable-flow-derivative","source_observable":self.summary(),
                        "evaluation_time":float(time)},validation_state=self.validation_state,
            value_type=self.value_type)

    @classmethod
    def identity(cls, state_dimension: int, *, name: str = "state") -> "Observable":
        dimension = int(state_dimension)
        return cls(
            name,
            dimension,
            dimension,
            lambda x: np.asarray(x, dtype=float),
            lambda x: np.eye(dimension, dtype=float),
            {"kind": "identity", "state_indices": list(range(dimension))},
        )

    @classmethod
    def component(
        cls, state_dimension: int, index: int, *, name: str | None = None
    ) -> "Observable":
        dimension = int(state_dimension); component_index = int(index)
        if not 0 <= component_index < dimension:
            raise ValueError("Observable component index is outside the state dimension.")
        row = np.zeros((1, dimension), dtype=float); row[0, component_index] = 1.0
        return cls(
            name or f"x[{component_index}]",
            1,
            dimension,
            lambda x, i=component_index: np.asarray((x[i],), dtype=float),
            lambda x, value=row: value,
            {"kind": "state-component", "state_index": component_index},
        )

    @classmethod
    def stack(cls, observables: Sequence["Observable"], *, name: str = "stack") -> "Observable":
        items = tuple(observables)
        if not items:
            raise ValueError("Observable.stack requires at least one observable.")
        state_dimension = items[0].state_dimension
        if any(item.state_dimension != state_dimension for item in items):
            raise ValueError("Stacked observables must share one state dimension.")
        output_dimension = sum(item.output_dimension for item in items)

        def evaluate(x):
            return np.concatenate([np.asarray(item.evaluate(x)).reshape(-1) for item in items])

        jacobian = None
        if all(item.has_analytic_jacobian for item in items):
            jacobian = lambda x: np.vstack([item.jacobian(x) for item in items])
        derivative_rows=None
        if all(item.derivative_observables is not None for item in items):
            derivative_rows=tuple(row for item in items for row in item.derivative_observables)
        return cls(
            name,
            output_dimension,
            state_dimension,
            evaluate,
            jacobian,
            {"kind": "stack", "members": [item.summary() for item in items]},
            validation_state=items[0].validation_state,
            value_type="complex" if any(item.value_type=="complex" for item in items) else "real",
            derivative_observables=derivative_rows,
        )


class ObservableCatalog(Mapping[str, Observable]):
    """Immutable, alias-aware observables owned by one dynamical system.

    Canonical keys are stable API names suitable for configuration files and
    Viewer requests. Human labels, units, provenance, and coordinate-system
    expectations live in each observable's descriptor.
    """

    def __init__(self, observables: Mapping[str, Observable], *, aliases: Mapping[str, str] | None=None):
        items={str(key):value for key,value in observables.items()}
        if not items:
            raise ValueError("ObservableCatalog requires at least one observable.")
        if any(not key.strip() for key in items):
            raise ValueError("Observable catalog keys must be nonempty.")
        if not all(isinstance(value,Observable) for value in items.values()):
            raise TypeError("ObservableCatalog values must be Observable objects.")
        dimensions={value.state_dimension for value in items.values()}
        if len(dimensions)!=1:
            raise ValueError("All catalog observables must share one state dimension.")
        normalized={str(alias):str(target) for alias,target in dict(aliases or {}).items()}
        missing={target for target in normalized.values() if target not in items}
        if missing:
            raise KeyError(f"Observable aliases reference missing keys: {sorted(missing)}")
        self._observables=items
        self._aliases=normalized

    def __getitem__(self, key: str) -> Observable:
        requested=str(key)
        return self._observables[self._aliases.get(requested,requested)]

    def __iter__(self) -> Iterator[str]:
        return iter(self._observables)

    def __len__(self) -> int:
        return len(self._observables)

    @property
    def aliases(self) -> dict[str,str]:
        return dict(self._aliases)

    def summary(self) -> dict[str,Any]:
        return {"encoding":"observable-catalog/v1","observables":{
            key:value.summary() for key,value in self._observables.items()},"aliases":dict(self._aliases)}


def state_observable_catalog(
    state_dimension: int,
    *,
    component_names: Sequence[str] | None=None,
    component_labels: Sequence[str] | None=None,
    coordinate_system: str="native",
    extras: Mapping[str,Observable] | None=None,
    aliases: Mapping[str,str] | None=None,
) -> ObservableCatalog:
    """Build the guaranteed full-state and state-component catalog entries."""
    dimension=int(state_dimension)
    names=tuple(component_names or (f"x{index}" for index in range(dimension)))
    labels=tuple(component_labels or names)
    if len(names)!=dimension or len(labels)!=dimension:
        raise ValueError("State observable names and labels must match state_dimension.")
    items: dict[str,Observable]={}
    components=[]
    for index,(name,label) in enumerate(zip(names,labels)):
        base=Observable.component(dimension,index,name=name)
        item=Observable(base.name,base.output_dimension,base.state_dimension,base.function,
                        base.jacobian_function,{**base.descriptor,"label":label,
                        "state_coordinate_system":coordinate_system})
        items[name]=item; components.append(item)
    items["state"]=Observable.stack(components,name="state")
    items["state"]=Observable(items["state"].name,items["state"].output_dimension,
        items["state"].state_dimension,items["state"].function,items["state"].jacobian_function,
        {**items["state"].descriptor,"label":"Full state","state_coordinate_system":coordinate_system})
    for key,value in dict(extras or {}).items():
        if key in items:
            raise KeyError(f"Observable extra {key!r} collides with a state observable.")
        items[str(key)]=value
    built_aliases={f"x[{index}]":name for index,name in enumerate(names)}
    built_aliases.update({f"state_{index}":name for index,name in enumerate(names)})
    built_aliases.update(dict(aliases or {}))
    return ObservableCatalog(items,aliases=built_aliases)


def coerce_observable(
    value: Observable | int | str | None,
    state_dimension: int,
    *,
    system: Any | None = None,
    coordinate_system: str | None = None,
) -> Observable:
    """Normalize an object, component index, or named system observable.

    String selection is deliberately resolved by the owning dynamical system,
    so names such as ``far_field_acoustic_echo`` and ``total_energy`` never
    become global aliases with ambiguous physical meaning.
    """
    if value is None:
        return Observable.identity(state_dimension)
    if isinstance(value, Observable):
        if value.state_dimension != state_dimension:
            raise ValueError("Observable state dimension does not match the dynamical system.")
        return value
    if isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_)):
        return Observable.component(state_dimension, int(value))
    if isinstance(value, str):
        if system is None or not hasattr(system, "observable"):
            raise ValueError(
                "Named observables require their dynamical system; pass system=... "
                "or archive metadata from which DynSysLab can reconstruct it."
            )
        resolved = system.observable(value, coordinate_system=coordinate_system)
        if resolved.state_dimension != state_dimension:
            raise ValueError("Named observable state dimension does not match the H-DMD input.")
        return resolved
    raise TypeError(
        "observable must be an Observable, a system observable name, a "
        "state-component index, or None."
    )


def coerce_coordinate_observable(
    value: Observable | int | str | Sequence[Observable | int | str] | None,
    state_dimension: int,
    *,
    system: Any | None = None,
    coordinate_system: str | None = None,
    name: str | None = None,
) -> Observable:
    """Normalize an eigenfunction coordinate map ``y=h(x)``.

    Unlike the scalar H-DMD boundary, SINDy coordinates may be a named system
    preset or an ordered collection of arbitrary scalar/vector observables.
    """
    if isinstance(value,str) and system is not None and hasattr(system,"observable_coordinate_presets"):
        preset=dict(system.observable_coordinate_presets(coordinate_system)).get(value)
        if preset is not None:
            value=preset
            if name is None: name=str(value if isinstance(value,str) else "_".join(map(str,preset)))
    if isinstance(value,(Observable,str,int,np.integer)) or value is None:
        return coerce_observable(value,state_dimension,system=system,coordinate_system=coordinate_system)
    if isinstance(value,Sequence) and not isinstance(value,(bytes,bytearray)):
        items=tuple(coerce_observable(item,state_dimension,system=system,
            coordinate_system=coordinate_system) for item in value)
        return Observable.stack(items,name=name or "observable_coordinates["+",".join(item.name for item in items)+"]")
    raise TypeError("coordinate_observable must be one observable, a named coordinate preset, or an ordered observable sequence.")


def observable_from_summary(summary: dict[str, Any]) -> Observable:
    """Reconstruct a serializable built-in observable descriptor."""
    payload=dict(summary or {}); descriptor=dict(payload.get("descriptor") or payload)
    kind=descriptor.get("kind")
    state_dimension=int(payload.get("state_dimension") or descriptor.get("state_dimension") or 0)
    name=str(payload.get("name") or "observable")
    if kind=="identity": return Observable.identity(state_dimension,name=name)
    if kind in {"state-component","state_index"}:
        index=int(descriptor.get("state_index",descriptor.get("index")))
        return Observable.component(state_dimension,index,name=name)
    if kind=="stack":
        members=[observable_from_summary(item) for item in descriptor.get("members",())]
        return Observable.stack(members,name=name)
    if kind=="analytic-expression":
        from .analytic_expressions import AnalyticExpression
        expressions=tuple(AnalyticExpression.from_summary(item) for item in descriptor.get("expressions",()))
        if not expressions: raise ValueError("Archived analytic observable has no expressions.")
        axis_names=descriptor.get("derivative_axis_names")
        retained={key:value for key,value in descriptor.items() if key not in {"kind","expressions"}}
        validation=payload.get("validation_state")
        return Observable.from_expressions(name,expressions,state_dimension=state_dimension,
            descriptor=retained,validation_state=None if validation is None else tuple(validation),
            derivative_axis_names=axis_names)
    raise ValueError(f"Observable descriptor kind {kind!r} is not archive-reconstructable.")


__all__ = ["Observable", "ObservableCatalog", "state_observable_catalog",
           "coerce_observable", "coerce_coordinate_observable", "observable_from_summary"]
