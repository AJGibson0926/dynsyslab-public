"""Export a fitted field in the signed coordinates used by replay and Viewer."""
import numpy as np

from .field_codec import encode
from .integration import affine_equations, polynomial_arrays
from .models import SparseField


def modal_equations(model):
    """Return JSON-compatible equations, without refitting or sampling the field."""
    chart = model.chart
    context = dict(coordinates="signed", lags=chart.lags.tolist(), provenance=model.metadata)
    if model.space == "direct" and not model.feature_asinh and isinstance(model.field, SparseField):
        if model.field.spec["family"] == "polynomial":
            powers, coefficients = polynomial_arrays(model.field)
            degree = powers.sum(axis=1)
            if not np.any(coefficients[:, degree > 1]):
                A, bias = affine_equations(model.field)
                return dict(schema="dynsyslab.miscellaneous-computes/affine/v1", **context,
                            A=A.tolist(), bias=(bias / chart.gain).tolist())
            # The regression uses b = gain * xi and fits db/dt. Convert each active monomial.
            rows = []
            for row in coefficients:
                rows.append([dict(exponents=powers[j].tolist(),
                                  coefficient=float(row[j] * chart.gain ** (int(degree[j]) - 1)))
                             for j in np.flatnonzero(row)])
            return dict(schema="dynsyslab.miscellaneous-computes/polynomial-equations/v1",
                        **context, rows=rows)

    # The same evaluator specification used by ModalModel.save describes composite and local laws.
    arrays = {}
    tree = encode(model.field, arrays)
    if model.space == "decoded":
        rhs = "xi_dot = field(xi @ decoder) @ encoder"
    elif model.feature_asinh:
        rhs = "xi_dot = field(asinh(gain * xi / feature_asinh)) / gain"
    else:
        rhs = "xi_dot = field(gain * xi) / gain"
    return dict(schema="dynsyslab.miscellaneous-computes/field-equations/v1", **context,
                row_vector_equation=rhs, field_space=model.space, gain=chart.gain,
                feature_asinh=model.feature_asinh, encoder=chart.encoder.tolist(),
                decoder=chart.decoder.tolist(), channel_scales=chart.scales.tolist(),
                field=tree, arrays={k: v.tolist() for k, v in arrays.items()})
