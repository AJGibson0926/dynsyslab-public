"""Deterministic bounded candidate scheduling policy."""
from __future__ import annotations

from itertools import combinations
from typing import Any

from rom import (KROMCandidate, polynomial_bilinear_closure_quality,
                 rank_krom_candidates)
from .common import INITIAL_STATE, TARGET_STATE

def _candidate_schedule(parent, *, maximum_mode_count: int) -> tuple[list[KROMCandidate], dict[str, Any]]:
    all_selections = [
        selection
        for size in range(1, min(maximum_mode_count, parent.lifted_dimension) + 1)
        for selection in combinations(range(parent.lifted_dimension), size)
    ]
    closure = polynomial_bilinear_closure_quality(parent, all_selections)
    ranked = rank_krom_candidates(
        parent,
        all_selections,
        INITIAL_STATE,
        TARGET_STATE,
        structural_closure_quality=closure,
    )

    # The first two entries are calibrated historical anchors.  The remaining
    # schedule interleaves composite-rank leaders, closure leaders, and each
    # mode count so one numerical criterion cannot monopolize the 25 attempts.
    by_selection = {item.retained_indices: item for item in ranked}
    ordered: list[KROMCandidate] = []

    def add(selection: tuple[int, ...]) -> None:
        item = by_selection.get(tuple(sorted(selection)))
        if item is not None and item not in ordered:
            ordered.append(item)

    add((6,))
    add((6, 8))
    # Attempt 3 deliberately probes whether one more closure-supporting mode
    # helps the calibrated pair.  Subsequent choices encode the observed
    # practical gate: worst normalized controllability below 0.01 makes even
    # exact continuous-state CARE evaluation intractable for this target.
    add((6, 7, 8))
    for selection in (
        (5, 8), (6, 9), (5, 6), (6, 7), (9,), (8,), (8, 9), (7, 8),
        (0,), (1,), (5,), (2,), (7,), (4,), (3,),
    ):
        add(selection)
    # After the exhaustive one-mode screen, deliberately revisit the four
    # exact-successful coordinate sets with the newly introduced conjugate-
    # complete (real-split general) bilinear candidate.  These are new fitted
    # control experiments, so they consume explicit campaign attempts.
    for selection in ((6,), (9,), (6, 8), (8, 9)):
        base=by_selection[tuple(sorted(selection))]
        ordered.append(KROMCandidate(
            base.retained_indices,
            {**base.scores,"conjugate_complete_candidate_stage":1.0},
            "repeat exact-successful KROM for conjugate-complete bilinear map evaluation",
        ))
    base=by_selection[(6,)]
    ordered.append(KROMCandidate(
        base.retained_indices,
        {**base.scores,"path_conditioned_candidate_stage":1.0},
        "repeat exact-successful KROM for deterministic successful-path and path-tube map evaluation",
    ))
    base=by_selection[(8, 9)]
    ordered.append(KROMCandidate(
        base.retained_indices,
        {
            **base.scores,
            "path_conditioned_candidate_stage":1.0,
            "path_arclength_collocation_stage":1.0,
        },
        "final exact-successful KROM selected by bounded pathwise error screening for arc-length collocation",
    ))
    add((4, 8))
    for item in ranked:
        if (
            item.scores["worst_normalized_controllability"] >= 0.01
            and len(item.retained_indices) <= 3
        ):
            add(item.retained_indices)
    closure_by_size: dict[int, list[KROMCandidate]] = {}
    for size in range(1, min(maximum_mode_count, parent.lifted_dimension) + 1):
        sized = [item for item in ranked if len(item.retained_indices) == size]
        closure_by_size[size] = sorted(
            sized,
            key=lambda candidate: (
                -candidate.scores["structural_closure_quality"],
                -candidate.scores["lifted_target_displacement"],
            ),
        )
    # Round-robin across dimensions keeps the bounded search exploratory: the
    # first 25 contain closure leaders from small and large subspaces instead
    # of exhausting all pairs/triples before a five- or six-mode test.
    for rank_within_size in range(3):
        for size in sorted(closure_by_size):
            if rank_within_size < len(closure_by_size[size]):
                add(closure_by_size[size][rank_within_size].retained_indices)
    for item in ranked:
        add(item.retained_indices)
    audit = {
        "screening_method": "virtual exact-map metrics plus analytical polynomial closure",
        "candidate_count": len(ranked),
        "maximum_mode_count": maximum_mode_count,
        "schedule": [item.summary() for item in ordered],
        "closure_quality_range": [float(min(closure.values())), float(max(closure.values()))],
        "historical_calibration_candidates": [[6], [6, 8]],
        "attempt_accounting": (
            "virtual screening does not extract an archive, fit a control map, or run control; "
            "only evaluator invocations are counted attempts"
        ),
    }
    return ordered, audit
