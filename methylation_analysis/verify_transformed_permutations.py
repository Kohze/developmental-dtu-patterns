"""Independently verify exhaustive stage-randomisation p-values.

The legacy output column and filename use ``exact`` to preserve the released
machine schema. Exactness is conditional on exchangeable stage labels; this
script does not establish valid time-series inference for ordered development.

This check deliberately does not import ``replicate_robustness.py``.  It rebuilds
stage means from the released replicate tables and uses a residual-maker matrix
or a first-difference matrix to transform every one of the 8! stage
permutations.  A non-zero exit status indicates disagreement with the released
full-family audit.
"""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
STAGES = ["10.5", "11.5", "12.5", "13.5", "14.5", "15.5", "16.5", "0"]


def stage_profiles(frame: pd.DataFrame, keys: list[str], value: str) -> dict[tuple, np.ndarray]:
    table = frame.pivot_table(index=keys, columns="stage", values=value, aggfunc="mean")
    table = table.reindex(columns=STAGES)
    return {
        key if isinstance(key, tuple) else (key,): row.to_numpy(dtype=float)
        for key, row in table.iterrows()
    }


def transform_matrix(kind: str, length: int) -> np.ndarray:
    if kind == "first_difference":
        matrix = np.zeros((length, length - 1), dtype=float)
        for column in range(length - 1):
            matrix[column, column] = -1.0
            matrix[column + 1, column] = 1.0
        return matrix
    if kind == "detrended":
        design = np.column_stack([np.ones(length), np.arange(length, dtype=float)])
        return np.eye(length) - design @ np.linalg.inv(design.T @ design) @ design.T
    raise ValueError(f"Unknown transform: {kind}")


def row_correlations(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    left = left - left.mean(axis=1, keepdims=True)
    right = right - right.mean()
    denominator = np.sqrt(np.sum(left**2, axis=1) * np.sum(right**2))
    return np.divide(
        left @ right,
        denominator,
        out=np.full(left.shape[0], np.nan, dtype=float),
        where=denominator > 0,
    )


def exact_p(x: np.ndarray, y: np.ndarray, permutations: np.ndarray, kind: str) -> float:
    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        raise ValueError("Verification set unexpectedly contains an incomplete profile")
    operator = transform_matrix(kind, len(x))
    transformed_y = y @ operator
    observed = abs(row_correlations((x @ operator)[None, :], transformed_y)[0])
    permuted = x[permutations] @ operator
    correlations = row_correlations(permuted, transformed_y)
    return float(np.mean(np.abs(correlations) >= observed - 1e-12))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=HERE / "results")
    parser.add_argument("--tolerance", type=float, default=1e-12)
    args = parser.parse_args()
    results = args.results_dir.resolve()

    audit = pd.read_csv(results / "replicate_robustness_all.csv")
    methylation = pd.read_csv(results / "replicate_level_methylation.csv", dtype={"stage": str})
    usage = pd.read_csv(results / "replicate_level_isoform_usage.csv", dtype={"stage": str})
    methyl_profiles = stage_profiles(
        methylation, ["tissue", "gene_id", "region"], "mean_weighted"
    )
    usage_profiles = stage_profiles(usage, ["tissue", "isoform_id"], "rep_if")
    selected = audit.loc[
        (audit["detrended_global_q"] < 0.05) | (audit["difference_global_q"] < 0.05)
    ].copy()
    permutations = np.asarray(list(itertools.permutations(range(8))), dtype=np.int8)

    records: list[dict[str, object]] = []
    for row in selected.itertuples(index=False):
        x = methyl_profiles[(row.tissue, row.gene_id, row.region)]
        y = usage_profiles[(row.tissue, row.isoform_id)]
        for kind, reported_column in (
            ("detrended", "detrended_exact_permutation_p"),
            ("first_difference", "difference_exact_permutation_p"),
        ):
            independent = exact_p(x, y, permutations, kind)
            reported = float(getattr(row, reported_column))
            discrepancy = abs(independent - reported)
            records.append(
                {
                    "gene_id": row.gene_id,
                    "isoform_id": row.isoform_id,
                    "region": row.region,
                    "tissue": row.tissue,
                    "transform": kind,
                    "reported_exact_p": reported,
                    "independent_exact_p": independent,
                    "absolute_discrepancy": discrepancy,
                    "within_tolerance": discrepancy <= args.tolerance,
                }
            )

    report = pd.DataFrame(records)
    output = results / "transformed_exact_independent_verification.csv"
    report.to_csv(output, index=False)
    maximum = float(report["absolute_discrepancy"].max()) if len(report) else np.nan
    print(f"Verified {len(selected)} selected rows ({len(report)} transformed tests).")
    print(f"Maximum absolute discrepancy: {maximum:.3g}")
    print(f"Wrote {output}")
    if not report["within_tolerance"].all():
        raise SystemExit("Independent exact-p verification failed")


if __name__ == "__main__":
    main()
