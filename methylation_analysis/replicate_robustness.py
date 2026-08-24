"""Hostile robustness audit of archived methylation–isoform associations.

The WGBS and RNA replicates are independent samples, not paired observations.
Accordingly, the developmental stage is the inferential unit. Replicates are
used to estimate stage means, profile reliability and resampling uncertainty;
they are never treated as 16 paired methylation/RNA observations.
"""

from __future__ import annotations

import argparse
import itertools
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / "results"
DEFAULT_INPUT_DIR = Path(os.environ.get("METHYLATION_INPUT_DIR", ROOT)).expanduser()
STAGES = ["10.5", "11.5", "12.5", "13.5", "14.5", "15.5", "16.5", "0"]
# Ordinal stage positions, not elapsed embryonic days. In particular, P0 is
# treated as the next ordered observation after E16.5 for the post hoc
# detrending diagnostic; adjacent differences are likewise not time-scaled.
STAGE_ORDER = np.arange(8, dtype=float)
RNG = np.random.default_rng(20260730)


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    adjusted = np.full(p.shape, np.nan)
    valid = np.isfinite(p)
    pv = p[valid]
    if not len(pv):
        return adjusted
    order = np.argsort(pv)
    ranked = pv[order]
    correction = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    correction = np.minimum.accumulate(correction[::-1])[::-1]
    restored = np.empty_like(correction)
    restored[order] = np.minimum(correction, 1.0)
    adjusted[valid] = restored
    return adjusted


def safe_pearson(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 4 or np.nanstd(x[valid]) == 0 or np.nanstd(y[valid]) == 0:
        return np.nan, np.nan
    x_valid = x[valid]
    y_valid = y[valid]
    x_centred = x_valid - x_valid.mean()
    y_centred = y_valid - y_valid.mean()
    denominator = np.sqrt(np.sum(x_centred**2) * np.sum(y_centred**2))
    correlation = float(np.sum(x_centred * y_centred) / denominator)
    correlation = float(np.clip(correlation, -1, 1))
    if abs(correlation) == 1:
        return correlation, 0.0
    degrees_freedom = len(x_valid) - 2
    statistic = correlation * np.sqrt(
        degrees_freedom / max(1 - correlation**2, np.finfo(float).tiny)
    )
    p_value = float(2 * stats.t.sf(abs(statistic), degrees_freedom))
    return correlation, p_value


def safe_correlation(x: np.ndarray, y: np.ndarray) -> float:
    """Fast Pearson correlation for diagnostics that do not require a p-value."""
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 4:
        return np.nan
    x_valid = x[valid]
    y_valid = y[valid]
    x_centred = x_valid - x_valid.mean()
    y_centred = y_valid - y_valid.mean()
    denominator = np.sqrt(np.sum(x_centred**2) * np.sum(y_centred**2))
    if denominator == 0:
        return np.nan
    return float(np.clip(np.sum(x_centred * y_centred) / denominator, -1, 1))


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 4:
        return np.nan
    return safe_correlation(
        stats.rankdata(x[valid], method="average"),
        stats.rankdata(y[valid], method="average"),
    )


def linear_residual(values: np.ndarray) -> np.ndarray:
    valid = np.isfinite(values)
    output = np.full_like(values, np.nan, dtype=float)
    if valid.sum() < 4:
        return output
    design = np.column_stack([np.ones(valid.sum()), STAGE_ORDER[valid]])
    coefficients = np.linalg.lstsq(design, values[valid], rcond=None)[0]
    output[valid] = values[valid] - design @ coefficients
    return output


def leave_one_out(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    values = []
    for index in range(len(x)):
        keep = np.ones(len(x), dtype=bool)
        keep[index] = False
        values.append(safe_correlation(x[keep], y[keep]))
    return np.asarray(values)


def profile_from(frame: pd.DataFrame, value: str, keys: list[str]) -> dict[tuple, np.ndarray]:
    pivot = frame.pivot_table(index=keys, columns="stage", values=value, aggfunc="mean")
    pivot = pivot.reindex(columns=STAGES)
    return {index if isinstance(index, tuple) else (index,): row.to_numpy(float)
            for index, row in pivot.iterrows()}


def replicate_profile_from(
    frame: pd.DataFrame, value: str, keys: list[str]
) -> dict[tuple, np.ndarray]:
    pivot = frame.pivot_table(
        index=keys, columns=["stage", "replicate"], values=value, aggfunc="mean"
    )
    desired = pd.MultiIndex.from_product([STAGES, [1, 2]], names=["stage", "replicate"])
    pivot = pivot.reindex(columns=desired)
    return {
        index if isinstance(index, tuple) else (index,): row.to_numpy(float).reshape(8, 2)
        for index, row in pivot.iterrows()
    }


def exhaustive_stage_permutation_p(
    x_profiles: np.ndarray, y_profiles: np.ndarray, batch_size: int = 32
) -> np.ndarray:
    """Two-sided exhaustive stage-label-randomisation p-values.

    Developmental stage, rather than assay replicate, is the cross-omic unit.
    For complete eight-stage profiles this enumerates all 8! assignments in
    memory-bounded batches. Profiles with fewer complete stages are enumerated
    separately over their available stages. The observed assignment is part of
    the randomisation space, so no Monte Carlo correction is required. The
    enumeration is exact only under exchangeable stage labels; it is not a
    general time-series independence test.
    """
    if x_profiles.shape != y_profiles.shape or x_profiles.ndim != 2:
        raise ValueError("x_profiles and y_profiles must be equal 2D matrices")

    output = np.full(x_profiles.shape[0], np.nan, dtype=float)
    complete_count = np.sum(np.isfinite(x_profiles) & np.isfinite(y_profiles), axis=1)

    x_scale = np.full(x_profiles.shape[0], np.nan, dtype=float)
    y_scale = np.full(y_profiles.shape[0], np.nan, dtype=float)
    fully_observed = complete_count == x_profiles.shape[1]
    x_scale[fully_observed] = np.std(x_profiles[fully_observed], axis=1)
    y_scale[fully_observed] = np.std(y_profiles[fully_observed], axis=1)
    complete_rows = np.flatnonzero(
        (complete_count == x_profiles.shape[1])
        & np.isfinite(x_scale)
        & np.isfinite(y_scale)
        & (x_scale > 0)
        & (y_scale > 0)
    )
    permutations = np.asarray(
        list(itertools.permutations(range(x_profiles.shape[1]))), dtype=np.int8
    )
    for start in range(0, len(complete_rows), batch_size):
        rows = complete_rows[start : start + batch_size]
        x = x_profiles[rows].astype(float, copy=True)
        y = y_profiles[rows].astype(float, copy=True)
        x -= x.mean(axis=1, keepdims=True)
        y -= y.mean(axis=1, keepdims=True)
        x /= np.sqrt(np.sum(x**2, axis=1, keepdims=True))
        y /= np.sqrt(np.sum(y**2, axis=1, keepdims=True))
        observed = np.abs(np.sum(x * y, axis=1))
        permuted_x = x[:, permutations]
        correlations = np.einsum("bpj,bj->bp", permuted_x, y, optimize=True)
        output[rows] = np.mean(
            np.abs(correlations) >= observed[:, None] - 1e-12, axis=1
        )

    for row in np.flatnonzero((complete_count >= 4) & (complete_count < 8)):
        valid = np.isfinite(x_profiles[row]) & np.isfinite(y_profiles[row])
        x = x_profiles[row, valid].astype(float, copy=True)
        y = y_profiles[row, valid].astype(float, copy=True)
        x -= x.mean()
        y -= y.mean()
        x_norm = np.sqrt(np.sum(x**2))
        y_norm = np.sqrt(np.sum(y**2))
        if x_norm == 0 or y_norm == 0:
            continue
        x /= x_norm
        y /= y_norm
        observed = abs(float(x @ y))
        local_permutations = np.asarray(
            list(itertools.permutations(range(len(x)))), dtype=np.int8
        )
        correlations = x[local_permutations] @ y
        output[row] = np.mean(np.abs(correlations) >= observed - 1e-12)

    return output


def exhaustive_transformed_stage_permutation_p(
    x_profiles: np.ndarray,
    y_profiles: np.ndarray,
    transform: str,
    batch_size: int = 16,
) -> np.ndarray:
    """Exhaustive stage-label-randomisation p-values after a profile transform.

    The transform is recomputed after every stage permutation, rather than
    permuting already transformed residuals or differences. This preserves the
    ordinal ordering used by the detrending/first-difference specification; it
    does not model unequal elapsed time between E16.5 and P0.
    """
    if transform not in {"detrended", "first_difference"}:
        raise ValueError(f"Unsupported transform: {transform}")
    if x_profiles.shape != y_profiles.shape or x_profiles.ndim != 2:
        raise ValueError("x_profiles and y_profiles must be equal 2D matrices")

    def apply(values: np.ndarray) -> np.ndarray:
        if transform == "first_difference":
            return np.diff(values, axis=-1)
        time = np.arange(values.shape[-1], dtype=float)
        centred_time = time - time.mean()
        centred = values - values.mean(axis=-1, keepdims=True)
        slope = np.sum(centred * centred_time, axis=-1, keepdims=True) / np.sum(
            centred_time**2
        )
        return centred - slope * centred_time

    def normalize(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        centred = values - values.mean(axis=-1, keepdims=True)
        norm = np.sqrt(np.sum(centred**2, axis=-1, keepdims=True))
        normalized = np.divide(
            centred,
            norm,
            out=np.full_like(centred, np.nan, dtype=float),
            where=norm > 0,
        )
        return normalized, norm[..., 0]

    output = np.full(x_profiles.shape[0], np.nan, dtype=float)
    complete_count = np.sum(np.isfinite(x_profiles) & np.isfinite(y_profiles), axis=1)
    complete_rows = np.flatnonzero(complete_count == x_profiles.shape[1])
    permutations = np.asarray(
        list(itertools.permutations(range(x_profiles.shape[1]))), dtype=np.int8
    )

    for start in range(0, len(complete_rows), batch_size):
        rows = complete_rows[start : start + batch_size]
        x_raw = x_profiles[rows].astype(float, copy=False)
        y_raw = y_profiles[rows].astype(float, copy=False)
        observed_x, observed_x_norm = normalize(apply(x_raw))
        observed_y, observed_y_norm = normalize(apply(y_raw))
        observed = np.abs(np.sum(observed_x * observed_y, axis=1))

        permuted_raw = x_raw[:, permutations]
        permuted_x, permuted_norm = normalize(apply(permuted_raw))
        correlations = np.einsum("bpj,bj->bp", permuted_x, observed_y, optimize=True)
        valid = (
            (observed_x_norm > 0)
            & (observed_y_norm > 0)
            & np.isfinite(observed)
        )
        counts = np.mean(
            np.abs(correlations) >= observed[:, None] - 1e-12, axis=1
        )
        output[rows[valid]] = counts[valid]

    for row in np.flatnonzero((complete_count >= 5) & (complete_count < 8)):
        valid = np.isfinite(x_profiles[row]) & np.isfinite(y_profiles[row])
        x_raw = x_profiles[row, valid].astype(float, copy=False)
        y_raw = y_profiles[row, valid].astype(float, copy=False)
        local_permutations = np.asarray(
            list(itertools.permutations(range(len(x_raw)))), dtype=np.int8
        )
        observed_x, observed_x_norm = normalize(apply(x_raw[None, :]))
        observed_y, observed_y_norm = normalize(apply(y_raw[None, :]))
        if observed_x_norm[0] <= 0 or observed_y_norm[0] <= 0:
            continue
        permuted_x, _ = normalize(apply(x_raw[local_permutations]))
        observed = abs(float(observed_x[0] @ observed_y[0]))
        correlations = permuted_x @ observed_y[0]
        output[row] = np.mean(np.abs(correlations) >= observed - 1e-12)

    return output


def bootstrap_correlation(
    methylation: np.ndarray, usage: np.ndarray, draws: int = 5000
) -> tuple[float, float, float, float]:
    if methylation.shape != (8, 2) or usage.shape != (8, 2):
        return np.nan, np.nan, np.nan, np.nan
    choice_m = RNG.integers(0, 2, size=(draws, 8))
    choice_u = RNG.integers(0, 2, size=(draws, 8))
    rows = np.arange(8)[None, :]
    x = methylation[rows, choice_m]
    y = usage[rows, choice_u]
    valid = np.isfinite(x) & np.isfinite(y)
    enough = valid.sum(axis=1) >= 6
    x = np.where(valid, x, np.nan)
    y = np.where(valid, y, np.nan)
    x_centered = x - np.nanmean(x, axis=1, keepdims=True)
    y_centered = y - np.nanmean(y, axis=1, keepdims=True)
    numerator = np.nansum(x_centered * y_centered, axis=1)
    denominator = np.sqrt(
        np.nansum(x_centered**2, axis=1) * np.nansum(y_centered**2, axis=1)
    )
    correlations = np.divide(
        numerator,
        denominator,
        out=np.full(draws, np.nan),
        where=(denominator > 0) & enough,
    )
    correlations = correlations[np.isfinite(correlations)]
    if not len(correlations):
        return np.nan, np.nan, np.nan, np.nan
    median, low, high = np.quantile(correlations, [0.5, 0.025, 0.975])
    positive = np.mean(correlations > 0)
    return float(median), float(low), float(high), float(positive)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Reproduce the complete methylation--isoform robustness audit."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        help=(
            "Directory containing methylation_expression_correlation_data.csv. "
            "Defaults to METHYLATION_INPUT_DIR or the historical workspace root."
        ),
    )
    args = parser.parse_args()
    input_dir = (args.input_dir or DEFAULT_INPUT_DIR).resolve()
    archived_path = input_dir / "methylation_expression_correlation_data.csv"
    if not archived_path.is_file():
        raise FileNotFoundError(f"Missing archived correlation table: {archived_path}")

    methylation = pd.read_csv(
        RESULTS / "replicate_level_methylation.csv", dtype={"stage": str}
    )
    usage = pd.read_csv(
        RESULTS / "replicate_level_isoform_usage.csv", dtype={"stage": str}
    )
    archived = pd.read_csv(archived_path)
    archived["stage_dummy"] = 1
    isoform_gene = archived[["isoform_id", "gene_id"]].drop_duplicates()
    if isoform_gene["isoform_id"].duplicated().any():
        raise ValueError("An isoform maps to more than one gene in the archived table")
    usage_with_gene = usage.merge(isoform_gene, on="isoform_id", how="left")

    methyl_keys = ["tissue", "gene_id", "region"]
    usage_keys = ["tissue", "isoform_id"]
    methyl_mean = profile_from(methylation, "mean_weighted", methyl_keys)
    methyl_unweighted = profile_from(methylation, "mean_unweighted", methyl_keys)
    methyl_strand = profile_from(
        methylation, "mean_weighted_strand_sensitive", methyl_keys
    )
    usage_mean = profile_from(usage, "rep_if", usage_keys)
    expression_mean = profile_from(usage, "expression", usage_keys)
    gene_expression_mean = profile_from(
        usage_with_gene.dropna(subset=["gene_id"])
        .groupby(["tissue", "gene_id", "stage", "replicate"], as_index=False)[
            "expression"
        ]
        .sum(),
        "expression",
        ["tissue", "gene_id"],
    )
    methyl_replicates = replicate_profile_from(
        methylation, "mean_weighted", methyl_keys
    )
    usage_replicates = replicate_profile_from(usage, "rep_if", usage_keys)

    quality = (
        methylation.groupby(methyl_keys)
        .agg(
            median_cpg=("n_cpg", "median"),
            min_cpg=("n_cpg", "min"),
            median_coverage=("total_coverage", "median"),
            strand_cpg_gain=("n_cpg", "median"),
            strand_sensitive_cpg=("n_cpg_strand_sensitive", "median"),
        )
        .reset_index()
    )
    quality["strand_cpg_gain"] = (
        quality["strand_cpg_gain"] - quality["strand_sensitive_cpg"]
    )

    records = []
    exact_x_profiles = []
    exact_y_profiles = []
    for row in archived.itertuples(index=False):
        m_key = (row.tissue, row.gene_id, row.region)
        u_key = (row.tissue, row.isoform_id)
        x = methyl_mean.get(m_key, np.full(8, np.nan))
        x_unweighted = methyl_unweighted.get(m_key, np.full(8, np.nan))
        x_strand = methyl_strand.get(m_key, np.full(8, np.nan))
        y = usage_mean.get(u_key, np.full(8, np.nan))
        y_expression = expression_mean.get(u_key, np.full(8, np.nan))
        y_gene_expression = gene_expression_mean.get(
            (row.tissue, row.gene_id), np.full(8, np.nan)
        )
        exact_x_profiles.append(x)
        exact_y_profiles.append(y)

        raw_r, raw_p = safe_pearson(x, y)
        unweighted_r, unweighted_p = safe_pearson(x_unweighted, y)
        strand_r, strand_p = safe_pearson(x_strand, y)
        detrended_r, detrended_p = safe_pearson(linear_residual(x), linear_residual(y))
        difference_r, difference_p = safe_pearson(np.diff(x), np.diff(y))
        isoform_expression_r, isoform_expression_p = safe_pearson(x, y_expression)
        gene_expression_r, gene_expression_p = safe_pearson(x, y_gene_expression)
        loo = leave_one_out(x, y)
        original_sign = np.sign(raw_r) if np.isfinite(raw_r) and raw_r != 0 else np.nan
        loo_same = (
            np.mean(np.sign(loo[np.isfinite(loo)]) == original_sign)
            if np.isfinite(original_sign)
            else np.nan
        )

        m_reps = methyl_replicates.get(m_key, np.full((8, 2), np.nan))
        u_reps = usage_replicates.get(u_key, np.full((8, 2), np.nan))
        meth_reliability = safe_correlation(m_reps[:, 0], m_reps[:, 1])
        usage_reliability = safe_correlation(u_reps[:, 0], u_reps[:, 1])

        records.append(
            {
                "gene_id": row.gene_id,
                "isoform_id": row.isoform_id,
                "region": row.region,
                "tissue": row.tissue,
                "archived_r": row.correlation,
                "archived_p": row.p_value,
                "archived_q": row.adj_p_value,
                "raw_weighted_r": raw_r,
                "raw_weighted_p": raw_p,
                "raw_unweighted_r": unweighted_r,
                "raw_unweighted_p": unweighted_p,
                "strand_sensitive_r": strand_r,
                "strand_sensitive_p": strand_p,
                "spearman_rho": safe_spearman(x, y),
                "linear_detrended_r": detrended_r,
                "linear_detrended_p": detrended_p,
                "first_difference_r": difference_r,
                "first_difference_p": difference_p,
                "isoform_absolute_expression_r": isoform_expression_r,
                "isoform_absolute_expression_p": isoform_expression_p,
                "gene_total_expression_r": gene_expression_r,
                "gene_total_expression_p": gene_expression_p,
                "loo_min_r": np.nanmin(loo) if np.isfinite(loo).any() else np.nan,
                "loo_max_r": np.nanmax(loo) if np.isfinite(loo).any() else np.nan,
                "loo_same_sign_fraction": loo_same,
                "methylation_replicate_reliability": meth_reliability,
                "usage_replicate_reliability": usage_reliability,
                "n_complete_stages": int(np.sum(np.isfinite(x) & np.isfinite(y))),
            }
        )

    audit = pd.DataFrame(records).merge(quality, on=methyl_keys, how="left")
    audit["raw_weighted_exact_permutation_p"] = exhaustive_stage_permutation_p(
        np.asarray(exact_x_profiles, dtype=float),
        np.asarray(exact_y_profiles, dtype=float),
    )
    audit["detrended_exact_permutation_p"] = exhaustive_transformed_stage_permutation_p(
        np.asarray(exact_x_profiles, dtype=float),
        np.asarray(exact_y_profiles, dtype=float),
        "detrended",
    )
    audit["difference_exact_permutation_p"] = exhaustive_transformed_stage_permutation_p(
        np.asarray(exact_x_profiles, dtype=float),
        np.asarray(exact_y_profiles, dtype=float),
        "first_difference",
    )
    audit["raw_global_q"] = bh_adjust(audit["raw_weighted_p"].to_numpy())
    audit["unweighted_global_q"] = bh_adjust(audit["raw_unweighted_p"].to_numpy())
    audit["strand_sensitive_global_q"] = bh_adjust(
        audit["strand_sensitive_p"].to_numpy()
    )
    audit["detrended_global_q"] = bh_adjust(audit["linear_detrended_p"].to_numpy())
    audit["difference_global_q"] = bh_adjust(audit["first_difference_p"].to_numpy())
    audit["raw_weighted_exact_global_q"] = bh_adjust(
        audit["raw_weighted_exact_permutation_p"].to_numpy()
    )
    audit["detrended_exact_global_q"] = bh_adjust(
        audit["detrended_exact_permutation_p"].to_numpy()
    )
    audit["difference_exact_global_q"] = bh_adjust(
        audit["difference_exact_permutation_p"].to_numpy()
    )
    audit["archived_sign_reproduced"] = (
        np.sign(audit["archived_r"]) == np.sign(audit["raw_weighted_r"])
    )
    audit["raw_nominal_reproduced"] = (
        audit["archived_sign_reproduced"] & (audit["raw_weighted_p"] < 0.05)
    )
    audit.to_csv(RESULTS / "replicate_robustness_all.csv", index=False)

    # Audit both the originally selected rows and any discoveries that emerge
    # after reconstruction.  The union prevents the counter-audit from being
    # conditioned solely on the original selection.
    hits = audit.loc[
        (audit["archived_q"] < 0.05) | (audit["raw_global_q"] < 0.05)
    ].copy()
    hits["selection_origin"] = np.select(
        [
            (hits["archived_q"] < 0.05) & (hits["raw_global_q"] < 0.05),
            hits["archived_q"] < 0.05,
            hits["raw_global_q"] < 0.05,
        ],
        ["both", "archived_only", "reconstructed_only"],
        default="neither",
    )
    permutations = np.asarray(list(itertools.permutations(range(8))), dtype=np.int8)
    bootstrap_rows = []
    for row in hits.itertuples(index=False):
        m_key = (row.tissue, row.gene_id, row.region)
        u_key = (row.tissue, row.isoform_id)
        x = methyl_mean[m_key]
        y = usage_mean[u_key]
        observed = abs(row.raw_weighted_r)
        complete = np.isfinite(x) & np.isfinite(y)
        if complete.all():
            # Vectorise the full 8! enumeration.  This is numerically equivalent
            # to calling Pearson's r once per permutation but avoids millions of
            # Python-level function calls during the hostile audit.
            permuted_x = x[permutations]
            permuted_x = permuted_x - permuted_x.mean(axis=1, keepdims=True)
            centred_y = y - y.mean()
            denominator = np.sqrt(
                np.sum(permuted_x**2, axis=1) * np.sum(centred_y**2)
            )
            permuted = np.divide(
                permuted_x @ centred_y,
                denominator,
                out=np.full(len(permutations), np.nan),
                where=denominator > 0,
            )
        else:
            # Archived hits are expected to have all eight stages, but retain a
            # exhaustive randomisation if a reconstructed profile contains a gap.
            x_complete = x[complete]
            y_complete = y[complete]
            complete_permutations = np.asarray(
                list(itertools.permutations(range(len(x_complete)))), dtype=np.int8
            )
            permuted_x = x_complete[complete_permutations]
            permuted_x = permuted_x - permuted_x.mean(axis=1, keepdims=True)
            centred_y = y_complete - y_complete.mean()
            denominator = np.sqrt(
                np.sum(permuted_x**2, axis=1) * np.sum(centred_y**2)
            )
            permuted = np.divide(
                permuted_x @ centred_y,
                denominator,
                out=np.full(len(complete_permutations), np.nan),
                where=denominator > 0,
            )
        # This is an exhaustive randomisation p-value conditional on treating
        # stage labels as exchangeable.  It is not exact time-series inference,
        # because ordered developmental stages are not generally exchangeable.
        exact_p = np.nanmean(np.abs(permuted) >= observed - 1e-12)
        circular = np.asarray(
            [safe_correlation(np.roll(x, shift), y) for shift in range(8)]
        )
        circular_rank = 1 + int(np.sum(np.abs(circular[1:]) >= observed - 1e-12))
        median, low, high, positive = bootstrap_correlation(
            methyl_replicates[m_key], usage_replicates[u_key]
        )
        archived_positive = row.archived_r > 0
        sign_probability = positive if archived_positive else 1 - positive
        bootstrap_rows.append(
            {
                "gene_id": row.gene_id,
                "isoform_id": row.isoform_id,
                "region": row.region,
                "tissue": row.tissue,
                "stage_permutation_p": exact_p,
                "circular_abs_rank": circular_rank,
                "bootstrap_median_r": median,
                "bootstrap_low_r": low,
                "bootstrap_high_r": high,
                "bootstrap_archived_sign_probability": sign_probability,
            }
        )
    hits = hits.merge(
        pd.DataFrame(bootstrap_rows),
        on=["gene_id", "isoform_id", "region", "tissue"],
        how="left",
    )
    # Do not BH-adjust only this selected subset and call it confirmatory. With
    # eight stages, the finite randomisation resolution is itself part of the
    # design audit; stage_permutation_p remains a per-profile sensitivity under
    # an exchangeability null, not exact time-series inference.
    hits["detrended_same_sign"] = (
        np.sign(hits["archived_r"]) == np.sign(hits["linear_detrended_r"])
    )
    hits["difference_same_sign"] = (
        np.sign(hits["archived_r"]) == np.sign(hits["first_difference_r"])
    )
    hits["bootstrap_interval_excludes_zero"] = (
        (hits["bootstrap_low_r"] > 0) | (hits["bootstrap_high_r"] < 0)
    )
    hits["hostile_support_score"] = (
        hits["archived_sign_reproduced"].astype(int)
        + hits["raw_nominal_reproduced"].astype(int)
        + (hits["loo_same_sign_fraction"] == 1).astype(int)
        + hits["detrended_same_sign"].astype(int)
        + hits["difference_same_sign"].astype(int)
        + (hits["circular_abs_rank"] == 1).astype(int)
        + hits["bootstrap_interval_excludes_zero"].astype(int)
        + (hits["bootstrap_archived_sign_probability"] >= 0.95).astype(int)
        + (hits["methylation_replicate_reliability"] > 0.5).astype(int)
        + (hits["usage_replicate_reliability"] > 0.5).astype(int)
    )
    hits = hits.sort_values(
        ["hostile_support_score", "stage_permutation_p", "archived_q"],
        ascending=[False, True, True],
    )
    hits.to_csv(RESULTS / "replicate_robustness_archived_hits.csv", index=False)
    (
        hits.sort_values(
            ["raw_global_q", "hostile_support_score", "stage_permutation_p"],
            ascending=[True, False, True],
            na_position="last",
        )
        .drop_duplicates("gene_id")
        .head(20)
        .to_csv(RESULTS / "hostile_candidate_shortlist.csv", index=False)
    )

    archived_hits = hits.loc[hits["selection_origin"].isin(["archived_only", "both"])]
    reconstructed_only_hits = hits.loc[hits["selection_origin"] == "reconstructed_only"]

    summary = pd.DataFrame(
        [
            ("all_tests", len(audit)),
            ("archived_fdr_rows", int((audit["archived_q"] < 0.05).sum())),
            ("reconstructed_parametric_global_fdr_rows", int((audit["raw_global_q"] < 0.05).sum())),
            ("overlap_between_fdr_sets", int(((audit["archived_q"] < 0.05) & (audit["raw_global_q"] < 0.05)).sum())),
            ("reconstructed_parametric_global_fdr_genes", int(audit.loc[audit["raw_global_q"] < 0.05, "gene_id"].nunique())),
            ("reconstructed_stage_randomisation_global_fdr_rows", int((audit["raw_weighted_exact_global_q"] < 0.05).sum())),
            ("reconstructed_stage_randomisation_nominal_rows", int((audit["raw_weighted_exact_permutation_p"] < 0.05).sum())),
            ("detrended_parametric_global_fdr_rows", int((audit["detrended_global_q"] < 0.05).sum())),
            ("detrended_stage_randomisation_nominal_rows", int((audit["detrended_exact_permutation_p"] < 0.05).sum())),
            ("detrended_stage_randomisation_global_fdr_rows", int((audit["detrended_exact_global_q"] < 0.05).sum())),
            ("difference_parametric_global_fdr_rows", int((audit["difference_global_q"] < 0.05).sum())),
            ("difference_stage_randomisation_nominal_rows", int((audit["difference_exact_permutation_p"] < 0.05).sum())),
            ("difference_stage_randomisation_global_fdr_rows", int((audit["difference_exact_global_q"] < 0.05).sum())),
            ("unweighted_global_fdr_rows", int((audit["unweighted_global_q"] < 0.05).sum())),
            ("strand_sensitive_global_fdr_rows", int((audit["strand_sensitive_global_q"] < 0.05).sum())),
            ("targeted_union_rows", len(hits)),
            ("targeted_archived_hit_rows", len(archived_hits)),
            ("targeted_reconstructed_only_rows", len(reconstructed_only_hits)),
            ("archived_hits_raw_sign_reproduced", int(archived_hits["archived_sign_reproduced"].sum())),
            ("archived_hits_raw_nominal_same_sign", int(archived_hits["raw_nominal_reproduced"].sum())),
            ("archived_hits_leave_one_stage_out_same_sign", int((archived_hits["loo_same_sign_fraction"] == 1).sum())),
            ("archived_hits_detrended_same_sign", int(archived_hits["detrended_same_sign"].sum())),
            ("archived_hits_first_difference_same_sign", int(archived_hits["difference_same_sign"].sum())),
            ("archived_hits_circular_shift_strongest", int((archived_hits["circular_abs_rank"] == 1).sum())),
            ("archived_hits_bootstrap_interval_excludes_zero", int(archived_hits["bootstrap_interval_excludes_zero"].sum())),
            ("archived_hits_stage_randomisation_p_lt_0.05", int((archived_hits["stage_permutation_p"] < 0.05).sum())),
            ("targeted_union_leave_one_stage_out_same_sign", int((hits["loo_same_sign_fraction"] == 1).sum())),
            ("targeted_union_detrended_same_sign", int(hits["detrended_same_sign"].sum())),
            ("targeted_union_first_difference_same_sign", int(hits["difference_same_sign"].sum())),
            ("targeted_union_circular_shift_strongest", int((hits["circular_abs_rank"] == 1).sum())),
            ("targeted_union_bootstrap_interval_excludes_zero", int(hits["bootstrap_interval_excludes_zero"].sum())),
            ("targeted_union_stage_randomisation_p_lt_0.05", int((hits["stage_permutation_p"] < 0.05).sum())),
        ],
        columns=["quantity", "value"],
    )
    summary.to_csv(RESULTS / "replicate_robustness_summary.csv", index=False)
    print(summary.to_string(index=False))
    print("\nLeading hostile-audit candidates:")
    print(
        hits[
            [
                "gene_id",
                "isoform_id",
                "region",
                "archived_r",
                "raw_weighted_r",
                "linear_detrended_r",
                "first_difference_r",
                "bootstrap_low_r",
                "bootstrap_high_r",
                "hostile_support_score",
            ]
        ]
        .head(20)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
