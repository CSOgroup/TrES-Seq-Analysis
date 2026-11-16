"""Clustering agreement metrics used in the resolution-sweep panels."""

from __future__ import annotations

import numpy as np
from sklearn.metrics import (
    adjusted_mutual_info_score,
    adjusted_rand_score,
    completeness_score,
    homogeneity_score,
    mutual_info_score,
    normalized_mutual_info_score,
    v_measure_score,
)


def _entropy(labels) -> float:
    _, counts = np.unique(np.asarray(labels), return_counts=True)
    p = counts / counts.sum()
    return float(-(p * np.log2(p + 1e-15)).sum())


def conditional_entropy(y_true, y_pred) -> float:
    """H(y_true | y_pred) in bits. sklearn mutual_info_score is in nats."""
    h_t = _entropy(y_true)
    mi_bits = mutual_info_score(y_true, y_pred) / np.log(2)
    return float(h_t - mi_bits)


def compare_clusterings(labels_a, labels_b, name_a: str = "A", name_b: str = "B") -> dict:
    """External metrics between two partitions of the same cells."""
    a = np.asarray(labels_a).astype(str)
    b = np.asarray(labels_b).astype(str)
    return {
        f"n_{name_a}": int(len(np.unique(a))),
        f"n_{name_b}": int(len(np.unique(b))),
        "ARI": float(adjusted_rand_score(a, b)),
        "AMI": float(adjusted_mutual_info_score(a, b)),
        "NMI": float(normalized_mutual_info_score(a, b)),
        f"homogeneity_{name_a}_given_{name_b}": float(homogeneity_score(b, a)),
        f"completeness_{name_a}_given_{name_b}": float(completeness_score(b, a)),
        f"v_measure_{name_a}_given_{name_b}": float(v_measure_score(b, a)),
        f"homogeneity_{name_b}_given_{name_a}": float(homogeneity_score(a, b)),
        f"completeness_{name_b}_given_{name_a}": float(completeness_score(a, b)),
        f"v_measure_{name_b}_given_{name_a}": float(v_measure_score(a, b)),
        f"H_{name_a}_given_{name_b}": conditional_entropy(a, b),
        f"H_{name_b}_given_{name_a}": conditional_entropy(b, a),
        f"H_{name_a}": _entropy(a),
        f"H_{name_b}": _entropy(b),
    }


def vs_celltype(labels, cell_types) -> dict:
    """Metrics for Leiden clusters versus cell-type labels."""
    metrics = compare_clusterings(labels, cell_types, "clusters", "celltype")
    return {
        "n_clusters": int(metrics["n_clusters"]),
        "n_types": int(metrics["n_celltype"]),
        "ARI": metrics["ARI"],
        "NMI": metrics["NMI"],
        "AMI": metrics["AMI"],
        "H_celltype_given_clusters": metrics["H_celltype_given_clusters"],
        "H_clusters_given_celltype": metrics["H_clusters_given_celltype"],
        "homogeneity_clusters_given_celltype": metrics[
            "homogeneity_clusters_given_celltype"
        ],
        "completeness_clusters_given_celltype": metrics[
            "completeness_clusters_given_celltype"
        ],
    }
