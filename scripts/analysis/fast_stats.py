"""Vectorised cluster bootstrap (same estimand as phageannoscope.analysis.cluster_bootstrap).

Each replicate resamples relatedness blocks with replacement (multinomial counts),
then recomputes per-level proportions from per-block sums.  Kept in scripts/ so
that the pipeline code hash (run identity) is not changed mid-study.
"""

import numpy as np


def cluster_bootstrap(df, value, by, cluster="block", reps=2000, seed=30102026):
    rng = np.random.default_rng(seed)
    levels = sorted(df[by].unique())
    codes, uniq = np.unique(df[cluster].to_numpy(), return_inverse=True)[1], None
    n_clusters = codes.max() + 1
    sums = np.zeros((n_clusters, len(levels)))
    counts = np.zeros((n_clusters, len(levels)))
    vals = df[value].to_numpy(dtype=float)
    lv = df[by].to_numpy()
    for j, level in enumerate(levels):
        mask = lv == level
        sums[:, j] = np.bincount(codes[mask], weights=vals[mask], minlength=n_clusters)
        counts[:, j] = np.bincount(codes[mask], minlength=n_clusters)
    draws = np.empty((reps, len(levels)))
    for r in range(reps):
        w = rng.multinomial(n_clusters, np.full(n_clusters, 1 / n_clusters))
        with np.errstate(invalid="ignore", divide="ignore"):
            draws[r] = (w @ sums) / (w @ counts)
    point = sums.sum(0) / counts.sum(0)
    out = {"levels": {}, "n_clusters": int(n_clusters), "reps": reps,
           "method": "vectorised multinomial cluster bootstrap"}
    for j, level in enumerate(levels):
        out["levels"][str(level)] = {"n": int(counts[:, j].sum()), "proportion": float(point[j]),
                                     "ci95": [float(np.nanpercentile(draws[:, j], 2.5)),
                                              float(np.nanpercentile(draws[:, j], 97.5))]}
    if len(levels) == 2:
        d = draws[:, 1] - draws[:, 0]
        out["difference_pp"] = {"estimate": float(100 * (point[1] - point[0])),
                                "ci95": [float(100 * np.nanpercentile(d, 2.5)),
                                         float(100 * np.nanpercentile(d, 97.5))],
                                "contrast": f"{levels[1]} minus {levels[0]}"}
    return out
