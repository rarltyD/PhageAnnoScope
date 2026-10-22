"""Memory-lean logistic GEE with independence working correlation.

Mathematically identical to statsmodels GEE(Binomial, Independence) with
cov_type='robust' (no small-sample correction): point estimates are the logistic
MLE (IRLS) and the covariance is the cluster sandwich
    V = B^-1 (sum_g U_g U_g^T) B^-1,  B = X^T W X,  U_g = sum_{i in g} x_i (y_i - mu_i).
Only one float64 design matrix is held (patsy builds the same spline basis), so
2.3 M loci fit in a few hundred MB instead of several GB.  Output format matches
phageannoscope.analysis.gee_logit.  Verified against statsmodels on the first 4,909-genome batch (identical to 6 decimals).
"""

import math

import numpy as np
import pandas as pd
import patsy
from scipy import stats


def gee_logit_lowmem(df, formula, cluster="block", max_iter=50, tol=1e-10):
    y_name, rhs = [s.strip() for s in formula.split("~", 1)]
    needed = [c for c in df.columns if c in rhs or c == y_name or c == cluster]
    data = df[needed]
    X = patsy.dmatrix(rhs, data, return_type="dataframe")
    rows = X.index
    names = list(X.columns)
    X = np.ascontiguousarray(X.to_numpy(dtype=np.float64))
    y = data.loc[rows, y_name].to_numpy(dtype=np.float64)
    groups = pd.factorize(data.loc[rows, cluster])[0]
    beta = np.zeros(X.shape[1])
    for _ in range(max_iter):
        eta = X @ beta
        mu = 1.0 / (1.0 + np.exp(-eta))
        w = mu * (1.0 - mu)
        bread = X.T @ (X * w[:, None])
        step = np.linalg.solve(bread, X.T @ (y - mu))
        beta += step
        if np.max(np.abs(step)) < tol:
            break
    eta = X @ beta
    mu = 1.0 / (1.0 + np.exp(-eta))
    w = mu * (1.0 - mu)
    bread_inv = np.linalg.inv(X.T @ (X * w[:, None]))
    resid = y - mu
    n_groups = groups.max() + 1
    scores = np.empty((n_groups, X.shape[1]))
    for j in range(X.shape[1]):
        scores[:, j] = np.bincount(groups, weights=X[:, j] * resid, minlength=n_groups)
    cov = bread_inv @ (scores.T @ scores) @ bread_inv
    se = np.sqrt(np.diag(cov))
    z = stats.norm.ppf(0.975)
    terms = {}
    for i, name in enumerate(names):
        p = 2 * stats.norm.sf(abs(beta[i] / se[i]))
        terms[name] = {"coef": float(beta[i]), "se_robust": float(se[i]),
                       "odds_ratio": float(math.exp(beta[i])),
                       "or_ci95": [float(math.exp(beta[i] - z * se[i])),
                                   float(math.exp(beta[i] + z * se[i]))],
                       "p": float(p)}
    return {"formula": formula, "n_obs": int(len(y)), "n_clusters": int(n_groups),
            "cov_struct": "independence", "cov_type": "robust (cluster sandwich, numpy IRLS)",
            "terms": terms}
