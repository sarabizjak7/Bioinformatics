import numpy as np

import scanpy as sc

def normalize(counts):
    """Regress out the effects of sequencing depth using analytic Pearson
    residuals for an NB model with fixed theta."""
    counts_sum0 = np.sum(counts, axis=0, keepdims=True)
    counts_sum1 = np.sum(counts, axis=1, keepdims=True)
    counts_sum = np.sum(counts)

    # Get residuals
    theta = 100
    mu = counts_sum1 @ counts_sum0 / counts_sum
    z = (counts - mu) / np.sqrt(mu + mu ** 2 / theta)

    # Clip to sqrt(n)
    n = counts.shape[0]
    z[z > np.sqrt(n)] = np.sqrt(n)
    z[z < -np.sqrt(n)] = -np.sqrt(n)

    return z
