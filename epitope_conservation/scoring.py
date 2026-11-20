import numpy as np
from scipy.spatial.distance import cdist
from scipy.optimize import least_squares


def fit_rmsd(coords1, coords2):
    """Return the Kabsch-aligned RMSD, or infinity if coordinate shapes differ."""
    c1 = np.asarray(coords1, dtype=float)
    c2 = np.asarray(coords2, dtype=float)
    if c1.shape != c2.shape or len(c1) == 0:
        return float("inf")
    center1 = c1.mean(axis=0)
    center2 = c2.mean(axis=0)
    x = c1 - center1
    y = c2 - center2
    u, _, vt = np.linalg.svd(x.T @ y)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    fitted = (rotation @ y.T).T + center1
    return float(np.sqrt(np.mean(np.sum((x - (fitted - center1)) ** 2, axis=1))))
