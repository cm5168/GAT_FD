import numpy as np

def charpath(D, diagonal_dist=0, infinite_dist=1):
    """
    Characteristic path length, global efficiency and related statistics

    Parameters
    ----------
    D : np.ndarray
        Distance matrix
    diagonal_dist : int, optional
        Include distances on the main diagonal (default: 0)
    infinite_dist : int, optional
        Include infinite distances in calculation (default: 1)

    Returns
    -------
    lambda_ : float
        Network characteristic path length
    efficiency : float
        Network global efficiency
    ecc : np.ndarray
        Nodal eccentricity
    radius : float
        Network radius
    diameter : float
        Network diameter
    """
    n = D.shape[0]
    if np.isnan(D).any():
        raise ValueError('The distance matrix must not contain NaN values')
    D = D.copy()
    if (not diagonal_dist) or (diagonal_dist is None):
        D[np.arange(n), np.arange(n)] = np.nan  # set diagonal distance to NaN
    if (infinite_dist is not None) and (not infinite_dist):
        D[np.isinf(D)] = np.nan  # ignore infinite path lengths

    Dv = D[~np.isnan(D)]  # get non-NaN indices of D

    # Mean of entries of D(G)
    lambda_ = np.mean(Dv)

    # Efficiency: mean of inverse entries of D(G)
    efficiency = np.mean(1.0 / Dv)

    # Eccentricity for each vertex
    ecc = np.nanmax(D, axis=1)

    # Radius of graph
    radius = np.min(ecc)

    # Diameter of graph
    diameter = np.max(ecc)

    return lambda_, efficiency, ecc, radius, diameter
