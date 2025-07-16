import numpy as np

def distance_inv(A):
    """
    Compute the inverse‑distance matrix D for binary adjacency A,
    following the MATLAB distance_inv subfunction.
    """
    n = A.shape[0]
    Lpath = A.copy().astype(int)     # paths of length l
    D = A.copy().astype(int)         # distance matrix
    l = 1

    # grow paths until no new connections appear
    while True:
        l += 1
        Lpath = Lpath @ A
        # new shortest paths of length l where D==0
        idx = (Lpath != 0) & (D == 0)
        if not idx.any():
            break
        D[idx] = l

    # assign inf to disconnected pairs and diagonal
    # note: np.eye gives floats, so cast mask
    mask_diag = np.eye(n, dtype=bool)
    D[(D == 0) | mask_diag] = np.inf

    return 1.0 / D  # invert

def efficiency_bin(A, local=False):
    """
    Global or local efficiency of binary undirected network A.

    Parameters
    ----------
    A : (n,n) array_like
        Binary adjacency matrix (0/1), directed or undirected.
    local : bool, optional
        If False (default), compute global efficiency (scalar).
        If True, compute local efficiency (n‑vector).

    Returns
    -------
    E : float or ndarray
        Global efficiency (float) if local=False, otherwise local efficiency
        as a length‑n numpy array.
    """
    A = np.array(A, copy=True)
    n = A.shape[0]

    # clear self‑loops & binarize
    np.fill_diagonal(A, 0)
    A = (A != 0).astype(int)

    if local:
        E = np.zeros(n, dtype=float)
        for u in range(n):
            # neighbors including outgoing or incoming
            neighbors = np.where((A[u, :] == 1) | (A[:, u] == 1))[0]
            if len(neighbors) < 2:
                continue

            # build the induced subgraph
            subA = A[np.ix_(neighbors, neighbors)]
            # symmetrized adjacency vector for denominator
            sa = A[u, neighbors] + A[neighbors, u]
            # inverse‑distance matrix on subgraph
            e = distance_inv(subA)
            se = e + e.T

            numer = 0.5 * np.sum((sa[:, None] * sa[None, :]) * se)
            denom = sa.sum()**2 - np.sum(sa**2)

            if denom > 0:
                E[u] = numer / denom
        return E

    else:
        e = distance_inv(A)
        # sum over all pairs, normalize by n*(n-1)
        return np.nansum(e) / (n**2 - n)
