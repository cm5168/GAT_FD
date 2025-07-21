import numpy as np

def efficiency_bin(A, local=False):
    """
    EFFICIENCY_BIN     Global efficiency, local efficiency.

    Eglob = efficiency_bin(A)
    Eloc = efficiency_bin(A, local=True)

    The global efficiency is the average of inverse shortest path length,
    and is inversely related to the characteristic path length.

    The local efficiency is the global efficiency computed on the
    neighborhood of the node, and is related to the clustering coefficient.

    Inputs:     A,              binary undirected or directed connection matrix
                local,          optional argument
                                    local=0 computes global efficiency (default)
                                    local=1 computes local efficiency

    Output:     Eglob,          global efficiency (scalar)
                Eloc,           local efficiency (vector)
    """
    n = len(A)
    A = np.array(A, dtype=float)
    np.fill_diagonal(A, 0)  # clear diagonal
    A = (A != 0).astype(float)  # enforce double precision/binary

    if local:
        E = np.zeros(n)
        for u in range(n):
            V = np.where((A[u, :] != 0) | (A[:, u] != 0))[0]  # neighbors
            if V.size == 0:
                continue
            sa = A[u, V] + A[V, u]  # symmetrized adjacency vector
            e = distance_inv(A[np.ix_(V, V)])  # inverse distance matrix
            se = e + e.T  # symmetrized inverse distance matrix
            numer = np.sum((sa[:, None] @ sa[None, :]) * se) / 2  # numerator
            if numer != 0:
                denom = np.sum(sa) ** 2 - np.sum(sa ** 2)  # denominator
                if denom != 0:
                    E[u] = numer / denom  # local efficiency
        return E
    else:
        e = distance_inv(A)
        return np.sum(e) / (n ** 2 - n)


def distance_inv(A_):
    l = 1
    A_ = np.array(A_, dtype=float)
    n_ = len(A_)
    Lpath = A_.copy()
    D = A_.copy()

    Idx = (Lpath != 0) & (D == 0)
    while np.any(Idx):
        l += 1
        Lpath = Lpath @ A_
        Idx = (Lpath != 0) & (D == 0)
        D[Idx] = l

    D[(D == 0) | np.eye(n_, dtype=bool)] = np.inf  # assign inf to disconnected and diagonal
    D = 1.0 / D  # invert distance
    return D
