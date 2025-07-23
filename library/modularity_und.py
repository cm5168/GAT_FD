import numpy as np

def modularity_und(A, gamma=1):
    """
    Optimal community structure and modularity (identical to MATLAB modularity_und)

    Inputs:
        A: undirected weighted/binary connection matrix (numpy array)
        gamma: resolution parameter (default=1)
    Outputs:
        Ci: optimal community structure (numpy array, 1-based labels)
        Q: maximized modularity (float)
    """
    A = np.array(A, dtype=float)
    N = A.shape[0]
    K = np.sum(A, axis=0)
    m = np.sum(K)
    B = A - gamma * np.outer(K, K) / m
    Ci = np.ones(N, dtype=int)
    cn = 1
    U = [1]
    ind = np.arange(N)
    Bg = B.copy()
    Ng = N
    while U:
        # Eigen decomposition
        D, V = np.linalg.eig(Bg)
        i1 = np.argmax(np.real(D))
        v1 = V[:, i1]
        S = np.ones(Ng)
        S[v1 < 0] = -1
        q = S.T @ Bg @ S
        if q > 1e-10:
            qmax = q
            np.fill_diagonal(Bg, 0)
            indg = np.ones(Ng)
            Sit = S.copy()
            while np.any(indg == 1):
                Qit = qmax - 4 * Sit * (Bg @ Sit)
                Qit_masked = Qit * indg
                imax = np.nanargmax(Qit_masked)
                qmax = Qit_masked[imax]
                Sit[imax] = -Sit[imax]
                indg[imax] = np.nan
                if qmax > q:
                    q = qmax
                    S = Sit.copy()
            if np.abs(np.sum(S)) == Ng:
                U.pop(0)
            else:
                cn += 1
                Ci[ind[S == 1]] = U[0]
                Ci[ind[S == -1]] = cn
                U = [cn] + U
        else:
            U.pop(0)
        if U:
            ind = np.where(Ci == U[0])[0]
            bg = B[np.ix_(ind, ind)]
            Bg = bg - np.diag(np.sum(bg, axis=0))
            Ng = len(ind)
    s = np.tile(Ci, (N, 1)).T
    Qmat = (s == s.T) * B / m
    Q = np.sum(Qmat)
    return Ci, Q