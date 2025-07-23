import numpy as np

def makerandCIJ_und(N, K):
    """
    Synthetic undirected random network (identical to MATLAB makerandCIJ_und)

    Inputs:
        N: number of vertices
        K: number of edges
    Output:
        CIJ: undirected random connection matrix (N x N)
    Note: no connections are placed on the main diagonal.
    """
    ind = np.triu(~np.eye(N, dtype=bool), 1)  # upper triangle, no diagonal
    i = np.flatnonzero(ind)
    rp = np.random.permutation(len(i))
    irp = i[rp]
    CIJ = np.zeros((N, N), dtype=int)
    CIJ.flat[irp[:K]] = 1
    CIJ = CIJ + CIJ.T  # symmetrize
    return CIJ