import numpy as np

def degrees_und(CIJ):
    """
    Degree of nodes in an undirected graph.
    
    Parameters:
    - CIJ: 2D numpy array
        Undirected (binary or weighted) connection matrix.

    Returns:
    - deg: 1D numpy array
        Node degree. Weight information is discarded.
    """
    # Ensure CIJ is binary
    CIJ = (CIJ != 0).astype(float)

    # Degree is the sum of connections per node (row or column)
    deg = np.sum(CIJ, axis=0)  # or axis=1, since it's undirected and symmetric

    return deg
