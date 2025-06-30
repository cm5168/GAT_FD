import numpy as np
from scipy.special import gamma
def spm_Gpdf(u, h, l):
    """
    Gamma probability density function as used in SPM.
    u : array-like
        Time vector (non-negative).
    h : float
        Shape parameter.
    l : float
        Scale parameter.
    Returns
    -------
    pdf : ndarray
        Gamma PDF evaluated at u.
    """
    u = np.array(u, dtype=float)
    pdf = np.where(
        u >= 0,
        (u ** (h - 1) * np.exp(-u / l)) / (gamma(h) * (l ** h)),
        0.0
    )
    return pdf
def spm_hrf(RT, P=None, T=None):
    """
    Haemodynamic response function
    Parameters
    ----------
    RT : float
        Scan repeat time.
    P : array-like, optional
        Parameters of the response function [p1..p7].
        Defaults to [6, 16, 1, 1, 6, 0, 32].
    T : int, optional
        Microtime resolution. Defaults to 16.
    Returns
    -------
    hrf : ndarray
        Haemodynamic response function.
    p : ndarray
        Parameters used for the response function.
    """
    # Default parameters
    p = np.array([6, 16, 1, 1, 6, 0, 32], dtype=float)
    if P is not None:
        p[:len(P)] = P
    # Microtime resolution
    fMRI_T = T if T is not None else 16

    # Time step
    dt = RT / fMRI_T
    # Time vector (MATLAB: u = [0:ceil(p(7)/dt)] - p(6)/dt;)
    u = np.arange(0, np.ceil(p[6] / dt) + 1) - (p[5] / dt)
    # Mixture of Gammas (use correct parameter order for spm_Gpdf)
    response = spm_Gpdf(u, p[0] / p[2], 1.0 / (dt / p[2]))
    undershoot = spm_Gpdf(u, p[1] / p[3], 1.0 / (dt / p[3])) / p[4]
    hrf_full = response - undershoot
    # Downsample to 1:RT resolution (MATLAB: hrf([0:floor(p(7)/RT)]*fMRI_T + 1))
    idx = (np.arange(0, np.floor(p[6] / RT) + 1) * fMRI_T).astype(int)
    hrf = hrf_full[idx]
    # Normalize
    hrf = hrf / np.sum(hrf)
    return hrf, p
