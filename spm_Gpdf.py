import numpy as np
from scipy.special import gammaln

def spm_Gpdf(x, h, l):
    """
    Gamma PDF as implemented in SPM.

    Parameters
    ----------
    x : array-like
        Gamma-variate (must be >= 0).
    h : array-like
        Shape parameter (must be > 0).
    l : array-like
        Scale parameter (must be > 0).

    Returns
    -------
    f : ndarray
        Probability density function of the Gamma distribution.
    """
    x = np.asarray(x, dtype=float)
    h = np.asarray(h, dtype=float)
    l = np.asarray(l, dtype=float)

    # Broadcast all to same shape
    x, h, l = np.broadcast_arrays(x, h, l)
    f = np.zeros_like(x)

    # Valid domain mask
    md = (h > 0) & (l > 0)
    if np.any(~md):
        f[~md] = np.nan
        import warnings
        warnings.warn("Returning NaN for out of range arguments")

    # Special cases at x == 0
    at_zero = (x == 0) & md

    # Case h < 1 => Inf
    f[at_zero & (h < 1)] = np.inf

    # Case h == 1 => l
    f[at_zero & (h == 1)] = l[at_zero & (h == 1)]

    # Case h > 1 => 0 (default remains zero)

    # General case where x > 0 and h,l valid
    Q = (x > 0) & md
    if np.any(Q):
        log_f = (h[Q] * np.log(l[Q]) + (h[Q] - 1) * np.log(x[Q]) - l[Q] * x[Q] - gammaln(h[Q]))
        f[Q] = np.exp(log_f)

    return f
