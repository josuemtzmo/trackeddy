import numpy as np


def approx_RD_lat(x) -> float:
    """
    approx_RD_lat

    Function to approximate the meridional average of the Rossby Radius of deformation

    Parameters
    ----------
    x : float
        latitude of interest

    Returns
    -------
    float
        Approximation of Rossby radius of deformation
    """
    b = 0
    c = 0.15
    y = 200 * np.exp(-((x - b) ** 2) / 2 * c**2) + 60 * np.exp(
        -((x - b) ** 2) / 40 * c**2
    )
    return y
