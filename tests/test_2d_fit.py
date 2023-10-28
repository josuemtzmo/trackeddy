#################################
##      Import packages        ##
#################################

import numpy as np
import pytest

from trackeddy.eddy import *
from trackeddy.geometry import *
from trackeddy.trackeddy import *

#################################
## Test 1: Check the detection ##
##    of a steady gaussian     ##
#################################


def gaussfit(loc):
    """
    Test the gaussian fitting.
    """
    # Domain:
    x, y = np.linspace(-10, 10, 50), np.linspace(-10, 10, 50)
    X, Y = np.meshgrid(x, y)
    # Gaussian:
    gauss = gaussian((X, Y, 1), loc, loc, 3, 3, 0)
    gauss = gauss.reshape(50, 50)
    # Fake eddy
    eddy = Eddy(None, None)
    eddy.eddy_maxima = [1, loc, loc]
    eddy.ellipse_params = (3, 3, 0)

    # Checking fitting:
    fitting = Fit_Surface(eddy, gauss, X, Y)

    # fitting.construct_bounds()
    fitted_gauss, fitdict = fitting._fitting()
    return np.sum(abs(gauss - fitted_gauss))


@pytest.mark.trackeddy
def test_gaussfit():
    gausf = gaussfit(1)
    assert gausf <= 1e-5


@pytest.mark.trackeddy
def test_gaussfit_origin():
    gausf = gaussfit(0)
    assert gausf <= 1e-5


# TODO add test when gaussian is not symmetric
