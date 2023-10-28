import numpy as np
import pytest

from trackeddy.generator.gaussian_field_functions import *
from trackeddy.geometry import *
from trackeddy.trackeddy import *

#################################
## Test 1: Check the detection ##
##   of a 2 gaussian (+ and +) ##
##   moving zonal direction    ##
#################################


def two_positive_gaussian_track():
    """
    Test the tracking during certain timesteps of a simple gaussian.
    """
    # Number of timesteps:
    time = 40
    # Generation of a gaussian moving on the zonal direction:
    data = moveGaussian(
        600, 100, np.array([[x, x * 0 + 450] for x in np.linspace(100, 500, 40)]), time
    ) - moveGaussian(
        600, 100, np.array([[x, x * 0 + 250] for x in np.linspace(100, 500, 40)]), time
    )
    # Replacing the coordinates of the Domain:
    y = np.linspace(0, 10, 600)
    x = np.linspace(0, 10, 600)

    # Tracking the eddy over the level 0.2 over 40 timesteps:

    dataset = xr.Dataset(
        {"ssh": (["time", "lat", "lon"], data)}, coords={"lon": x, "lat": y}
    )

    TEddy = TrackEddy(dataset=dataset, variable="ssh")

    lin_levels = np.arange(-0.5, 0.6, 0.5)

    TEddy.filter = None

    track_in_time = TEddy.time_tracking(t0=0, tf=time, lin_levels=lin_levels, ntimes=5)

    eddy_0 = track_in_time.loc[0].xs(0, level="index").index
    eddy_1 = track_in_time.loc[1].xs(0, level="index").index

    return eddy_0, eddy_1


@pytest.mark.trackeddy
def test_2eddy_detection():
    eddy_0, eddy_1 = two_positive_gaussian_track()
    assert len(eddy_0) == len(eddy_1)


#################################
## Test 2: Check the detection ##
##  randomly moving gaussians  ##
##                             ##
#################################
# TODO add more tests here


#################################
## Test 3: Check discontinuous ##
##  track of randomly moving   ##
##           gaussians         ##
#################################
