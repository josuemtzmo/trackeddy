#################################
##      Import packages        ##
#################################
import numpy as np
import pytest
import xarray as xr

import trackeddy.generator.field_generator as fg
from trackeddy.generator.gaussian_field_functions import *
from trackeddy.trackeddy import *

#################################
##   Import tools to create    ##
##     syntetic fields         ##
#################################


#################################
## Test 1: Check the detection ##
##    of a steady gaussian     ##
#################################


def gauss_n_fit(n, size, res, gap, maxlat=20):
    """
    Test the number of eddies identified during 40 timesteps in a random walker gaussian field.
    """
    a = size
    b = size
    t = 1

    xx = np.linspace(10, maxlat, res)
    yy = np.linspace(10, maxlat, res)

    gf = fg.Generate_field(a, b, n, xx, yy, "Nint")

    data = gf.assemble_field(t, gap)

    x = np.linspace(10, maxlat, res + gap * 2)
    y = np.linspace(10, maxlat, res + gap * 2)

    dataset = xr.Dataset(
        {"ssh": (["time", "lat", "lon"], data)}, coords={"lon": x, "lat": y}
    )
    TEddy = TrackEddy(dataset=dataset, variable="ssh")

    lin_levels = np.arange(-1, 1.5, 0.5)

    TEddy.filter = None

    track_in_time = TEddy.time_tracking(t0=0, tf=t, lin_levels=lin_levels, ntimes=5)

    return track_in_time


@pytest.mark.parametrize(
    ("n", "size", "res", "gap"),
    [
        (100, 0.1, 400, 50),
        (7, 0.4, 2000, 300),
    ],
)
@pytest.mark.trackeddy
def test_gauss_large_n_fit(n, size, res, gap):
    track_in_time = gauss_n_fit(n, size, res, gap)
    identified_eddies = len(track_in_time.index.levels[0])
    assert identified_eddies == n


@pytest.mark.trackeddy
def test_gauss_extra_large():
    track_in_time = gauss_n_fit(2, 10, 300, 100, 80)
    assert track_in_time.empty


def gauss_mult_n_fit(n, t):
    """
    Test the number of eddies identified during 40 timesteps in a random walker gaussian field.
    """
    a = 0.07
    b = 0.07
    t0 = 0

    xx = np.linspace(10, 12.5, 300)
    yy = np.linspace(10, 12.5, 300)

    gf = fg.Generate_field(a, b, n, xx, yy, "Nint")

    data = gf.assemble_field(1)
    data = np.repeat(data, t)
    data = data.reshape(400, 400, t).T

    x = np.linspace(10, 13.5, 400)
    y = np.linspace(10, 13.5, 400)

    dataset = xr.Dataset(
        {"ssh": (["time", "lat", "lon"], data)}, coords={"lon": x, "lat": y}
    )
    TEddy = TrackEddy(dataset=dataset, variable="ssh")

    lin_levels = np.arange(-1, 1.5, 0.5)

    TEddy.filter = None

    track_in_time = TEddy.time_tracking(t0=0, tf=t, lin_levels=lin_levels, ntimes=5)

    track_in_time.index.get_level_values(level=0), track_in_time.index.get_level_values(
        level=1
    )

    return len(track_in_time.index.levels[0]) * len(track_in_time.index.levels[1])


@pytest.mark.parametrize(
    ("n", "t"),
    [
        (3, 10),
        (4, 9),
        (6, 7),
        (9, 4),
        (10, 3),
        (12, 1),
        (7, 1),
    ],
)
@pytest.mark.trackeddy
def test_gauss_mult_n_fit(n, t):
    assert gauss_mult_n_fit(n, t) == n * t
