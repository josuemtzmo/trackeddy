#################################
##      Import packages        ##
#################################
import numpy as np
import pytest

from trackeddy.trackeddy import *


def orientation(angles):
    sigma_x = 1.5
    sigma_y = 1
    amplitude = 1
    xo = 0
    yo = 0
    X = np.linspace(-5, 5, 300)
    Y = np.linspace(-5, 5, 300)
    x, y = np.meshgrid(X, Y)

    t = len(angles)
    data = np.zeros((t, 300, 300))

    count = 0
    for theta in angles:
        cos_phi = np.cos(theta)
        sin_phi = np.sin(theta)
        a = (cos_phi**2) / (2 * sigma_x**2) + (sin_phi**2) / (2 * sigma_y**2)
        b = (np.sin(2 * theta)) / (4 * sigma_x**2) - (np.sin(2 * theta)) / (
            4 * sigma_y**2
        )
        c = (sin_phi**2) / (2 * sigma_x**2) + (cos_phi**2) / (2 * sigma_y**2)
        data[count, :, :] = amplitude * np.exp(
            -(a * (x - xo) ** 2 + 2 * b * (x - xo) * (y - yo) + c * (y - yo) ** 2)
        )
        count = count + 1

    x = np.linspace(10, 12, 300)
    y = np.linspace(10, 12, 300)

    dataset = xr.Dataset(
        {"ssh": (["time", "lat", "lon"], data)}, coords={"lon": x, "lat": y}
    )

    TEddy = TrackEddy(dataset=dataset, variable="ssh")
    lin_levels = np.arange(-1, 1.5, 0.5)

    TEddy.filter = None

    track_in_time = TEddy.time_tracking(t0=0, tf=t, lin_levels=lin_levels, ntimes=5)

    return track_in_time


angle = [np.radians(31)]
allowed_error = np.radians(1)


@pytest.mark.parametrize(("angle", "allowed_error"), [(angle, allowed_error)])
@pytest.mark.trackeddy
def test_ellipse_orientation(angle, allowed_error):
    track_in_time = orientation(angle)
    identified_angle = track_in_time.xs(0, level="index").ellipse_params_theta.values
    assert abs(angle[0] - identified_angle) <= allowed_error


@pytest.mark.parametrize(("angle", "allowed_error"), [(angle, allowed_error)])
@pytest.mark.trackeddy
def test_gaussian_orientation(angle, allowed_error):
    track_in_time = orientation(angle)
    identified_angle = track_in_time.xs(0, level="index").gaussian_params_theta.values
    assert abs(angle[0] - identified_angle) <= allowed_error


# Test a diversity of angles
angle = np.linspace(0, np.pi - 0.01, 10)


@pytest.mark.parametrize(("angle", "allowed_error"), [(angle, allowed_error)])
@pytest.mark.trackeddy
def test_mult_orientation(angle, allowed_error):
    track_in_time = orientation(angle)
    identified_angle = track_in_time.xs(0, level="index").ellipse_params_theta
    identified_angles = np.round(identified_angle.values, 3)

    identified_angles = [
        id_angle if id_angle >= 0 else id_angle + np.pi
        for id_angle in identified_angles
    ]

    print(np.max(identified_angles - angle))

    assert np.max(identified_angles - angle) < allowed_error
