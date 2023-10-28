import math

import numpy as np
from scipy.optimize import least_squares

from trackeddy import _cntr as cntr


def extract_contours(X, Y, data, level) -> set:
    """
    extract_contours Extract contours using the legacy matplotlib function.
    For more info look at:
        - https://github.com/matplotlib/legacycontour

    Parameters
    ----------
    X : array
        2D X coordinate of data
    Y : array
        2D Y coordinate of data
    data : array
        Data to extract contour
    level : float
        Level of contour

    Returns
    -------
    set
        Arrays containing the contour and type of segment of the contour
    """
    # Convert to arrays to use the cntr library
    c = cntr.Cntr(X, Y, data)
    # Extract contours in the level
    res = c.trace(level)
    # result is a list of arrays of vertices and path codes
    # (see docs for matplotlib.path.Path)

    nseg = len(res) // 2
    segments = res[:nseg]
    codes = res[nseg:]

    return segments, codes


class Fit_Surface:
    def __init__(self, eddy, eddy_data, X, Y, mode="gaussian") -> None:
        self.X = X
        self.Y = Y
        self.eddy_max = eddy.eddy_maxima[0]
        self.x_eddy = eddy.eddy_maxima[2]
        self.y_eddy = eddy.eddy_maxima[1]
        self.level = eddy.level

        # Initial guess to start fitting feature
        self.initial_guess = np.hstack((self.x_eddy, self.y_eddy, eddy.ellipse_params))
        # TODO is the best to pass the coords from the eddy_data or
        # should I provide the self.X and self.Y cropped? In fact, it may not matter.

        self.data = eddy_data  # - eddy.level
        self.mode = mode

    def _fitting(self):
        coords = (self.X, self.Y, self.eddy_max)

        fitdict = self.fit_curve(coords)

        fitted_curve = gaussian(coords, *fitdict)
        fitted_data = fitted_curve.reshape(*self.X.shape)

        return fitted_data, fitdict

    def construct_bounds(self):
        xdelta = np.mean(np.diff(self.X, axis=1))
        ydelta = np.mean(np.diff(self.Y, axis=0))

        LCoordbounds = np.array(
            [[self.initial_guess[0] - xdelta], [self.initial_guess[1] - ydelta]]
        )
        UCoordbounds = np.array(
            [[self.initial_guess[0] + xdelta], [self.initial_guess[1] + ydelta]]
        )

        # Limit the bounds for the initial guess to 70% of their ellipse value.

        Lellipsebounds = np.array(
            [
                [init - 0.7 * init] if init != 0 else [-1e-4]
                for init in self.initial_guess[2:]
            ]
        )
        Uellipsebounds = np.array(
            [
                [init + 0.7 * init] if init != 0 else [1e-4]
                for init in self.initial_guess[2:]
            ]
        )

        Lbounds = np.vstack((LCoordbounds, Lellipsebounds))
        Ubounds = np.vstack((UCoordbounds, Uellipsebounds))

        bounds = np.hstack((Lbounds, Ubounds)).T
        bounds.sort(axis=0)

        return bounds

    def fit_curve(self, coords):
        bounds = self.construct_bounds()

        res = least_squares(
            gaussian_residual,
            self.initial_guess,
            args=(coords, self.data),
            bounds=bounds,
        )

        fitdict = np.array(res.x)

        return fitdict


def gaussian(coords, xo, yo, sigma_x, sigma_y, theta):
    """
    *************** gaussian *******************
    Build a 2D gaussian.
    Notes:
        Remmember to do g.ravel().reshape(len(x),len(y)) for plotting purposes.
    Args:
        coords [x,y] (list|array): Coordinates in x and y.
        amplitude (float): Amplitud of gaussian.
        x0 , yo (float): Center of Gausian.
        sigma_x,sigma_y (float): Deviation.
        theta (Float): Orientation.
        offset (Float): Gaussian Offset.
    Returns:
        g.ravel() (list|array) - Gaussian surface in a list.
    Usage:
        Check scan_eddym function.
    """
    x = coords[0]
    y = coords[1]
    amplitude = coords[2]

    cos_phi = np.cos(theta)
    sin_phi = np.sin(theta)
    a = (cos_phi**2) / (2 * sigma_x**2) + (sin_phi**2) / (2 * sigma_y**2)
    b = (np.sin(2 * theta)) / (4 * sigma_x**2) - (np.sin(2 * theta)) / (
        4 * sigma_y**2
    )
    c = (sin_phi**2) / (2 * sigma_x**2) + (cos_phi**2) / (2 * sigma_y**2)
    g = amplitude * np.exp(
        -(a * (x - xo) ** 2 + 2 * b * (x - xo) * (y - yo) + c * (y - yo) ** 2)
    )

    return g.ravel()


def gaussian_residual(popt, coords, data2fit):
    gauss = gaussian(coords, *popt).reshape(np.shape(data2fit))

    # residual = abs(data2fit - gauss)

    # residual = np.nansum(abs(data2fit - gauss))
    # residual = np.nanmean(abs(data2fit - gauss))
    # This seems to outperform the exp()
    residual = np.nanstd(abs(data2fit - gauss)) * np.nanmean(abs(data2fit - gauss))

    # residual = np.exp(np.float64(np.nanmean(abs(data2fit - gauss)))) - 1
    return residual


def fit_ellipse(x, y):
    """
    **************** fit_ellipse *****************
    Fitting of an ellipse to an array of positions.

    Function translated form Matlab to python by Josue Martinez Moreno,
    the original source:
    Copyright (c) 2003, Ohad Gal
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in
    the documentation and/or other materials provided with the distribution
    For more information go to the main source:
    https://www.mathworks.com/matlabcentral/fileexchange/3215-fit-ellipse?requestedDomain=www.mathworks.com
    Notes:

    Args:
        x,y (array): Coordinates of the datapoints to fit an ellipse.
        diagnostics (boolean): Used to display all the statistics and
        plots to identify bugs.
    Returns:
        ellipse_t (dict) - This dictionary contains useful parameters
        describing completly the ellipsoid ajusted.
        status (boolean) - This value will be true if and only if the
        fit corresponds to a ellipse.
    Usage:
    R = np.arange(0,2*pi, 0.01)
    x = 1.5*np.cos(R) + 2 + 0.1*np.random.rand(len(R))
    y = np.sin(R) + 1. + 0.1*np.random.rand(len(R))
    ellipse,status=fit_ellipse(x,y,diagnostics=False)
    """
    orientation_tolerance = 1e-3

    x = x[:]
    y = y[:]

    mean_x = np.mean(x)
    mean_y = np.mean(y)

    xp = x - mean_x
    yp = y - mean_y

    X = np.array([xp**2, xp * yp, yp**2, xp, yp]).T

    a = np.sum(X, axis=0)
    b = np.dot(X.T, X)

    x2 = np.linalg.solve(b.T, a.T)

    a, b, c, d, e = x2

    if b == 0 and a < c:
        anglexaxis_rad = 0
    elif b == 0 and c < a:
        anglexaxis_rad = np.pi / 2
    else:
        anglexaxis_rad = np.arctan((c - a - np.sqrt((a - c) ** 2 + b**2)) / b)

    if min(abs(b / a), abs(b / c)) > orientation_tolerance:
        # TODO: Replace this non sign definite orientation_rad for anglexaxis
        # which is a sign definite.
        orientation_rad = 1 / 2 * np.arctan(b / (c - a))
        cos_phi = np.cos(orientation_rad)
        sin_phi = np.sin(orientation_rad)
        a, b, c, d, e = [
            a * cos_phi**2 - b * cos_phi * sin_phi + c * sin_phi**2,
            0,
            a * sin_phi**2 + b * cos_phi * sin_phi + c * cos_phi**2,
            d * cos_phi - e * sin_phi,
            d * sin_phi + e * cos_phi,
        ]
        mean_x, mean_y = (
            cos_phi * mean_x - sin_phi * mean_y,
            sin_phi * mean_x + cos_phi * mean_y,
        )
    else:
        orientation_rad = 0
        cos_phi = np.cos(orientation_rad)
        sin_phi = np.sin(orientation_rad)

    # final ellipse parameters
    X0 = mean_x - (d / 2) / a
    Y0 = mean_y - (e / 2) / c
    F = 1 + (d**2) / (4 * a) + (e**2) / (4 * c)
    a = np.sqrt(abs(F / a))
    b = np.sqrt(abs(F / c))

    # rotate the axes backwards to find the center point of the original TILTED ellipse
    R = np.array([[cos_phi, sin_phi], [-sin_phi, cos_phi]])

    theta_r = np.linspace(0, 2 * np.pi, len(y))
    ellipse_x_r = X0 + a * np.cos(theta_r)
    ellipse_y_r = Y0 + b * np.sin(theta_r)
    rotated_ellipse = np.dot(R, np.array([ellipse_x_r, ellipse_y_r]))

    # Ensure that a is always larger than b
    b, a = np.sort([a, b])

    ellipse_params = np.array((a, b, anglexaxis_rad))

    return rotated_ellipse.T, ellipse_params


def arc_length(curve):
    """
    Args:
    points: type arrays two values [[x, y], [x, y]]
    Returns:
    acc_length: curve length
    Descriptions:
    Calculate the length of the curve
    """

    acc_length = 0
    for i in range(0, len(curve) - 1):
        acc_length += math.dist(curve[i], curve[i + 1])

    return acc_length


def compute_similarity(curve1, curve2):
    geo_avg_curve_len = math.sqrt(arc_length(curve1) * arc_length(curve2))

    freshet_dist = frdist(curve1, curve2)

    result = max(1 - freshet_dist / (geo_avg_curve_len / math.sqrt(2)), 0)
    return np.array(round(result, 4))


def _c(ca, i, j, p, q):
    if ca[i, j] > -1:
        return ca[i, j]
    elif i == 0 and j == 0:
        ca[i, j] = np.linalg.norm(p[i] - q[j])
    elif i > 0 and j == 0:
        ca[i, j] = max(_c(ca, i - 1, 0, p, q), np.linalg.norm(p[i] - q[j]))
    elif i == 0 and j > 0:
        ca[i, j] = max(_c(ca, 0, j - 1, p, q), np.linalg.norm(p[i] - q[j]))
    elif i > 0 and j > 0:
        ca[i, j] = max(
            min(
                _c(ca, i - 1, j, p, q),
                _c(ca, i - 1, j - 1, p, q),
                _c(ca, i, j - 1, p, q),
            ),
            np.linalg.norm(p[i] - q[j]),
        )
    else:
        ca[i, j] = float("inf")

    return ca[i, j]


def frdist(P, Q):
    """
    Computes the discrete Fréchet distance between
    two curves. The Fréchet distance between two curves in a
    metric space is a measure of the similarity between the curves.
    The discrete Fréchet distance may be used for approximately computing
    the Fréchet distance between two arbitrary curves,
    as an alternative to using the exact Fréchet distance between a polygonal
    approximation of the curves or an approximation of this value.

    This is a Python 3.* implementation of the algorithm produced
    in Eiter, T. and Mannila, H., 1994. Computing discrete Fréchet distance.
    Tech. Report CD-TR 94/64, Information Systems Department, Technical
    University of Vienna.
    http://www.kr.tuwien.ac.at/staff/eiter/et-archive/cdtr9464.pdf

    Function dF(P, Q): real;
        input: polygonal curves P = (u1, . . . , up) and Q = (v1, . . . , vq).
        return: δdF (P, Q)
        ca : array [1..p, 1..q] of real;
        function c(i, j): real;
            begin
                if ca(i, j) > −1 then return ca(i, j)
                elsif i = 1 and j = 1 then ca(i, j) := d(u1, v1)
                elsif i > 1 and j = 1 then ca(i, j) := max{ c(i − 1, 1), d(ui, v1) }
                elsif i = 1 and j > 1 then ca(i, j) := max{ c(1, j − 1), d(u1, vj) }
                elsif i > 1 and j > 1 then ca(i, j) :=
                max{ min(c(i − 1, j), c(i − 1, j − 1), c(i, j − 1)), d(ui, vj ) }
                else ca(i, j) = ∞
                return ca(i, j);
            end; /* function c */

        begin
            for i = 1 to p do for j = 1 to q do ca(i, j) := −1.0;
            return c(p, q);
        end.

    Parameters
    ----------
    P : Input curve - two dimensional array of points
    Q : Input curve - two dimensional array of points

    Returns
    -------
    dist: float64
        The discrete Fréchet distance between curves `P` and `Q`.

    Examples
    --------
    >>> from frechetdist import frdist
    >>> P=[[1,1], [2,1], [2,2]]
    >>> Q=[[2,2], [0,1], [2,4]]
    >>> frdist(P,Q)
    >>> 2.0
    >>> P=[[1,1], [2,1], [2,2]]
    >>> Q=[[1,1], [2,1], [2,2]]
    >>> frdist(P,Q)
    >>> 0
    """

    p = P if len(P) >= len(Q) else Q
    q = Q if len(P) >= len(Q) else P

    p = np.array(p, np.float64)
    q = np.array(q, np.float64)

    len_p = len(p)
    len_q = len(q)

    # Force comparison between curves smaller than 300 points to avoid
    # RecursionError in _c
    c = 2
    while len_p > 300 and len_q > 300:
        p = p[0::c]
        q = q[0::c]
        len_p = len(p)
        len_q = len(q)
        c += 1

    ca = np.ones((len_p, len_q), dtype=np.float64) * -1

    dist = _c(ca, len_p - 1, len_q - 1, p, q)

    return dist


def eccentricity(a, b):
    """
    *************** eccentricity *******************
    This function calculate the eccentricity of a ellipse.
    Notes:

    Args:
        a (float): Mayor axis of an ellipse
        b (float): Minor axis of an ellipse
    Returns:
        eccen (float) - Eccentricity of the ellipsoid with parameters a,b.
    Usage:
        a=0.5
        b=0.3
        eccen=eccentricity(a,b)
    """
    a = abs(a)
    b = abs(b)

    eccen = np.sqrt(1 - (abs(b) ** 2 / abs(a) ** 2))
    return np.array(eccen)


def area_latlon_polygon(lonlat):
    # https://gis.stackexchange.com/questions/413349/calculating-area-of-lat-lon-polygons-without-transformation-using-geopandas

    contour_lon = lonlat[:, 0]
    contour_lat = lonlat[:, 1]

    lon_rad = np.deg2rad(contour_lon)
    lat_rad = np.deg2rad(contour_lat)

    a = np.sin(lat_rad / 2) ** 2 + np.cos(lat_rad) * np.sin(lon_rad / 2) ** 2
    colat = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    az = np.arctan2(np.cos(lat_rad) * np.sin(lon_rad), np.sin(lat_rad)) % (2 * np.pi)

    # Calculate diffs
    # daz = np.diff(az) % (2*pi)
    daz = np.diff(az)
    daz = (daz + np.pi) % (2 * np.pi) - np.pi

    deltas = np.diff(colat) / 2
    colat = colat[0:-1] + deltas

    # Perform integral
    integrands = (1 - np.cos(colat)) * daz

    # Integrate
    area = abs(sum(integrands)) / (4 * np.pi)

    radius = 6378137
    area_km = (area * 4 * np.pi * radius**2) * 1e-6  # km^2
    return np.array(area_km)
