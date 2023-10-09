################################
##  Set diagnostics to True   ##
##  If you want display the   ##
##      Tracking process.     ##
################################

diagnostics = False

import os
import random

#################################
##      Import packages        ##
#################################
import sys
import time

import pytest
from numpy import *
from pylab import *

import trackeddy.generator.field_generator as fg
from trackeddy.generator.gaussian_field_functions import *
from trackeddy.savedata import *
from trackeddy.tracking import *

#################################
##   Import tools to create    ##
##     syntetic fields         ##
#################################


n = 2
a = 0.1
b = 0.1
t0 = 0
t = 1

xx = linspace(10, 12, 200)
yy = linspace(10, 12, 200)

gf = fg.Generate_field(a, b, n, xx, yy, "Nint")

data = abs(gf.assemble_field(t))

x = linspace(10, 13, 300)
y = linspace(10, 13, 300)

preferences = {"ellipse": 0.85, "eccentricity": 0.85, "gaussian": 0.8}
eddytd = {}
eddytdn = {}

levels = {"max": data.max(), "min": 0.1, "step": 0.1}
eddytd = analyseddyzt(
    data,
    x,
    y,
    t0,
    t,
    1,
    levels,
    preferences=preferences,
    areamap="",
    mask="",
    maskopt="forcefit",
    destdir="",
    physics="",
    diagnostics=False,
    plotdata=False,
    pprint=True,
    debug=False,
)


@pytest.mark.ttrackeddy_data
def test_data2npy():
    save_data("./test.npy", eddytd)
    assert os.path.isfile("./test.npy")


@pytest.mark.ttrackeddy_data
def test_tracknpy2nc():
    track2nc = Trackeddy2dataset("./test.npy", "./", "nc")
    track2nc.file2dict()
    track2nc.trackeddy2nc()
    assert os.path.isfile("./output_000.nc")


@pytest.mark.ttrackeddy_data
def test_trackdata2nc():
    track2nc = Trackeddy2dataset(eddytd, "./", "nc")
    track2nc.trackeddy2nc()
    assert os.path.isfile("./output_001.nc")


@pytest.mark.ttrackeddy_data
def test_rm_files():
    os.remove("./test.npy")
    os.remove("./output_000.nc")
    os.remove("./output_001.nc")
