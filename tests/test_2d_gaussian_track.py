################################
##  Set diagnostics to True   ##
##  If you want display the   ##
##      Tracking process.     ##
################################

diagnostics=False

#################################
##      Import packages        ##
#################################
import sys

from trackeddy.tracking import *
from trackeddy.datastruct import *
from trackeddy.geometryfunc import *
from trackeddy.init import *
from trackeddy.physics import *
from numpy import *
from pylab import *
import random
import pytest

#################################
##   Import tools to create    ##
##     syntetic fields         ##
#################################

from gaussian_field_functions import *

#################################
## Test 1: Check the detection ##
##    of a steady gaussian     ##
#################################

def gaussfit():
    '''
    Test the gaussian fitting.
    '''
    # Domain:
    x,y=linspace(-10,10,50),linspace(-10,10,50)
    X,Y=meshgrid(x,y)
    # Gaussian:
    gaussian=twoD_Gaussian((X,Y,1,0,0),  2.5, 2.5, 0, 0, 0 , 0)
    zz=gaussian.reshape(50,50)
    # Checking fitting:
    gausssianfitp,R2=fit2Dcurve(zz,(x,y,zz.max(),0,0),0,initial_guess='',\
                                date='',diagnostics=diagnostics)
    gaussianfit=twoD_Gaussian((X,Y,1,0,0),*gausssianfitp)
    return R2

@pytest.mark.test_trackeddy
def test_gaussfit():
    gausf=gaussfit()
    assert gausf >= 0.99

#################################
## Test 2: Check the detection ##
##   of a 2 gaussian (+ and +) ##
##   moving zonal direction    ##
#################################

def two_positive_gaussian_track():
    '''
    Test the tracking during certain timesteps of a simple gaussian.
    '''
    #Number of timesteps:
    time=40
    # Generation of a gaussian moving on the zonal direction:
    zz=moveGaussian(600,100,array([[x,x*0+450] for x in linspace(100,500,40)]),time)+\
       moveGaussian(600,100,array([[x,x*0+250] for x in linspace(100,500,40)]),time)
    # Replacing the coordinates of the Domain:
    lat=linspace(0,10,600)
    lon=linspace(0,10,600)

    # Tracking the eddy over the level 0.2 over 40 timesteps:

    gaussian=analyseddyzt(zz[:,:,:],lon,lat,0,40,1,0.2,areamap='',mask=''
                    ,maskopt='contour',destdir='',physics='',diagnostics=False,pprint=False)
    
    syntetic_gaussian=reconstruct_syntetic(shape(zz),lon,lat,gaussian)
     
    positive=len(gaussian['eddyn_0']['time'])
    return positive,time

@pytest.mark.test_trackeddy
def test_2eddydetection():
    positive,time=two_positive_gaussian_track()
    assert positive == time

#################################
## Test 3: Check the detection ##
##  of 4 gaussians (2+ and 2-) ##
##   moving zonal direction    ##
#################################

def multiple_gaussian_track():
    '''
    Test the tracking multiple gaussians. 
    '''
    #Number of timesteps:
    time=40
    # Generation of a gaussian moving on the zonal direction:
    zz=moveGaussian(600,50,array([[x,x*0+150] for x in linspace(100,500,40)]),40)+\
        moveGaussian(600,50,array([[x,x*0+350] for x in linspace(100,500,40)]),40)-\
        moveGaussian(600,50,array([[600-x,x*0+250] for x in linspace(100,500,40)]),40)-\
        moveGaussian(600,50,array([[600-x,x*0+450] for x in linspace(100,500,40)]),40)
    # Replacing the coordinates of the Domain:
    lat=linspace(0,10,600)
    lon=linspace(0,10,600)

    # Tracking the eddy over the level 0.2 over 40 timesteps:
    levels = {'max':zz[:,:,:].max(),'min':0.1,'step':0.25}
    gaussianp=analyseddyzt(zz[:,:,:],lon,lat,0,40,1,levels,areamap='',mask=''
                    ,maskopt='contour',destdir='',physics='',diagnostics=False,pprint=False)
    levels = {'max':zz[:,:,:].min(),'min':-0.1,'step':-0.25}
    gaussiann=analyseddyzt(zz[:,:,:],lon,lat,0,40,1,levels,areamap='',mask=''
                    ,maskopt='contour',destdir='',physics='',diagnostics=False,pprint=False)
    syntetic_gaussian=reconstruct_syntetic(shape(zz),lon,lat,gaussianp)+reconstruct_syntetic(shape(zz),lon,lat,gaussiann)
    
    positive=len(gaussianp['eddyn_0']['time'])
    negative=len(gaussiann['eddyn_0']['time'])
    return positive,negative,time

@pytest.mark.test_trackeddy
def test_3eddydetection():
    positive,negative,time=multiple_gaussian_track()
    assert positive == time and negative == time