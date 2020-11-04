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
from trackeddy.physics import *
from numpy import *
from pylab import *
import cmocean as cm
import random
import pytest


def orientation(angles):
    sigma_x=1.5
    sigma_y=1
    amplitude=1
    xo=0
    yo=0
    X=linspace(-5,5,300)
    Y=linspace(-5,5,300)
    x,y=meshgrid(X,Y)

    t=len(angles)
    data = zeros((t,300,300))

    count=0
    for theta in angles:
        cos_phi = np.cos(theta)
        sin_phi = np.sin(theta)
        a = (cos_phi**2)/(2*sigma_x**2) + (sin_phi**2)/(2*sigma_y**2)
        b = (np.sin(2*theta))/(4*sigma_x**2) - (np.sin(2*theta))/(4*sigma_y**2)
        c = (sin_phi**2)/(2*sigma_x**2) + (cos_phi**2)/(2*sigma_y**2)
        data[count,:,:] = amplitude*exp(-(a*(x-xo)**2 + 2*b*(x-xo)*(y-yo) + c*(y-yo)**2))
        count=count+1

    x = linspace(10,12,300)
    y = linspace(10,12,300)

    preferences={'ellipse':0.85,'eccentricity':0.85,'gaussian':0.8}
    eddytd={}
    eddytdn={}

    t0 = 0

    levels = {'max':0.05,'min':0.05,'step':0.05}
    eddytd = analyseddyzt(data,x,y,t0,t,1,levels,preferences=preferences,areamap='',mask='',maskopt='forcefit'\
                        ,destdir='',physics='',diagnostics=False,plotdata=False,pprint=False,debug=False)
    return eddytd

angle=np.array([31])*np.pi/180
allowed_error=1*np.pi/180
@pytest.mark.parametrize(('angle', 'allowed_error'), [(angle,allowed_error)])

@pytest.mark.ttrackeddy
def test_ellipse_orientation(angle,allowed_error):
    eddytd=orientation(angle)
    assert abs(angle[0] - eddytd['eddyn_0']['angle']) <= allowed_error

@pytest.mark.parametrize(('angle', 'allowed_error'), [(angle,allowed_error)])

@pytest.mark.ttrackeddy
def test_gaussian_orientation(angle,allowed_error):
    eddytd=orientation(angle)
    assert abs(angle[0] - eddytd['eddyn_0']['2dgaussianfit'][0][2]) <= allowed_error

# angle=np.linspace(0,2*np.pi,37)
# @pytest.mark.parametrize(('angle', 'allowed_error'), [
#     zip(angle,allowed_error+angle*0)])

# @pytest.mark.ttrackeddy
# def test_mult_orientation(n,t):
#     if 
#     assert  gauss_mult_n_fit(n,t) == n * t


# count=0
# for theta in angles:
#     cos_phi = np.cos(theta)
#     sin_phi = np.sin(theta)
#     a = (cos_phi**2)/(2*sigma_x**2) + (sin_phi**2)/(2*sigma_y**2)
#     b = (np.sin(2*theta))/(4*sigma_x**2) - (np.sin(2*theta))/(4*sigma_y**2)
#     c = (sin_phi**2)/(2*sigma_x**2) + (cos_phi**2)/(2*sigma_y**2)
#     data[count,:,:] = amplitude*exp(-(a*(x-xo)**2 + 2*b*(x-xo)*(y-yo) + c*(y-yo)**2))
#     plt.figure()
#     plt.title(round(math.degrees(theta),2))
#     plt.pcolormesh(X,Y,data[count,:,:])
#     plt.gca().set_aspect('equal', adjustable='box')
#     count=count+1

# for ii in range(t):
#     plt.figure(figsize=(6,3))
#     ax1=plt.subplot(1, 2, 1)
#     ax1.pcolormesh(x,y,data[ii,:,:])
#     ax1.plot(eddytd['eddyn_0']['ellipse'][ii][0],eddytd['eddyn_0']['ellipse'][ii][1])
#     ax1.set_title(round(math.degrees(angles[ii]),2))
#     plt.gca().set_aspect('equal', adjustable='box')
#     ax2=plt.subplot(1, 2, 2)
#     ax2.pcolormesh(x,y,pos_f[ii,:,:])
#     ax2.plot(eddytd['eddyn_0']['ellipse'][ii][0],eddytd['eddyn_0']['ellipse'][ii][1])
#     ax2.set_title(str(round(math.degrees(eddytd['eddyn_0']['angle'][ii][0]),2))+'/'+str(round(math.degrees(eddytd['eddyn_0']['2dgaussianfit'][ii][2]),2)))
#     ax2.set_xlabel(eddytd['eddyn_0']['level'][ii][0])
#     plt.gca().set_aspect('equal', adjustable='box')
    

# for ii in range(20):
#     print(eddytd['eddyn_0']['angle'][ii][0])
#     print(eddytd['eddyn_0']['2dgaussianfit'][ii][2])
# print(gf.eddies['eddy_n0']['angle'],eddytd['eddyn_0']['angle'])
# print(gf.eddies['eddy_n0']['angle'],eddytdn['eddyn_0']['angle'])