import numpy as np
from sympy.physics.vector import curl
import numpy.ma as ma
#from trackeddy.init import *
import gsw as gs
import seawater as sw
from netCDF4 import Dataset
import os
import xarray

def okuboweissparm(u,v,lat,lon,z):
    if z==0:
        du=np.gradient(u[:,:],axis=1)
        dv=np.gradient(v[:,:],axis=1)
        du=np.gradient(u[:,:],axis=0)
        dv=np.gradient(v[:,:],axis=0)
    else:
        du=np.gradient(u[z,:,:],axis=1)
        dv=np.gradient(v[z,:,:],axis=1)
        du=np.gradient(u[z,:,:],axis=0)
        dv=np.gradient(v[z,:,:],axis=0)
    distmlon=sw.dist(0,lon,'km')[0][:]*1000
    distmlat=sw.dist(lat,0,'km')[0][:]*1000
    mlon=np.cumsum(distmlon)
    mlat=np.cumsum(distmlat)
    mlon = np.hstack((mlon,mlon[-1]))
    mlat = np.hstack((mlat,mlat[-1]))
    dy=np.gradient(mlat)
    dx=np.gradient(mlon)
    dX,dY = np.meshgrid(dx,dy)
    
    du_dx=du/dX
    du_dy=du/dY
    dv_dx=dv/dX
    dv_dy=dv/dY
    sn=du_dx-dv_dy
    ss=dv_dx+du_dy
    w=vorticity2D(u,v,lon,lat)
    owparm=sn**2+ss**2-w**2
    return owparm
    
def vorticity2D(u,v,lon,lat):
    dv=np.gradient(v[:,:],axis=1)
    du=np.gradient(u[:,:],axis=0)
    distmlon=sw.dist(0,lon,'km')[0][:]*1000
    distmlat=sw.dist(lat,0,'km')[0][:]*1000
    mlon=np.cumsum(distmlon)
    mlat=np.cumsum(distmlat)
    mlon = np.hstack((mlon,mlon[-1]))
    mlat = np.hstack((mlat,mlat[-1]))
    dy=np.gradient(mlat)
    dx=np.gradient(mlon)
    dX,dY = np.meshgrid(dx,dy)
    dv_dx=dv/dX
    du_dy=du/dY
    w=dv_dx-du_dy    
    return w

def geovelfield(ssha,lon,lat,mask='',anomval=100):
    try:
        ma.filled(ssha,np.nan)
    except:
        pass
    Lon,Lat=np.meshgrid(lon,lat)
    u=np.zeros(np.shape(ssha))
    v=np.zeros(np.shape(ssha))
    print(np.shape(ssha))
    for ii in range(np.shape(ssha)[0]):
        v[ii,1:]=gs.geostrophic_velocity(ssha[ii,:], Lon[ii,:], Lat[ii,:], 0, axis=0)[0]*9.81
    for jj in range(np.shape(ssha)[1]):
        u[1:,jj]=gs.geostrophic_velocity(ssha[:,jj], Lon[:,jj], Lat[:,jj], 0, axis=1)[0]*9.81

    if mask != '':
        u= np.ma.masked_array(u, mask)
        v= np.ma.masked_array(v, mask)
    return u,v


    
def KE(u,v):
    ke=(1/2)*(u**2+v**2)
    return ke
    
def PVort(S,T,P,U,V):
    theta=np.zeros(np.shape(S))
    for ii in range(0,np.shape(S)[1]):
        for jj in range(0,np.shape(S)[2]):
            theta[:,ii,jj]=sw.ptmp(S[:,ii,jj],T[:,ii,jj],P)
    rho=sw.dens(S,T,P)
    omega=7.2921159*10**(-5)
    zeta=2*omega+curl(U,V)
    Q=(1/rho)*(zeta)*np.gradient(theta)

def coriolis(lat):
    return gs.f(lat)
    
def rossbyR(lon,lat):
    '''
    **************rossbyR***************
    Barotropic Rossby radius
    '''
    try:
        path=os.path.expanduser(os.path.dirname(os.path.realpath(__file__)))
        RrD_file=xarray.open_mfdataset(path+'/../input/rossby_g.nc')
        lon=round(lon,2)
        lat=round(lat,2)
        if lon<0:
            lon=360+lon
        RrD = RrD_file.RrD.sel(lon=[lon],lat=[lat],method='nearest').values
    except:
        RrD=[[np.inf]]
    return RrD[0][0]*1000

def ssh2ke(data,x,y,mask,anomval=100):
    u,v = geovelfield(data,x,y,mask,anomval)
    return KE(u,v)

def checkmesoscalearea(checkarea,lon,lat,ellipsex,ellipsey,contourx='',contoury=''):
    if checkarea==True:
        areachecker=(2*np.pi*(rossbyR(np.mean(lon),np.mean(lat))))**2
        ellipsarea=gs.distance([[ellipsex.max()],[ellipsex.min()]],[[ellipsey.mean()],[ellipsey.mean()]],axis=0)[0][0]*\
                   gs.distance([[ellipsey.mean()],[ellipsey.mean()]],[[ellipsey.max()],[ellipsey.min()]],axis=0)[0][0]
        if contourx!='' or contoury!='':
            contarea=gs.distance([[contoury.max()],[contoury.min()]],[[contoury.mean()],[contoury.mean()]],axis=0)[0][0]*\
                   gs.distance([[contoury.mean()],[contoury.mean()]],[[contoury.max()],[contoury.min()]],axis=0)[0][0]
        else:
            contarea=None
    else:
        areachecker=np.inf
        ellipsarea=gs.distance([[ellipsex.max()],[ellipsex.min()]],[[ellipsey.mean()],[ellipsey.mean()]],axis=0)[0][0]*\
                   gs.distance([[ellipsey.mean()],[ellipsey.mean()]],[[ellipsey.max()],[ellipsey.min()]],axis=0)[0][0]
        if contourx!='' or contoury!='':
            contarea=gs.distance([[contoury.max()],[contoury.min()]],[[contoury.mean()],[contoury.mean()]],axis=0)[0][0]*\
                   gs.distance([[contoury.mean()],[contoury.mean()]],[[contoury.max()],[contoury.min()]],axis=0)[0][0]
        else:
            contarea=None
    return areachecker,ellipsarea,contarea