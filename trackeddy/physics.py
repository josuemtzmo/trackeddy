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

def coriolis(lat,dist):
    omega = 7.292115e-5
    if lat < 5 and lat >-5:
        f = 2*omega*sin(lat)
    else:
        beta = 2*omega*cos(lat) / 6.3781e6 
        y = dist
        f = 2*omega*sin(lat) + beta*y
    return f
    
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

def checkscalearea(checkarea,lon,lat,ellipsex,ellipsey,contourx='',contoury=''):
    if checkarea==False:
        ellipsarea=None
        contarea=None
        areachecker=None
        
    elif checkarea==None or type(checkarea)==dict:
        if checkarea==None:
            checkarea={'mesoscale':2*np.pi}
            
        ellipsarea=gs.distance([[ellipsex.max()],[ellipsex.min()]],[[ellipsey.mean()],[ellipsey.mean()]],axis=0)[0][0]*\
                   gs.distance([[ellipsex.mean()],[ellipsex.mean()]],[[ellipsey.max()],[ellipsey.min()]],axis=0)[0][0]
        if contourx!='' or contoury!='':
            contarea=gs.distance([[contourx.max()],[contourx.min()]],[[contoury.mean()],[contoury.mean()]],axis=0)[0][0]*\
                   gs.distance([[contourx.mean()],[contourx.mean()]],[[contoury.max()],[contoury.min()]],axis=0)[0][0]
        else:
            contarea=None
        if len(checkarea.keys())==1:
            if 'mesoscale' in checkarea.keys():
                areachecker=(checkarea['mesoscale']*(rossbyR(np.mean(lon),np.mean(lat))))**2
            elif 'field' in checkarea.keys():
                try:
                    path=os.path.expanduser(os.path.dirname(os.path.realpath(__file__)))
                    area_file=xarray.open_mfdataset(path+checkarea['mesoscale']['path'])
                    areachecker = (RrD_file.RrD.sel(lon=[lon],lat=[lat],method='nearest').values) * checkarea['mesoscale']['factor']
                except:
                    areachecker=False
            elif 'constant' in checkarea.keys():
                if checkarea['constant']==np.inf or checkarea['constant']==None:
                    areachecker=None
                else:
                    areachecker = checkarea['constant']
            else:
                raise Exception('The Area Check dictionary should have one option: mesoscale, field or constant.')
        else:
            raise Exception('The Area Check dictionary should have only one option: mesoscale, field or constant.')
        
        if areachecker == None:
            areastatus={'status':True,'check':areachecker,'ellipse':ellipsarea,'contour':contarea}
        elif areachecker == False:
            areastatus={'status':False,'check':None,'ellipse':None,'contour':None}
        elif ellipsarea <= areachecker and contarea != None:
            if contarea <= areachecker:
                areastatus={'status':True,'check':areachecker,'ellipse':ellipsarea,'contour':contarea}
            else:
                areastatus={'status':False,'check':None,'ellipse':None,'contour':None}
        elif ellipsarea <= areachecker and contarea == None:
            areastatus={'status':True,'check':areachecker,'ellipse':ellipsarea,'contour':None}
        else:
            areastatus={'status':False,'check':None,'ellipse':None,'contour':None}
    else:
        raise Exception('Unexpected dictionary format. Check the Area Check documentation.')
    return areastatus


def phase_angle(v,u,lon,lat):
    phi=arctan(v,u)