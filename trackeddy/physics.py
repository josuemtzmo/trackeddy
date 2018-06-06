import numpy as np
import seawater as sw
from sympy.physics.vector import curl
import numpy.ma as ma
from trackeddy.init import *

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
    distmlat=sw.dist(0,lat,'km')[0][:]*1000
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
    distmlat=sw.dist(0,lat,'km')[0][:]*1000
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


    
#def vorticity3D(u,v,w,z):
#    dv_x=np.gradient(v[z,:,:],axis=1)
#    du_y=np.gradient(u[z,:,:],axis=0)
#    w=dv_x-du_y
#    return w

def geovelfield(ssha,lon,lat,mask=None,anomval=100):
    try:
        ma.filled(ssha,np.nan)
    except:
        pass
    distmlon=sw.dist(0,lon,'km')[0][:]*1000
    distmlat=sw.dist(0,lat,'km')[0][:]*1000
    mlon=np.cumsum(distmlon)
    mlat=np.cumsum(distmlat)
    dy=np.gradient(mlat)
    dx=np.gradient(mlon)
    detay,detax=np.gradient(ssha)
    omega = 7.2921e-5
    g=9.81
    f=2*omega*np.sin(np.deg2rad(lat))
    u=np.zeros(np.shape(ssha))
    v=np.zeros(np.shape(ssha))
    for ii in range(np.shape(ssha)[1]-1):
        detaxdy=detax[:,ii]/dx[ii]
        v[:,ii]=(g/f)*(detaxdy)
    for jj in range(np.shape(ssha)[0]-1):
        detaydx=detay[jj,:]/dy[jj]
        u[jj,:]=-(g/f[jj])*(detaydx)
    u[u>anomval]=np.nan
    v[v>anomval]=np.nan
    u[u<-anomval]=np.nan
    v[v<-anomval]=np.nan
    if mask != None:
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
    omega=7.2921159e-5
    f=2*omega*np.sin(np.deg2rad(lat))
    return f
    
def rossbyR(lat,g=9.81, D=3688):
    '''
    **************rossbyR***************
    Barotropic Rossby radius
    '''
    f = coriolis(lat)
    Lr = np.sqrt(g*D)/f
    return abs(Lr)

def ssh2ke(data,x,y,mask,anomval=100):
    u,v = geovelfield(data,x,y,mask,anomval)
    return KE(u,v)

def checkmesoscalearea(checkarea,lat,ellipsex,ellipsey,contourx,contoury):
    if checkarea==True:
        areachecker=rossbyR(np.mean(lat),g=9.81, D=3688)**2
        ellipsarea=sw.dist(0,[ellipsex.max(),ellipsex.min()],'km')[0][:]*1000*\
                   sw.dist(0,[ellipsey.max(),ellipsey.min()],'km')[0][:]*1000
        contarea=sw.dist(0,[contourx.max(),contourx.min()],'km')[0][:]*1000*\
                 sw.dist(0,[contoury.max(),contoury.min()],'km')[0][:]*1000
    else:
        areachecker=np.inf
        ellipsarea=sw.dist(0,[ellipsex.max(),ellipsex.min()],'km')[0][:]*1000*\
                   sw.dist(0,[ellipsey.max(),ellipsey.min()],'km')[0][:]*1000
        contarea=sw.dist(0,[ellipsex.max(),ellipsex.min()],'km')[0][:]*1000*\
                   sw.dist(0,[ellipsey.max(),ellipsey.min()],'km')[0][:]*1000
    return areachecker,ellipsarea,contarea