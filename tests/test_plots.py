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
from trackeddy.plotfunc import *
import numpy as np
import random
import pytest

@pytest.mark.ttrackeddy
def test_mapplot():
    x=np.linspace(0,360,360)
    y=np.linspace(-90,90,180)

    X,Y=np.meshgrid(x,y)

    jets = np.zeros((len(x),len(y)))
    r = Y
    k_y = random.uniform(2, 3)
    phase = random.uniform(0, 1)
    k_x = random.uniform(1, 2)
    amp = 0.3
    jets = amp*np.cos((k_y*(k_y*Y+phase+np.sin(k_x*X))))

    title='Test'

    basemap_mplot(x,y,jets,title,projection='ortho',lat_0=-90,\
        lon_0=-100,boundinglat=-30,scale='Lin',vmin='',vmax='',\
        xan=1,yan=1,figsize=(3,3),fontsize=15,dpi=50)

    basemap_mplot(x,y,jets,title,projection='ortho',lat_0=-90,\
        lon_0=-100,boundinglat=-30,scale='SymLog',vmin='',vmax='',\
        xan=1,yan=1,figsize=(3,3),fontsize=15,dpi=50)

    basemap_mplot(x,y,abs(jets),title,projection='ortho',lat_0=-90,\
        lon_0=-100,boundinglat=-30,scale='Log',vmin=0,vmax=10,\
        xan=1,yan=1,figsize=(3,3),fontsize=15,dpi=50)

    basemap_mplot(x,y,abs(jets)+1,title,projection='ortho',lat_0=-90,\
        lon_0=-100,boundinglat=-30,scale='Log2+',vmin=1,vmax=10,\
        xan=1,yan=1,figsize=(3,3),fontsize=15,dpi=50)

    basemap_mplot(x,y,[jets,jets],[title,title],projection='ortho',lat_0=-90,\
        lon_0=-100,boundinglat=-30,scale='Linc',vmin=-1,vmax=1,\
        xan=2,yan=1,figsize=(3,3),fontsize=15,dpi=50)
    
    basemap_mplot(x,y,[jets,jets],[title,title],projection='ortho',lat_0=-90,\
        lon_0=-100,boundinglat=-30,scale='SymLog',vmin=-1,vmax=1,\
        xan=1,yan=2,figsize=(3,3),fontsize=15,dpi=50)
    
    basemap_mplot(x,y,[jets,jets,jets,jets],[title,title,title,title],projection='ortho',lat_0=-90,\
        lon_0=-100,boundinglat=-30,scale='Lin',vmin=-1,vmax=1,\
        xan=2,yan=2,figsize=(3,3),fontsize=15,dpi=50)
    