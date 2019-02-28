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
import cmocean as cm
import random
import pytest

#################################
##   Import tools to create    ##
##     syntetic fields         ##
#################################

from trackeddy.utils.gaussian_field_functions import *
import trackeddy.utils.field_generator as fg

#################################
## Test 1: Check the detection ##
##    of a steady gaussian     ##
#################################

def gauss_n_fit(n):
    '''
    Test the number of eddies identified during 40 timesteps in a random walker gaussian field.
    '''
    a  = 0.1
    b  = 0.1
    t0 = 0
    t  = 1

    xx=linspace(10,12,200)
    yy=linspace(10,12,200)

    gf=fg.Generate_field(a,b,n,xx,yy,'Nint')   

    data = gf.assemble_field(t)

    x=linspace(10,13,300)
    y=linspace(10,13,300)
            
    preferences={'ellipse':0.85,'eccentricity':0.85,'gaussian':0.8}
    eddytd={}
    eddytdn={}

    try:
        levels = {'max':data.max(),'min':0.1,'step':0.1}
        eddytd = analyseddyzt(data,x,y,t0,t,1,levels,preferences=preferences,areamap='',mask='',maskopt='forcefit'\
                            ,destdir='',physics='',diagnostics=False,plotdata=False,pprint=True,debug=False)
    except:
        print("No positive")

    try:
        levels  = {'max':data.min(),'min':-0.1,'step':-0.1}
        eddytdn = analyseddyzt(data,x,y,t0,t,1,levels,preferences=preferences,areamap='',mask='',maskopt='forcefit'\
                    ,destdir='',physics='',diagnostics=False,plotdata=False,pprint=True,debug=False)
    except:
        print("No negative")
    
    return len(eddytd.keys())+len(eddytdn.keys())

@pytest.mark.ttrackeddy
def test_gauss_n_fit():
    assert gauss_n_fit(7) == 7 

def gauss_mult_n_fit(n,t):
    '''
    Test the number of eddies identified during 40 timesteps in a random walker gaussian field.
    '''
    a  = 0.07
    b  = 0.07
    t0 = 0

    xx=linspace(10,12.5,300)
    yy=linspace(10,12.5,300)

    gf=fg.Generate_field(a,b,n,xx,yy,'Nint')   

    data = gf.assemble_field(t)

    x=linspace(10,13.5,400)
    y=linspace(10,13.5,400)
            
    preferences={'ellipse':0.85,'eccentricity':0.85,'gaussian':0.8}
    eddytd={}
    eddytdn={}

    
    try:
        levels = {'max':data.max(),'min':0.1,'step':0.1}
        eddytd = analyseddyzt(data,x,y,t0,t,1,levels,preferences=preferences,areamap='',mask='',maskopt='forcefit'\
                            ,destdir='',physics='',diagnostics=False,plotdata=False,pprint=True,debug=False)
    except:
        print("No positive")

    try:
        levels  = {'max':data.min(),'min':-0.1,'step':-0.1}
        eddytdn = analyseddyzt(data,x,y,t0,t,1,levels,preferences=preferences,areamap='',mask='',maskopt='forcefit'\
                    ,destdir='',physics='',diagnostics=False,plotdata=False,pprint=True,debug=False)
    except:
        print("No negative")

    posn=sum([len(item['time']) for keys,item in eddytd.items()])
    negn=sum([len(item['time']) for keys,item in eddytdn.items()])

    # print(posn,negn)
    # print([item['time'] for keys,item in eddytd.items()])

    # for ii in range(shape(data)[0]):
    #     contourf(x,y,data[ii,:,:],vmin=-1,vmac=1,cmap=cm.cm.balance)
    #     colorbar()
    #     savefig("deleteme_%03d_n%d.png" % (ii,n))
    #     close()

    return posn+negn


@pytest.mark.parametrize(('n', 't'), [
    (3, 10),
    (4, 9),
    (5, 8),
    (6, 7),
    (7, 6),
    (8, 5),
    (9, 4),
    (10, 3),
    (11, 2),
    (12, 1),
])

@pytest.mark.ttrackeddy
def test_gauss_mult_n_fit(n,t):
    assert  gauss_mult_n_fit(n,t) == n * t