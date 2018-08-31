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

def dist(loc1,loc2):
    return sqrt((loc1[0]-loc2[0])**2 + (loc2[1]-loc1[1])**2)

def checkposition(eddies,x,y,):
    for key,item in eddies.items():
        for key1,item1 in eddies.items():
            xc=item['loc'][0][0]
            yc=item['loc'][0][1]
            xc1=item1['loc'][0][0]
            yc1=item1['loc'][0][1]
            distance=dist([x[xc],y[yc]],[x[xc1],y[yc1]])
            #print(distance)
            checker = (distance < 5*a or distance < 5*b ) and key1!=key
            #print(checker)
            while checker:
                newx=randint(0,xlen)
                newy=randint(0,ylen)
                eddies[key1]={'loc':[[newx, newy]],'grow':True,'radius':a,\
                             'amp':random.choice([-1,1])}
                xc1=newx
                yc1=newy
                distance=dist([x[xc],y[yc]],[x[xc1],y[yc1]])
                checker = (distance < 5*a or distance < 5*b ) and key1!=key

def make_random_walk(indexs, steps):
    move_dict = {1: go_up,
                 2: go_right,
                 3: go_left,
                 4: go_down,
                 5: go_downleft,
                 6: go_downright,
                 7: go_upleft,
                 8: go_upright,
                 }
    for _ in range(steps):
        for ii in indexs:
            move_in_a_direction = move_dict[random.randint(1, 8)]
            movcood=move_in_a_direction(ii,5)
    return indexs[0]+movcood[0],indexs[1]+movcood[1]

def gaussfit():
    '''
    Test the number of eddies identified during 40 timesteps in a random walker gaussian field.
    '''
    
    n=7
    steps=1
    x=linspace(0,1,100)
    y=linspace(0,1,100)

    xlen=len(x)
    ylen=len(y)

    count=0

    X,Y=meshgrid(x,y)
    a=b=0.1
    count=a
    data=zeros(shape(X))
    diagnostics=True
    rmlist=[]

    eddies={'eddy_n%s' % ii:{'loc':[[randint(0,xlen), randint(0,ylen)]],'grow':True,'radius':a,\
                             'amp':random.choice([-1,1])} for ii in range(n)}
    
    checkposition(eddies,x,y)
    
    time=1
    moving=True
    count=0
    while time <= 40:
        [item['loc'].append(list(make_random_walk([item['loc'][time-1][0],item['loc'][time-1][1]],1))) for key,item in eddies.items()]
        time=time+1
    #Next lines are useful to remove boundary problems.
    newx=linspace(-1,2,300)
    newy=linspace(-1,2,300)
    data=zeros((40,len(newx),len(newy)))
    X,Y=meshgrid(newx,newy)
    for t in range(40):
        for keys, item in eddies.items():
            gauss=twoD_Gaussian((X,Y,item['amp'],newx[item['loc'][t][0]+100],newy[item['loc'][t][1]+100]),\
                        item['radius'], item['radius'], 0, slopex=0, slopey=0, offset=0).reshape(shape(X))
            data[t,:,:]=data[t,:,:]+gauss
            
    eddytd=analyseddyzt(data,x,y,0,40,1,0.5,0.1,0.1\
                    ,data_meant='',areamap='',mask=''\
                    ,destdir='',physics='',diagnostics=False\
                    ,plotdata=False,pprint=False)
    
    eddytdn=analyseddyzt(data,x,y,0,40,1,-0.5,-0.1,-0.1\
                    ,data_meant='',areamap='',mask=''\
                    ,destdir='',physics='',diagnostics=False\
                    ,plotdata=False,pprint=False)
    
    
    return R2

@pytest.mark.trackeddy
def test_gaussfit():
    assert gaussfit() <= 0.99