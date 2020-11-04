from __future__ import print_function

import matplotlib
matplotlib.use('Agg')
import numpy as np
import numpy.ma as ma
import pylab as plt
from trackeddy.datastruct import *
from trackeddy.geometryfunc import *
from trackeddy.physics import *
from trackeddy.savedata import *
import seawater as sw
from scipy import ndimage
from astropy import convolution
import sys
import time
import pdb 

def scan_eddym(data,lon,lat,levels,date,areamap,mask='',destdir='',physics='',eddycenter='masscenter',maskopt='contour',preferences=None,mode='gaussian',basemap=False,checkgauss=True,areaparms=None,usefullfit=False,diagnostics=False,plotdata=False,debug=False):
    '''scan_eddym wraps the identification and Gaussian fit functions.

    Function to identify each eddy using closed contours,
    also this function checks if the elipse adjusted have
    a consistent eccentricity, vorticty and other parameters.

    Parameters
    ----------
    args:
        data: array
            Sea Surface Height in cm.
        lon: array|list
            Longitude of data field.
        lat: array|list 
            Latitude of data field.
        levels: array|list
            Discrete levels where the code will find the closed contours.
        date: int | datetime 
            Date in julian days.
        areamap: list
            Section of interest [lon0,lon1,lat0,lat1].
        mask: np.mask 
            Continent mask.
    
    Returns
    -------
    eddys: dict
        OrderedDict of identified eddies.
    check: boolean
        Status of identification.
    total_contours: int
        Total number of contours analysed.

    Example
    -------
    data: 
        Read netCDF4 data with mask or create a mask for your data.
    lon: 
        Import your longitude coordinates, if your grid is regular, you can use a linspace instead.
    lat: 
        Import your latitude coordinates (same as above).
    levels: 
        List of the levels in which ones you want to find the eddies.
    date: 
        Date as Julian Days.
    areamap: 
        array([[0,len(lon)],[0,len(lat)]]) Array with the index 
        of your area of interest.
    
    I used some auxilar functions, each one has their respective author.
    Author: Josue Martinez Moreno, 2017
    '''

    # Defining lists inside dictionary 
    ellipse_path=[]
    contour_path=[]
    mayoraxis_eddy=[]
    minoraxis_eddy=[]
    gaussianfitdict=[]
    gaussfit2d=[]
    # Data shape
    shapedata=np.shape(data)
    # Diagnostics to list, which allows to print multiple diagnostics at the same time. 
    # (Be carefull because it uses a lot of memory)
    if type(diagnostics) != list:
        diagnostics=[diagnostics]
    # Check if data is masked
    elif data is ma.masked:
        print('Invalid data data, must be masked')
        return
    # Make sure the shape of the data is identical to the grid.
    elif shapedata == [len(lat), len(lon)]:
        print('Invalid  data size, should be [length(lat) length(lon]')
        return
    #Check that area to analyze is valid
    elif np.shape(areamap) == shapedata:
        if np.shape(areamap) == [1, 1] | len(areamap) != len(lat):
            print('Invalid areamap, using NaN for eddy surface area')
        return
    #Check the number of levels is valid.
    elif len(levels)!= 2:
        print('Invalid len of levels, please use the function for multiple levels or use len(levels)==2')
        return
    #Saving mask for future post-processing.  
    if mask!='':
        data=np.ma.masked_array(data, mask)
    datanan=data.filled(np.nan)
    # Obtain the contours of a surface (contourf), this aproach is better than the contour.
    if len(np.shape(lon))== 1 and len(np.shape(lat)) == 1:
        Lon,Lat=np.meshgrid(lon,lat)
    else:
        Lon,Lat=lon,lat
    # Extract min value and max value from the coordinates.
    min_x=Lon[0,0]
    min_y=Lat[0,0]
    max_x=Lon[-1,-1]
    max_y=Lat[-1,-1]
    # Plot contours according to the data.   
    if len(shapedata)==3:
        CS=plt.contourf(Lon[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],\
                        Lat[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],\
               datanan[date,areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],\
                        levels=levels)
    else:
        CS=plt.contourf(Lon[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],\
                        Lat[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],\
                datanan[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],\
                        levels=levels)
    if preferences==None:
        preferences={'ellipse':0.85,'eccentricity':0.85,'gaussian':0.8}
    # Close the contour plot.
    plt.close()
    # Extracting detected contours.
    CONTS=CS.allsegs[:][:]
    # Define variables used in the main loop.
    total_contours=0
    eddyn=0
    threshold=7
    threshold2D=20
    numverlevels=np.shape(CONTS)[0]
    fiteccen=1
    areachecker=np.inf
    ellipsarea=np.inf
    center_eddy=[np.nan,np.nan]
    center_extrem=[np.nan,np.nan]
    contarea=np.inf
    checkm=False
    checkM=False
    gaussarea=[False,0]
    gausssianfitp=[0,0,0,0,0,0]
    xx=np.nan
    yy=np.nan
    areastatus={'check':None,'contour':None,'ellipse':None}
    # Loop in contours of the levels defined.
    for ii in range(0,numverlevels):
        if debug==True:
            print("\n *******EDDY******")
            pdb.set_trace()
        CONTSlvls=CONTS[ii]
        numbereddies=np.shape(CONTSlvls)[0]
        # Loop over all the close contours.
        for jj in range(0,numbereddies):
            if debug==True:
                print("\n ******* EDDY N ******")
                pdb.set_trace()
            check=False
            CONTeach=CONTSlvls[jj]
            
            #Relevant contour values 
            xidmin,xidmax=find2l(lon,lon,CONTeach[:,0].min(),CONTeach[:,0].max())
            yidmin,yidmax=find2l(lat,lat,CONTeach[:,1].min(),CONTeach[:,1].max())
            
            if xidmin<=threshold-1:
                xidmin=+threshold-1
            elif xidmax>=len(lon)-threshold:
                xidmax=len(lon)-threshold
            if yidmin<=threshold-1:
                yidmin=threshold-1
            elif yidmax>=len(lat)-threshold:
                yidmax=len(lat)-threshold
            lon_contour=lon[xidmin-threshold+1:xidmax+threshold]
            lat_contour=lat[yidmin-threshold+1:yidmax+threshold]
            
            cmindex=find(CONTeach[:,1],CONTeach[:,1].max())
            xmindex,ymindex=find2l(lon,lat,CONTeach[cmindex,0],CONTeach[cmindex,1])
            centertop=[ymindex-yidmin+threshold-2,xmindex-xidmin+threshold-1]
            
            if len(shapedata)==3:
                data4gauss=datanan[date,yidmin-threshold+1:yidmax+threshold,xidmin-threshold+1:xidmax+threshold].copy()
                data_in_contour=insideness_contour(data4gauss,centertop,levels,maskopt=maskopt,diagnostics=diagnostics)
            else:
                data4gauss=datanan[yidmin-threshold+1:yidmax+threshold,xidmin-threshold+1:xidmax+threshold].copy()
                data_in_contour=insideness_contour(data4gauss,centertop,levels,maskopt=maskopt,diagnostics=diagnostics)
            
            checkcontour = check_closecontour(CONTeach,lon_contour,lat_contour,data4gauss)
            
            if checkcontour==False:
                xx=np.nan
                yy=np.nan
                center=[np.nan,np.nan]
            else:
                checke=False
                ellipse,status,r2=fit_ellipse(CONTeach[:,0],CONTeach[:,1],diagnostics=diagnostics)
                if status==True and r2 >=preferences['ellipse'] and preferences['ellipse'] < 1:
                    checke=True
                if status==True:
                    ellipseadjust,checke=ellipsefit(CONTeach[:,1],ellipse['ellipse'][1],\
                                              ellipsrsquarefit=preferences['ellipse'],\
                                              diagnostics=diagnostics)
                if checke==True:
                    
                    center = [ellipse['X0_in'],ellipse['Y0_in']]
                    phi = ellipse['phi']
                    axes = [ellipse['a'],ellipse['b']]
                    R = np.arange(0,2.1*np.pi, 0.1)
                    a,b = axes
                    eccen=eccentricity(a,b)
                    
                    if eccen<=preferences['eccentricity']:
                        #Ellipse coordinates.
                        xx = ellipse['ellipse'][0]
                        yy = ellipse['ellipse'][1]

                        mayoraxis = ellipse['majoraxis']
                        minoraxis = ellipse['minoraxis']

                        areastatus=checkscalearea(areaparms,lon_contour,lat_contour,\
                                                    xx,yy,CONTeach[:,0],CONTeach[:,1])
                        
                        if eddycenter == 'maximum':
                            center_eddy=contourmaxvalue(data_in_contour,lon_contour,\
                                                 lat_contour,levels,date,threshold)
                            center_eddy[3]=center_eddy[3]+xidmin-threshold+1
                            center_eddy[4]=center_eddy[4]+yidmin-threshold+1

                            center_extrem=center_eddy

                        elif eddycenter == 'masscenter':
                            center_eddy=centroidvalue(CONTeach[:,0],CONTeach[:,1],\
                                                data_in_contour,lon_contour,\
                                                lat_contour,levels,date,threshold)
                            center_extrem=contourmaxvalue(data_in_contour,lon_contour,\
                                                        lat_contour,levels,date)
                            center_extrem[3]=center_extrem[3]+xidmin-threshold+1
                            center_extrem[4]=center_extrem[4]+yidmin-threshold+1
                        checkM=False
                        checkm=False 
                        
                        if areastatus['status']:
                            if checkgauss==True:
                                if len(shapedata)==3:
                                    profile,checkM=extractprofeddy(mayoraxis,\
                                                   data4gauss,lon_contour,lat_contour,50,\
                                                   gaus='One',kind='linear',\
                                                   gaussrsquarefit=preferences['gaussian'],\
                                                   diagnostics=diagnostics)
                                    if checkM==True:
                                        profile,checkm=extractprofeddy(minoraxis,\
                                                       data4gauss,lon_contour,\
                                                       lat_contour,50,\
                                                       gaus='One',kind='linear',\
                                                       gaussrsquarefit=preferences['gaussian'],\
                                                       diagnostics=diagnostics)
                                else:
                                    profile,checkM=extractprofeddy(mayoraxis,\
                                                   data4gauss,lon_contour,lat_contour,50,\
                                                   gaus='One',kind='linear',\
                                                   gaussrsquarefit=preferences['gaussian'],\
                                                   diagnostics=diagnostics)
                                    if checkM==True:
                                        profile,checkm=extractprofeddy(minoraxis,\
                                                       data4gauss,lon_contour,\
                                                       lat_contour,50,\
                                                       gaus='One',kind='linear',\
                                                       gaussrsquarefit=preferences['gaussian'],\
                                                       diagnostics=diagnostics)
                                #print(checkM,checkm)
                                if checkM==True and checkm==True: 
                                    if levels[0] > 0:
                                        level=levels[0]
                                        extremvalue=np.nanmax(data_in_contour)
                                    else:
                                        level=levels[1]
                                        extremvalue=np.nanmin(data_in_contour)

                                    initial_guess=[a,b,phi,0,0,0]
                                    #print("\n ---pre fit---")
                                    #pdb.set_trace()
                                    fixvalues=[lon_contour,lat_contour,extremvalue,\
                                               center_extrem[0],center_extrem[1]]
                                    gausssianfitp,R2=fit2Dcurve(data_in_contour,\
                                                  fixvalues,\
                                                  level,initial_guess=initial_guess,date='',\
                                                  mode=mode,diagnostics=diagnostics)
                                    #gausssianfitp=initial_guess
                                    # Buf fix for anomalous big Gaussians
                                    #print('++++++++++',abs(gausssianfitp[0]) < 2*np.pi*(xx.max()-xx.min()))
                                    if abs(gausssianfitp[0]) < 2*np.pi*(xx.max()-xx.min()) or abs(gausssianfitp[1]) < 2*np.pi*(xx.max()-xx.min()):
                                        fiteccen=eccentricity(gausssianfitp[0],gausssianfitp[1])
                                        gausscheck2D = checkgaussaxis2D(a,b,\
                                                                        gausssianfitp[0],\
                                                                        gausssianfitp[1])
                                        
                                        #print('=======',gausscheck2D,fiteccen)
                                        if fiteccen <= preferences['eccentricity'] and gausscheck2D==True:
                                            if xidmin <= threshold2D:
                                                xidmin= threshold2D
                                            elif xidmax>=len(lon)-threshold2D:
                                                xidmax=len(lon)-threshold2D
                                            if yidmin <= threshold2D:
                                                xidmin= threshold2D
                                            elif yidmax>=len(lat)-threshold2D:
                                                yidmax=len(lat)-threshold2D
                                            fixvalues[0]=lon[xidmin-threshold2D+1:xidmax+threshold2D]
                                            fixvalues[1]=lat[yidmin-threshold2D+1:yidmax+threshold2D]
                                            gaussarea= gaussareacheck(fixvalues,level,\
                                                                      areaparms,gausssianfitp,\
                                                                      areastatus['contour'])
                                            if  gaussarea[0]==True: 
                                                check=True
                            else:
                                print('Checkgauss need to be True to reconstruct the field.')           
                    if check==True:
                        if usefullfit==False and mode=='gaussian':
                            gausssianfitp[-1]=0
                            gausssianfitp[-2]=0
                            gausssianfitp[-3]=0
                        ellipse_path.append([xx,yy])
                        contour_path.append([CONTeach[:,0],CONTeach[:,1]])
                        mayoraxis_eddy.append([mayoraxis[0],mayoraxis[1]])
                        minoraxis_eddy.append([minoraxis[0],minoraxis[1]])
                        gaussianfitdict.append([gausssianfitp])
                        gaussfit2d.append([R2])
                        #Switch from the ellipse center to the position of 
                        #the maximum value inside de contour
                        if eddyn==0:
                            position_selected=[center_eddy]
                            position_max=[center_extrem]
                            position_ellipse=[center]
                            total_eddy=[[eddyn]]
                            area=[[areastatus['contour'],areastatus['check'],gaussarea[1]]]
                            angle=[phi]
                            if CS.levels[0] > 0:
                                level=CS.levels[0]
                            else:
                                level=CS.levels[1]
                            levelm=[[level]]
                        else:
                            position_selected=np.vstack((position_selected,center_eddy))
                            position_max=np.vstack((position_max,center_extrem))
                            position_ellipse=np.vstack((position_ellipse,center))
                            total_eddy=np.vstack((total_eddy,eddyn))
                            area=np.vstack((area,[areastatus['contour'],areastatus['check'],gaussarea[1]]))
                            angle=np.vstack((angle,phi))
                            
                            if CS.levels[0] > 0:
                                levelprnt=CS.levels[0]
                                levelm=np.vstack((levelm,levelprnt))
                            else:
                                levelprnt=CS.levels[1]
                                levelm=np.vstack((levelm,levelprnt))
                        eddyn=eddyn+1
                    #diagnostics=True
                    if ("ellipse" in diagnostics) or ("all" in diagnostics) or (True in diagnostics):# and  check == True:

                        print("Eddy Number (No time tracking):", eddyn)
                        print("Ellipse parameters")
                        print("Ellipse center = ",  center)
                        print("Mass center = ",  center_eddy)
                        print("angle of rotation = ",  phi)
                        print("axes (a,b) = ", axes)
                        print("Eccentricity (ellips,gauss) = ",eccen,fiteccen)
                        print("Area (rossby,cont,ellips,gauss) = ",
                              areastatus['check'],areastatus['contour'],
                              areastatus['ellipse'],gaussarea[1])
                        print("Ellipse adjust = ", ellipseadjust, checke)
                        print("Mayor Gauss fit = ", checkM)
                        print("Minor Gauss fit = ", checkm)
                        print("2D Gauss fit (Fitness, R^2)=",gausssianfitp,gaussfit2d)
                        print("Conditions | Area | Ellipse | Eccen | Gaussians ")
                        if areastatus['ellipse'] == None or areastatus['check'] == None or areastatus['contour'] == None:
                            print("           | ", False," | ", checke ,"| ",\
                                  eccen <= preferences['eccentricity'] and fiteccen <= preferences['eccentricity'] ,\
                                  " | ", checkM and checkm) 
                        else:
                            print("           | ", areastatus['ellipse'] < areastatus['check'] and areastatus['contour'] < areastatus['check'] and  gaussarea[0],\
                              " | ", checke ,"| ",
                                  eccen <= preferences['eccentricity'] and fiteccen <= preferences['eccentricity'] ,\
                              " | ", checkM and checkm)
                        
                    if ("all" in diagnostics) or (True in diagnostics): #and plotdata == True:
                        f, (ax1, ax2) = plt.subplots(1, 2,figsize=(13, 6))
                        if len(shapedata)==3:
                            ax1.contourf(lon[areamap[0,0]:areamap[0,1]],lat[areamap[1,0]:areamap[1,1]],\
                                data[date,areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]])
                            cc=ax2.pcolormesh(lon[areamap[0,0]:areamap[0,1]],\
                                              lat[areamap[1,0]:areamap[1,1]],\
                              data[date,areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],\
                              vmin=data[date,areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]].min(),\
                              vmax=data[date,areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]].max())
                            cca=ax2.contour(lon[areamap[0,0]:areamap[0,1]],\
                                            lat[areamap[1,0]:areamap[1,1]],\
                                            data[date,areamap[1,0]:areamap[1,1],\
                                            areamap[0,0]:areamap[0,1]],levels=levels,cmap='jet')
                            ax2.clabel(cca, fontsize=9, inline=1)
                        else:
                            cca=ax1.contourf(lon[areamap[0,0]:areamap[0,1]],\
                                             lat[areamap[1,0]:areamap[1,1]],\
                                             datanan[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],\
                                             levels=levels)
                            ax1.plot(CONTeach[:,0],CONTeach[:,1],'-r')
                            ax2.plot(CONTeach[:,0],CONTeach[:,1],'-r')
                            ax2.pcolormesh(lon[areamap[0,0]:areamap[0,1]],\
                                           lat[areamap[1,0]:areamap[1,1]],\
                                           datanan[areamap[1,0]:areamap[1,1],\
                                           areamap[0,0]:areamap[0,1]],vmin=-20,vmax=20)
                            plt.show()
                            f, (ax1, ax2) = plt.subplots(1, 2,figsize=(13, 6))
                            ax1.contourf(lon[areamap[0,0]:areamap[0,1]],\
                                         lat[areamap[1,0]:areamap[1,1]],\
                                         data[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]])
                            cc=ax2.pcolormesh(lon[areamap[0,0]:areamap[0,1]],\
                                              lat[areamap[1,0]:areamap[1,1]],\
                              data[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],\
                              vmin=data[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]].min(),\
                              vmax=data[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]].max())
                            cca=ax2.contour(lon[areamap[0,0]:areamap[0,1]],\
                                            lat[areamap[1,0]:areamap[1,1]],\
                                            data[areamap[1,0]:areamap[1,1],areamap[0,0]:areamap[0,1]],\
                                            levels=levels,cmap='jet')
                            ax2.clabel(cca, fontsize=9, inline=1)
                        ax1.plot(CONTeach[:,0],CONTeach[:,1],'*r')
                        ax1.plot(xx,yy,'-b')
                        ax1.plot(center[0],center[1],'ob')
                        f.subplots_adjust(right=0.8)
                        cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
                        f.colorbar(cc, cax=cbar_ax)
                        ax2.plot(CONTeach[:,0],CONTeach[:,1],'-r')
                        ax2.plot(xx,yy,'-b')
                        ax2.plot(center[0],center[1],'ob')
                        idxelipcheck,idyelipcheck=find2l(lon,lat,center[0],center[1])
                        ax2.plot(lon[idxelipcheck],lat[idyelipcheck],'om')
                        ax2.plot(center_eddy[0],center_eddy[1],'oc')
                        ax2.plot(center_extrem[0],center_extrem[1],'*g')
                        ax2.set_ylim([CONTeach[:,1].min(),CONTeach[:,1].max()])
                        ax2.set_xlim([CONTeach[:,0].min(),CONTeach[:,0].max()])
                        plt.show()
                        plt.close()
                        
                total_contours=total_contours+1
        try:
            position_selected=np.array(position_selected)
            position_max=np.array(position_max)
            position_ellipse=np.array(position_ellipse)
            area=np.array(area)
            levelm=np.array(levelm)
            mayoraxis_eddy=np.array(mayoraxis_eddy)
            minoraxis_eddy=np.array(minoraxis_eddy)
            gaussianfitdict=np.array(gaussianfitdict)
            gaussfit2d=np.array(gaussfit2d)
            eddys=dict_eddym(contour_path,ellipse_path,position_selected,\
                             position_max,position_ellipse,\
                             mayoraxis_eddy,minoraxis_eddy,\
                             area,angle,total_eddy,levelm,gaussianfitdict,gaussfit2d)
            check=True
            
        except:
            check=False
            eddys=0
        #if destdir!='':
        #    save_data(destdir+'day'+str(date)+'_one_step_cont'+str(total_contours)+'.dat', variable)
    return eddys,check,total_contours

def analyseddyzt(data,x,y,t0,t1,tstep,levels,areamap='',mask='',physics='',eddycenter='masscenter',preferences=None,checkgauss=True,areaparms=None,maskopt='contour',mode='gaussian',filters=None,timeanalysis='closest',destdir='',saveformat='nc',diagnostics=False,plotdata=False,pprint=False,debug=False):
    '''Identify each eddy using closed contours.

    Function to identify each eddy using closed contours, 
    moving in time and contour levels.
    
    Parameters
    ----------
    args:
        data: array
            Sea Surface Height in cm.
        x: array|list
            Longitude of data field.
        y: array|list 
            Latitude of data field.
        levels: array|list
            Discrete levels where the code will find the closed contours.
        areamap: list
            Section of interest [lon0,lon1,lat0,lat1].
        mask: np.mask 
            Continent mask.
    
    Returns
    -------
        

    Example
    -------
        data: 
            Read netCDF4 data with mask or create a mask for your data.
    I used some auxilar functions, each one has his respective author.
    Author: Josue Martinez Moreno, 2017
    '''

    if pprint == True:
        pp = Printer(); 
    if len(np.shape(data)) < 3:
        raise Exception('If you whant to analyze in time the data need to be 3d [i.e. data(t,x,y)]')
        #return
    if areamap=='':
        areamap = np.array([[0,len(x)],[0,len(y)]])
    if mask == '':
        if ma.is_masked(data):
            if len(np.shape(data))<3:
                mask = ma.getmask(data[:,:])
            else:
                mask = ma.getmask(data[0,:,:])
        else:
            if len(np.shape(data))<3:
                mask = np.zeros(np.shape(data[:,:]))
            else:
                mask = np.zeros(np.shape(data[0,:,:]))
    if preferences==None:
        preferences={'ellipse':0.85,'eccentricity':0.85,'gaussian':0.8}
    #Check area
    if areaparms==None:
        areaparms={'mesoscale':2*np.pi}
    #Define list of levels based on the user defined dictionary
    if type(levels) == dict:
        keycheck = ('max' in levels.keys() and 'min' in levels.keys() and 'step' in levels.keys())
        difftozero = (levels['min'] != 0 or levels['max'] != 0 or levels['step'] != 0)
        if keycheck and difftozero:
            del keycheck,difftozero
            levellist  = np.arange(levels['min'],levels['max']+levels['step'],levels['step'])
            farlevel   = levellist[0]
            if abs(levellist)[0] < abs(levellist)[-1]:
                levellist = np.flipud(levellist)
                farlevel  = levellist[0]
        elif levels==dict and ('max' in levels.keys() and 'min' in levels.keys() and 'step' in levels.keys()):
            raise ValueError("Not all the level parameters are defined, make sure the dictionary contains the keys: 'max','min','step'")
        else:
            raise ValueError("The level parameters set up is incorrect, it can be 0.")
    #Define list of levels based on the user defined list.
    elif type(levels) == list or type(levels) == np.ndarray:
        levellist  = levels
        farlevel   = levellist[0]
        if abs(levellist)[0] < abs(levellist)[-1]:
            levellist = np.flipud(levellist)
            farlevel  = levellist[0]
    elif type(levels) ==int or type(levels) == float: 
        farlevel=levels
        levellist=[levels]
    else: 
        raise ValueError("Levels should be a dictionary or a list. Read the documentation to know more about it.")
                    
    if pprint==True:
        pp.timepercentprint(0,1,1,0,"Init time")
    if pprint==True:
        pp = Printer(); 
    numbereddieslevels=0
    
    for ii in range(t0,t1,tstep):
        checkcount = 0 
        # Defining filters
        if filters == None or ('time' not in filters.keys() and 'spatial' not in filters.keys()):
            filters = {'time':{'type':None,'t':None,'t0':None,'value':None},'spatial':{'type':None,'window':None,'mode':None}}
            dataanomaly=data[ii,:,:]
        # Check if the expected filters are defined
        elif 'time' not in filters.keys():
            filters['time'] = {'type':None,'t':None,'t0':None,'value':None}
            dataanomaly=data[ii,:,:]
        elif 'spatial' not in filters.keys():
            filters['spatial'] = {'type':None,'window':None,'mode':None}
            dataanomaly=data[ii,:,:]
        elif len(filters.keys()) != 2:
            raise ValueError("Unexpected dictionary, documentation at: \n https://trackeddy.readthedocs.io/en/latest/pages/Methods.html")
        # Apply temporal filter
        if filters['time']['type'] == None and filters['time']['value'] == None and (filters['time']['t0'] == None or filters['time']['t'] == None):
            dataanomaly=data[ii,:,:]
            #print('No time filter apply')
        # Check if the user selects to remove a predefined or calculated historical filter.
        elif filters['time']['type'] == 'historical' and type(filters['time']['value']) != type(None):
            dataanomaly  = ma.masked_array(data[ii,:,:]-filters['time']['value'], mask)
        elif filters['time']['type'] == 'historical' and filters['time']['value'] == None:
            dataanomaly = ma.masked_array(data[ii,:,:]-np.nanmean(data,axis=0), mask)
        # Check if the user selects an time orthogonal filter.
        elif filters['time']['type'] != 'historical' and (filters['time']['t'] != None or  filters['time']['t0'] != None):
            dataanomaly = ma.masked_array(data[ii,:,:]-np.nanmean(data[ii-t0:ii+t],axis=0), mask)
        elif filters['time']['type'] != 'historical' and  filters['time']['type'] != None and (filters['time']['t'] == None or  filters['time']['t0'] == None):
            raise ValueError("T and T0 are undefined, include the definition in the dictionary.")
        else:
            raise ValueError("Define the filter argument like: /n filters={'time':{'type':'orthogonal','t':1,'t0':shape(data)[0]},'spatial':{'type':'moving','window':70,'mode':'uniform'}}")
        # Apply spatial filter.
        if filters['spatial']['type'] == None and filters['spatial']['window'] == None and filters['spatial']['mode'] == None:
            pass
            #print('No spatial filter apply')
        elif filters['spatial']['type'] == 'moving' and filters['spatial']['window'] != None:
            if filters['spatial']['mode'] == 'uniform':
                nofilterdata = data[ii,:,:]
                if filters['spatial']['window']%2 == 0:
                    ker=np.ones((filters['spatial']['window']+1,filters['spatial']['window']+1))
                else:
                    ker=np.ones((filters['spatial']['window'],filters['spatial']['window']))
                nofilterdata = nofilterdata - convolution.convolve(nofilterdata, kernel = ker)
                dataanomaly = ma.masked_array(nofilterdata, mask)
            if filters['spatial']['mode'] == 'gaussian':
                raise Warning('ndimage.gaussian_filter may create artefacts near nan values. Therefore, data is filled with zeros.')
                data = data.filled(fill_value=0)
                nofilterdata = data[ii,:,:]
                nofilterdata = nofilterdata - ndimage.gaussian_filter(nofilterdata, size = filters['spatial']['window'])
                dataanomaly = ma.masked_array(nofilterdata, mask)       
        # Check if the user selects an moving average meridional filter.
        elif filters['spatial']['type'] == 'meridional':
            dataanomaly = ma.masked_array(data[ii,:,:]-np.nanmean(np.squeeze(data[ii,:,:]),axis=0), mask)
        # Check if the user selects an moving average zonal filter.
        elif filters['spatial']['type'] == 'zonal':
            dataanomaly = ma.masked_array((data[ii,:,:].T-np.nanmean(np.squeeze(data[ii,:,:]),axis=1)).T, mask)
        else:
            raise ValueError("Define the filter argument like: /n filters={'time':{'type':'orthogonal','t':1,'t0':shape(data)[0]},'spatial':{'type':'moving','window':70,'mode':'uniform'}}")
            
        for ll in range(0,len(levellist)):
            if levellist[ll]<0:
                levels_scan=[-np.inf,levellist[ll]]
            else:
                levels_scan=[levellist[ll],np.inf]
            
            eddies,check,numbereddies=scan_eddym(dataanomaly,x,y,levels_scan,ii\
                          ,areamap,mask=mask,destdir=destdir\
                          ,physics=physics,eddycenter=eddycenter,maskopt=maskopt\
                          ,checkgauss=checkgauss,areaparms=areaparms\
                          ,preferences=preferences,mode=mode\
                          ,diagnostics=diagnostics,plotdata=plotdata,debug=debug)
            if check==True and checkcount==0:
                eddzcheck=True
                checkcount=1
            else:
                eddzcheck=False
            if eddies!=0 and check==True:
                if levellist[ll] == farlevel or eddzcheck==True:
                    eddz = dict_eddyz(dataanomaly,x,y,ii,ll,levellist,farlevel,eddies,diagnostics=diagnostics,debug=debug)
                else:
                    eddz = dict_eddyz(dataanomaly,x,y,ii,ll,levellist,farlevel,eddies,eddz,diagnostics=diagnostics,debug=debug)
            if pprint==True:
                numbereddieslevels=numbereddieslevels+numbereddies
                pp.timepercentprint(t0,t1,tstep,ii,'# of E '+ str(numbereddies),[0,len(levellist),ll])
        if ii==t0:
            eddytd=dict_eddyt(ii,eddz,analysis=timeanalysis,debug=debug)
        else:
            eddytd=dict_eddyt(ii,eddz,eddytd,data=dataanomaly,analysis=timeanalysis,x=x,y=y,debug=debug) 
        if pprint==True:
            pp.timepercentprint(t0,t1,tstep,ii,'# of E '+ str(numbereddieslevels))
    if destdir!='':
        if saveformat=='nc':
            eddync(destdir+str(level)+'.nc',eddytd)
        else:
            np.save(destdir+str(level)+'.npy',eddytd)
    return eddytd