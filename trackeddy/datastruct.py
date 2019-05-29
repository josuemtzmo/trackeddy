import numpy as np
import netCDF4 as nc4
from datetime import datetime
import pylab as plt
import seawater as sw
import warnings
from trackeddy.savedata import *
from trackeddy.geometryfunc import *
from trackeddy.physics import *
warnings.filterwarnings("ignore")
import pdb

def dict_eddym(contour, ellipse, position_selected,position_max,position_ellipse,majoraxis_eddy,minoraxis_eddy,area,angle,number,level,gaussianfitdict,gaussfit2d):
    '''
    ********************** dict_eddym **********************
    Create a dictionary with a basic structure of the eddies. Can contain the structure of one eddy or multple eddies.
    Notes: 
        
    Args:
        contour (array): Close contour used to track any gaussian (Eddy).
        ellipse (array): Ellipse ajusted to identify a gaussian shape. 
        position_selected (array): Positions (X, Y) describing the location of the feature. 
        position_ellipse (array): Positions (X, Y) describing the location of the ajusted ellipse. 
        majoraxis_eddy (array): Positions (X, Y) describing the mayor axis of the ajusted ellipse.
        minoraxis_eddy (array): Positions (X, Y) describing the minor axis of the ajusted ellipse.
        area (array): Area of the feature.
        angle (array): Angle of the ajusted ellipse starting on the 0 degree convention (counterclockwise).
        number (list): Number of the feature. Flag used to track features in time.
        level (array): Contour level used during the identification of the feature.
        GaussianFitdict(dict): Syntetic field of a 2D gaussian adjust.
    Returns:
        contour_data - Dictionary containing all the relevant information of the features.
        The new dictionary containing all the information has the following form:        
        {'Contour':contour,'Ellipse':ellipse,'Position':position_selected,\
        'PositionEllipse':position_ellipse,'MajorAxis':majoraxis_eddy,\
        'MinorAxis':minoraxis_eddy,'Area':area,'Angle':angle,'EddyN':number,'Level':level,'2DGaussianFit':gausssianfitp}
    Usage:
        eddys=dict_eddym(contour_path,ellipse_path,position_selected,position_ellipse,mayoraxis_eddy,minoraxis_eddy,\
                     area,angle,total_eddy,level)
    '''
    contour_data={'Contour':contour,'Ellipse':ellipse,'Position':position_selected,'PositionExtreme':position_max,'PositionEllipse':position_ellipse,'MajorAxis':majoraxis_eddy,'MinorAxis':minoraxis_eddy,'Area':area,'Angle':angle,'EddyN':number,'Level':level,'2DGaussianFit':gaussianfitdict,'2DGaussianR2':gaussfit2d}
    return contour_data

def datastruct_time(ts,eddys,eddydt):
    dictime={'eddyn_'+str(eddys['EddyN'][0][0]):{'neddy':eddys['EddyN'][0],'time':np.array([ts]),\
            'position_default':[eddys['Position'][0]],\
            'area':eddys['Area'][0],'ellipse':[eddys['Ellipse'][0]],'contour':[eddys['Contour'][0]],\
            'angle':eddys['Angle'][0],'position_maxvalue':[eddys['PositionExtreme'][0]],\
            'position_eddy':eddys['PositionEllipse'][0],'level':eddys['Level'][0],\
            'majoraxis':eddys['MajorAxis'][0],'minoraxis':eddys['MinorAxis'][0],\
            '2dgaussianfit':eddys['2DGaussianFit'][0],'timetracking':True}}

    for nn in range(1,len(eddys['EddyN'])):
        dictime['eddyn_'+str(eddys['EddyN'][nn][0])]={'neddy':eddys['EddyN'][nn],'time':np.array([ts]),\
                'position_default':[eddys['Position'][nn]],'area':eddys['Area'][nn],\
                'ellipse':[eddys['Ellipse'][nn]],'contour':[eddys['Contour'][nn]],\
                'angle':eddys['Angle'][nn],'position_maxvalue':[eddys['PositionExtreme'][nn]],\
                'position_eddy':eddys['PositionEllipse'][nn],'level':eddys['Level'][nn],\
                'majoraxis':eddys['MajorAxis'][nn],'minoraxis':eddys['MinorAxis'][nn],\
                '2dgaussianfit':eddys['2DGaussianFit'][nn],'timetracking':True}
    return dictime

def dict_eddyt(ts,eddys,eddydt='',data="",x="",y="",analysis="overlap",maxvalue='maxvalue',coords='latlon',diagnostics=False,debug=False):
    '''
    ********************** dict_eddyt **********************
    Create a dictionary with all the eddies and it's track on time. When 'ts==0' 
    or eddydt is not defined it start creating a new dictionary otherwise 
    it grows the dictionary.
    Notes:
        Check for some possible bugs where it's saving more than once some features.
    Args:
        
        ts (int): Step or counter used to move in time.
        eddys (dict): Dictionary containing all the relevant information of 
            the features. Check "dict_eddym".
        eddydt(dict): Output dictionary of this function used to grow the 
            features. Check "dict_eddyt".
        analysis (str): Multiple options to track features on time: 
            insideness: Check the coordenates of the var.max() at time t is inside 
                        close contour at time t+1.
            overlap: Check the overlap between close contours between time t and T+1. 
                     The saved contour is the one which contains the maximum overlap.
            closest: Calculates the distance between maximum values or the center of 
                     mass and correlates the closest features.
        maxvalue (str): Point inside the eddy which will be compared.
            maxvalue: Maximum value inside the contour.
            default: Eddy center of mass or maximum value.
        data (matrix): Field analysed, neccessary depending on the analysis used.
    Returns:
        eddydt - Dictionary containig each eddy and it's track on time.
        All the keys have the following form:
        {eddyn_0:{neddy,time,position,ellipse,contour,...},
        eddyn_1:{neddy,time,position,ellipse,contour,...},...}
    Usage:
        if tt==0:
            eddytd=dict_eddyt(tt,eddys)
        else:
            eddytd=dict_eddyt(tt,eddys,eddytd) 
        
    '''
    if debug==True:
        print("\n *******TIME******")
        pdb.set_trace()
    #neweddies=[]
    threshold=7
    if eddydt=='':
        eddydt=datastruct_time(ts,eddys,eddydt)
        return eddydt
    
    eddyt1={str(ii): '' for ii in range(0,len(eddys['EddyN']))}
    eddyt0={str(eddydt[key]['neddy'][0]): '' for key in eddydt}# if eddydt[key]['timetracking']==True}
    
    if analysis=='insideness':
        pass
            
    elif analysis=='overlap':
        for t0key in eddyt0.keys():
            t0contour=eddydt['eddyn_'+str(t0key)]['contour']
            levels=eddydt['eddyn_'+str(t0key)]['level']
            
            posxs=int(eddydt['eddyn_'+str(t0key)]['position_maxvalue'][-1][0])
            posys=int(eddydt['eddyn_'+str(t0key)]['position_maxvalue'][-1][1])
            
            
            idxmxs=int(eddydt['eddyn_'+str(t0key)]['position_maxvalue'][-1][3])
            idxmys=int(eddydt['eddyn_'+str(t0key)]['position_maxvalue'][-1][4])
            
            #print(t0contour[-1][0])
            xidmin,xidmax=find2l(x,x,t0contour[-1][0].min(),t0contour[-1][0].max())
            yidmin,yidmax=find2l(y,y,t0contour[-1][1].min(),t0contour[-1][1].max())
            
            #print(xidmin,xidmax,yidmin,yidmax)
            if yidmin<threshold:
                xcontour=x[xidmin-threshold+1:xidmax+threshold]
                ycontour=y[yidmin:yidmax+threshold]
                datacontours=data[yidmin:yidmax+threshold,xidmin-threshold+1:xidmax+threshold]*1
            elif xidmin<threshold:
                xcontour=x[xidmin:xidmax+threshold]
                ycontour=y[yidmin-threshold+1:yidmax+threshold]
                datacontours=data[yidmin-threshold+1:yidmax+threshold,xidmin:xidmax+threshold]*1
            else:
                xcontour=x[xidmin-threshold+1:xidmax+threshold]
                ycontour=y[yidmin-threshold+1:yidmax+threshold]
                datacontours=data[yidmin-threshold+1:yidmax+threshold,xidmin-threshold+1:xidmax+threshold]*1
            
            plt.pcolormesh(xcontour,ycontour,datacontours)
            plt.plot(posxs,posys,'or')
            plt.plot(t0contour[-1][0],t0contour[-1][1])
            plt.show()
            
            if levels>0:
                datacontours[datacontours<levels]=0
                datacontours[datacontours>=levels]=1
            else:
                datacontours[datacontours>levels]=0
                datacontours[datacontours<=levels]=1
            markers=ndimage.label(datacontours)[0]
            
            for t1key in range(0,len(eddys['EddyN'])):
                t1position=eddys['PositionExtreme'][t1key]
                
            
            print(markers.max())
            plt.pcolormesh(xcontour,ycontour,markers)
            plt.plot(posxs,posys,'or')
            plt.plot(t0contour[-1][0],t0contour[-1][1])
            plt.show()
            
    elif analysis=='none':
        neweddies=[int(key) for key in eddyt1 if eddyt1[key]=='']    
        eddydt=addtimetrack(ts,eddydt,eddys,neweddies,debug=debug)   
        
    else:
        for t0key in eddyt0.keys():
            eddydist=[]
            t0position=eddydt['eddyn_'+str(t0key)]['position_'+maxvalue][-1]
            for t1key in range(0,len(eddys['EddyN'])):
                t1position=eddys['PositionExtreme'][t1key]
                if coords=='latlon':
                    eddydist.append(abs(sw.dist((t1position[1],t0position[1]),\
                             (t1position[0],t0position[0]),'km')[0][0]*1000))

            eddydist=np.asarray(eddydist)
            t1track=str(find(eddydist,eddydist.min()))
            t0track=str(eddydt['eddyn_'+str(t0key)]['neddy'][0])
            L_R=rossbyR(np.mean((eddys['PositionExtreme'][int(t1track)][1],\
                        eddydt['eddyn_'+t0track]['position_'+maxvalue][-1][1])),\
                        np.mean((eddys['PositionExtreme'][int(t1track)][1],\
                        eddydt['eddyn_'+t0track]['position_'+maxvalue][-1][1])))
            majoraxdist= sw.dist(eddys['MajorAxis'][t1key][1],eddys['MajorAxis'][t1key][0],'km')[0][0]*1000 
            sw.dist((t1position[1],t0position[1]),\
                             (t1position[0],t0position[0]),'km')[0][0]*1000
            if eddydist.min()<L_R and eddydist.min()<majoraxdist and (eddyt0[t0track]!=True and eddyt1[t1track]!=True):
                eddyt0[t0track]=True
                eddyt1[t1track]=True
                eddydt=jointimetrack(ts,eddydt,eddys,int(t0track),int(t1track))
        
        ## Add new eddies to the dictionary 
        neweddies=[int(key) for key in eddyt1 if eddyt1[key]=='']    
        eddydt=addtimetrack(ts,eddydt,eddys,neweddies,debug=debug)   
        ## Count time steps without appearance 
        oldeddies=[int(key) for key in eddyt0 if eddyt0[key]=='']
        for oldeddy in oldeddies:
            if type(eddydt['eddyn_'+str(oldeddy)]['timetracking'])==5:
                eddydt['eddyn_'+str(oldeddy)]['timetracking']=False
            elif eddydt['eddyn_'+str(oldeddy)]['timetracking']==True:
                eddydt['eddyn_'+str(oldeddy)]['timetracking']=1
            elif type(eddydt['eddyn_'+str(oldeddy)]['timetracking'])==int:
                eddydt['eddyn_'+str(oldeddy)]['timetracking']=eddydt['eddyn_'+str(oldeddy)]['timetracking']+1
    return eddydt

def addtimetrack(ts,eddydt,eddys,neweddies,debug=False):
    if debug==True:
        print("\n *******New eddy TIME******")
        pdb.set_trace()
    for neweddy in neweddies:
        number=len(eddydt.keys())
        #print(neweddy)
        eddydt['eddyn_'+str(number)]={'neddy':[number],'time':np.array([ts]),\
                            'position_default':[eddys['Position'][neweddy]],\
                            'position_maxvalue':[eddys['PositionExtreme'][neweddy]],\
                            'area':eddys['Area'][neweddy],\
                            'angle':eddys['Angle'][neweddy],\
                            'position_eddy':eddys['PositionEllipse'][neweddy],\
                            'ellipse':[eddys['Ellipse'][neweddy]],\
                            'contour':[eddys['Contour'][neweddy]],\
                            'level':[eddys['Level'][neweddy]],\
                            'minoraxis':eddys['MinorAxis'][neweddy],\
                            'majoraxis':eddys['MajorAxis'][neweddy],\
                            '2dgaussianfit':eddys['2DGaussianFit'][neweddy],\
                            'timetracking':True}
    return eddydt
  

def jointimetrack(ts,eddydt,eddys,t0track,t1track):
    eddydt['eddyn_'+str(t0track)]['time']=np.vstack((eddydt['eddyn_'+str(t0track)]['time'],ts))
    eddydt['eddyn_'+str(t0track)]['position_default']=np.vstack((eddydt['eddyn_'+str(t0track)]['position_default'],eddys['Position'][t1track]))
    eddydt['eddyn_'+str(t0track)]['position_maxvalue']=np.vstack((eddydt['eddyn_'+str(t0track)]['position_maxvalue'],eddys['PositionExtreme'][t1track]))
    eddydt['eddyn_'+str(t0track)]['area']=np.vstack((eddydt['eddyn_'+str(t0track)]['area'],eddys['Area'][t1track]))
    eddydt['eddyn_'+str(t0track)]['angle']=np.vstack((eddydt['eddyn_'+str(t0track)]['angle'],eddys['Angle'][t1track]))
    
    eddydt['eddyn_'+str(t0track)]['position_eddy']=np.vstack((eddydt['eddyn_'+str(t0track)]['position_eddy'],eddys['PositionEllipse'][t1track]))
    eddydt['eddyn_'+str(t0track)]['level']=np.vstack((eddydt['eddyn_'+str(t0track)]['level'],eddys['Level'][t1track]))
    eddydt['eddyn_'+str(t0track)]['minoraxis']=np.vstack((eddydt['eddyn_'+str(t0track)]['minoraxis'],eddys['MinorAxis'][t1track]))
    eddydt['eddyn_'+str(t0track)]['majoraxis']=np.vstack((eddydt['eddyn_'+str(t0track)]['majoraxis'],eddys['MajorAxis'][t1track]))
    eddydt['eddyn_'+str(t0track)]['2dgaussianfit']=np.vstack((eddydt['eddyn_'+str(t0track)]['2dgaussianfit'],eddys['2DGaussianFit'][t1track]))
    
    eddydt['eddyn_'+str(t0track)]['ellipse']=eddydt['eddyn_'+str(t0track)]['ellipse']+[eddys['Ellipse'][t1track]]
    eddydt['eddyn_'+str(t0track)]['contour']=eddydt['eddyn_'+str(t0track)]['contour']+[eddys['Contour'][t1track]]
    
    eddydt['eddyn_'+str(t0track)]['timetracking']=True
    return eddydt


def dict_eddyz(data,x,y,ts,levelindex,levellist,maxlevel,eddys,eddz='',threshold=1.5,diagnostics=False,debug=False):
    '''
    ********************** dict_eddyz **********************
    Create a dictionary with all the eddies and it's develop in delta eta, 
    where the biggest contour is assigned as the contour that will be tracked.
    Notes:
        Check for some possible bugs where it's saving more than once some features.
    Args:
        ts (int): Step or counter used to move in time.
        ll (int): Level analysed used to move in delta eta, where delta eta 
        is the changes in elevation of the surface.
        maxlevel (int): Max level defined by the user. 
        eddys (dict): Dictionary containing all the relevant information 
        of the features. Check "dict_eddym".
        eddz(dict): Output dictionary of this function used to grow the 
        features. Check "dict_eddyz".
        threshold (float): Used to grow the detection; R=R0+threshold, 
        where R0 is the radius of the contour.
        diagnostics (boolean): Used to display all the statistics 
        and plots to identify bugs.
    Returns:
        eddz - Dictionary containing the largest contour and parameters 
        of each features.
        All the keys have the following form:
        {eddyn_0:{neddy,time,position,ellipse,contour},
        eddyn_1:{neddy,time,position,ellipse,contour}}
    Usage:
        if ll == 0:
            eddz = dict_eddyz(ii,ll,farlevel,eddies,diagnostics=diagnostics)
        else:
            eddz = dict_eddyz(ii,ll,farlevel,eddies,eddz,diagnostics=diagnostics)
    '''
    athresh=1000
    threshold=7
    if type(diagnostics) != list:
        diagnostics=[diagnostics]
    if eddz=='' or maxlevel==levellist[levelindex]:
        eddz=eddys
    else:         
        #Always check in the next one because probable we will have more contours if we get closer to 0.
        contour=eddz['Contour']
        ellipse=eddz['Ellipse']
        position=eddz['Position']
        position_max=eddz['PositionExtreme']
        position_ellipse=eddz['PositionEllipse']
        area=eddz['Area']
        angle=eddz['Angle']
        majoraxis=eddz['MajorAxis']
        minoraxis=eddz['MinorAxis']
        number=eddz['EddyN']
        level=eddz['Level']
        gauss2d=np.array(eddz['2DGaussianFit'])
        gaussfit2d=np.array(eddz['2DGaussianR2'])
        
        eddysmaskdata={str(ii): '' for ii in range(0,len(eddys['EddyN']))}
        eddyzmaskdata={str(ii): '' for ii in range(0,len(eddz['EddyN']))}

        neweddies=[]
        for nn0 in range(0,len(eddys['EddyN'])):
            if debug==True:
                print("*******LEVELS******")
                pdb.set_trace()
            xscontour=eddys['Contour'][nn0][0]
            yscontour=eddys['Contour'][nn0][1]
            levels=eddys['Level'][nn0]
            areas=eddys['Area'][nn0]
            gauss2dfits=eddys['2DGaussianR2'][nn0]
            coordmxs=eddys['PositionExtreme'][nn0][0]
            coordmys=eddys['PositionExtreme'][nn0][1]
            magms=eddys['PositionExtreme'][nn0][2]
            idxmxs=int(eddys['PositionExtreme'][nn0][3])
            idxmys=int(eddys['PositionExtreme'][nn0][4])

            xidmin,xidmax=find2l(x,x,xscontour.min(),xscontour.max())
            yidmin,yidmax=find2l(y,y,yscontour.min(),yscontour.max())

            replacecontour=False
            n=0
            while n != -1 :
                if n>=len(eddz['EddyN']) or all([val=='Done' for val in eddysmaskdata.values()]):
                    neweddies.append([len(neweddies)+1,nn0])
                    n=-1
                elif eddysmaskdata[str(nn0)]!='Done' and eddyzmaskdata[str(n)] != 'Done':
                    xzcontour=eddz['Contour'][n][0]
                    yzcontour=eddz['Contour'][n][1]
                    levelz=eddz['Level'][n]
                    areaz=eddz['Area'][n]
                    gauss2dfitz=eddz['2DGaussianR2'][n]
                    coordmxz=eddz['PositionExtreme'][n][0]
                    coordmyz=eddz['PositionExtreme'][n][1]
                    magmz=eddz['PositionExtreme'][n][2]
                    if coordmxz==coordmxs and coordmyz==coordmys and magmz==magms:
                        eddysindex=n
                        eddysmaskdata[str(nn0)]='Done'
                        eddyzmaskdata[str(n)]='Done'
                        ### Change and to or, but first solve the bug 
                        ### with the R2 of the 2d fit.
                        ### gauss2dfitz < gauss2dfits and 
                        if (areas[0] - areas[2]) < (areas[0]-areaz[2]):
                            replacecontour=True
                        else:
                            replacecontour=False
                        n=-1
                    else:
                        replacecontour=False
                        n=n+1
                else:    
                    n=n+1
            if replacecontour==True:
                if ("levels" in diagnostics) or ("all" in diagnostics) or (True in diagnostics):
                    if xidmin<7 and yidmin>7:
                        plt.pcolormesh(x[xidmin:xidmax+threshold],y[yidmin-threshold+1:yidmax+threshold],data[yidmin-threshold+1:yidmax+threshold,xidmin:xidmax+threshold])
                    elif yidmin<7 and xidmin>7:
                        plt.pcolormesh(x[xidmin-threshold+1:xidmax+threshold],y[yidmin:yidmax+threshold],data[yidmin:yidmax+threshold,xidmin-threshold+1:xidmax+threshold])
                    else:
                        plt.pcolormesh(x[xidmin-threshold+1:xidmax+threshold],y[yidmin-threshold+1:yidmax+threshold],data[yidmin-threshold+1:yidmax+threshold,xidmin-threshold+1:xidmax+threshold])
                    plt.colorbar()
                    plt.plot(eddys['Contour'][nn0][0],eddys['Contour'][nn0][1],'-r')
                    plt.plot(contour[eddysindex][0],contour[eddysindex][1],'-b')

                    plt.show()
                    print('Area:',areas/areaz)

                contour[eddysindex]=eddys['Contour'][nn0]
                ellipse[eddysindex]=eddys['Ellipse'][nn0]
                position[eddysindex]=eddys['Position'][nn0]
                position_max[eddysindex]=np.array(eddys['PositionExtreme'][nn0])
                position_ellipse[eddysindex]=eddys['PositionEllipse'][nn0]
                area[eddysindex]=eddys['Area'][nn0]
                angle[eddysindex]=eddys['Angle'][nn0]
                majoraxis[eddysindex]=eddys['MajorAxis'][nn0]
                minoraxis[eddysindex]=eddys['MinorAxis'][nn0]
                level[eddysindex]=eddys['Level'][nn0]
                gauss2d[eddysindex]=np.array(eddys['2DGaussianFit'][nn0])
                gaussfit2d[eddysindex]=np.array(eddys['2DGaussianR2'][nn0])
                
        for ii,jj in neweddies:
            contour=list(contour)
            contour.append(eddys['Contour'][jj])
            ellipse=list(ellipse)
            ellipse.append(eddys['Ellipse'][jj])
            majoraxis=list(majoraxis)
            majoraxis.append(eddys['MajorAxis'][jj])
            minoraxis=list(minoraxis)
            minoraxis.append(eddys['MinorAxis'][jj])
            position=np.vstack((position,[eddys['Position'][jj]]))
            position_max=np.vstack((position_max,[eddys['PositionExtreme'][jj]]))
            position_ellipse=np.vstack((position_ellipse,[eddys['PositionEllipse'][jj]]))
            area=np.vstack((area,eddys['Area'][jj]))
            angle=np.vstack((angle,eddys['Angle'][jj]))
            level=np.vstack((level,eddys['Level'][jj]))
            gauss2d=np.vstack((gauss2d,[eddys['2DGaussianFit'][jj]]))
            gaussfit2d=np.vstack((gaussfit2d,[eddys['2DGaussianR2'][jj]]))
            number=np.vstack((number,[len(number)]))
                
        eddz={'Contour':contour,'Ellipse':ellipse,'Position':position,'PositionExtreme':position_max,\
              'PositionEllipse':position_ellipse,'Area':area,'MajorAxis':majoraxis,\
              'MinorAxis':minoraxis,'Angle':angle,'EddyN':number,'Level':level,\
              '2DGaussianFit':gauss2d,'2DGaussianR2':gaussfit2d}
    return eddz

def joindict(dict1,dict2):
    '''
    ********************** joindict **********************
    Join two dictionaries containing the features parameters in time an delta eta. Useful when datasets are in multiple files.
    Notes:
        Check for some possible bugs where it's saving more than once some features.
    Args:
        dict1 (dict): Dictionary containing all the relevant information of the features. Check "dict_eddym".
        dict2 (dict): Dictionary containing all the relevant information of the features. Check "dict_eddym".
    Returns:
        dict1 - Dictionary containing both input dictionaries.
        Same keys that the initial dicionaries.
    Usage:
        dictjoin=joindict(eddydt1,eddydt2)
    '''
    checklist=[]
    checklist1=[]
    checklist2=[]
    for key, value in list(dict1.items()):
        check=False
        if type(value['time'])==int or len(value['time'])==1:
            eddyxt0=value['position'][0]
            eddyyt0=value['position'][1]
            if value['time']>=89-5:
                check=True
                timee=value['time']
        else:
            eddyxt0=value['position'][-1][0]
            eddyyt0=value['position'][-1][1]
            if value['time'][-1]>=89-5:
                check=True
                timee=value['time'][-1]
        if check==True:
            if type(value['time'])==int:
                lonmi0=value['contour'][0][0].min()
                lonma0=value['contour'][0][0].max()
                latmi0=value['contour'][0][1].min()
                latma0=value['contour'][0][1].max()
            else:
                lonmi0=value['contour'][-1][0].min()
                lonma0=value['contour'][-1][0].max()
                latmi0=value['contour'][-1][1].min()
                latma0=value['contour'][-1][1].max()
                
            for key1, value1 in list(dict2.items()):
                if type(value1['time'])==int:
                    ts=value1['time']
                    eddyxt1=value1['position'][0]
                    eddyyt1=value1['position'][1]
                else:
                    ts=int(value1['time'][0])
                    eddyxt1=value1['position'][0][0]
                    eddyyt1=value1['position'][0][1]
                if ts<=5:
                    if (eddyxt1<=lonma0 and eddyxt1>=lonmi0 and eddyyt1<=latma0 and eddyyt1>=latmi0) and\
                        (eddyxt0<=lonma0 and eddyxt0>=lonmi0 and eddyyt0<=latma0 and eddyyt0>=latmi0):
                        dict1[key]={'neddy':int(value['neddy']),'time':np.vstack((value['time'],value1['time']+timee)),\
                                            'position':np.vstack((value['position'],value1['position'])),\
                                            'area':np.vstack((value['area'],value1['area'])),\
                                            'angle':np.vstack((value['angle'],value1['angle'])),\
                                            'ellipse':value['ellipse']+value1['ellipse'],\
                                            'contour':value['contour']+value1['contour'],\
                                            'position_eddy':np.vstack((value['position_eddy'],value1['position_eddy'])),\
                                            'level':np.vstack((value['level'],value1['level'])),\
                                            'minoraxis':np.vstack((value['minoraxis'],value1['minoraxis'])),\
                                            'majoraxis':np.vstack((value['majoraxis'],value1['majoraxis']))}
                        checklist.append(int(value['neddy']))
                        checklist1.append(int(value1['neddy']))
                        checklist2.append(key)
    neweddycount=1
    for key1, value1 in list(dict2.items()):
        if type(value1['neddy']!=checklist1) is np.ndarray:
            check=(value1['neddy']!=checklist1).all()
        else:
            check=(value1['neddy']!=checklist1)
        if check:
            dict1['eddyn_'+str(len(dict1)+neweddycount)]={'neddy':len(dict1)+neweddycount,\
                                            'time':value1['time']+timee,\
                                            'position':value1['position'],\
                                            'area':value1['area'],\
                                            'angle':value1['angle'],\
                                            'ellipse':value1['ellipse'],\
                                            'contour':value1['contour'],\
                                            'position_eddy':value1['position_eddy'],\
                                            'level':value1['level'],\
                                            'minoraxis':value1['minoraxis'],\
                                            'majoraxis':value1['majoraxis']}
            neweddycount=neweddycount+1
    return dict1