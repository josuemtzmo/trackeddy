import numpy as np
import numpy.ma as ma
import pylab as plt
import cmocean as cm
from scipy import linalg
from scipy.interpolate import interp2d,interp1d
from scipy.optimize import curve_fit,leastsq,least_squares

from scipy.optimize import minimize

from scipy import ndimage
from scipy.stats import pearsonr
from scipy import stats
import cartopy.crs as ccrs
from trackeddy.physics import *
from trackeddy.printfunc import *
import pdb

def fit_ellipse(x,y,diagnostics=False):
    '''
    **************** fit_ellipse *****************
    Fitting of an ellipse to an array of positions.
    
    Function translated form Matlab to python by Josue Martinez Moreno,
    the original source:
    Copyright (c) 2003, Ohad Gal 
    All rights reserved.

    Redistribution and use in source and binary forms, with or without 
    modification, are permitted provided that the following conditions are 
    met:

    * Redistributions of source code must retain the above copyright 
    notice, this list of conditions and the following disclaimer. 
    * Redistributions in binary form must reproduce the above copyright 
    notice, this list of conditions and the following disclaimer in 
    the documentation and/or other materials provided with the distribution
    For more information go to the main source:
    https://www.mathworks.com/matlabcentral/fileexchange/3215-fit-ellipse?requestedDomain=www.mathworks.com
    Notes:
    
    Args:
        x,y (array): Coordinates of the datapoints to fit an ellipse.
        diagnostics (boolean): Used to display all the statistics and plots to identify bugs. 
    Returns:
        ellipse_t (dict) - This dictionary contains useful parameters describing completly the ellipsoid ajusted.
        status (boolean) - This value will be true if and only if the the fit corresponds to a ellipse.
    Usage:
    R = np.arange(0,2*pi, 0.01)
    x = 1.5*np.cos(R) + 2 + 0.1*np.random.rand(len(R))
    y = np.sin(R) + 1. + 0.1*np.random.rand(len(R))
    ellipse,status=fit_ellipse(x,y,diagnostics=False)
    '''
    orientation_tolerance = 1e-3
    x=x[:]
    y=y[:]
    mean_x=np.mean(x)
    mean_y=np.mean(y)
    xp=x-mean_x
    yp=y-mean_y
    X=np.array([xp**2,xp*yp,yp**2,xp,yp]).T
    
    a=np.sum(X,axis=0)
    b=np.dot(X.T,X)
    
    x2=np.linalg.lstsq(b.T, a.T)
    
    res=[np.linalg.norm(ii) for ii in b-a*x2[0]]
    
    r2 = np.mean(1 - res / (a.T.size * a.T.var()))
    
    a,b,c,d,e=x2[0]

    if b == 0 and a<c:
        anglexaxis_rad=0
    elif b==0 and c<a:
        anglexaxis_rad=np.pi/2
    else:
        anglexaxis_rad = np.arctan((c-a-np.sqrt((a-c)**2+b**2))/b)

    if ( min(abs(b/a),abs(b/c)) > orientation_tolerance ):
        #TODO: Replace this non sign definite orientation_rad for anglexaxis 
        # which is a sign definite.
        orientation_rad = 1/2 * np.arctan(b/(c-a))
        cos_phi = np.cos( orientation_rad )
        sin_phi = np.sin( orientation_rad )
        a,b,c,d,e = [a*cos_phi**2 - b*cos_phi*sin_phi + c*sin_phi**2,0,a*sin_phi**2 + b*cos_phi*sin_phi + \
                     c*cos_phi**2,d*cos_phi - e*sin_phi,d*sin_phi + e*cos_phi]
        mean_x,mean_y=cos_phi*mean_x - sin_phi*mean_y,sin_phi*mean_x + cos_phi*mean_y       
    else:
        orientation_rad = 0
        cos_phi = np.cos(orientation_rad)
        sin_phi = np.sin(orientation_rad)
    test = a*c
    if test>0 :#and r2>0.8 and r2<1:
        detect='Ellipse'
        status=True
    elif test==0:
        detect='Parabola'
        status=False
    else:
        detect='Hyperbola'
        status=False
    if status==True:        
        # final ellipse parameters
        X0          = mean_x - (d/2)/a
        Y0          = mean_y - (e/2)/c
        F           = 1 + (d**2)/(4*a) + (e**2)/(4*c)
        a           = np.sqrt(abs(F/a))
        b           = np.sqrt(abs(F/c))

        long_axis   = 2*max(a,b)
        short_axis  = 2*min(a,b)

        # rotate the axes backwards to find the center point of the original TILTED ellipse
        R           = np.array([[ cos_phi,sin_phi],[-sin_phi,cos_phi ]])
        P_in        = np.dot(R, np.array([X0,Y0]))
        X0_in       = P_in[0]
        Y0_in       = P_in[1]
        
        ver_line        = np.array([ [X0,X0], [Y0-1*b, Y0+1*b]])
        horz_line       = np.array([ [X0-1*a,X0+1*a], [Y0,Y0] ])
        
        new_ver_line    = np.dot(R,ver_line)
        new_horz_line   = np.dot(R,horz_line)
        
        ### TODO Rotate the axis when the slope is positive.

        theta_r         = np.linspace(0,2*np.pi,len(y));
        ellipse_x_r     = X0 + a*np.cos(theta_r)
        ellipse_y_r     = Y0 + b*np.sin(theta_r)
        rotated_ellipse =  np.dot(R, np.array([ellipse_x_r,ellipse_y_r]))

        # pack ellipse into a structure
        ellipse_t = {'a':a,'b':b,'phi':anglexaxis_rad,'X0':X0,'Y0':Y0,\
                     'X0_in':X0_in,'Y0_in':Y0_in,'long_axis':long_axis,\
                     'short_axis':short_axis,'minoraxis':new_horz_line,\
                     'majoraxis':new_ver_line,'ellipse':rotated_ellipse\
                     ,'status':'Cool'}
        
    else:
        # report an empty structure
        ellipse_t = {'a':'','b':'','phi':'','X0':'','Y0':'',\
                     'X0_in':'','Y0_in':'','long_axis':'',\
                     'short_axis':'','minoraxis':'',\
                     'majoraxis':'','ellipse':'',\
                     'status':detect}
    if (("ellipse" in diagnostics) or ("all" in diagnostics) or (True in diagnostics)) and status==True:
        # draw
        plt.plot( x,y,'b',label='data');
        plt.plot( new_ver_line[0],new_ver_line[1],'k',label='minor axis' )
        plt.plot( new_horz_line[0],new_horz_line[1],'b',label='major axis')
        plt.plot( rotated_ellipse[0],rotated_ellipse[1],'r',label='Fitted ellipse $R^2$=%f' %r2 )
        plt.legend(loc=1)
        plt.show()
        plt.plot(ellipse_x_r,ellipse_y_r,'m')
        plt.plot(horz_line[0],horz_line[1])
        plt.plot(ver_line[0],ver_line[1])
        plt.show()
    return ellipse_t,status,r2

def PolyArea(x,y):
    '''
    *************** PolyArea *******************
    Calculate the area of a poligon.
    Notes:
    
    Args:
        x,y (array): Coordinates of the datapoints to fit an ellipse.
    Returns:
        area (float) -  Area contained by the poligon.
    Usage:
        R = np.arange(0,2*pi, 0.01)
        x = 1.5*np.cos(R) + 2 + 0.1*np.random.rand(len(R))
        y = np.sin(R) + 1. + 0.1*np.random.rand(len(R))
        area=PolyArea(x,y)
    '''
    area=0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
    return area

def eccentricity(a,b):
    '''
    *************** eccentricity *******************
    This function calculate the eccentricity of a ellipse.
    Notes:
        
    Args:
        a (float): Mayor axis of an ellipse
        b (float): Minor axis of an ellipse
    Returns:
        eccen (float) - Eccentricity of the ellipsoid with parameters a,b.
    Usage:
        a=0.5
        b=0.3
        eccen=eccentricity(a,b)
    '''
    a=abs(a)
    b=abs(b)
    if b>a:
        b1=a
        a=b
        b=b1
    eccen=np.sqrt(1-(abs(b)**2/abs(a)**2))
    return eccen

def find2l(arrayx,arrayy,valuex,valuey):
    '''
    *************** find2l *******************
    Find values in two list of values.
    Notes:
        
    Args:
        arrayx (list|array): Array where it will look the closer index to the valuex
        arrayy (list|array): Array where it will look the closer index to the valuey
        valuex (int|float): Value to look for in arrayx.
        valuey (int|float): Value to look for in arrayy.
    Returns:
        idx,idy (int) - Index of closer values.
    Usage:
        arrayx=[0,1,2,3]
        arrayy=[4,5,6,7]
        valuex=2.2
        valuey=6.6
        indexes=find2l(arrayx,arrayy,valuex,valuey)
    '''
    idx=(np.abs(arrayx-valuex)).argmin()
    idy=(np.abs(arrayy-valuey)).argmin()
    return idx,idy

def find(array,value):
    '''
    *************** find *******************
    Find values in a list of values.
    Notes:
        
    Args:
        array (list|array): Array where it will look the closer index to the value
        value (int|float): Value to look for in array.
    Returns:
        idx - Index of closer values.
    Usage:
        array=[0,1,2,3]
        value=2.2
        idx=find(array,value)
    '''
    idx=int(np.mean(np.where(array==value)))
    return idx

def find2D(array,value):
    '''
    *************** find2D *******************
    Find values in a 2D array of values.
    Notes:
        
    Args:
        array (list|array): 2D Array where it will look the closer index to the value
        value (int|float): Value to look for in array.
    Returns:
        yp,xp - Index of closer values.
    Usage:
        array=[[0,1,2,3],[0,1,4,3]]
        value=2
        idx,idy=find2D(array,value)
    '''
    yp,xp=np.where(array==value)
    return yp,xp
    
def contourmaxvalue(var,x,y,levels,date=''):
    '''
    *************** contourmaxvalue *******************
    Find the maximum value inside an specific contour.
    Notes:
        
    Args:
        contcoordx (list|array): Contains the coordinates in X of the contour of the field var. 
        contcoordy (list|array): Contains the coordinates in Y of the contour of the field var.
        var (array): 3D Matrix representing a surface (np.shape(var)=(date,len(x),len(y))).
        x (list|arrat): Contains the coordinate X of the grid of var.
        y (list|arrat): Contains the coordinate Y of the grid of var.
        levels (list): Level of the extracted contour.
        date (int): Used if len(var)==3 (i.e. Var contains a time dimension).
    Returns:
        coord (list) - Location of the max value in the grid.
    Usage:
        center_eddy=contourmaxvalue(contcoordx,contcoordx,sshnan,lon,lat,levels,date)
    '''
    if len(np.shape(var))==3 or date=='':
        if levels[0]>0:
            sshextrem=np.nanmax(var[date,:,:])
        else:
            sshextrem=np.nanmin(var[date,:,:])
        indexes=find2D(var[date,:,:],sshextrem)
    else:
        if levels[0]>0:
            sshextrem=np.nanmax(var)
        else:
            sshextrem=np.nanmin(var)
        indexes=find2D(var[:,:],sshextrem)
    coord=[x[indexes[1][0]],y[indexes[0][0]],sshextrem,indexes[1][0],indexes[0][0]]
    return coord

def centroidvalue(contcoordx,contcoordy,var,x,y,levels,date,threshold=1):
    '''
    *************** centroidvalue *******************
    Find the centroid inside an specific contour.
    Notes:
        
    Args:
        contcoordx (list|array): Contains the coordinates in X of the contour of the field var. 
        contcoordy (list|array): Contains the coordinates in Y of the contour of the field var.
        var (array): 3D Matrix representing a surface (np.shape(var)=(date,len(x),len(y))).
        x (list|arrat): Contains the coordinate X of the grid of var.
        y (list|arrat): Contains the coordinate Y of the grid of var.
        levels (list): Level of the extracted contour.
        date (int): Used if len(var)==3 (i.e. Var contains a time dimension).
    Returns:
        coord (list) - Location of the centroid in the grid.
    Usage:
        center_eddy=centroidvalue(contcoordx,contcoordx,sshnan,lon,lat,levels,date)
    '''
    if len(np.shape(var))==3:
        var=var.filled(0)
        var=np.abs(var)
        sum_T=np.nansum(var[date,:,:])
        sum_X=np.nansum(var[date,:,:],axis=0)
        sum_Y=np.nansum(var[date,:,:],axis=1)
        XM=0
        for ii in range(len(sum_X)):
            XM=XM+(sum_X[ii]*x[ii])
        YM=0
        for ii in range(len(sum_Y)):
            YM=YM+(sum_Y[ii]*y[ii])
        xcpos=XM/sum_T
        ycpos=YM/sum_T
    else:
        var=var.filled(0)
        var=np.abs(var)
        sum_T=np.nansum(var)
        sum_X=np.nansum(var,axis=0)
        sum_Y=np.nansum(var,axis=1)
        XM=0
        for ii in range(len(sum_X)):
            XM=XM+(sum_X[ii]*x[ii])
        YM=0
        for ii in range(len(sum_Y)):
            YM=YM+(sum_Y[ii]*y[ii])
        xcpos=XM/sum_T
        ycpos=YM/sum_T
    coord=np.asarray([xcpos,ycpos])
    return coord

def adjust1Gaus(x,y):
    '''
    *************** adjust1Gaus *******************
    Fit one gaussian in a curve curve.
    Notes:
    
    Args:
        x(list|array): Coordinates in x of data.
        y(list|array): Data to be ajusted with one gaussian. 
    Returns:
        gausfit(list|array) - Data ajusted.
    Usage:
        x=np.arange(-5,5,0.1)
        x0=0
        a=3
        sigma=2
        gaussian=gaus(x,a,x0,sigma)
        gaussianfit=adjust1Gaus(x,gaussian)
    '''
    gauss_fit = lambda p, x: p[0]*(1/np.sqrt(2*np.pi*(p[2]**2)))*np.exp(-(x-p[1])**2/(2*p[2]**2)) #1d Gaussian func
    e_gauss_fit = lambda p, x, y: (gauss_fit(p,x) -y) #1d Gaussian fit

    v0= [1,10,1,1,30,1] #inital guesses for Gaussian Fit. - just do it around the peaks
    out = leastsq(e_gauss_fit, v0[:], args=(x, y), maxfev=2000, full_output=1) #Gauss Fit
    v = out[0] #fit parameters out
    covar = out[1] #covariance matrix output

    gausfit = gauss_fit(v,x) 
    return gausfit

def twoD_Paraboloid(coords, amplitude, xo, yo, a, b,offset):
    '''
    *************** twoD_Gaussian *******************
    Build a 2D gaussian.
    Notes:
        Remmember to do g.ravel().reshape(len(x),len(y)) for plotting purposes. 
    Args:
        coords [x,y] (list|array): Coordinates in x and y.
        amplitude (float): Amplitud of gaussian.
        x0 , yo (float): Center of Gausian.
        sigma_x,sigma_y (float): Deviation.
        theta (Float): Orientation.
        offset (Float): Gaussian Offset.
    Returns:
        g.ravel() (list|array) - Gaussian surface in a list.
    Usage:
        Check scan_eddym function.
    '''
    x=coords[0]
    y=coords[1]
    amplitude = coords[2]
    
    xo = float(coords[3])
    yo = float(coords[4])
    
    #a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    #b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    #c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    
    #g =offset - amplitude * (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2))
    
    g = -amplitude*(((x-x0)/a)**2+((y-y0)/b)**2) + offset
    
    return g.ravel()

#def twoD_Gaussian(coords,amplitude,xo,yo, sigma_x, sigma_y, theta, slopex=0, slopey=0, offset=0):
def twoD_Gaussian(coords, sigma_x, sigma_y, theta, slopex=0, slopey=0, offset=0):
    '''
    *************** twoD_Gaussian *******************
    Build a 2D gaussian.
    Notes:
        Remmember to do g.ravel().reshape(len(x),len(y)) for plotting purposes. 
    Args:
        coords [x,y] (list|array): Coordinates in x and y.
        amplitude (float): Amplitud of gaussian.
        x0 , yo (float): Center of Gausian.
        sigma_x,sigma_y (float): Deviation.
        theta (Float): Orientation.
        offset (Float): Gaussian Offset.
    Returns:
        g.ravel() (list|array) - Gaussian surface in a list.
    Usage:
        Check scan_eddym function.
    '''
    x=coords[0]
    y=coords[1]
    amplitude = coords[2]
    
    xo = float(coords[3])
    yo = float(coords[4])

    #print(sigma_y,sigma_x,sigma_y/sigma_x)
    if sigma_y or sigma_x != 0:
        #a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
        #b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
        #c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
        #g = amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
        cos_phi = np.cos(theta)
        sin_phi = np.sin(theta)
        a = (cos_phi**2)/(2*sigma_x**2) + (sin_phi**2)/(2*sigma_y**2)
        b = (np.sin(2*theta))/(4*sigma_x**2) - (np.sin(2*theta))/(4*sigma_y**2)
        c = (sin_phi**2)/(2*sigma_x**2) + (cos_phi**2)/(2*sigma_y**2)
        g = amplitude*np.exp(-(a*(x-xo)**2 + 2*b*(x-xo)*(y-yo) + c*(y-yo)**2))
    else: 
        g = (x-xo)*0 + (y-yo)*0
    return g.ravel()

def gaussian2Dresidual(popt, coords, varm):
    g=twoD_Gaussian(coords,*popt).reshape(np.shape(varm))
    residual = np.exp(np.abs(np.float128(np.nanmean(abs(varm - g))))) - 1
    #print('Residual:',np.nanmean(residual))
    return residual

def paraboloid2Dresidual(popt,coords,varm):
    residual = np.exp(np.abs(varm - twoD_Gaussian(coords,*popt))) - 1
    return residual

def correlation_coefficient(data, data1):
    product = np.mean((data - data.mean()) * (data1 - data1.mean()))
    stds = data1.std() * data1.std()
    if stds == 0:
        return 0
    else:
        product = product/stds
        return product


def fit2Dcurve(var,values,level,initial_guess='',date='',mode='gaussian',diagnostics=False):
    '''
    *************** fit2Dgaussian *******************
    Fit a surface to the data.
    Notes:
        
    Args:
        
    Returns:
        
    Usage:
        
    '''

    if type(diagnostics) != list:
        diagnostics=[diagnostics]
    Lon, Lat = np.meshgrid(values[0], values[1])
    coords=(Lon,Lat,values[2],values[3],values[4])
    if date!='':
        varm=var[date,:,:]*1
    else:
        varm=var*1
    mask=ma.getmask(varm[:,:])
    
    if initial_guess=='':
        initial_guess = [1,1,0,0,0,0]   
        
    if mode == 'parabolic':
        popt, pcov, infodict,mesg,ier = leastsq(paraboloid2Dresidual, initial_guess,\
                                                args=(coords, varm.ravel()),full_output=True)#,\
        fitdict = popt
        fitted_curve = twoD_Paraboloid(coords, *fitdict)
    elif mode == 'best':
        popt, pcov = leastsq(paraboloid2Dresidual, initial_guess, args=(coords, varm.ravel()))
        popt, pcov = leastsq(gaussian2Dresidual, initial_guess, args=(coords, varm.ravel())) 
        fitdict = popt
    elif mode == 'both':
        popt, pcov = leastsq(gaussian2Dresidual, initial_guess, args=(coords, varm.ravel())) 
        popt, pcov = leastsq(paraboloid2Dresidual, initial_guess, args=(coords, varm.ravel())) 
        fitdict = popt
    else:
        #print("\n ----Fit----")
        #pdb.set_trace()
        res = minimize(gaussian2Dresidual, initial_guess,args=(coords,varm),
               method='SLSQP',options={'xtol': 1e-12, 'disp': False})
        fitdict = res.x
    
        #popt, pcov, infodict,mesg,ier = leastsq(gaussian2Dresidual, initial_guess,\
        #                                        args=(coords, varm.ravel()),\
        #                                        xtol=1e-10,maxfev=1000000,
        #                                        full_output=True) 
        #fitdict = popt
    
    fitted_curve = twoD_Gaussian(coords, *fitdict)
    fittedata=fitted_curve.reshape(len(values[1]), len(values[0]))
    
    try:
        R2=correlation_coefficient(varm.ravel(),fittedata.ravel())
    except:
        R2=0
    
    if type(varm) == ma.core.MaskedArray:
        o=sum(sum(varm.filled(0)))
        g=sum(sum(ma.masked_array(fittedata, varm.mask).filled(0)))
        if o*0.9 <= g and g <= o*1.1:
            R2=0
            
    if ("2dfit" in diagnostics) or ("all" in diagnostics) or (True in diagnostics):
        print('R^2 2D fitting:',R2)
        print('OPT steps:',res.nfev)
        print("             |sigmaX|sigmaY|Theta|slopeX|slopeY|")
        print("initial guess|" + ''.join(str(e)+'|' for e in initial_guess))
        print("Fit.         |" + ''.join(str(e)+'|' for e in fitdict))
        f, (ax1, ax2,ax3) = plt.subplots(1, 3,figsize=(25,7),sharey=True)
        p=ax1.pcolormesh(values[0],values[1],varm,vmin=varm.min(),vmax=varm.max())
        ax1.set_title('Original Field')
        ax1.axis('equal')
        plt.colorbar(p,ax=ax1)
        p=ax2.pcolormesh(values[0], values[1],fittedata,vmin=varm.min(),vmax=varm.max())
        ax2.set_title('2D Gauss Fit')
        ax2.axis('equal')
        plt.colorbar(p,ax=ax2)
        p=ax3.pcolormesh(values[0], values[1],varm-fittedata,\
                       cmap=cm.cm.balance)
        ax3.axis('equal')
        ax3.set_title('Difference between Fit & Original')
        plt.colorbar(p,ax=ax3)
        plt.show()
        plt.close()
    return fitdict,R2

def rsquard(y,yfit):
    '''
    *************** rsquard *******************
    Calculate the Pearson Coefficient.
    Notes:
        Make sure the x grid coincide at least with indexes for y and y1.
    Args:
        y(list|array): Original data.
        yfit(list|array): Data ajusted. 
    Returns:
        R2 (float): Pearson Coefficient.
    Usage:
        x=np.arange(-5,5,0.1)
        x0=0
        a=3
        sigma=2
        gaussian=gaus(x,a,x0,sigma)+gaus(x,a-2,x0+2,sigma-1)
        gaussianfit=adjustMGaus(x,gaussian)
        R2=rsquard(gaussian,gaussianfit)
    '''
    #yhat=yfit            # or [p(z) for z in x]
    #ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    #ssreg = np.sum((y-yhat)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    #sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    #R2 = 1 - ssreg / sstot
    R2=pearsonr(y,yfit)[0]**2
    return R2 

def ellipsefit(y,yfit,ellipsrsquarefit=0.85,diagnostics=False):
    '''
    *************** ellipsefit *******************
    Check the fitness of an ellipse in a curve.
    Notes:
        
    Args:
        y(list|array): Original data.
        yfit(list|array): Ellipse ajusted to contour. 
        ellipsrsquarefit (float): [0 > ellipsrsquarefit < 1]
            Pearson Coefficient to validate an ellipse.
        diagnostics (boolean): Used to display all the 
            statistics and plots to identify bugs.
    Returns:
        Rsquard (float) - Fitness of the ellipse.
        check (boolean) - True if gaussian adjust is greater than gaussrsquarefit.
    Usage:
        Check scan_eddym function.
    '''
    x=range(0,len(y))
    f=interp1d(x, y)
    xnew=np.linspace(0,len(y)-1,len(yfit))
    eddy2fit=f(xnew)
    
    indxed=find(yfit,yfit.max())
    indxrd=find(eddy2fit,eddy2fit.max())
    
    if yfit[indxed]==yfit[indxed+1] and yfit[indxed]==yfit[indxed-1]:
        indxed=find(yfit,yfit.min())
        indxrd=find(eddy2fit,eddy2fit.min())
    
    eddy2fit=list(eddy2fit)*2
    eddyfitdisplace=np.zeros(len(yfit))
    for ii in range(len(yfit)):
        eddyfitdisplace[ii]=eddy2fit[indxrd-indxed+ii]
    Rsquard=rsquard(eddyfitdisplace,yfit)
    if ("ellipse" in diagnostics) or ("all" in diagnostics) or (True in diagnostics):
        plt.figure()
        plt.title('Ellipse Fit')
        plt.plot(yfit,'-b')
        plt.plot(eddyfitdisplace,'-r')
        plt.text(0, np.mean(yfit), str(round(Rsquard,2)))
        plt.show()
        print(sum(yfit),sum(eddyfitdisplace))
    if Rsquard>=ellipsrsquarefit:# and Rsquard < 1:
        check=True
    else:
        check=False
    return Rsquard,check
    
def extractprofeddy(axis,field,lon,lat,n,gaus='One',kind='linear',gaussrsquarefit=0.65,varname='',diagnostics=False,plotprofile=False):
    '''
    *************** extractprofeddy *******************
    Extracts the profile inside a segment.
    Notes:
        Primary used to extract the mayor or minor axis of an ellipse.
    Args:
        axis (list): Coordinates of an segment.
        field (array): Surface where the profile will be extracted.
        lon,lat (array|list): Coordinates fo the field.
        n (int): Number of desired divisions in the segment.
        gaus (Default:One|None|Multiple): Ajustment to segment.
        kind (Default:linear|cubic|etc): Type of interpolation inside
            segment (For more information check scipy.interpolate interp2d)
        gaussrsquarefit (float): [0 > ellipsrsquarefit < 1]
            Pearson Coefficient to validate an gaussian.
        varname (str): Name of variable, just used for plots.
        diagnostics (boolean): Used to display all the 
            statistics and plots to identify bugs.
    Returns:
        axisdata (array) - Data extracted in the segment.
        check (boolean) - True if gaussian adjust is greater than gaussrsquarefit.
        field_interp (array) -
    Usage:
        Check scan_eddym function.
    '''
    try:
        fieldnan=field.filled(np.nan)
    except:
        field[~np.isfinite(field)]=np.nan
        fieldnan=field
        
    if type(diagnostics) != list:
        diagnostics=[diagnostics]

    #fieldnan=field
    #print(nanmax(field),nanmin(field))
    ycoord=np.linspace(axis[1,0],axis[1,1],n)
    xcoord=np.linspace(axis[0,0],axis[0,1],n)
    
    field2interp=interp2d(lon, lat, fieldnan[:,:], kind=kind)
    field_interp = field2interp(xcoord,ycoord)
    axisdata=np.zeros([n])
    for ii in range(n):
        axisdata[ii]=field_interp[ii,ii]

    n = len(axisdata)  #the number of data
    x = np.array(range(n))
    y = axisdata
    
    p = np.poly1d(np.polyfit(x, y, 1))
    slope, intercept, r_value, p_value, std_err = stats.linregress(y,p(x))
    linearfit = rsquard(y,p(x))
    
    if gaus=='None':
        Rsquared=1
    else:
        if gaus=='One':
            gausfit=adjust1Gaus(x,y)
        elif gaus=='Multiple':
            gausfit=adjustMGaus(x,y)
        else:
            print('Select a gausian method to adjust.')
            return
        
        Rsquared = rsquard(y,gausfit)
        #print('R2_gauss: ',Rsquared)
    #print('test 746 geo:',Rsquared)
    if Rsquared >= gaussrsquarefit and linearfit < 0.5:
        check=True
    else:
        check=False
        
    if ("all" in diagnostics) or ("gauss" in diagnostics) or (True in diagnostics):
        plt.plot(xcoord,ycoord)
        plt.title('Gauss_fit')
        plt.pcolormesh(xcoord,ycoord,field_interp)
        plt.show()
        print('std',varname,' vs fit',Rsquared)
        plt.plot(x,y,'b+:',label='data')
        plt.plot(x,gausfit,'ro:',label='fit')
        plt.legend()
        plt.title('Fit for Time Constant')
        plt.xlabel('Position (n)')
        plt.ylabel(varname)
        plt.show()
    if plotprofile==True:
        return y,gausfit,check
    else:
        return y,check

def eddylandcheck(contour,lon,lat,var,diagnostics=False):
    '''
    *************** eddylandcheck *******************
    Check if the contour is surrounded by land.
    Notes:
        
    Args:
        
    Returns:
        
    Usage:
        
    '''
    if type(diagnostics) != list:
        diagnostics=[diagnostics]
    land = {'True':{},'False':{}}
    lcount = 0
    fcount = 0
    checkland = True
    for ii in range(0,len(contour[:,0])):
        idxcheck,idycheck = find2l(lon,lat,contour[ii,0],contour[ii,1])
        if idxcheck == len(lon):
            ar=0
        elif idycheck == len(lat):
            ar=0
        else:
            ar=1
        checkarea=var[idycheck-ar:idycheck+ar,idxcheck-ar:idxcheck+ar]
        
        if np.isnan(checkarea).any() and type(var)!=ma.core.MaskedArray:
            land['True'][str(lcount)] = {'x':idxcheck,'y':idycheck}
            lcount=lcount+1
        elif type(var)==ma.core.MaskedArray and checkarea.mask.any() :
            land['True'][str(lcount)] = {'x':idxcheck,'y':idycheck}
            lcount=lcount+1
        else:
            land['False'][str(fcount)] = {'x':idxcheck,'y':idycheck}
            fcount=fcount+1
        
    if len(land['True']) >= len(contour[:,0])/2:
        checkland=False
        
    if ("landcheck" in diagnostics) or ("all" in diagnostics) or (True in diagnostics):
        plt.figure()
        plt.pcolormesh(lon,lat,var,cmap=cm.cm.deep)
        plt.plot(contour[:,0],contour[:,1])
        for key,item in land.items():
            for key1,item1 in item.items():
                if key == 'True':
                    plt.plot(lon[item1['x']],lat[item1['y']],'or',label='Invalid land grid points ($n_i = %d$)' %lcount if key1 == str(len(item.keys())-1) else "")                
                else:
                    plt.plot(lon[item1['x']],lat[item1['y']],'og',label='Valid grid points ($n_v = %d$)' %(len(contour[:,0])-lcount) if key1 == str(len(item.keys())-1)  else "")
        plt.title('Land check $n_g = %d$' % len(contour[:,0]) )
        plt.legend()
    return checkland

def reconstruct_syntetic(varshape,lon,lat,eddytd,mode='gaussian',rmbfit=False,usefullfit=False,diagnostics=False,one_time=None,debug=False):
    '''
    *************** reconstruct_syntetic *******************
    Recunstruct the syntetic field using the gaussian 
    parameters saved in the dictionary of eddies.
    Notes:
        
    Args:
        
    Returns:
        
    Usage:
    
    '''
    if debug==True:
        print("\n ******* Reconstruct ******")
        pdb.set_trace()
    Lon,Lat=np.meshgrid(lon,lat)
    fieldfit=np.zeros(varshape)
    if type(diagnostics) != list:
        diagnostics=[diagnostics]
    pp =  Printer(); 
    keys=tuple(eddytd.keys())
    loop_len=len(keys)
    for xx in range(0,loop_len):
        key=keys[xx]
        counter=0
        if one_time!=None and type(one_time)==int:
            timeloop=[one_time]
        else:
            timeloop=range(0,len(eddytd[key]['time']))
        for tt in timeloop:
            ttt=eddytd[key]['time'][tt]
            level=eddytd[key]['level'][tt]
            maxposition=[Lon,Lat,eddytd[key]['position_maxvalue'][counter][2],\
                         eddytd[key]['position_maxvalue'][counter][0],\
                         eddytd[key]['position_maxvalue'][counter][1]]
            curvefit=eddytd[key]['2dgaussianfit'][counter]
            if isinstance(curvefit, np.float64):
                curvefit=eddytd[key]['2dgaussianfit']
            #Remove the slope and constant in the reconstruction of the eddy.
            if mode == 'parabolic':
                fittedcurve=twoD_Paraboloid(maxposition, *curvefit[:])
                #print(level)
                if level>0:
                    fittedcurve[fittedcurve<0]=0
                else:
                    fittedcurve[fittedcurve>0]=0
            elif mode == 'best':
                print('Work in progress')            
            elif mode == 'both':
                print('Work in progress')
            else:
                #if usefullfit==False:
                #    curvefit[-1]=0
                #    curvefit[-2]=0
                #    curvefit[-3]=0
                #print(gaussfit)
                fittedcurve=twoD_Gaussian(maxposition, *curvefit)
            #print(fittedcurve[0])
            if np.isnan(fittedcurve[0]):
            #or (curvefit[0]/curvefit[1]+curvefit[1]/curvefit[0])/2>1.7:
            #or curvefit[0]/curvefit[1]>1.7 or curvefit[1]/curvefit[0]>1.7:
                fittedcurve=np.zeros(np.shape(fittedcurve))
            else:
                fieldfit[ttt,:,:]=fieldfit[ttt,:,:]+fittedcurve.reshape(len(lat),len(lon))
            #print(fieldfit[ttt,:,:])
            #plt.pcolormesh(fieldfit[ttt,:,:])
            #plt.show()
            counter=counter+1
        if ("reconstruct" in diagnostics) or ("all" in diagnostics) or (True in diagnostics):
            ax = plt.axes(projection=ccrs.PlateCarree())
            print('key: ',key,'Level: ',level)
            plt.pcolormesh(lon[::5],lat[::5],fieldfit[0,::5,::5])
            ax.coastlines()
            plt.colorbar()
            plt.show()
        pp.timepercentprint(0,loop_len,1,xx,key)
    return fieldfit
    
def insideness_contour(data,center,levels,mask=False,maskopt=None,diagnostics=False):
    '''
    
    '''
    if type(diagnostics) != list:
        diagnostics=[diagnostics]
    data_rmove=np.array(np.zeros(np.shape(data)))
    if (type(levels)==int or len(levels)==1 ) and levels < 0:
        data_rmove[data<levels]=1
    elif (type(levels)==int or len(levels)==1 ) and levels > 0:
        data_rmove[data>levels]=1
    elif levels[1]<0:
        data_rmove[data<levels[1]]=1
    else:
        data_rmove[data>levels[0]]=1
        
    markers,features=np.asarray(ndimage.label(data_rmove))
    
    if markers.max()!=1 and maskopt==None:
        markers=markers*0
        returnmasked=True    
        
    elif markers.max()!=1 and maskopt=='contour' or maskopt=='forcefit':
        if center[1]==np.shape(markers)[1]:
            markers[markers!=markers[center[0],center[1]-1]]=0
        elif center[0]==np.shape(markers)[0]:
            markers[markers!=markers[center[0]-1,center[1]]]=0
        else:
            markers[markers!=markers[center[0],center[1]]]=0
        markers=markers.max()-markers
        returnmasked=True 
            
    elif maskopt=='contour' or maskopt=='forcefit':
        markers=1-markers
    else:
        markers=markers*0  
    
    if levels[0]<0 and maskopt!='forcefit':
        markers[data>0]=1
    elif levels[0]>0 and maskopt!='forcefit':
        markers[data<0]=1
    elif levels[0]<0 and maskopt=='forcefit':
        markers[np.multiply(data<levels[1]/2, data>levels[1]-levels[0]/4)]=0
        data[np.multiply(data<levels[1]/2, data>levels[1]-levels[0]/4)]=0
    elif levels[0]>0 and maskopt=='forcefit':
        markers[np.multiply(data>levels[0]/2, data<levels[0]-levels[0]/4)]=0
        data[np.multiply(data>levels[0]/2, data<levels[0]-levels[0]/4)]=0
        
    maskeddata=ma.masked_array(data, markers)
    if ("contours" in diagnostics) or ("all" in diagnostics) or (True in diagnostics):
        print('markerceter:',markers[center[0],center[1]])
        f, (ax1, ax2,ax3) = plt.subplots(1, 3,figsize=(15,7), sharey=True)
        ax1.pcolormesh(data)
        ax1.set_title('original data')
        m=ax2.pcolormesh(markers)
        plt.colorbar(m,ax=ax2)
        msk=ax2.plot(center[1],center[0],'or')
        ax2.set_title('identified eddy mask')
        ax3.pcolormesh(maskeddata)
        ax3.set_title('masked data')
        plt.show()
    if mask==True:
        return maskeddata,markers 
    else:
        return maskeddata
              
    
def gaussareacheck(values,level,areaparms,gauss2dfit,contour_area,contour_x=None,contour_y=None):
    #print('gauss check')
    Lon, Lat = np.meshgrid(values[0], values[1])
    coords=(Lon,Lat,values[2],values[3],values[4])
    fitted_curve = twoD_Gaussian(coords, *gauss2dfit)
    fittedata = fitted_curve.reshape(len(values[1]),len(values[0]))
    #fittedata = ma.masked_array(fittedata, mask)
    try:
        if level>0:
            CS=plt.contour(values[0],values[1],fittedata,levels=[level,np.inf])
        else:
            CS=plt.contour(values[0],values[1],fittedata,levels=[level,np.inf])
        plt.close()
    
        CONTS=CS.allsegs[0][0]
        areastatus = checkscalearea(areaparms,np.mean(CONTS[:,0]),np.mean(CONTS[:,1]),CONTS[:,0],CONTS[:,1])
    except:
        return False,0
        
    if areastatus['ellipse'] == None:
        test=False
    elif (contour_area*1.5 > areastatus['ellipse']) and areastatus['status']: #and area[1]!=0:
        test=True
    else:
        test=False
    #print('---------',test)
    return test,areastatus['ellipse']

def checkgaussaxis2D(a,b,a_g,b_g):
    #print('a',a,a_g,'b',b,b_g)
    if a*2<a_g or b*2<b_g:
        return False
    else:
        return True
    
def check_closecontour(contour,lon_contour,lat_contour,var):
    if (len(contour[:,0]) | len(contour[:,1])) <= 8:
        return False
    elif contour[0,0] != contour[-1,0] and contour[0,1] !=contour[-1,1]:
        return False
    else:
        checkland=eddylandcheck(contour,lon_contour,lat_contour,var)
        if checkland == False:
            return False
        else: 
            return True
