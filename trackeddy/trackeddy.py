from __future__ import print_function

import numpy as np
import numpy.ma as ma
import pylab as plt

from trackeddy.physics import *
from trackeddy.geometryfunc import *
from trackeddy.decorators import *

import seawater as sw
from scipy import ndimage
from astropy import convolution
import sys
import time
import pdb 
from skimage.measure import find_contours, EllipseModel
import scipy.ndimage.filters as filters
import scipy.special as sspecial
import pandas as pd

import copy

earth_radius = 6371e3 #meters

np.seterr(all='ignore')

# New version: 31.5 s ± 1.02 s per loop (mean ± std. dev. of 7 runs, 1 loop each)
# Old version: 48.9 s ± 1.68 s per loop (mean ± std. dev. of 7 runs, 1 loop each)

# New version memory: peak memory: 255.98 MiB, increment: 0.03 MiB
# Old version memory: peak memory: 232.28 MiB, increment: -0.34 MiB


class TrackEddy(object):
    """
    
    """
    @check_input
    def __init__(self,dataset,coords=None,variable=None,**kwargs):
        self.DataArray = dataset[variable].squeeze()
        self.coords = coords
        self.identified_eddies={}
        
        dx,dy = self.coords['X'].diff('X').mean() , self.coords['Y'].diff('Y').mean()
        
        self.preferences = {'ellipse_error':[2*dy,2*dx],'eccentricity':0.5,'gaussian':0.8}
        self.spatial_filter = {'type':'convolution','kernel':10}
        self.mask_contour_opt = "contour"
        self.fit_gaussian=True

    @check_single_level
    def _scan_eddy_single_level(self,level,polarity='both',geo=True):
        # We are using skimage.measure.find_contours, which returns the contour using indexes, thus for each contour we will convert to xy space.
        self.identified_eddies_level = pd.DataFrame({})
        close_contours_ji = extract_contours(self.DataArray,self.level,self.coords)
        number_close_contours = len(close_contours_ji)
        coords_range = coordinates_range(self.coords)        
        
        contours=[]
        contours_rossby=[]
        n = 1
        for n_contour in range(0,number_close_contours):
            # Check if contour is close.
            if not is_contour_close(close_contours_ji[n_contour]):
                continue
            # Convert ij contour to xy contour.
            contour_yx_deg = contour_ji_2_yx(close_contours_ji[n_contour],
                                coords_range,
                                self.DataArray.shape)
            
            # contours.append(contour_yx_deg)
            # Compute area of polygon.
            area = poly_area(close_contours_ji[n_contour],
                                coords_range,
                                self.DataArray.shape,
                                geo)
            # Compute equivalent radius.
            radius = np.sqrt(area)/np.pi
            Rossby_radius = rossbyR(*np.flipud(contour_yx_deg.mean(axis=0)))
            # Check area of closed contour.
            if radius >= Rossby_radius:
                continue
            # Fit ellipse
            ellipse_dict = ellipse_fit_leastsq(contour_yx_deg)
            # Check if fitting a ellipse was succesful.
            if not ellipse_dict:
                continue
            # Compute residuals, sum of residual distances must be smaller than 
            # twice the grid cell.
            error_y = self.preferences['ellipse_error'][0]
            error_x = self.preferences['ellipse_error'][1]
            sum_error_distances = np.sum(ellipse_dict['residual']) 
            if sum_error_distances > error_y or sum_error_distances > error_x:
                continue
            # Compute eccentricity using fitted ellipse parameters.
            eccen = eccentricity(ellipse_dict['a'],ellipse_dict['b'])
            # 
            if eccen <= self.preferences['eccentricity']:
                continue

            # Get centre location of contour 
            # (i.e. top center location of close contour)
            xslice, yslice, centertop = get_centre_location(self.coords,contour_yx_deg)

            # Extract sliced data as np.array
            data4gauss = self.DataArray.isel({'X':xslice,'Y':yslice}).copy()
            # Check that subset of data is 2D.
            if 0 in data4gauss.shape:
                continue
            # Eddy centre is the centroid of the close contour.
            eddy_c_location = find_centre_and_max(contour_yx_deg,data4gauss,centertop)
            
            sign = eddy_c_location[3]
            #TODO add gaussian profile test?
            if not self.fit_gaussian:
                # TODO: finish coding single level dict
                # TODO: add residual field 
                data_dict = {
                    'n': n,
                    'amp': eddy_c_location[2],
                    'x': eddy_c_location[1], 
                    'y': eddy_c_location[0],
                    'a': ellipse_dict['a'],
                    'b': ellipse_dict['b'],
                    'theta': ellipse_dict['theta'],
                    'level': sign*abs(level),
                    'radius': radius,
                    'area': area,
                    'residual': 0,
                    #'region_indexes':np.hstack((xslice, yslice))
                    'polarity':sign
                }
                identified_eddies_level =  self.single_level_dict(**data_dict)
                n+=1
                continue

            # Extract data inside contour and mask regions of the data.

            data_inside_contour = extract_inside_contour(data4gauss.values,centertop, sign, sign*abs(level), self.mask_contour_opt)
            # FIT gaussian (Keep consistency with ji indexes)
            fixvalues = [self.coords['Y'].isel(Y=yslice).values,self.coords['X'].isel(X=xslice).values]
            # Initial guess follows the convention:
            # [Amplitud,  x0 , y0, ellipse_mayor_axis, ellipse_minor_axis, ellipse_angle, theta]
            initial_guess = [ eddy_c_location[2], eddy_c_location[1], eddy_c_location[0],\
                              ellipse_dict['a'],ellipse_dict['b'],ellipse_dict['theta'] ]

            gausssianfitp  = minimize_surface_fitting(data_inside_contour,fixvalues,initial_guess=initial_guess)

            gaussian_params = gausssianfitp.x
            
            data_dict = {
                'n': n,
                'amp': eddy_c_location[2],
                'x': eddy_c_location[1], 
                'y': eddy_c_location[0],
                'a': ellipse_dict['a'],
                'b': ellipse_dict['b'],
                'theta': ellipse_dict['theta'],
                'level': sign*abs(level),
                'gausian_amp': gaussian_params[0],
                'gausian_x': gaussian_params[1],
                'gausian_y': gaussian_params[2],
                'gausian_a': gaussian_params[3],
                'gausian_b': gaussian_params[4],
                'gausian_theta': gaussian_params[5],
                'radius': radius,
                'area': area,
                'residual': gausssianfitp.fun,
                'polarity':sign
                #'region_indexes':np.hstack((xslice, yslice))
            }

            self.identified_eddies_level = self.identified_eddies_level.append(self.single_level_dict(**data_dict))
            
            n+=1
            #print(gausssianfitp)
            # gausssianfitp, R2 = fit2Dcurve(data_inside_contour,\
            #                 fixvalues,\
            #                 level,initial_guess=initial_guess,date='',\
            #                 mode=mode,diagnostics=diagnostics)

            
            # theta_r = np.linspace(0, 2 * np.pi, len(contour_yx_deg))
            # yx = ellipse_dict['generator'].predict_xy(theta_r) + contour_yx_deg.mean(axis=0)

            #contours_rossby.append(contour_yx_deg)

            # import cmocean as cm
            # plt.figure(figsize=(6,3))

            # Lon,Lat = np.meshgrid(fixvalues[1],fixvalues[0])
            # init_gaussian = construct_gaussian((Lat,Lon),*initial_guess).reshape(np.shape(data4gauss))

            # opt_gaussian = construct_gaussian((Lat,Lon),*gausssianfitp.x).reshape(np.shape(data4gauss))

            # plt.subplot(121)
            # plt.plot(contour_yx_deg[:,1],contour_yx_deg[:,0])
            # plt.contourf(Lon,Lat,data4gauss,vmin=-0.1,vmax=0.1,cmap=cm.cm.balance)
            # plt.contour(Lon,Lat,init_gaussian,colors = 'k')
            # plt.contour(Lon,Lat,opt_gaussian,vmin=-0.1,vmax=0.1,cmap=cm.cm.balance)

            # test = np.where(data4gauss==data4gauss.max())
            # plt.plot(Lon[test[0],test[1]],Lat[test[0],test[1]],'*c')

            # plt.plot(eddy_c_location[1],eddy_c_location[0],'.m')

            # plt.subplot(122)
            # plt.contourf(Lon,Lat,data4gauss-opt_gaussian,vmin=-0.1,vmax=0.1,cmap=cm.cm.balance)

            # plt.colorbar()
            # plt.savefig('./figures/tmp_{0:05}.png'.format(n_contour))
            # plt.close()

        #return identified_eddies_level #contours, contours_rossby

    @check_multiple_levels
    def _scan_eddy_multiple_levels(self,levels,polarity='both',fit_gaussian=True,geo=True):
        """
        
        """
        #Counter of level analysis
        n_level_counter=0
        for level in self.levels:
            
            self._scan_eddy_single_level(level, polarity,geo=True)
            self._select_optimal_fit(n_level_counter)
            
            n_level_counter+=1

    
    def _select_optimal_fit(self,n_level_counter):
        if n_level_counter == 0:
            self.multilevel_identified_eddies = self.identified_eddies_level
        else:
            loc_residual_gausian = self.multilevel_identified_eddies[['x','y','residual']]
            

        
        

    def _scan_eddy_in_time(self):
        pass


    # TODO: export same data structure with fitted gaussian and without.
    def single_level_dict(self,**kwargs):
        return pd.DataFrame(kwargs,[kwargs['n']])

    def filter_data_spatially(self):

        if self.spatial_filter['kernel']%2 == 0 :
            ker=np.ones((self.spatial_filter['kernel']+1,self.spatial_filter['kernel']+1))
        else :
            ker=np.ones((self.spatial_filter['kernel'],self.spatial_filter['kernel']))
            
        self.DataArray = self.DataArray - convolution.convolve(self.DataArray , kernel = ker, preserve_nan=True)

    
def construct_ellipse(ellipse_dict):
    yx = ellipse_dict['generator'].predict_xy(theta_r) + contour_yx_deg.mean(axis=0)


def find_centre_and_max(contour_yx_deg,data_inside_contour,center):
    """
    
    """
    centre_contour = np.nanmean(contour_yx_deg,axis=0)
    sign = np.sign(data_inside_contour[center[0],center[1]]).values

    data_max = filters.maximum_filter(sign*data_inside_contour, 3)
    maxima_indx = np.where(sign*data_inside_contour == data_max)

    maxima_loc = np.array([ 
        data_inside_contour.Y.isel(Y=maxima_indx[0]).values, 
        data_inside_contour.X.isel(X=maxima_indx[1]).values 
        ], dtype = float)
    # If multiple local maxima, then select the one with the smalles 
    # distance to the centre_contour
    if np.shape(maxima_indx)[1] != 1:
        dist = np.sqrt(abs(np.sum(centre_contour[:,np.newaxis] - maxima_loc)))
        index = np.argmin(dist)
        data_max = data_inside_contour.isel({'X':maxima_indx[1][index],'Y':maxima_indx[0][index]})
    else:
        data_max = ((sign*data_inside_contour).max())*sign

    return np.array([ *centre_contour, data_max.values, sign], dtype = float)

#plt.plot(contours_xy_deg[:,1],contours_xy_deg[:,0],'--g')
#plt.savefig('./tmp_{0:05}.png'.format(n_contour))
#plt.close()

# plt.plot(yx[:,1],yx[:,0])
# plt.plot(contours_xy_deg[:,1],contours_xy_deg[:,0],'--g')
# plt.show()

def get_centre_location(coords,contours,halo = 3):
    """
    """
    # Find max and min of coordinates
    yidmin, yidmax = find2l(coords['Y'],coords['Y'],contours[:,0].min(),contours[:,0].max())
    xidmin, xidmax = find2l(coords['X'],coords['X'],contours[:,1].min(),contours[:,1].max())
    # Find centre of contour. 
    cmindex = find(contours[:,1], contours[:,1].max())
    # Locate position of top centre of contour.
    ymindex,xmindex = find2l(coords['Y'],coords['X'],contours[cmindex,0],contours[cmindex,1])
    # Index location of centre contour.
    centertop = [ ymindex - yidmin + halo, xmindex - xidmin + halo ]
    # Slice of coordinates
    xslice = slice(xidmin.values-halo,xidmax.values+halo)
    yslice = slice(yidmin.values-halo,yidmax.values+halo)
    return xslice, yslice, centertop

def ellipse_fit_leastsq(contour):
    """
    """
    # Fit ellipse using skimage.measure
    ellipse_m = EllipseModel()
    mean_loc_contour = contour.mean(axis=0)
    if ellipse_m.estimate( contour - mean_loc_contour ):
        # Invert b and a due to ji convention
        xc, yc, b, a, theta = ellipse_m.params
        residual = ellipse_m.residuals( contour - mean_loc_contour )
        return {'x0' : xc, 'y0' : yc, 'a' : a, 'b': b, 'theta': theta, 'generator': ellipse_m, 'residual': residual}

def is_contour_close(contour):
    """
    """
    return contour[0,0] == contour[-1,0] and contour[0,1] == contour[-1,1]

def rossby_area(contour):
    """
    """
    # Flip contour from y,x to x,y
    rossby_radius = rossbyR(*np.flipud(contour.mean(axis=0)))
    rossby_area = np.pi*(rossby_radius**2)
    return rossby_area


def poly_area(contour_ji,coords_range,shape,geo=True):
    """
    """
    if geo:
        verts= contour_ji_2_meters(contour_ji,coords_range,shape)
    else: 
        verts= contour_ji_2_yx(contour_ji,coords_range,shape)
    verts_roll = np.roll(verts, 1, axis=0)
    area_elements = ((verts_roll[:,0] + verts[:,0]) *
                     (verts_roll[:,1] - verts[:,1]))
    return abs(area_elements.sum())/2.0


def contour_ji_2_meters(contour_ji,coords_range,shape):
    """
    """
    # Convert to meters, to compute area of polygon.
    yx_arc = np.pi*earth_radius*coords_range['range']/180
    yx_arc_min = np.pi*earth_radius*coords_range['mins']/180
    return ((contour_ji/shape)*yx_arc) + yx_arc_min

def extract_contours(data,level,coords):
    """
    """
    if isinstance(level,int) or isinstance(level,float):
        contours_ji = find_contours(data.values, level=level)
    elif len(level)==2:
        contours_ji = np.hstack((find_contours(data.values, level=level[0]),find_contours(data.values, level=level[1])))
    else:
        raise ValueError('Levels must be an integer, float, or array (len(level)==2)')
    return contours_ji

def contour_ji_2_yx(contour_ji,coords_range,shape,geo=True):
    """
    Convert to yx coordinates
    """
    return ((contour_ji/(np.array(shape,dtype=int)-1))*coords_range['range']) + coords_range['mins']

def coordinates_range(coords):
    """
    """
    xmin = min(coords['X']).values
    ymin = min(coords['Y']).values
    x_range = ( max(coords['X']) - xmin).values
    y_range = ( max(coords['Y']) - ymin).values
    
    return {'range':np.array([y_range,x_range],dtype=float),'mins':np.array([ymin,xmin],dtype=float)}

def construct_gaussian(coords,amplitude,xo,yo, sigma_x, sigma_y, theta, offset=0):
    '''
    *************** construct_gaussian *******************
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
    x=np.asarray(coords[1],dtype=np.float64)
    y=np.asarray(coords[0],dtype=np.float64)
 
    cos_phi = np.cos(theta)
    sin_phi = np.sin(theta)
    a = (cos_phi**2)/(2*sigma_x**2) + (sin_phi**2)/(2*sigma_y**2)
    b = (np.sin(2*theta))/(4*sigma_x**2) - (np.sin(2*theta))/(4*sigma_y**2)
    c = (sin_phi**2)/(2*sigma_x**2) + (cos_phi**2)/(2*sigma_y**2)
    g = amplitude * np.exp(-(a*(x-xo)**2 + 2*b*(x-xo)*(y-yo) + c*(y-yo)**2))

    return g.ravel()

def _gaussian_compute_residual(gaussian_params, coords, original_data):
    # Construct data.
    fitted_gaussian = construct_gaussian(coords,*gaussian_params).reshape(np.shape(original_data))
    # Cost function.
    residual = np.exp(np.nanmean(abs(original_data - fitted_gaussian))) - 1
    return residual


def minimize_surface_fitting(data,values,initial_guess='',etol=1e-2):
    '''
    *************** minimize_surface_fitting *******************
    Fit a surface to the data.
    Notes:
        
    Args:
        
    Returns:
        
    Usage:
        
    '''
    Lon, Lat = np.meshgrid(values[1], values[0])
    coords=(Lat,Lon)
    
    if initial_guess=='':
        initial_guess = [1,1,1,0,0,0]   

    #TODO constrains of initial guess i.e. abs(initial_guess[0])
    error = 1
    n = 0
    # Initial bound percentage for a and b. 
    a_b_percent = 4 
    # Decrease bounds until we get a result smaller than etol.
    while error > etol:    
        # Construct bounds
        gaussian_bounds = (range_bounds(initial_guess[0],0.01), # Amplitud bounds
                        range_bounds(initial_guess[1],0.01), # x0 bounds
                        range_bounds(initial_guess[2],0.01), # y0 bounds
                        [1e-5 if ii <= 0 else ii if ii < np.pi/2.1 else np.sign(ii)*np.pi/4 for ii in range_bounds(initial_guess[3],a_b_percent)], # a bounds between 1e-5 and the upper bound
                        [1e-5 if ii <= 0 else ii if ii < np.pi/2.1 else np.sign(ii)*np.pi/4 for ii in range_bounds(initial_guess[4],a_b_percent)], # b bounds between 1e-5 and the upper bound
                        range_bounds(initial_guess[5],1)) # theta bounds
        # Construct constrains
        # gaussian_constraints = ({'type':'eq', 'fun': positive_constrain})
        # Minimize function
        minimized_function = minimize(_gaussian_compute_residual, initial_guess,
                                args=(coords,data), method='SLSQP',
                                tol = 1e-5,
                                options={'maxiter': 100, 'disp': False},
                                bounds = gaussian_bounds)#,
                                #constraints = gaussian_constraints)
        # Double percentage
        a_b_percent = a_b_percent/2
        
        
        # Stop while if number of iterations is larger than 4.
        if  minimized_function.fun < error and n < 4:
            gausssianfitp = copy.deepcopy(minimized_function)
        elif n >= 4 or minimized_function.fun >= gausssianfitp.fun:
            break
        else:
            gausssianfitp = copy.deepcopy(minimized_function)
        
        n+=1
        error = minimized_function.fun
        #print(n,error,minimized_function.nfev)
    return gausssianfitp

def positive_constrain(x):
    return float(x[3]>=0) + float(x[4]>=0)

def range_bounds(value,percent=0.1):
    # If value is smaller than 0.1, then explore range of values dependent to 
    # one order of magnitude larger than its current value.
    if np.log10(abs(value)) < -1 and np.isfinite(np.log10(abs(value))):
        return (value - 10**np.log10(abs(value)) * percent, value + 10**np.log10(abs(value)) * percent)
    # if its exactly zero then return a range of 0.01.
    elif not np.isfinite(np.log10(abs(value))): 
        return (-0.01,0.01)
    # else return a range within a given percent of the value.
    else:
        return (value - abs(value) * percent, value + abs(value) * percent)


def extract_inside_contour(data,center,sign,levels,mask=False,maskopt=None,diagnostics=False):
    '''
    
    '''
    if type(diagnostics) != list:
        diagnostics=[diagnostics]
    data_rmove=np.zeros(np.shape(data))

    data_rmove[sign * data > levels]=1
        
    markers,features=ndimage.label(data_rmove)

    if markers.max()!=1 and maskopt==None:
        markers=markers*0
        returnmasked=True

    elif markers.max()!=1 and maskopt=='maskomax': 
        for ii in range(features):
            if ii != markers[center[0],center[1]-1]:
                markers[markers==ii]==1
        #markers=markers-markers[center[0],center[1]-1]
        returnmasked=True

    elif markers.max() != 1 and (maskopt == 'contour' or maskopt == 'forcefit'):
        if center[1] == np.shape(markers)[1]:
            markers[markers!=markers[center[0],center[1]-1]] = 0
        elif center[0] == np.shape(markers)[0]:
            markers[markers!=markers[center[0]-1,center[1]]] = 0
        else:
            markers[markers!=markers[center[0],center[1]]] = 0
        markers=markers.max()-markers
        returnmasked=True 
            
    elif maskopt=='contour' or maskopt=='forcefit':
        markers=1-markers
        returnmasked=True
    else:
        markers=markers*0  
        returnmasked=True
    
    if levels < 0 and maskopt!='forcefit':
        markers[data>0]=1
    elif levels > 0 and maskopt!='forcefit':
        markers[data<0]=1
    elif levels < 0 and maskopt=='forcefit':
        markers[ np.multiply( data < levels/2, data > levels-levels/4) ]=0
        data[ np.multiply( data < levels/2, data > levels-levels/4) ]=0
    elif levels > 0 and maskopt=='forcefit':
        markers[ np.multiply(data > levels/2, data < levels-levels/4) ]=0
        data[ np.multiply(data > levels/2, data < levels-levels/4) ]=0
        
    maskeddata = ma.masked_array(data, markers)

    if mask==True:
        return maskeddata.copy(),markers 
    else:
        return maskeddata.copy()