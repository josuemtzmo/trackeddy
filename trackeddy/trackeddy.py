import numpy as np
import xarray as xr
import pandas as pd
from scipy import ndimage
from astropy import convolution
from sklearn.neighbors import BallTree
from scipy.optimize import least_squares, Bounds


from trackeddy import _cntr as cntr
from functools import wraps

# from frechetdist import frdist

from trackeddy.geometryfunc import area_latlon_polygon, eccentricity


# class eddy_identifier():
#     """
#     eddy_identifier A decorator

#     Parameters
#     ----------
#     f : _type_
#         _description_
#     """
#     def __init__(self, inline_func, *args, **kwargs):
#         self.calls = 0
#         self.df_store = pd.DataFrame({'' : []})
#         self.inline_func = inline_func

#     def __call__(self, *args, **kwargs):
#         table = self.inline_func(*args, id = self.calls, **kwargs)
#         if self.df_store.empty:
#             self.df_store = table
#         else:
#             self.df_store = pd.concat((self.df_store, table))
#         self.calls += 1
#         return self.df_store
        

def eddy_identifier(f):
    def wrapped(*args, **kwargs):
        if 'exit' in kwargs:
            wrapped.calls = 0
            wrapped.df_store = pd.DataFrame({'' : []})
            return wrapped
            
        table = f(*args, id = wrapped.calls, **kwargs)
        if wrapped.df_store.empty:
            wrapped.df_store = table
        else:
            wrapped.df_store = pd.concat((wrapped.df_store, table))
        wrapped.calls += 1
        return wrapped.df_store
    wrapped.calls = 0
    wrapped.df_store = pd.DataFrame({'' : []})
    return wrapped


class Eddy():
    def __init__(self,contour,level) -> None:
        self.contour = contour
        self.level = np.array(level)
        self.area_eddy = 0
        self.radius_eddy = 0
        self.ellipse = None
        self.ellipse_params=None
        self.contour_ellipse_error = 0
        self.ellipse_eccen = 0
        self.contour_center=None
        self.eddy_sign = 0
        self.eddy_maxima = 0
        self.gaussian_params=None
        self.contour_gaussian_error = 0

        # TODO remove these.
        self.contour_mask=None
        self.gaussian=None
        self.data_near_contour=None

    def to_dict(self,out_properties=None):

        if not out_properties:
            out_properties = ['identifier','contour','level','area_eddy','radius_eddy','ellipse_params','contour_ellipse_error','ellipse_eccen','contour_center','eddy_sign','eddy_maxima','gaussian_params','contour_gaussian_error']

        eddy_dict = {}
        for property in out_properties:
            attr = getattr(self,property)

            if not isinstance(attr,np.ndarray):
                attr = np.array(attr)
            #     raise ValueError('The convertion to dict is expecting ndarray objects, check the eddy.{0} variable.'.format(property))

            attr = attr.squeeze()
            if len(attr.shape) > 2:
                pass
                # TODO Ravel data to export matrixes. 
                # eddy_dict['contour_x']=attr[:,0]
            elif len(attr.shape) == 2 and 2 in attr.shape:
                eddy_dict['contour_path_x']=attr[:,0]
                eddy_dict['contour_path_y']=attr[:,1]
            elif not attr.shape:
                eddy_dict[property]=attr
            else:
                match property:
                    case 'ellipse_params':
                        eddy_dict[property+'_a']=attr[0]
                        eddy_dict[property+'_b']=attr[1]
                        eddy_dict[property+'_theta']=attr[2]
                    case 'contour_center':
                        eddy_dict['contour_x']=attr[0]
                        eddy_dict['contour_y']=attr[1]
                    case 'eddy_maxima':
                        eddy_dict['maxima']=attr[0]
                        eddy_dict['maxima_x']=attr[1]
                        eddy_dict['maxima_y']=attr[2]
                    case 'gaussian_params':
                        eddy_dict[property+'_x']=attr[0]
                        eddy_dict[property+'_y']=attr[1]
                        eddy_dict[property+'_sigma_x']=attr[2]
                        eddy_dict[property+'_sigma_y']=attr[3]
                        eddy_dict[property+'_theta']=attr[4]
                    case _:
                        raise ValueError("{0} doesn't exist, make sure the property is a attribute of the eddy class.".format(property))
                    
        return eddy_dict            
        
    def to_table(self,out_properties=None,**kwards):

        eddy_dict = self.to_dict(out_properties)
        table = pd.DataFrame.from_dict(eddy_dict,**kwards).reset_index()
        
        return table

    def diagnose_eddy(self):
        # Print or plot the eddy
        pass
    
    @eddy_identifier
    def store(self,id):
        self.identifier=id
        table = self.to_table()
        return table

    @eddy_identifier
    def discarded(self,reason,id,exit=False):
        # self.identifier = id
        self.reason = np.array(reason)
        out_properties = ['reason','level','area_eddy','radius_eddy','contour_ellipse_error','contour_center','ellipse_eccen','eddy_sign','eddy_maxima','contour_gaussian_error'] 
        dict = self.to_dict(out_properties)
        table = pd.DataFrame(dict,index=[id])

        #TODO concatenate all the discarded eddies.
        return table
    
    def exit(self):
        # Reset decorator counter in both functions.
        self.store(exit=True)
        self.discarded(exit=True)

        

    
def extract_contours(X, Y, data, level):
    # Convert to arrays to use the cntr library
    c = cntr.Cntr(X, Y, data)
    # Extract contours in the level  
    res = c.trace(level)
    # result is a list of arrays of vertices and path codes
    # (see docs for matplotlib.path.Path)

    nseg = len(res) // 2
    segments = res[:nseg] 
    codes = res[nseg:]

    return segments,codes


def filter_data(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        filtered = self.treat_nan(*args)
        if kwargs.get('filter') == 'both':
            filtered = func(self, filtered, 'space')
            self.data2track = func(self, filtered, 'time')
        elif  kwargs.get('filter') == 'space' or  kwargs.get('filter') == 'time':
            self.data2track = func(self, *args, kwargs.get('filter'))
        else:
            self.data2track = filtered

    return wrapper



class TrackEddy():
    def __init__(self,path, variable, nan_treatment=False, **xrkwargs) -> None:
        self.Dataset = xr.open_mfdataset(path,**xrkwargs)
        self.rawdata = self.Dataset[variable]
        self.nan_treatment = nan_treatment
        self.filter='space'


        self.identify_coords()
        self.check_coords()

        self.identification_criteria = {'ellipse_fit':0.85,'eccentricity':0.85,'gaussian_fit':0.8, 'max_area': 2*np.pi }

        self.S_filter_setup = {'filter_type':'convolution','mode':'uniform','kernel':10}

        self.T_filter_setup = {}

    def check_coords(self):
        if ( len(self.Dataset[self.coords['x']].shape) == 1 and len(self.Dataset[self.coords['y']].shape) ==1 ):
            self.X, self.Y = np.meshgrid(self.Dataset[self.coords['x']],self.Dataset[self.coords['y']])
        
        elif ( len(self.Dataset[self.coords['x']].shape) == 2 and len(self.Dataset[self.coords['y']].shape) ==2 ):
            self.X = self.Dataset[self.coords['x']]
            self.Y = self.Dataset[self.coords['y']]
        
        else: 
            raise ValueError("The coordinates of the dataset should be a 1D or 2D array, check the vars: '{0}', {1}.".format(self.coords['x'],self.coords['y']))


    # TODO add a wrapper to check if the coords where properly identified, if not rise error.
    def identify_coords(self) -> None:
        self.coords={'time':'','x':'','y':''}
        for coord_name in self.rawdata.coords.keys():
            if 'time' in coord_name or coord_name == 't':
                self.coords['time']=coord_name
            if 'lon' in coord_name or coord_name == 'x':
                self.coords['x']=coord_name
            if 'lat' in coord_name or coord_name == 'y':
                self.coords['y']=coord_name


    def setup(self) -> None:
        pass

    @filter_data
    def _filter_data_(self, data_snapshot, filter=None) ->  xr.DataArray:

        if filter=='time':
            data2track = self._time_filter(data_snapshot,**self.T_filter_setup)
        elif filter=='space':
            data2track = self._space_filter(data_snapshot,**self.S_filter_setup)
        else:
            raise ValueError("Select a filter to extract field anomaly. The current options are: 'time' and 'space'.")
        
        return data2track


    def _time_filter(self,filter):
        pass

    def _space_filter(self, data2track, filter_type, mode='uniform', kernel=10) -> xr.DataArray:
        
        # Apply spatial filter.
        if filter_type == 'convolution':
            if mode == 'uniform':
                tmp_data = data2track.squeeze().copy()
                if kernel%2 == 0:
                    ker=np.ones((kernel+1,kernel+1))
                else:
                    ker=np.ones((kernel,kernel))
                convolved_data = convolution.convolve(tmp_data, kernel = ker)
                tmp_data = tmp_data - convolved_data

                data2track = tmp_data
            if mode == 'gaussian':
                raise Warning('ndimage.gaussian_filter may create artefacts near nan values. Therefore, data is filled with zeros.')
                tmp_data = data2track.squeeze().copy()
                tmp_data = tmp_data - ndimage.gaussian_filter(tmp_data, size = kernel)
                data2track = tmp_data, mask
        # Check if the user selects an moving average meridional filter.
        elif filter_type == 'meridional':
            data2track = data2track - data2track.mean(self.coords['x'])
        # Check if the user selects an moving average zonal filter.
        elif filter_type == 'zonal':
            data2track = data2track - data2track.mean(self.coords['y'])
        else:
            raise ValueError("Define the filter_type argument between the options: 'convolution', 'meridional', and 'zonal'")
        return data2track

    def _detect_snapshot(self, time, levels):
        
        if isinstance(time, int) and self.coords['time']:
            data_snapshot = self.rawdata.isel({self.coords['time']:time}).squeeze()
        elif isinstance(time, str) and self.coords['time']:
            data_snapshot = self.rawdata.isel({self.coords['time']:time}).squeeze()
        else:
            data_snapshot = self.rawdata.squeeze()

        self._filter_data_(data_snapshot,filter=self.filter)

        joint_eddies = pd.DataFrame({'' : []})

        for level in levels:

            eddies_current, discarded = self._detect_one_level(level)
            
            # If no eddies identified, then continue
            if eddies_current.empty:
            # print(eddies_current.empty)
                continue
            
            if joint_eddies.empty : 
                joint_eddies = eddies_current
                # Continue, since no check needs to be done.
                continue
            
            # TODO pass the previous and current eddy to the _detect_nearest function
            index, nearest, d = self._detect_nearest(joint_eddies,eddies_current)

            joint_eddies, eddies_current = self._update_better_eddy(joint_eddies,eddies_current,index,nearest,d)
            #TODO function with decorator that handles the merging of eddy tracks.
            
            joint_eddies = self._merge_reindex_eddies(joint_eddies, eddies_current)

            # try:
            #     eddies_current.loc[22].plot.line(x='contour_path_x',y='contour_path_y',color='k')
            #     joint_eddies.loc[22].plot.line(x='contour_path_x',y='contour_path_y',color='r')
            #     print(level)
            # except:
            #     pass

        return joint_eddies
        # return eddies_previous, eddies_current


    def _detect_nearest(self,p_eddy,n_eddy):

        # Extract from table the first element in the index 'index' 
        eddy_info_p_level = p_eddy.xs(0,level=('index'))
        eddy_info_n_level = n_eddy.xs(0,level=('index'))
        # Extract coordinates of maxima within the contour. 
        eddy_coords_p = eddy_info_p_level[['maxima_y','maxima_x']].values
        eddy_coords_n = eddy_info_n_level[['maxima_y','maxima_x']].values
        # Extract eddy radius of the new eddies.
        eddy_radius_n = eddy_info_n_level['radius_eddy'].values
        # Search nearest points between the previous eddies and the new eddies. 
        # Important to pass ((lat, lot)). 
        ball_tree = BallTree( np.deg2rad(eddy_coords_p), metric='haversine')
        distance, index = ball_tree.query(np.deg2rad(eddy_coords_n))
        # If distance is smaller than radius of eddy, then mark them as within the eddy radius.
        distance_between_nearest = ( distance * 6371 ) - np.expand_dims(eddy_radius_n,1)
        # If the distance is negative, then the current eddy is within the previous eddy.
        nearest, _ = np.where(distance_between_nearest < 0)
        return index, nearest, ball_tree
    
    def _update_better_eddy(self,p_eddy,n_eddy,index,nearest,d):
        # Get unique identifier for each eddy identified at the previous and current level
        previous_index = p_eddy.index.get_level_values(level=0)
        current_index = n_eddy.index.get_level_values(level=0)
        # Maximum index in the previous table
        previous_eddies_count = previous_index.unique().max()
        # Shift the current index by the largest unique identifier + 1 in the previous eddy table. As if all the new eddies are unique.
        new_identifier = (current_index+previous_eddies_count+1).unique()
        new_index = n_eddy.index.set_levels(new_identifier, level=0)
        n_eddy.index = new_index

        replaced_eddies = []
        for duplicated in nearest:
            
            shifted_counter = duplicated+previous_eddies_count+1

            previous_eddy = p_eddy.xs((index[duplicated][0], 0),level=['identifier','index'])
            current_eddy = n_eddy.xs((shifted_counter, 0),level=['identifier','index'])

            p_ellipse_error = previous_eddy.contour_ellipse_error.values
            c_ellipse_error = current_eddy.contour_ellipse_error.values

            p_gauss_error = previous_eddy.contour_gaussian_error.values
            c_gauss_error = current_eddy.contour_gaussian_error.values
            
            #TODO check if this is the best condition. 
            if p_ellipse_error > c_ellipse_error and p_gauss_error > c_gauss_error:
                n_eddy = n_eddy.rename(index={shifted_counter:index[duplicated][0]},level=0)
                replaced_eddies.append(index[duplicated][0]) 
                print('better')
            else:
                # Delete the input that does a worst job compared to the previous level
                print('worst')
                n_eddy = n_eddy.drop(shifted_counter,level=0)
        
        if replaced_eddies:
            p_eddy = p_eddy.drop(replaced_eddies,level=0)

        return p_eddy, n_eddy

    def _merge_reindex_eddies(self,p_eddy,n_eddy):
        
        # Merge eddy tables
        new_eddy_table = pd.concat([p_eddy,n_eddy])
        # Get new length of unique identifiers 
        new_eddy_count = len(new_eddy_table.index.get_level_values(level=0).unique())

        # Create new monotonously increasing index
        new_identifier = np.arange(0,new_eddy_count,dtype=int)
        
        # Reset table to avoid issues with codes.
        reset_table = new_eddy_table.reset_index().set_index(['identifier','index'])

        # Create new index for the table
        new_index = reset_table.index.set_levels(new_identifier, level=0)
        # Assign the new index to the table
        reset_table.index = new_index
        
        #Important that index is sorted, so when computing distances, the values will be assigned to an monotonously increasing index.
        return reset_table.sort_index(level=['identifier','index'])

    def treat_nan(self,data_snapshot, nan_value=0):
        if not self.nan_treatment:
            data2track = data_snapshot.copy()
        elif self.nan_treatment:
            data2track = data_snapshot.fillna(nan_value).copy()
        else:
            raise ValueError("The nan_treatment can only be True or False")
        return data2track
    

    def _detect_one_level(self,level):

        contours, _ = extract_contours( self.X, self.Y, self.data2track.values, level)
        
        for contour in contours:
            
            eddy = Eddy(contour,level)

            # Brute force the removal of contours that are too large
            eddy.area_eddy = area_latlon_polygon(eddy.contour)
            eddy.radius_eddy = np.array(np.sqrt(eddy.area_eddy/np.pi))

            eddy.contour_center = np.mean(eddy.contour,0)

            # Continue to next contour if eddy size is too large.
            # TODO add criteria selected by user here.
            if eddy.radius_eddy > 200:
                discarded = eddy.discarded('area')
                continue
                
            
            # TODO clean the fit_ellipse function and input only eddy.contour
            eddy.ellipse, eddy.ellipse_params = fit_ellipse(eddy.contour[:,0],eddy.contour[:,1])

            # Check if area is similar between fitted ellipse and contour
            # Discarding by area first improves the time. 
            area_ellipse = area_latlon_polygon(eddy.ellipse)
            if not np.isclose(eddy.area_eddy,area_ellipse, rtol=1-self.identification_criteria['ellipse_fit']): 
            # Continue to next contour if ellipse fit is bad
                discarded = eddy.discarded('ellipse_area')
                continue
            
            # Check eccentricity of ellipse.
            eddy.ellipse_eccen = eccentricity(*eddy.ellipse_params[0:2])
            if eddy.ellipse_eccen > self.identification_criteria['eccentricity']:
                discarded = eddy.discarded(reason='eccentricity')
                continue

            eddy.contour_ellipse_error = compute_similarity(eddy.contour,eddy.ellipse)
            
            # Check if the similarity between the contours matches the user criteria.
            if eddy.contour_ellipse_error > self.identification_criteria['ellipse_fit']:
                discarded = eddy.discarded(reason='similarity')
                continue

            # TODO Extract the area from the contour.
            data_near_contour  = self._data_in_contour(eddy)
            
            if eddy.eddy_maxima[0] == 0:
                discarded = eddy.discarded(reason='eddy_is_island')
                continue

            # TODO Extract center of mass of contour.

            #TODO FIT GAUSSIAN USING THE eddy_sign
            X, Y = np.meshgrid(data_near_contour[self.coords['x']].values, data_near_contour[self.coords['y']].values)

            # Create object to handle fitting of surfaces.
            F_surface = Fit_Surface(eddy, data_near_contour, X, Y)
            
            # Fit the chosen feature.
            eddy.gaussian, eddy.gaussian_params = F_surface._fitting()

            # Extract the contour of the gaussian fitted + the eddy contour, since all the gaussians decay to zero, for a fair comparison.
            # print(np.where(eddy.gaussian < eddy.level, eddy.contour_mask,1))
            # tmp_mask = eddy.contour_mask * np.where(eddy.gaussian < eddy.level, eddy.contour_mask,1)
            # print(tmp_mask)
            # masked_gauss = np.where(tmp_mask, eddy.gaussian, eddy.level)
            # eddy.gaussian[ ~ tmp_mask] = eddy.level
            # print(masked_gauss)
            #TODO check when eddy is negative in positive level, it may fail, but they may be identified when tracking in negative levels.
            gauss_contour, _ = extract_contours(X, Y, eddy.gaussian, eddy.level)

            # No contour was extracted. This is an issue with the gaussian fitting optimization.
            if not gauss_contour:
                discarded = eddy.discarded(reason='gaussian_fit_failed')
                continue

            # TODO Move to a decorator of extract_contours.
            # Fix in case the mask cuts the contour
            if len(gauss_contour) > 1 : 
                gauss_contour = np.vstack(gauss_contour)
            else: 
                gauss_contour=gauss_contour[0]
            
            # Compute similarity between the eddy contour and the gaussian contour at the same level.
            eddy.contour_gaussian_error = compute_similarity(eddy.contour,gauss_contour)

            # If the similarity between contours doesn't match criteria the eddy is descarde.
            if eddy.contour_gaussian_error < self.identification_criteria['gaussian_fit']:
                discarded = eddy.discarded(reason='gaussian_check')
                continue
            
            # Copy identifier into the eddy object.            
            eddy_info = eddy.store()


        #TODO move this to a decorator.
        # In case no eddies are identified
        try:
            eddy_info = eddy_info.set_index(['identifier','index'])
        except UnboundLocalError:
            eddy_info = pd.DataFrame({'' : []})
        try:
            eddy.exit()
        except UnboundLocalError:
            discarded=[]
        
        return eddy_info, discarded
    


    def _data_in_contour(self,eddy):

        pt = np.expand_dims(np.mean(eddy.contour,axis=0),0)

        # To speed up the process of finding the nearest point, we mask values near the point that are close by 10% of its value, with a largest cap of 3 degrees if coordinates are geographical. The atol fixes the issue in x or y coordinates close to 0.
        Xmask = np.isclose(self.X,pt[0][0],rtol=0.1, atol=3)
        Ymask = np.isclose(self.Y,pt[0][1],rtol=0.1, atol=3)

        sift_x = np.argmax(np.sum(Xmask*Ymask,0))
        sift_y = np.argmax(np.sum(Xmask*Ymask,1))

        mask_shape = np.unique(np.sum(Xmask*Ymask,0))[1], np.unique(np.sum(Xmask*Ymask,1))[1]

        x = self.X[Xmask*Ymask]
        y = self.Y[Xmask*Ymask]

        # Another option is to search through all the X and Y space, but it's 100 times slower to compute. 
        # TODO add an option to bruteforce in case it's needed, but I can't see a case right now. 
        # x = self.X.ravel()
        # y = self.Y.ravel()

        coords = np.vstack((x,y)).T

        # This allows support for non structured grids, but it's expensive to compute.
        ball_tree = BallTree(coords)
        dist, ind  = ball_tree.query(eddy.contour, k=1)

        # Extract coordinates from the mask shape
        coord_slice = np.hstack( np.unravel_index(ind, mask_shape) ).squeeze()

        # Shift back to the original coordinates
        coord_ind = coord_slice + [sift_y,sift_x]
        
        eddy.data_near_contour,contour_coords_grid = self._get_data_in_contour( coord_ind )

        # data_near_contour is loaded here, this should speed up the gaussian fit, since we don't need to reload the data
        masked_data_near_contour, eddy_sign, eddy_maxima, eddy_contour_mask = self._mask_data_in_contour(contour_coords_grid,eddy.data_near_contour,eddy.level)

        eddy.eddy_sign = eddy_sign
        eddy.eddy_maxima = eddy_maxima
        eddy.contour_mask = eddy_contour_mask
        return masked_data_near_contour

    def _get_data_in_contour(self, coord_ind,threshold = 'auto'):
        
        # Get corners
        TR_corner = np.max(coord_ind,axis=0)
        BL_corner = np.min(coord_ind,axis=0)
        if threshold == 'auto':
            diag_gridpoints= TR_corner - BL_corner
            thresh_value = int(0.5*np.max(diag_gridpoints))
            if thresh_value <=3:
                thresh_value = 3
        elif isinstance(threshold, int):
            thresh_value = threshold
        else: 
            thresh_value = 0

        TR_corner = TR_corner + [ thresh_value+1, thresh_value+1]
        BL_corner = BL_corner - [ thresh_value, thresh_value ]

        #Make sure that corners are never larger or smaller than the dataset dimensions
        while (TR_corner > self.X.shape).any():
            index_max = np.argmax(TR_corner - self.X.shape)
            TR_corner[index_max]=TR_corner[index_max]-1
        
        while (BL_corner < (0,0)).any():
            index_min = np.argmin(BL_corner - self.X.shape)
            BL_corner[index_min] = 0

        data_near_contour = self.data2track.isel({self.coords['x']: slice(BL_corner[1],TR_corner[1]),self.coords['y']: slice(BL_corner[0],TR_corner[0])})
        
        contour_coords_grid = coord_ind - BL_corner

        return data_near_contour, contour_coords_grid
    
    def _mask_data_in_contour(self,coord_ind_contour, data_near_contour,level,mask='contour'):
        
        #TODO Define if eddy is cyclonic or anticyclonic depending on the values inside the contour or not.

        inside_contour =  np.zeros(data_near_contour.shape)
        inside_contour[coord_ind_contour[:,0],coord_ind_contour[:,1]]=1
        inside_contour = ndimage.binary_fill_holes(inside_contour)

        data_inside_contour = data_near_contour.where(inside_contour,np.nan).load()

        # TODO is there a faster way to check if the sign of the contour interior.
        mean_value_inside_contour = data_inside_contour.mean()

        # This is used to differentiate between cyclonic and anticyclonic eddies
        eddy_sign = np.sign(mean_value_inside_contour).values
        # This computation allow to extract the eddy maxima independent to their sign.
        eddy_maxima = eddy_sign * ( eddy_sign * data_inside_contour ).max()

        eddy_argmaxima = ( eddy_sign * data_inside_contour ).argmax()
        eddy_idxmaxima = np.unravel_index(eddy_argmaxima, data_inside_contour.shape)
        
        eddy_xlocmax = data_near_contour[self.coords['x']].isel({self.coords['x']: eddy_idxmaxima[1]})
        eddy_ylocmax = data_near_contour[self.coords['y']].isel({self.coords['y']: eddy_idxmaxima[0]})

        eddy_maxima_loc = np.vstack((eddy_maxima, eddy_xlocmax, eddy_ylocmax))

        if mask == 'forcefit':
            # TODO reimplement the forcefit option, make a masked ring and then set to zero.
            pass

        # After some testing, leaving the nans outside the contour in data_inside_contour allows for better fitting of the gaussian and eddy detection.
        return data_inside_contour, eddy_sign, eddy_maxima_loc, inside_contour
    
    def plot_eddy_detection_multilevel(self,df_eddy_multilevel_store):
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt
        import cmocean as cm

        n_colors =list(mcolors.XKCD_COLORS.keys())

        levels = df_eddy_multilevel_store.level.unique()
        color_index = np.linspace(0,len(n_colors)-1,len(levels),dtype=int)

        plt.figure(figsize=(10,6),dpi=300)
        plt.pcolormesh(self.X, self.Y, self.data2track,cmap=cm.cm.balance,vmin=-0.5,vmax=0.5)

        for eddy in df_eddy_multilevel_store.index.get_level_values(level=0):
            c_path_x = df_eddy_multilevel_store.loc[eddy].contour_path_x
            c_path_y = df_eddy_multilevel_store.loc[eddy].contour_path_y

            inx = np.where(df_eddy_multilevel_store.loc[eddy].level[0]==levels)[0][0]
            
            plt.plot(c_path_x,c_path_y,color=n_colors[color_index[inx]])
        
        if len(levels) < 5:
            for level in range(len(levels)):
                plt.plot(c_path_x[0],c_path_y[0],color=n_colors[level],label=levels[level])
            plt.legend()


    def plot_eddy_detection_in_level(self, eddies, discarded, plot_args):
        if not plot_args:
            plot_args = {'alpha':0.5,'markersize':1}

        import matplotlib.pyplot as plt
        import cmocean as cm
        fig,ax = plt.subplots(1,2,figsize=(10,5),dpi=300)

        ax[0].pcolormesh(self.X, self.Y, self.data2track.squeeze().values, cmap = cm.cm.balance, vmin=-0.5, vmax=0.5)

        for idx in eddies.index.get_level_values(0).unique():
            Xcontour = eddies.loc[idx]['contour_path_x']
            Ycontour = eddies.loc[idx]['contour_path_y']

            x_cont = eddies.loc[idx]['contour_x']
            y_cont = eddies.loc[idx]['contour_y']

            ax[0].plot(Xcontour,Ycontour,'-m')
            ax[0].plot(x_cont,y_cont,'.k',markersize=1)

        ax[1].pcolormesh(self.X,self.Y , self.data2track.squeeze().values,cmap = cm.cm.balance,vmin=-0.5,vmax=0.5)

        for idx in discarded.index:
            x_cont = discarded.loc[idx]['contour_x']
            y_cont = discarded.loc[idx]['contour_y']
            
            match discarded.loc[idx].reason:
                case 'eccentricity':
                    ax[1].plot(x_cont,y_cont,'.r',**plot_args)
                case 'area':
                    ax[1].plot(x_cont,y_cont,'.',color='gray',**plot_args)
                case 'similarity':
                    ax[1].plot(x_cont,y_cont,'.c',**plot_args)
                case 'ellipse_area':
                    ax[1].plot(x_cont,y_cont,'.m',**plot_args)
                case 'eddy_is_island':
                    ax[1].plot(x_cont,y_cont,'.',color='orange',**plot_args)
                case 'gaussian_fit_failed':
                    ax[1].plot(x_cont,y_cont,'.',color='pink',**plot_args)
                case _:
                    ax[1].plot(x_cont,y_cont,'.b',**plot_args)

        label = ['Fail eccentricity','Area fail','Contour != ellipse', 'Ellipse area fail', 'Eddy is island', 'Gaussian fit failed','Contour != Gaussian']

        colors=['r','gray','c','m','orange','pink','b']
        [ax[1].plot(self.X[0,0],self.Y[0,0], '.',color=colors[c], label = label[c], alpha=0.5 )for c in range(len(colors))]

        ax[0].set_title('Identified')
        ax[1].set_title('Discarded')

        ax[1].legend(loc='upper center',  bbox_to_anchor=(0, -0.05),
          ncol=3)

    
class Fit_Surface():

    def __init__(self,eddy,eddy_data,X,Y,mode='gaussian') -> None:

        self.X = X
        self.Y = Y
        self.eddy_max = eddy.eddy_maxima[0]
        self.x_eddy = eddy.eddy_maxima[1]
        self.y_eddy = eddy.eddy_maxima[2]
        self.level = eddy.level
        
        # Initial guess to start fitting feature
        self.initial_guess = np.hstack((self.x_eddy,self.y_eddy,eddy.ellipse_params))
        #TODO is the best to pass the coords from the eddy_data or should I provide the self.X and self.Y cropped? In fact, it may not matter. 
        
        self.data = eddy_data - eddy.level
        self.mode = mode

        
    def _fitting(self):

        coords=(self.X,self.Y,self.eddy_max)

        fitdict = self.fit_curve(coords)
        
        fitted_curve = gaussian(coords,  *fitdict)
        fitted_data=fitted_curve.reshape(*self.X.shape)

        return fitted_data, fitdict
    
    def construct_bounds(self):

        xdelta = np.mean(np.diff(self.X,axis=1))
        ydelta = np.mean(np.diff(self.Y,axis=0))

        LCoordbounds = np.array([[ self.initial_guess[0] - xdelta ], [self.initial_guess[1] - ydelta ]])
        UCoordbounds = np.array([[ self.initial_guess[0] + xdelta ], [self.initial_guess[1] + ydelta ]])

        # Limit the bounds for the initial guess to 50% of their ellipse value. 

        Lellipsebounds = np.array([[ init-0.5*init ] for init in self.initial_guess[2:]])
        Uellipsebounds = np.array([[ init+0.5*init ] for init in self.initial_guess[2:]])

        Lbounds = np.vstack((LCoordbounds,Lellipsebounds))
        Ubounds = np.vstack((UCoordbounds,Uellipsebounds))

        bounds = np.hstack((Lbounds,Ubounds)).T
        bounds.sort(axis=0)

        return bounds
        

    def fit_curve(self, coords):

        # res = minimize(gaussian_residual, self.initial_guess, args=(coords,self.data),method='SLSQP')
        bounds = self.construct_bounds()
        
        res = least_squares(gaussian_residual, self.initial_guess, args=(coords,self.data), bounds = bounds)
        
        fitdict = np.array(res.x)
        
        return fitdict
        
def gaussian(coords, xo, yo, sigma_x, sigma_y, theta):  
    '''
    *************** gaussian *******************
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

    cos_phi = np.cos(theta)
    sin_phi = np.sin(theta)
    a = (cos_phi**2)/(2*sigma_x**2) + (sin_phi**2)/(2*sigma_y**2)
    b = (np.sin(2*theta))/(4*sigma_x**2) - (np.sin(2*theta))/(4*sigma_y**2)
    c = (sin_phi**2)/(2*sigma_x**2) + (cos_phi**2)/(2*sigma_y**2)
    g = amplitude*np.exp(-(a*(x-xo)**2 + 2*b*(x-xo)*(y-yo) + c*(y-yo)**2))
    
    return g.ravel()

def gaussian_residual(popt, coords, data2fit):

    gauss = gaussian(coords,*popt).reshape(np.shape(data2fit))
    
    # residual = abs(data2fit - gauss)
    
    # residual = np.nansum(abs(data2fit - gauss))
    # residual = np.nanmean(abs(data2fit - gauss))
    #This seems to outperform the exp()
    residual = np.nanstd(abs(data2fit - gauss)) * np.nanmean(abs(data2fit - gauss))

    # residual = np.exp(np.float64(np.nanmean(abs(data2fit - gauss)))) - 1
    return residual


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
    
    xp = x-mean_x
    yp = y-mean_y
    
    # D1 = np.vstack([xp**2, xp*yp, yp**2]).T
    # D2 = np.vstack([xp, yp, np.ones(len(xp))]).T
    # S1 = D1.T @ D1
    # S2 = D1.T @ D2
    # S3 = D2.T @ D2
    # T = -np.linalg.solve(S3,S2.T)
    # M = S1 + S2 @ T
    # C = np.array(((0, 0, 2), (0, -1, 0), (2, 0, 0)), dtype=float)
    # M = np.linalg.solve(C, M)
    # eigval, eigvec = np.linalg.eig(M)
    # con = 4 * eigvec[0] * eigvec[2] - eigvec[1]**2
    # ak = eigvec[:, np.nonzero(con > 0)[0]]
    
    # H2 = np.dot(eigvec,np.diag(eigval),np.linalg.inv(eigvec))
    # Error = np.linalg.norm(M - H2)

    # # print(Error)
    # # print(np.concatenate((ak, T @ ak)).ravel())
    # try:
    #     a,b,c,d,e,f = np.concatenate((ak, T @ ak)).ravel()
    # except:
    #     a,b,c,d,e,f = np.zeros(6)

    X = np.array([xp**2,xp*yp,yp**2,xp,yp]).T
    
    a = np.sum(X,axis=0)
    b = np.dot(X.T,X)
    
    x2  = np.linalg.solve(b.T, a.T)
    a,b,c,d,e=x2

    if b == 0 and a<c:
        anglexaxis_rad=0
    elif b == 0 and c<a:
        anglexaxis_rad=np.pi/2
    else:
        anglexaxis_rad = np.arctan((c-a-np.sqrt((a-c)**2+b**2))/b)

    if ( min(abs(b/a),abs(b/c)) > orientation_tolerance ):
        # TODO: Replace this non sign definite orientation_rad for anglexaxis 
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
    
    theta_r         = np.linspace(0,2*np.pi,len(y));
    ellipse_x_r     = X0 + a*np.cos(theta_r)
    ellipse_y_r     = Y0 + b*np.sin(theta_r)
    rotated_ellipse = np.dot(R, np.array([ellipse_x_r,ellipse_y_r]))

    # Ensure that a is always larger than b
    b, a = np.sort([a,b])   
    
    ellipse_params = np.array((a,b,anglexaxis_rad))

    return rotated_ellipse.T, ellipse_params


import math 
def arc_length(curve):
    '''
    Args:
    points: type arrays two values [[x, y], [x, y]]
    Returns:
    acc_length: curve length
    Descriptions:
    Calculate the length of the curve
    '''

    acc_length = 0
    for i in range(0, len(curve)-1):
        acc_length += math.dist(curve[i], curve[i+1])

    return acc_length


def compute_similarity(curve1,curve2):

    geo_avg_curve_len = math.sqrt(
    arc_length(curve1) *  arc_length(curve2))

    freshet_dist = frdist(curve1,curve2)

    result = max(1 - freshet_dist / (geo_avg_curve_len / math.sqrt(2)), 0)
    return np.array(round(result, 4))

def _c(ca, i, j, p, q):

    if ca[i, j] > -1:
        return ca[i, j]
    elif i == 0 and j == 0:
        ca[i, j] = np.linalg.norm(p[i]-q[j])
    elif i > 0 and j == 0:
        ca[i, j] = max(_c(ca, i-1, 0, p, q), np.linalg.norm(p[i]-q[j]))
    elif i == 0 and j > 0:
        ca[i, j] = max(_c(ca, 0, j-1, p, q), np.linalg.norm(p[i]-q[j]))
    elif i > 0 and j > 0:
        ca[i, j] = max(
            min(
                _c(ca, i-1, j, p, q),
                _c(ca, i-1, j-1, p, q),
                _c(ca, i, j-1, p, q)
            ),
            np.linalg.norm(p[i]-q[j])
            )
    else:
        ca[i, j] = float('inf')

    return ca[i, j]


def frdist(P, Q):
    """
    Computes the discrete Fréchet distance between
    two curves. The Fréchet distance between two curves in a
    metric space is a measure of the similarity between the curves.
    The discrete Fréchet distance may be used for approximately computing
    the Fréchet distance between two arbitrary curves,
    as an alternative to using the exact Fréchet distance between a polygonal
    approximation of the curves or an approximation of this value.

    This is a Python 3.* implementation of the algorithm produced
    in Eiter, T. and Mannila, H., 1994. Computing discrete Fréchet distance.
    Tech. Report CD-TR 94/64, Information Systems Department, Technical
    University of Vienna.
    http://www.kr.tuwien.ac.at/staff/eiter/et-archive/cdtr9464.pdf

    Function dF(P, Q): real;
        input: polygonal curves P = (u1, . . . , up) and Q = (v1, . . . , vq).
        return: δdF (P, Q)
        ca : array [1..p, 1..q] of real;
        function c(i, j): real;
            begin
                if ca(i, j) > −1 then return ca(i, j)
                elsif i = 1 and j = 1 then ca(i, j) := d(u1, v1)
                elsif i > 1 and j = 1 then ca(i, j) := max{ c(i − 1, 1), d(ui, v1) }
                elsif i = 1 and j > 1 then ca(i, j) := max{ c(1, j − 1), d(u1, vj) }
                elsif i > 1 and j > 1 then ca(i, j) :=
                max{ min(c(i − 1, j), c(i − 1, j − 1), c(i, j − 1)), d(ui, vj ) }
                else ca(i, j) = ∞
                return ca(i, j);
            end; /* function c */

        begin
            for i = 1 to p do for j = 1 to q do ca(i, j) := −1.0;
            return c(p, q);
        end.

    Parameters
    ----------
    P : Input curve - two dimensional array of points
    Q : Input curve - two dimensional array of points

    Returns
    -------
    dist: float64
        The discrete Fréchet distance between curves `P` and `Q`.

    Examples
    --------
    >>> from frechetdist import frdist
    >>> P=[[1,1], [2,1], [2,2]]
    >>> Q=[[2,2], [0,1], [2,4]]
    >>> frdist(P,Q)
    >>> 2.0
    >>> P=[[1,1], [2,1], [2,2]]
    >>> Q=[[1,1], [2,1], [2,2]]
    >>> frdist(P,Q)
    >>> 0
    """

    p = P if len(P) >= len(Q) else Q
    q = Q if len(P) >= len(Q) else P

    p = np.array(p, np.float64)
    q = np.array(q, np.float64)

    len_p = len(p)
    len_q = len(q)

    # if len_p == 0 or len_q == 0:
    #     raise ValueError('Input curves are empty.')

    # if len_p != len_q or len(p[0]) != len(q[0]):
    #     raise ValueError('Input curves do not have the same dimensions.')



    ca = (np.ones((len_p, len_q), dtype=np.float64) * -1)

    dist = _c(ca, len_p-1, len_q-1, p, q)
    return dist