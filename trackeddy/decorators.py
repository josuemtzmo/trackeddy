import functools
import numpy as np
import xarray as xr
import inspect
import warnings

def check_input(init): 
    """
    Decorator to check init input of class TrackEddy
    """ 
    arg_names = inspect.getfullargspec(init)[0]
    
    @functools.wraps(init)
    def test_input(self, *args, **kwargs):  
        for name, value in zip(arg_names[1:], args):
            setattr(self, name, value)
        
        # Check dataset time:
        input_type = type(self.dataset)
        # Load dataset if string is provided
        if input_type ==str:
            if '*' in self.dataset:
                dataset = xr.open_mfdataset(self.dataset)
            else:
                dataset = xr.open_dataset(self.dataset)
            var_names = dataset.data_vars.keys()
        # Is input a Dataset 
        elif input_type == xr.core.dataset.Dataset:
            var_names = dataset.data_vars.keys()
        # Is input a Dataarray
        elif input_type == xr.core.dataset.DataArray:
            # Convert to dataset
            dataset = dataset.to_dataset(name = dataset.name)
        else:
            raise ValueError("Argument 'dataset' must be a string, xr.Dataset or xr.DataArray")

        # Check dimensions
        dims = dataset.dims

        # Dataset must be 2D or 3D.
        if len(dims) == 1:
            raise ValueError("xarrayMannKendall requires at least a 2D dataarray (x,y)")
        elif len(dims) > 3:
            raise ValueError("Currently xarrayMannKendall only supports 2D (x,t) and 3D dataarray (x,y,t)")
        
        # Rename coordinates to X,Y,T
        rename_dims = rename_dict(dims)
        transpose_order = sorted(rename_dims.values())[::-1] # i.e. [Y,X,T]
        dataset = dataset.rename(rename_dims).transpose(*transpose_order)

        # Build coordinates unless provided.
        # TODO: ADD CONDITION WHEN shape(COORDS) != shape(Dataset)
        if 'coords' in dir(self) or 'coords' in kwargs:
            if 'coords' in dir(self):
                coords = self.coords
            elif 'coords' in kwargs:
                coords = kwargs['coords']
            if type(coords)!= dict:
                raise Warning('Coords argument is not a dictionary, we will use the dimensions {0},{1}.'.format(flip_dict(rename_dims)['X'],flip_dict(rename_dims)['Y']))
                coords = {'X':dataset.X, 'Y':dataset.Y}
            elif 'X' not in coords.keys() or 'Y' not in coords.keys() : 
                rename_coords = rename_dict(coords)
                coords = {rename_coords[key]:item for key, item in coords.items()}
        else:
            coords = {'X':dataset.X, 'Y':dataset.Y}

        # 
        if 'variable' in dir(self) or 'variable' in kwargs:
            if 'variable' in dir(self):
                variable = self.variable
            elif 'variable' in kwargs:
                variable = kwargs['variable']
            if variable not in var_names:
                raise ValueError("Selected variable does not exist, make sure the variable is: {0} ".format(var_names))
        else:
            if len(var_names)==1:
                variable = list(var_names)[0]
            else:
                raise ValueError("Provide variable argument, Dataset contains multiple variables: {0} ".format(var_names))

        init(self,dataset, coords, variable, **kwargs)

    return test_input

def check_single_level(func):
    """
    docstring
    """
    @functools.wraps(func)
    def test_single_level(self, *args, **kwargs):  
        arg_names = inspect.getfullargspec(func)[0]
        for name, value in zip(arg_names[1:], args):
            setattr(self, name, value)

        for name, value in kwargs.items():
            setattr(self, name, value)

        if len(self.DataArray.dims)==3:
            raise ValueError("Dataset is 3D, use _scan_eddy_in_time instead")

        if not isinstance(self.level, float) and not isinstance(self.level, int):
            raise ValueError("_scan_eddy_single_level only supports one level at the time, use _scan_eddy_multiple_level instead")
        
        if self.spatial_filter['type'] == 'convolution':
            self.filter_data_spatially()
        
        if self.polarity == 'both':
            self.level = np.array([-abs(self.level),abs(self.level)])
        elif self.polarity == 'pos' and np.sign(self.level)==1:
            self.level = abs(self.level)
        elif self.polarity == 'pos' and np.sign(self.level)==-1:
            warnings.warn("Polarity and level sign are inconsistent, Polarity will replace level to: {0}".format(self.level), SyntaxWarning)
            self.level = abs(self.level)
        elif self.polarity == 'neg' and np.sign(self.level)==-1:
            self.level = -abs(self.level)
        elif self.polarity == 'neg' and np.sign(self.level)==1:
            warnings.warn("Polarity and level sign are inconsistent, Polarity will replace level to: {0}".format(-abs(self.level)), SyntaxWarning)
            self.level = -abs(self.level)
        else: 
            raise ValueError('polarity and level argument must be provided.')

        #contours, contours_rossby = func(self,*args, **kwargs)
        func(self,*args, **kwargs)

        #return contours, contours_rossby

    return test_single_level


def check_multiple_levels(func):
    """
    docstring
    """
    @functools.wraps(func)
    def test_multiple_level(self, *args, **kwargs):  
        arg_names = inspect.getfullargspec(func)[0]
        for name, value in zip(arg_names[1:], args):
            setattr(self, name, value)

        for name, value in kwargs.items():
            setattr(self, name, value)
        
        self.levels = np.array(self.levels,dtype=float)
        func(self,*args, **kwargs)

        #return contours, contours_rossby

    return test_multiple_level

def rename_dict(dims):
    rename_dims = {}
    for key in dims.keys():
        if "lon" in key.lower()  or key=="x":
            rename_dims[key] = 'X'
        elif "lat" in key.lower()  or key=="y":
            rename_dims[key] = 'Y'
        elif "time" in key.lower()  or key=="t":
            rename_dims[key] = 'T'
    return rename_dims

def flip_dict(dictionary):
    return {v: k for k, v in dictionary.items()}