import numpy as np
import netCDF4 as nc4
from datetime import datetime, timedelta
import pylab as plt
import warnings
warnings.filterwarnings("ignore")



def vargeonc(filename,lat,lon,var,tt,varname,init_time=datetime(1993, 1, 1),nc_description='',units='',dt='',dim='2D',format='NETCDF4'):
    '''
    *************Save Variable to netCDF ***********
    Function to save a single variable to netCDF file,
    Usage:
    
    Example:

    Author: Josue Martinez Moreno, 2017
    '''
    f = nc4.Dataset(filename,'w', format=format)
    f.createDimension('lon', len(lon))
    f.createDimension('lat', len(lat))
    f.createDimension('time', tt)
    longitude = f.createVariable('lon', 'f4', 'lon')
    latitude = f.createVariable('lat', 'f4', 'lat')
    
    time = f.createVariable('time', 'i4', 'time')
    
    longitude[:] = lon
    latitude[:] = lat
    
    if tt==0:
        time[:]=tt
    else:
        calendar = 'standard'
        t_units = 'days since 1970-01-01 00:00'
        time.units=t_units
        time.calendar=calendar
        dates=[init_time + timedelta(days=i) for i in range(0,tt)]
        time[:] = nc4.date2num(dates, units=units, calendar=calendar)
   
    time.cartesian_axis = "Time"
        
    if dim == '3D':
        f.createDimension('z', len(z))
        levels = f.createVariable('z', 'i4', 'z')
        varnc = f.createVariable(varname, 'f4', ('time', 'z', 'lat', 'lon'))
        varnc[:,:,:,:] = var
        levels[:] = z
        levels.units = 'meters [m]'
    else:
        #f.createDimension('z', 1)
        #levels = f.createVariable('z', 'i4', 'z')
        #varnc = f.createVariable(varname, 'f4', ('time','z', 'lat', 'lon'))
        varnc = f.createVariable(varname, 'f4', ('time', 'lat', 'lon'))
        #levels[:] = 0
        #for ttt in range(0,tt):
        #    for yy in range(0,np.shape(var)[1]):
        #        for xx in range(0,np.shape(var)[2]):
        if tt==0:
            varnc[tt,:,:] = var
        else:
            varnc[:,:,:] = var
        #if tt==0:
        #    varnc[tt,0,:,:] = var
        #else:
        #    varnc[:,0,:,:] = var
    
    today = datetime.today()
    time_num = today.toordinal()
    
    #Add global attributes
    f.description = nc_description
    f.history = "Created " + today.strftime("%d/%m/%y")
    
    #Add local attributes to variable instances
    latitude.units = 'degrees_N'
    latitude.cartesian_axis = "Y"
    longitude.units = 'degrees_E'
    longitude.cartesian_axis = "X"
    
    varnc.units = units
    f.close()


def trackeddy2nc(path,eddytd):
    from netCDF4 import Dataset
    ncfile = Dataset(path,'w')

    for key,item in eddytd.items():
        grp=ncfile.createGroup(key)
    
def trackeddy2hdf5(path,eddytd):
    import h5py
    f = h5py.File(path,'w')

    for key,item in eddytd.items():
        grp=f.create_group(key)
        for key1,item1 in item.items():
            if key1!='ellipse' and key1!='contour':
                print(key1)
                grp.create_dataset(key1,data=item1)

                


def save_data(path, variable):
    '''
    *********Function to save data as txt ***********
    Usage:
    
    Example:

    Author: Josue Martinez Moreno, 2017
    '''
    with file(path, 'w') as variable_file:
        np.savetxt(variable_file, variable)
    variable_file.close() 
