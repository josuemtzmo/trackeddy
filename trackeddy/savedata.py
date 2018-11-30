import numpy as np
import netCDF4 as nc4
from datetime import datetime
import pylab as plt
import warnings
warnings.filterwarnings("ignore")



def vargeonc(filename,lat,lon,var,tt,varname,init_time=datetime.datetime(1993, 1, 1),nc_description='',units='',dt='',dim='2D',format='NETCDF4'):
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
        time[:]=[init_time + datetime.timedelta(days=i) for i in range(0,tt)]
        
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
    
    time.units = 'days since Jan 01, 0001'
    time.cartesian_axis = "T"
    varnc.units = units
    f.close()

def dimensionseddy(eddystructure):
    '''
    *********Function to get dimensions of eddies ***********
    Function to obtain dimensions of the eddy structure, to 
    save it as netCDF.
    Usage:
    
    Example:

    Author: Josue Martinez Moreno, 2017
    '''
    num_eddy=len(eddystructure)
    time=0
    for key, value in list(eddystructure.items()):
        if type(value['time'])!=int:
            if len(value['time'])>time:
                time=len(value['time'])
    i_contour=0
    for key, value in list(eddystructure.items()):
        if type(value['time'])==int:
            timelist=[value['time']]
        else:
            timelist=value['time']
        iitime=0
        for ii in timelist:
            if max(shape(value['contour'][iitime]))>i_contour:
                i_contour=max(shape(value['contour'][iitime]))
            iitime=iitime+1
    i_ellipse=50
    dimensionsdict={'time':time,'eddies':num_eddy,'i_contour':i_contour,'i_ellipse':i_ellipse,'xy':2}
    return dimensionsdict

###

def eddy2ncformat(eddystructure,dimensionsdict):
    '''
    *********Function to restructure of eddy structure ***********
    Function to arrange eddy structure, to netCDF structure.
    Usage:
    
    Example:

    Author: Josue Martinez Moreno, 2017
    '''
    ellipse_ncvar=zeros([dimensionsdict['time'],dimensionsdict['eddies'],dimensionsdict['xy']\
                         ,dimensionsdict['i_ellipse']])
    contour_ncvar=zeros([dimensionsdict['time'],dimensionsdict['eddies'],dimensionsdict['xy']
                         ,dimensionsdict['i_contour']])
    angle_ncvar=zeros([dimensionsdict['time'],dimensionsdict['eddies']])
    area_ncvar=zeros([dimensionsdict['time'],dimensionsdict['eddies']])
    position_ncvar=zeros([dimensionsdict['time'],dimensionsdict['eddies'],dimensionsdict['xy']])
    position_ellipse_ncvar=zeros([dimensionsdict['time'],dimensionsdict['eddies'],dimensionsdict['xy']])
    eddycount=0
    for key, value in list(eddystructure.items()):
        iitime=0
        if type(value['time'])==int:
            timelist=[value['time']]
        else:
            timelist=value['time']
        for ii in timelist:
            ellipse_ncvar[ii,eddycount,:,:]=value['ellipse'][iitime]
            contour_ncvar[ii,eddycount,:,0:shape(value['contour'][iitime])[1]]=value['contour'][iitime]
            angle_ncvar[ii,eddycount]=value['angle'][iitime]
            area_ncvar[ii,eddycount]=value['area'][iitime]
            position_ncvar[ii,eddycount,:]=value['position'][iitime]
            position_ellipse_ncvar[ii,eddycount,:]=value['position_eddy'][iitime]
            iitime=iitime+1
        eddycount=eddycount+1
    
    contour_ncvar[contour_ncvar==0]=9999
    eddynetCDFstructure={'ellipse_ncvar':ellipse_ncvar,'contour_ncvar':contour_ncvar,'angle_ncvar':angle_ncvar\
                         ,'area_ncvar':area_ncvar,'position_ncvar':position_ncvar\
                         ,'position_ellipse_ncvar':position_ellipse_ncvar}
    return eddynetCDFstructure

def eddync(filename,eddystructure):
    '''
    *********Function to get dimensions of eddies ***********
    Function to obtain dimensions of the eddy structure, to 
    save it as netCDF.
    Usage:
    
    Example:

    Author: Josue Martinez Moreno, 2017
    '''
    dimensionsdict=dimensionseddy(eddystructure)
    f = nc4.Dataset(filename,'w', format='NETCDF4')
    f.createDimension('time', dimensionsdict['time'])
    f.createDimension('eddies',dimensionsdict['eddies'])
    f.createDimension('i_contour', dimensionsdict['i_contour'])
    f.createDimension('i_ellipse', dimensionsdict['i_ellipse'])
    f.createDimension('xy', dimensionsdict['xy'])
    eddynetCDFstructure=eddy2ncformat(eddystructure,dimensionsdict)
    ellipse=f.createVariable('Ellipsoid', np.float64, ('time','eddies','xy','i_ellipse'))
    contour=f.createVariable('Contour', np.float64, ('time','eddies','xy','i_contour'), fill_value=9999, zlib=True)
    area=f.createVariable('Area', np.float64, ('time','eddies'))
    angle=f.createVariable('Angle', np.float64, ('time','eddies'))
    position=f.createVariable('Position', np.float64, ('time','eddies','xy'))
    position_ellipse=f.createVariable('Position_Ellipse', np.float64, ('time','eddies','xy'))
        
    ellipse[:] = eddynetCDFstructure['ellipse_ncvar']
    ellipse.units = 'Lon,Lat (degrees)'
    contour[:] = eddynetCDFstructure['contour_ncvar']
    contour.units = 'Lon,Lat (degrees)'
    area[:] = eddynetCDFstructure['area_ncvar']
    angle[:] = eddynetCDFstructure['angle_ncvar']
    angle.units = 'Degrees^2'
    position[:] = eddynetCDFstructure['position_ncvar']
    position.units = 'Lon,Lat (degrees)'
    position_ellipse[:] = eddynetCDFstructure['position_ellipse_ncvar']
    position_ellipse.units = 'Lon,Lat (degrees)'
    
    #Add global attributes
    today = datetime.today()
    f.description = 'Eddy tracking file in netCDF format.'
    f.history = "Created " + today.strftime("%d/%m/%y")
    
    time.units = 'days since Jan 01, 0001'
    
    f.close()
    
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