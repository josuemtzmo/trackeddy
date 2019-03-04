import numpy as np
import netCDF4 as nc4
from datetime import datetime, timedelta
import pylab as plt
import warnings
import os
import time
import subprocess
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


class Trackeddy2dataset():
    """
    infield =  Input field or File
    """
    def __init__(self,infield,outpath,outformat,source='',reference='',time_coverage=['','']):
        self.outpath = outpath
        self.format  = outformat 
        self.source  = source
        self.reference = reference
        self.time_coverage = time_coverage
        if type(infield) == str:
            self.file = infield
            self.eddydict = None
        elif type(infield) == dict:
            self.eddydict = infield
        else:
            raise ValueError("File should be a string or a dictionary.")

    def file2dict(self):
        if os.path.isfile(self.file):
            npyfile=np.load(self.file)
            self.eddydict=npyfile.item()
        elif os.path.isdir(self.file):
            file=[f for f in os.listdir(self.file) if f.split('.')[-1] == 'npy']
            for f in file:
                npyfile = np.load(f)
                if self.eddydict == None:
                    self.eddydict = npyfile.item()
                else: 
                    self.join2dicts(npyfile)
        else:
            raise ValueError("File or path doesn't exist.")

    def join2dicts(self,npyfile):
        pass
            
    def trackeddy2nc(self):
        if 'nc' not in self.outpath.split('/')[-1]:
            outf=[f for f in os.listdir(self.outpath) if 'nc' in f and 'output_' in f]
            count = len(outf)
            self.fpath=self.outpath+'output_{:03}.nc'.format(count)
        else:
            self.fpath=self.outpath
        self.ncfile = nc4.Dataset(self.fpath,'w')
        self.metadata()
        self.dict2nc()
        self.ncfile.close()

    def metadata(self):
        if self.format == 'nc' or self.format == 'netCDF' or self.format == 'netcdf':
            if self.source == None:
                self.source = 'Satellite Altimeter Observations of SLA'
            self.ncfile.title = 'Mesoscale Eddies from' + self.source
            self.ncfile.Metadata_Conventions = 'Unidata Dataset Discovery v1.0'
            self.ncfile.institution = 'Australian National University - National Computational Facility'
            self.ncfile.project = 'TrackEddy'
            self.ncfile.creator_url = 'https://trackeddy.readthedocs.io/en/latest/'
            self.ncfile.creator_email = 'josue.martinezmoreno@anu.edu.au'
            self.ncfile.creator_name = 'Josue Martinez Moreno'
            self.ncfile.license = ''
            self.ncfile.summary = 'This dataset contains eddy locations and parameters over global ocean from ' + self.source
            self.ncfile.comment = 'Surface product; mesoscale eddies; eddy field reconstruction'
            self.ncfile.input_product_reference = self.reference
            self.ncfile.date_created = time.strftime('%Y-%B-%d T %H:%M:%S GMT', time.gmtime())
            self.ncfile.time_coverage_duration = ''
            self.ncfile.time_coverage_start = self.time_coverage[0]
            self.ncfile.time_coverage_end = self.time_coverage[1]
            self.ncfile.timstandard_name_vocabulary = 'NetCDF Climate and Forecast (CF) Metadata Convention Standard Name Table v37'
            self.ncfile.product_version = "1.0."+self.source.split(' ')[0]
            process = subprocess.Popen(['git', 'rev-parse', '--short', 'HEAD'], shell=False, stdout=subprocess.PIPE)
            git_head_hash = process.communicate()[0].strip()
            self.ncfile.trackeddy_version = 'Version: 1.0 - Git Hash: {}'.format(git_head_hash.decode("utf-8"))
            self.ncfile.loaddata='How to load the data can be found at: '

    def dict2nc(self):
        for k,i in self.eddydict.items():
            grp=self.ncfile.createGroup(k.split('_')[1])
            time = grp.createDimension('time', None)
            for k1,i1 in i.items():
                if k1 == 'time':
                    times=grp.createVariable('time','f8',('time',),zlib=True)
                    if len(i1)==1:
                        times[:]=i1[0]
                    else:
                        times[:]=i1[:,0]
                elif k1 == 'level':
                    level=grp.createVariable('level','f8',('time',),zlib=True)
                    if len(i1)==1:
                        level[:]=i1[0]
                    else:
                        level[:]=i1[:,0]
                elif k1 == 'angle':
                    angle=grp.createVariable('angle','f8',('time',),zlib=True)
                    if len(i1)==1:
                        angle[:]=i1[0]
                    else:
                        angle[:]=i1[:,0]
                elif k1 == 'majoraxis':
                    majoraxis_lon=grp.createVariable('majoraxis_lon','f8',('time',),zlib=True)
                    majoraxis_lon[:]=i1[:,0]
                    majoraxis_lat=grp.createVariable('majoraxis_lat','f8',('time',),zlib=True)
                    majoraxis_lat[:]=i1[:,1]
                elif k1 == 'minoraxis':
                    minoraxis_lon=grp.createVariable('minoraxis_lon','f8',('time',),zlib=True)
                    minoraxis_lon[:]=i1[:,0]
                    minoraxis_lat=grp.createVariable('minoraxis_lat','f8',('time',),zlib=True)
                    minoraxis_lat[:]=i1[:,1]
                elif k1 == 'position_maxvalue':
                    lon=grp.createVariable('lon','f8',('time',),zlib=True)
                    if len(i1)==1:
                        lon[:]=i1[0][0]
                    else:
                        lon[:]=i1[:,0]
                    lat=grp.createVariable('lat','f8',('time',),zlib=True)
                    if len(i1)==1:
                        lat[:]=i1[0][1]
                    else:
                        lat[:]=i1[:,1]
                    amp=grp.createVariable('amp','f8',('time',),zlib=True)
                    if len(i1)==1:
                        amp[:]=i1[0][2]
                    else:
                        amp[:]=i1[:,2]
                elif k1 == 'position_default':
                    print(k)
                    CM_lon=grp.createVariable('CM_lon','f8',('time',),zlib=True)
                    if len(i1)==1:
                        CM_lon[:]=i1[0][1]
                    else:
                        CM_lon[:]=i1[:,0]
                    CM_lat=grp.createVariable('CM_lat','f8',('time',),zlib=True)
                    if len(i1)==1:
                        CM_lat[:]=i1[0][1]
                    else:
                        CM_lat[:]=i1[:,1]
                elif k1 == '2dgaussianfit':
                    sigma_x=grp.createVariable('sigma_x','f8',('time',),zlib=True)
                    sigma_x[:]=i1[:,0]
                    sigma_y=grp.createVariable('sigma_y','f8',('time',),zlib=True)
                    sigma_y[:]=i1[:,1]
                    theta=grp.createVariable('theta','f8',('time',),zlib=True)
                    theta[:]=i1[:,2]
                else:
                    pass
    
    def trackeddy2hdf5(self):
        import h5py
        f = h5py.File(path,'w')

        for key,item in self.eddytd.items():
            grp=f.create_group(key)
            for key1,item1 in item.items():
                if key1!='ellipse' and key1!='contour':
                    print(key1)
                    grp.create_dataset(key1,data=item1)    
    

def save_data(path, variable):
    '''
    *********Function to save data as npy***********
    Usage:
    
    Example:

    '''
    np.save(path, variable)
