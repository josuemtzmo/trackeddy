# Importing all libraries.
from pylab import *
from netCDF4 import Dataset
import os
import cmocean as cm
from trackeddy.tracking import *
import xarray as xr
from trackeddy.trackeddy import *

import time 

# Load Data
filepath = './input/dt_global_allsat_phy_l4_20160901.nc'

tic = time.time()
track_class = TrackEddy(filepath)


# One level analysis
# track_class._scan_eddy_single_level(level = -0.1, polarity='neg', geo=True)

# Multiple level analysis
levels = np.arange(0.1 , 0.3, 0.1)
track_class._scan_eddy_multiple_levels(levels = levels, polarity='neg', geo=True)

toc = time.time()

print('New code: ', toc-tic, 's, Identified eddies:', *np.shape(track_class.identified_eddies_level))


# tic = time.time()

# data = xr.open_dataset(filepath)

# filters = {'time':{'type':'historical','t':None,'t0':None,'value':None},
#            'spatial':{'type':'moving','window':120,'mode':'uniform'}}

# preferences={'ellipse':0.7,'eccentricity':0.95,'gaussian':0.7}

# levels = 0.1

# eddytd=analyseddyzt(data.sla.values,data.longitude.values,data.latitude.values,0,np.shape(data.sla)[0],1,levels,mask='',maskopt='contour',timeanalysis='none'\
#                     ,preferences=preferences,filters=filters,destdir='',physics='',diagnostics=False,pprint=False)

# toc = time.time()

# print('Old code: ', toc-tic, 's, Identified eddies:', len(eddytd))


# data = xr.open_dataset(filepath).squeeze()

# plt.figure( figsize = (5, 2), dpi = 300)
# # data.sla.plot.contourf(x='longitude',y='latitude')
# # data.sla.plot.contour(x='longitude',y='latitude',levels=[0.1])

# for ii in range( np.shape( closed_contours )[0]):
#     plt.plot( closed_contours[ii][:,1], closed_contours[ii][:,0], '-r' )

# for jj in range( np.shape( contours_rossby )[0]):
#     plt.plot( contours_rossby[jj][:,1], contours_rossby[jj][:,0], '-b',linewidth=0.5 )
# print(np.shape( closed_contours), np.shape(contours_rossby) )

# plt.show()

# # Open netcdf Dataset.
# ncfile     = Dataset(filepath)
# # Load data into memory
# sla        = ncfile.variables['sla'][:]
# lon        = ncfile.variables['longitude'][:]
# lat        = ncfile.variables['latitude'][:]

# # Define area of study
# areamap = array([[0,len(lon)],[0,len(lat)]]) # Global option

# # Time and spatial filter
# filters = {'time':{'type':None,'t':None,'t0':None,'value':None},
#            'spatial':{'type':'moving','window':50,'mode':'uniform'}}

# # Mesoscale scaling 
# checkarea={'mesoscale':2*np.pi}

# # Eddy definition criteria
# preferences={'ellipse':0.85,'eccentricity':0.85,'gaussian':0.8}

# # Levels to be analysed and to extract positive eddies from anomaly
# levels = {'max':sla[0,:,:].max(),'min':0.01,'step':0.03}

# positive_eddies=analyseddyzt(sla,lon,lat,0,1,1,levels,preferences=preferences
#              ,areamap=areamap,areaparms=checkarea,filters=filters
#              ,maskopt='contour',diagnostics='contour',pprint=True,debug=False)


# import cmocean as cm

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
# plt.show()

# data4gauss.plot.contourf(levels=np.linspace(0,level,10))
# plt.plot(contour_yx_deg[:,1],contour_yx_deg[:,0])
# plt.plot(eddy_c_location[0],eddy_c_location[1],'.m')

# contour_centroid = np.mean(contour_yx_deg,axis=0)
# plt.plot(contour_centroid[1],contour_centroid[0],'*c')

# plt.show()