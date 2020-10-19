# Importing all libraries.
from pylab import *
from netCDF4 import Dataset
import os
import cmocean as cm
from trackeddy.tracking import *
from trackeddy.datastruct import *
from trackeddy.geometryfunc import *
from trackeddy.physics import *

# Load Data
filepath = '../input/dt_global_allsat_phy_l4_20160901.nc'

track_class = TrackEddy(filepath)

closed_contours, contours_rossby = track_class._scan_eddy_single_level(level = 0.1,geo=True)

data = xr.open_dataset(filepath).squeeze()

plt.figure( figsize = (10, 5), dpi = 500)
# data.sla.plot.contourf(x='longitude',y='latitude')
# data.sla.plot.contour(x='longitude',y='latitude',levels=[0.1])

for ii in range( np.shape( closed_contours )[0]):
    plt.plot( closed_contours[ii][:,1], closed_contours[ii][:,0], '-r' )

for jj in range( np.shape( contours_rossby )[0]):
    plt.plot( contours_rossby[jj][:,1], contours_rossby[jj][:,0], '-b',linewidth=0.5 )
print(np.shape( closed_contours), np.shape(contours_rossby) )

plt.show()

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


