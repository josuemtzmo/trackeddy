import warnings
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use('Agg')
import pylab as plt
import cmocean as cm
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation

def basemap_mplot(x,y,data,title,projection='ortho',lat_0=-90,lon_0=-100,boundinglat=-30,resolution='c',scale='Lin',vmin='',vmax='',cmap=cm.cm.thermal,xan=1,yan=1,figsize=(5,5),fontsize=15,dpi=300,freq=[30,30],drawlinewidth=1):
    fig, ax = plt.subplots(yan, xan, figsize=figsize,dpi=dpi)
    
    X,Y=np.meshgrid(x,y)
    count=0
    ims=[]
    for ii in range(0,xan):
        for jj in range(0,yan):
            if xan==1 and yan==1:
                data=[data]
                title=[title]
                ttl=plt.title(title[count], fontsize=fontsize)
                ttl.set_position([.5, 1.05+fontsize*0.001])
                if projection == 'ortho' or projection=='mbtfpq':
                    map = Basemap(projection=projection,lat_0=lat_0,lon_0=lon_0,boundinglat=boundinglat,resolution=resolution,ax=ax)
                    labels=[[0,0,0,0],[1,0,0,0]]
                else:
                    map = Basemap(projection=projection,lat_0=lat_0,lon_0=lon_0,llcrnrlon=x.min(),llcrnrlat=y.min(),urcrnrlon=x.max(),urcrnrlat=y.max(),boundinglat=boundinglat,resolution=resolution,ax=ax)
                    labels=[[0,0,0,1],[1,0,0,0]]
                lonm,latm=map(X,Y)
                map.drawmeridians(np.arange(0,360,freq[1]),labels=labels[0],fontsize=int(fontsize*0.6),linewidth=drawlinewidth)
                map.drawparallels(np.arange(-90,90,freq[0]),labels=labels[1],fontsize=int(fontsize*0.6),linewidth=drawlinewidth)
                map.fillcontinents(color='gray',lake_color='gray')
                #map.drawlsmask(land_color = "gray",resolution = resolution,zorder=0)
                #map.drawcoastlines(linewidth=drawlinewidth)
                m=plt
            elif xan!=1 and yan==1:
                ttl=ax[ii].set_title(title[count], fontsize=fontsize)
                ttl.set_position([.5, 1.05+fontsize*0.001])
                if projection == 'ortho' or projection=='mbtfpq':
                    map = Basemap(projection=projection,lat_0=lat_0,lon_0=lon_0,boundinglat=boundinglat,resolution=resolution,ax=ax[ii])
                    labels=[[0,0,0,0],[1,0,0,0]]
                else:
                    map = Basemap(projection=projection,lat_0=lat_0,lon_0=lon_0,llcrnrlon=x.min(),llcrnrlat=y.min(),urcrnrlon=x.max(),urcrnrlat=y.max(),boundinglat=boundinglat,resolution=resolution,ax=ax[ii])
                    labels=[[0,0,0,1],[1,0,0,0]]
                lonm,latm=map(X,Y)
                map.drawmeridians(np.arange(0,360,freq[1]),labels=labels[0],fontsize=int(fontsize*0.6),linewidth=drawlinewidth)
                map.drawparallels(np.arange(-90,90,freq[0]),labels=labels[1],fontsize=int(fontsize*0.6),linewidth=drawlinewidth)
                map.fillcontinents(color='gray',lake_color='gray')
                #map.drawlsmask(land_color = "gray",resolution = resolution)
                #map.drawcoastlines(linewidth=drawlinewidth)
                m=ax[ii]
            elif xan==1 and yan!=1:
                ttl=ax[jj].set_title(title[count], fontsize=fontsize)
                ttl.set_position([.5, 1.01+fontsize*0.001])
                if projection == 'ortho' or projection=='mbtfpq':
                    map = Basemap(projection=projection,lat_0=lat_0,lon_0=lon_0,boundinglat=boundinglat,resolution=resolution,ax=ax[jj])
                    labels=[[0,0,0,0],[1,0,0,0]]
                else:
                    map = Basemap(projection=projection,lat_0=lat_0,lon_0=lon_0,llcrnrlon=x.min(),llcrnrlat=y.min(),urcrnrlon=x.max(),urcrnrlat=y.max(),boundinglat=boundinglat,resolution=resolution,ax=ax[jj])
                    labels=[[0,0,0,1],[1,0,0,0]]
                lonm,latm=map(X,Y)
                map.drawmeridians(np.arange(0,360,freq[1]),labels=labels[0],fontsize=int(fontsize*0.6),linewidth=drawlinewidth)
                map.drawparallels(np.arange(-90,90,freq[0]),labels=labels[1],fontsize=int(fontsize*0.6),linewidth=drawlinewidth)
                map.fillcontinents(color='gray',lake_color='gray')
                #map.drawlsmask(land_color = "gray",resolution = resolution)
                #map.drawcoastlines(linewidth=drawlinewidth)
                
                m=ax[jj]
            else:
                ttl=ax[jj,ii].set_title(title[count], fontsize=fontsize)
                ttl.set_position([.5, 1.05+fontsize*0.001])
                if projection == 'ortho' or projection=='mbtfpq':
                    map = Basemap(projection=projection,lat_0=lat_0,lon_0=lon_0,boundinglat=boundinglat,resolution=resolution,ax=ax[jj,ii])
                else:
                    map = Basemap(projection=projection,lat_0=lat_0,lon_0=lon_0,llcrnrlon=x.min(),llcrnrlat=y.min(),urcrnrlon=x.max(),urcrnrlat=y.max(),boundinglat=boundinglat,resolution=resolution,ax=ax[jj,ii])
                lonm,latm=map(X,Y)
                map.drawmeridians(np.arange(0,360,freq[1]),labels=[0,0,0,1],fontsize=int(fontsize*0.6),linewidth=drawlinewidth)
                map.drawparallels(np.arange(-90,90,freq[0]),labels=[1,0,0,0],fontsize=int(fontsize*0.6),linewidth=drawlinewidth)
                map.fillcontinents(color='gray',lake_color='gray')
                #map.drawlsmask(land_color = "gray",resolution = resolution)
                #map.drawcoastlines(linewidth=drawlinewidth)
                m=ax[jj,ii]
                
            lonm,latm=map(X,Y)
            if scale =='Lin' and vmin == '' and vmax == '':
                vmintemp=data[count].min()
                vmaxtemp=data[count].max()
                im=m.pcolormesh(lonm,latm,data[count],cmap=cmap,vmin=vmintemp,vmax=vmaxtemp,rasterized=True)
            elif scale == 'SymLog' and vmin == '' and vmax == '':
                vmintemp=data[count].min()
                vmaxtemp=data[count].max()
                im=m.pcolormesh(lonm,latm,data[count],cmap=cmap,norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,vmin=vmintemp, vmax=vmaxtemp),rasterized=True)
            elif scale == 'Lin' and type(vmin)==list and type(vmax)==list:
                im=m.pcolormesh(lonm,latm,data[count],cmap=cmap[count],vmin=vmin[count],vmax=vmax[count],rasterized=True)
            elif scale == 'Lin' and vmin != '' and vmax != '':
                im=m.pcolormesh(lonm,latm,data[count],cmap=cmap,vmin=vmin, vmax=vmax,rasterized=True)
            elif scale == 'Linc' and vmin != '' and vmax != '':
                im=m.contourf(lonm,latm,data[count],cmap=cmap,vmin=vmin, vmax=vmax,levels=np.linspace(vmin,vmax,(vmax+1)-vmin))
            elif scale == 'SymLog' and vmin != '' and vmax != '':
                im=m.pcolormesh(lonm,latm,data[count],cmap=cmap,norm=colors.SymLogNorm(linthresh=0.0001, linscale=0.0001,vmin=vmin, vmax=vmax),rasterized=True)
            elif scale == 'Log' and vmin != '' and vmax != '':
                im=m.pcolormesh(lonm,latm,abs(data[count]),cmap=cmap,norm=colors.LogNorm(vmin=vmin, vmax=vmax),rasterized=True)
            elif scale == 'Log2+' and vmin != '' and vmax != '':
                lev_exp = np.arange(np.floor(np.log10(0.001)-1),np.ceil(np.log10(0.3)+1))
                levs = np.power(10, lev_exp)
                im=m.contourf(lonm,latm,abs(data[count]),levs,cmap=cmap,norm=colors.LogNorm(vmin=vmin, vmax=vmax),alpha=0.5,rasterized=True)
            count=count+1
            ims.append(im)
    return fig,ims,ax,map