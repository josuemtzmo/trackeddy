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

def basemap_mplot(x,y,data,title,projection='ortho',lat_0=-90,lon_0=-100,boundinglat=-30,resolution='c',scale='Lin',vmin='',vmax='',cmap=cm.cm.thermal,xan=1,yan=1,figsize=(5,5),fontsize=15,dpi=300):
    fig, ax = plt.subplots(yan, xan, figsize=figsize,dpi=dpi)
    
    X,Y=np.meshgrid(x,y)
    count=0
    for ii in range(0,xan):
        for jj in range(0,yan):
            if xan==1 and yan==1:
                data=[data]
                title=[title]
                ttl=plt.title(title[count], fontsize=fontsize)
                ttl.set_position([.5, 1.05+fontsize*0.001])
                map = Basemap(projection=projection,lat_0=lat_0,lon_0=lon_0,boundinglat=boundinglat,resolution=resolution)
                lonm,latm=map(X,Y)
                map.drawmeridians(np.arange(0,360,30),labels=[1,1,0,0],fontsize=int(fontsize*0.6))
                map.drawparallels(np.arange(-90,90,30),labels=[1,1,0,0],fontsize=int(fontsize*0.6))
                map.fillcontinents(color='black',lake_color='aqua')
                map.drawcoastlines()
                map.drawcoastlines()
                m=plt
            elif xan!=1 and yan==1:
                ttl=ax[ii].set_title(title[count], fontsize=fontsize)
                ttl.set_position([.5, 1.05+fontsize*0.001])
                map = Basemap(projection=projection,lat_0=lat_0,lon_0=lon_0,boundinglat=boundinglat,resolution=resolution,ax=ax[ii])
                lonm,latm=map(X,Y)
                map.drawmeridians(np.arange(0,360,30),labels=[1,1,0,0],fontsize=int(fontsize*0.6))
                map.drawparallels(np.arange(-90,90,30),labels=[1,1,0,0],fontsize=int(fontsize*0.6))
                map.fillcontinents(color='black',lake_color='aqua')
                map.drawcoastlines()
                map.drawcoastlines()
                m=ax[ii]
            elif xan==1 and yan!=1:
                ttl=ax[jj].set_title(title[count], fontsize=fontsize)
                ttl.set_position([.5, 1.05+fontsize*0.001])
                map = Basemap(projection=projection,lat_0=lat_0,lon_0=lon_0,boundinglat=boundinglat,resolution=resolution,ax=ax[jj])
                lonm,latm=map(X,Y)
                map.drawmeridians(np.arange(0,360,30),labels=[1,1,0,0],fontsize=int(fontsize*0.6))
                map.drawparallels(np.arange(-90,90,30),labels=[1,1,0,0],fontsize=int(fontsize*0.6))
                map.fillcontinents(color='black',lake_color='aqua')
                map.drawcoastlines()
                map.drawcoastlines()
                m=ax[jj]
            else:
                ttl=ax[jj,ii].set_title(title[count], fontsize=fontsize)
                ttl.set_position([.5, 1.05+fontsize*0.001])
                map = Basemap(projection=projection,lat_0=lat_0,lon_0=lon_0,boundinglat=boundinglat,resolution=resolution,ax=ax[jj,ii])
                lonm,latm=map(X,Y)
                map.drawmeridians(np.arange(0,360,30),labels=[1,1,0,0],fontsize=int(fontsize*0.6))
                map.drawparallels(np.arange(-90,90,30),labels=[1,1,0,0],fontsize=int(fontsize*0.6))
                map.fillcontinents(color='black',lake_color='aqua')
                map.drawcoastlines()
                map.drawcoastlines()
                m=ax[jj,ii]
                
            lonm,latm=map(X,Y)
            if scale =='Lin' and vmin == '' and vmax == '':
                vmintemp=data[count].min()
                vmaxtemp=data[count].max()
                im=m.pcolormesh(lonm,latm,data[count],cmap=cmap,vmin=vmintemp,vmax=vmaxtemp)
            elif scale == 'SymLog' and vmin == '' and vmax == '':
                vmintemp=data[count].min()
                vmaxtemp=data[count].max()
                im=m.pcolormesh(lonm,latm,data[count],cmap=cmap,norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,vmin=vmintemp, vmax=vmaxtemp))
            elif scale == 'Lin' and vmin != '' and vmax != '':
                im=m.pcolormesh(lonm,latm,data[count],cmap=cmap,vmin=vmin, vmax=vmax)
            elif scale == 'SymLog' and vmin != '' and vmax != '':
                im=m.pcolormesh(lonm,latm,data[count],cmap=cmap,norm=colors.SymLogNorm(linthresh=0.01, linscale=0.01,vmin=vmin, vmax=vmax))
            elif scale == 'Log' and vmin != '' and vmax != '':
                im=m.pcolormesh(lonm,latm,abs(data[count]),cmap=cmap,norm=colors.LogNorm(vmin=vmin, vmax=vmax))
            elif scale == 'Log2+' and vmin != '' and vmax != '':
                lev_exp = np.arange(np.floor(np.log10(0.001)-1),np.ceil(np.log10(0.3)+1))
                levs = np.power(10, lev_exp)
                im=m.contourf(lonm,latm,abs(data[count]),levs,cmap=cmap,norm=colors.LogNorm(vmin=vmin, vmax=vmax),alpha=0.5)
            count=count+1
    return fig,im,ax 



def animation_plot(x,y,data,title,projection='ortho',lat_0=-90,lon_0=-100,resolution='c',scale='Lin',vmin='',vmax='',cmap=cm.cm.thermal,xan=1,yan=1,figsize=(5,5),fontsize=15):
    fig, ax = plt.subplots(xan, yan, figsize=figsize)
    X,Y=np.meshgrid(x,y)
    count=0
    for ii in range(0,xan):
        for jj in range(0,yan):
            if xan==1 and yan==1:
                data=[data]
                title=[title]
                ttl=plt.title(title[count], fontsize=fontsize)
                ttl.set_position([.5, 1.05])
                map = Basemap(projection=projection,lat_0=lat_0,lon_0=lon_0,resolution=resolution)
                lonm,latm=map(X,Y)
                map.drawmeridians(np.arange(0,360,30),labels=[1,1,0,0],fontsize=10)
                map.drawparallels(np.arange(0,360,30),labels=[1,1,0,0],fontsize=10)
                map.fillcontinents(color='black',lake_color='aqua')
                map.drawcoastlines()
                map.drawcoastlines()
                m=plt
            elif xan!=1 and yan==1:
                ttl=ax[ii].set_title(title[count], fontsize=fontsize)
                ttl.set_position([.5, 1.05])
                map = Basemap(projection=projection,lat_0=lat_0,lon_0=lon_0,resolution=resolution,ax=ax[ii])
                lonm,latm=map(X,Y)
                map.drawmeridians(np.arange(0,360,30),labels=[1,1,0,0],fontsize=10)
                map.drawparallels(np.arange(0,360,30),labels=[1,1,0,0],fontsize=10)
                map.fillcontinents(color='black',lake_color='aqua')
                map.drawcoastlines()
                map.drawcoastlines()
                m=ax[ii]
            elif xan==1 and yan!=1:
                ttl=ax[jj].set_title(title[count], fontsize=fontsize)
                ttl.set_position([.5, 1.05])
                map = Basemap(projection=projection,lat_0=lat_0,lon_0=lon_0,resolution=resolution,ax=ax[jj])
                lonm,latm=map(X,Y)
                map.drawmeridians(np.arange(0,360,30),labels=[1,1,0,0],fontsize=10)
                map.drawparallels(np.arange(0,360,30),labels=[1,1,0,0],fontsize=10)
                map.fillcontinents(color='black',lake_color='aqua')
                map.drawcoastlines()
                map.drawcoastlines()
                m=ax[jj]
            else:
                ttl=ax[ii,jj].set_title(title[count], fontsize=fontsize)
                ttl.set_position([.5, 1.05])
                map = Basemap(projection=projection,lat_0=lat_0,lon_0=lon_0,resolution=resolution,ax=ax[ii,jj])
                lonm,latm=map(X,Y)
                map.drawmeridians(np.arange(0,360,30),labels=[1,1,0,0],fontsize=10)
                map.drawparallels(np.arange(0,360,30),labels=[1,1,0,0],fontsize=10)
                map.fillcontinents(color='black',lake_color='aqua')
                map.drawcoastlines()
                map.drawcoastlines()
                m=ax[ii,jj]
                
            lonm,latm=map(X,Y)
            if scale =='Lin' and vmin == '' and vmax == '':
                vmintemp=data[count].min()
                vmaxtemp=data[count].max()
                im=m.pcolormesh(lonm,latm,data[count],cmap=cmap,vmin=vmintemp,vmax=vmaxtemp)
                cb1 = fig.colorbar(im,ax=m)
            elif scale == 'SymLog' and vmin == '' and vmax == '':
                vmintemp=data[count].min()
                vmaxtemp=data[count].max()
                im=m.pcolormesh(lonm,latm,data[count],cmap=cmap,norm=colors.SymLogNorm(linthresh=0.03,\
                                                                     linscale=0.03,vmin=vmintemp, vmax=vmaxtemp))
                cb1 = fig.colorbar(im,ax=m)
            elif scale == 'Lin' and vmin != '' and vmax != '':
                im=m.pcolormesh(lonm,latm,data[count],cmap=cmap,vmin=vmin, vmax=vmax)
                cb1 = fig.colorbar(im,ax=m)
            elif scale == 'SymLog' and vmin != '' and vmax != '':
                im=m.pcolormesh(lonm,latm,data[count],cmap=cmap,norm=colors.SymLogNorm(linthresh=0.01, linscale=0.01,vmin=vmin,\
                                                                                       vmax=vmax))
                cb1 = fig.colorbar(im,ax=m)
            elif scale == 'Log' and vmin != '' and vmax != '':
                im=m.pcolormesh(lonm,latm,abs(data[count]),cmap=cmap,norm=colors.LogNorm(vmin=vmin, vmax=vmax))
                cb1 = fig.colorbar(im,ax=m)
            elif scale == 'Log2+' and vmin != '' and vmax != '':
                lev_exp = np.arange(np.floor(np.log10(0.001)-1),np.ceil(np.log10(0.3)+1))
                levs = np.power(10, lev_exp)
                im=m.contourf(lonm,latm,abs(data[count]),levs,cmap=cmap,norm=colors.LogNorm(vmin=vmin, vmax=vmax),alpha=0.5)
                cb1 = fig.colorbar(im,ax=m)
            count=count+1
            
            def init():
                m.set_array([])
            return m

            def animate(iter):
                m.set_array(sqrt(data[count][iter,:,:]).ravel())
            return m
    fig.tight_layout(fig)
    anim2 = animation.FuncAnimation(fig,animate,frames=90,interval=100,blit=False,repeat=False)
    return anim2
