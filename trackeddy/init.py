import numpy as np

def arkbinterpmom5(tracer,vel_u,vel_v,lont,latt,lonu,latu):
    '''
    *************arkbinterpmom5*****************
    Function to interpolate the Arakawa B grid of the
    model MOM5. It changes the velocity field grid to 
    the tracer grid. It means this function adjust the
    spatial diference between grids.
    Usage:
    tracer= Tracer as Temperature, Salinity, SSH, etc
    vel_u=Component U of velocity field
    vel_v=Component V of velocity field
    lont=Longitude coordinate of tracer grid
    latt=Latitude coordinate of tracer grid
    lonu=Longitude coordinate of velocity field grid
    latu=Latitude coordinate of velocity field grid
    '''
    print('***********Interp 4 Arakawa B************')