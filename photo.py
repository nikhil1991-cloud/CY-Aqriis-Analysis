import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits

def get_flux_sky(flux,positions,pixel_radius,sky_annulus_1,sky_annulus_2):
    '''
    Input:
    
          flux: Flux map => Size (M,N)
        
          positions: Pixel positions of the star center in (x,y)
    
          pixel_radius: Radius of the circular photometry aperture in pixels (float32)
          
          sky_annulus_1: Lower limit of  annular ring used to calculate SKY counts (float32)
          
          sky_annulus_2: Upper limit of  annular ring used to calculate SKY counts (float32)
          
    Returns:
          
          Photometric counts => float32
    
          SKY counts => float32
    
          Photometry mask => size (M,N); 0=bad,1=good
          
          SKY mask => size (M,N); 0=bad,1=good
    '''
    Distance = np.zeros(np.shape(flux))
    X0 = np.linspace(0,np.shape(flux)[1]-1,np.shape(flux)[1])
    Y0 = np.linspace(0,np.shape(flux)[0]-1,np.shape(flux)[0])
    G = np.meshgrid(X0,Y0)
    X = G[0]
    Y = G[1]
    xc,yc = positions[1],positions[0]
    i=0
    for i in range (0,np.shape(flux)[0]):
        j=0
        for j in range (0,np.shape(flux)[1]):
            Distance[i][j] = np.sqrt((X[i][j]-xc)**2 + (Y[i][j]-yc)**2)

    Photometry_mask = np.ones(np.shape(flux))
    sky_mask = np.ones(np.shape(flux))
    i=0
    for i in range (0,np.shape(flux)[0]):
        j=0
        for j in range (0,np.shape(flux)[1]):
            if Distance[i][j] < pixel_radius:
                Photometry_mask[i,j] = 0
    i=0
    for i in range (0,np.shape(flux)[0]):
        j=0
        for j in range (0,np.shape(flux)[1]):
            if Distance[i][j] > sky_annulus_1 and Distance[i][j] < sky_annulus_2:
                sky_mask[i,j] = 0
    
    Photo = np.sum(np.ma.array(flux,mask=Photometry_mask))
    SKY = np.ma.median(np.ma.array(flux,mask=sky_mask))
    return Photo,SKY,Photometry_mask,sky_mask
