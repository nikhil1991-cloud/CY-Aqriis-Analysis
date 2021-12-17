from subprocess import call
from matplotlib import pyplot as plt
from astropy.io import fits
import numpy as np
import os
os.chdir('/Users/nikhil/Data/CY_Aqr/Light_Frame')


with open('/Users/nikhil/Data/CY_Aqr/ALL_files.txt') as f:
    Line = [line.rstrip('\n') for line in open('/Users/nikhil/Data/CY_Aqr/ALL_files.txt')]

q=0
for q in range (0,np.shape(Line)[0]):
    call(['sex','/Users/nikhil/Data/CY_Aqr/Light_Frame/'+str(Line[q])])
    call(['mv','segmentation.fits','/Users/nikhil/Data/CY_Aqr/Seg/'+str(Line[q])])
    call(['rm','apertures.fits'])
    call(['rm','cat.fits'])
    call(['rm','objects.fits'])
    
    hdu = fits.open('/Users/nikhil/Data/CY_Aqr/Seg/'+str(Line[q]))
    seg = hdu[0].data
    
    fig = plt.figure()
    plt.imshow(seg)
    fig.set_size_inches(70,50)
    fig.savefig('/Users/nikhil/Data/CY_Aqr/Seg_images/'+str(Line[q])+'.png')
    plt.close()
