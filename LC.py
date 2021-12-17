import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
import os
os.chdir('/Users/nikhil/code/MacAdam Obs')
from photo import get_flux_sky #import aperture photometry function from photo.py

#set up fontsize and padding for plots
ft_12 =100
ft_22 =100
ft_1 = 100
lw=10
ltw =110
ltwm = 95
pd=60

#read all the light frames
with open('/Users/nikhil/Data/CY_Aqr/ALL_files.txt') as f:
    Line = [line.rstrip('\n') for line in open('/Users/nikhil/Data/CY_Aqr/ALL_files.txt')]

#read all dark frames
with open('/Users/nikhil/Data/CY_Aqr/Dark.txt') as f:
    Line_dark = [line.rstrip('\n') for line in open('/Users/nikhil/Data/CY_Aqr/Dark.txt')]
    
#read all flat frames
with open('/Users/nikhil/Data/CY_Aqr/flat.txt') as f:
    Line_flat = [line.rstrip('\n') for line in open('/Users/nikhil/Data/CY_Aqr/flat.txt')]
    
ALL_dark_frames = np.zeros((1023,1536,len(Line_dark)))
ALL_flat_frames = np.zeros((1023,1536,len(Line_flat)))
dark=0
for dark in range (0,len(Line_dark)):
    dark_hdu = fits.open('/Users/nikhil/Data/CY_Aqr/Dark_Frame/'+str(Line_dark[dark]))
    ALL_dark_frames[:,:,dark] = dark_hdu[0].data

flat=0
for flat in range (0,len(Line_flat)):
    flat_hdu = fits.open('/Users/nikhil/Data/CY_Aqr/Flat_Frame/'+str(Line_flat[dark]))
    ALL_flat_frames[:,:,flat] = flat_hdu[0].data

#take the median of all dark and flat frames along the rd axis
Median_Dark = np.median(ALL_dark_frames,axis=2)
Median_Flat = np.median(ALL_flat_frames,axis=2)

#initiate arrays for CY_Aqr and Control star magnitudes
Col_CY = np.zeros(np.shape(Line)) #CY_Aqr - Control_star color in magnitudes
CYmag = np.zeros(np.shape(Line)) #CY_Aqr magnitude
Conmag = np.zeros(np.shape(Line)) #Control star magnitude

q=0
for q in range (0,np.shape(Line)[0]):
    hdu = fits.open('/Users/nikhil/Data/CY_Aqr/Light_Frame/'+str(Line[q])) #read the light frame
    exposure_sec = hdu[0].header['EXPTIME'] #get exposure time in sec
    flux = hdu[0].data #read flux
    hdu_seg = fits.open('/Users/nikhil/Data/CY_Aqr/Seg/'+str(Line[q])) #read the segmented mask for the brightest star
    seg = hdu_seg[0].data
    K_prime = np.unique(seg)
    idx_sky = np.where(K_prime==0)
    sky_c = np.median(flux[idx_sky])
    K = np.delete(K_prime,0)
    sum_all = np.zeros(np.shape(K))
    for num in range (0,len(K)):
        idx = np.where(seg==K[num])
        sum_all[num] = np.sum(flux[idx])
        
    Bright_mask = K[np.argsort(sum_all)][-1]
    Bright2_mask = K[np.argsort(sum_all)][-2]
    KMASK = np.ones(np.shape(seg))
    Pos_Brightest = int(np.median(np.where(seg==Bright_mask)[0])),int(np.median(np.where(seg==Bright_mask)[1]))
    Pos_2nd_Brightest = int(np.median(np.where(seg==Bright2_mask)[0])),int(np.median(np.where(seg==Bright2_mask)[1]))
    if Pos_Brightest[1] >400 and Pos_Brightest[0]>400:
       Pos_ref = Pos_2nd_Brightest
    else:
       Pos_ref = Pos_Brightest

    #Locate CY_Aqr and Control star w.r.t to the brightest star
    Pos_CY = Pos_ref[0] + 299, Pos_ref[1] + 521
    Pos_Con = Pos_ref[0] + 273,Pos_ref[1] + 246
    
    Calibrated_flux = (flux-Median_Dark)/Median_Flat
    
    #calculate Photometric flux, Sky counts and masks for CY_Aqr and control star
    CY_flux,CY_sky,CY_mask,CY_Smask = get_flux_sky(Calibrated_flux,Pos_CY,10,20,50)
    Con_flux,Con_sky,Con_mask,Con_smask = get_flux_sky(Calibrated_flux,Pos_Con,10,20,50)
    mag_CY = -2.5*np.log10((CY_flux - CY_sky)/exposure_sec)
    mag_Con = -2.5*np.log10((Con_flux - Con_sky)/exposure_sec)
    CYmag[q] = mag_CY
    Conmag[q] = mag_Con
    Col_CY[q] = mag_CY - mag_Con
    print(Line[q])
    print(Col_CY[q])
    

    #Plot the masked CY_Aqr and Control star with their Sky annulus
    fig=plt.figure()
    plt.imshow(np.ma.array(flux,mask=Con_mask))
    plt.imshow(np.ma.array(flux,mask=Con_smask))
    plt.imshow(np.ma.array(flux,mask=CY_mask))
    plt.imshow(np.ma.array(flux,mask=CY_Smask))
    plt.minorticks_on()
    plt.text(Pos_CY[1],Pos_CY[0]+100,'CY Aqr',fontsize=ft_12)
    plt.text(Pos_Con[1],Pos_Con[0]+100,'Control Star',fontsize=ft_12)
    plt.xlabel('X',fontsize=ft_12)
    plt.ylabel('Y',fontsize=ft_12)
    plt.tick_params(which='major',length=ltw,width=lw,direction='in',labelsize=ft_1,pad=pd,axis='x',bottom=True,top=True,left=True,right=True)
    plt.tick_params(which='minor',length=ltwm,width=lw,direction='in',labelsize=ft_1,pad=pd,axis='x',bottom=True,top=True,left=True,right=True)
    plt.tick_params(which='major',length=ltw,width=lw,direction='in',labelsize=ft_1,pad=pd,axis='y',bottom=True,top=True,left=True,right=True)
    plt.tick_params(which='minor',length=ltwm,width=lw,direction='in',labelsize=ft_1,pad=pd,axis='y',bottom=True,top=True,left=True,right=True)
    for axis in ['top','bottom','left','right']:
        plt.gca().spines[axis].set_linewidth(lw)
    
    fig.set_size_inches(70,50)
    fig.tight_layout()
    fig.savefig('/Users/nikhil/Data/CY_Aqr/Seg_images/'+str(Line[q])+'.png')
    plt.close()
    
    #Plot the flux image for comparison
    fig1=plt.figure()
    plt.imshow(flux,vmax=5000)
    plt.minorticks_on()
    plt.text(Pos_CY[1],Pos_CY[0]+100,'CY Aqr',fontsize=ft_12,c='white')
    plt.text(Pos_Con[1],Pos_Con[0]+100,'Control Star',fontsize=ft_12,c='white')
    plt.xlabel('X',fontsize=ft_12)
    plt.ylabel('Y',fontsize=ft_12)
    plt.tick_params(which='major',length=ltw,width=lw,direction='in',labelsize=ft_1,pad=pd,axis='x',bottom=True,top=True,left=True,right=True)
    plt.tick_params(which='minor',length=ltwm,width=lw,direction='in',labelsize=ft_1,pad=pd,axis='x',bottom=True,top=True,left=True,right=True)
    plt.tick_params(which='major',length=ltw,width=lw,direction='in',labelsize=ft_1,pad=pd,axis='y',bottom=True,top=True,left=True,right=True)
    plt.tick_params(which='minor',length=ltwm,width=lw,direction='in',labelsize=ft_1,pad=pd,axis='y',bottom=True,top=True,left=True,right=True)
    for axis in ['top','bottom','left','right']:
        plt.gca().spines[axis].set_linewidth(lw)
    
    fig1.set_size_inches(70,50)
    fig1.tight_layout()
    fig1.savefig('/Users/nikhil/Data/CY_Aqr/Seg_images/'+str(Line[q])+'_flux.png')
    plt.close()


#Convert the Date time to julian date from this csv file
DAT = np.array(pd.read_csv('/Users/nikhil/Data/CY_Aqr/Datetime_JD.csv'))
JD = DAT[:,3].astype('float64')
#Store the lightcurve in fits format
hdu0 = fits.PrimaryHDU(Col_CY)
hdu1 = fits.ImageHDU(CYmag)
hdu2 = fits.ImageHDU(Conmag)
hdu3 = fits.ImageHDU(JD)
new_hdul = fits.HDUList([hdu0,hdu1,hdu2,hdu3])
new_hdul.writeto('/Users/nikhil/Data/CY_Aqr/CY_Aqr_lightcurve.fits')
