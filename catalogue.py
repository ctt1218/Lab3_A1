#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 19:41:03 2020

@author: Charlie
"""

from astropy.io import fits,ascii
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt

hdulist = fits.open("masked_image.fits")
image_pixels = hdulist[0].data
masking = image_pixels == 0
masked = np.ma.masked_array(image_pixels,masking)
ylim,xlim = masked.data.shape
maxp_lim = 3475     
bkg_cutoff = 3457   # 3 standard deviations from mean

galaxynos = []          # galaxy number 
xpos = []             # x coord
ypos = []             # y coord
apsizes = []          # aperture radius
galcounts = []      # cumulative galaxy counts
galcounts_errs = []      # error on galcounts
backgrounds = []      # background around galaxy
brightness = []       # magnitude of brightness for each galaxy
mag_errs = []         # error on brightness

def make_plot(brightness,mag_errs):
    logNs = []     
    Ns = []        
    mags = np.linspace(min(brightness) + 0.5 ,max(brightness) + 3,50) 
    for m in mags:
        logN = np.log10(sum(i < m for i in brightness))
        N = sum(i < m for i in brightness)
        logNs.append(logN)
        Ns.append(N)

    logN_errs = [1/((np.log(10)*(n**(0.5)))) for n in Ns]
    plt.errorbar(mags,logNs,yerr = logN_errs, fmt = '.', elinewidth = 0.5, capsize = 2)
    gradient = (logNs[20] - logNs[9])/(mags[20] - mags[9])
    print(gradient)
    plt.plot(mags[8:22], [(gradient*mag) - 1.7 for mag in mags[8:22]], 'g--')
    plt.xlabel('magnitude')
    plt.ylabel('$log_{10}$N(<m)')
    plt.show() #plot of magnitude-count relationship

def generate_catalogue(): #produces catalogue of galaxies in masked image
    new_data = masked[~masked.masking].data
    galaxyno = 0                                 
    maxpix = np.amax(new_data)                   # brightest pixel in image
    magzpt = 25.3                              
    magzrr = 0.02                              
    
    while maxpix > maxp_lim:
        print(maxpix)
        new_data = masked[~masked.masking].data 
        maxpix = np.amax(new_data)      # brightest pixel
        y0,x0 = np.where(masked.data == maxpix)[0][0], np.where(masked.data == maxpix)[1][0] # location of pixel
        print(galaxyno)
        aperturesize=1   
        upper_apsize = aperturesize + 10                 
        if aperturesize > 0:    
            if masked.masking[y0,x0] == True:                   # prevents redetection of masked galaxies 
                masked.data[y0,x0] = 0
            else:
                galaxyno += 1                                # add to galaxy count
                pixels = []                               
                background_vals = []                         
                masked.masking[y0,x0] = True                  # mask galaxy         
                xi,xf,yi,yf = aperturesize    
                for y in range(yi,yf): 
                    for x in range(xi,xf): 
                        r = ((x-x0)**2 + (y-y0)**2)**0.5      
                        if (aperturesize < r < upper_apsize):
                            background_vals.append(masked.data[y,x])  
                        if r < aperturesize:
                            pixels.append(masked.data[y,x])   
                            masked.masking[y,x] = True          
                            masked.data[y,x] = 0
                        
                
                background_vals = [b for b in background_vals if bkg_cutoff > b > 0]
                background = np.mean(background_vals)         # background at area of galaxy
                       
                pixels = [p - background for p in pixels]  # subtract background values 
                pixels = [p for p in pixels if p > 0]  
                counts = np.sum(pixels)               
                c_err = counts**0.5                  # error on galaxy counts
                m = magzpt - 2.5*np.log10(counts)    # brightness 
                log_c_err = c_err/(counts*np.log(10)) #error on brightness
                m_err = m*((magzrr/magzpt)**2 + (log_c_err/np.log10(counts))**2)**0.5 
                
                galaxynos.append(galaxyno)
                xpos.append(x0)
                ypos.append(y0)
                apsizes.append(aperturesize)
                galcounts.append(counts)
                galcounts_errs.append(c_err)
                backgrounds.append(background)
                brightness.append(m)
                mag_errs.append(m_err)
        else:
            masked.data[y0,x0] = 0
    
    hdulist[0].data = image_pixels 
    data = Table([galaxynos,xpos,ypos,galcounts,galcounts_errs,backgrounds,brightness,mag_errs], names = ['galaxy number', 'x', 'y', 'aperture radius', 'galaxy type', 'counts', 'counts error', 'background', 'magnitude', 'magnitude error'])
    ascii.write(data, 'catalogue1.dat') 
    hdulist.writeto('finalmaskedmosaic.fits')       
   
generate_catalogue()
make_plot(brightness,mag_errs)