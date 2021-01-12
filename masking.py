#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 22:54:07 2020

@author: Charlie
"""



#f = open(os.path.expanduser("~/Documents/mosaic.fits"))



from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import time

hdulist=fits.open(r"/Users/Charlie/Documents/original_mosaic.fits")
image_pixels = hdulist[0].data

start = time.time()        

def mask(image_pixels):    #used to remove unwanted details from image

    for y in range(image_pixels.shape[0] ): 
        for x in range(image_pixels.shape[1]):
            
            r1 = ((x-1430)**2 + (y-3200)**2)**0.5
            r2 = ((x-1315)**2 + (y-4397)**2)**0.5
            r3 = ((x-1365)**2 + (y-4332)**2)**0.5
            r4 = ((x-560)**2 + (y-4098)**2)**0.5
            r5 = ((x-1461)**2 + (y-4032)**2)**0.5
            r6 = ((x-2133)**2 + (y-3760)**2)**0.5
            r7 = ((x-2465)**2 + (y-3417)**2)**0.5
            r8 = ((x-779)**2 + (y-3321)**2)**0.5
            r9 = ((x-976)**2 + (y-2773)**2)**0.5
            r10 = ((x-905)**2 + (y-2289)**2)**0.5
            r11 = ((x-2131)**2 + (y-2309)**2)**0.5
            r12 = ((x-2089)**2 + (y-1424)**2)**0.5 #circles with coordinates of areas where blooming occurs
            
            if r1 < 320:
                image_pixels[y,x] = 0 
            elif r2 < 25:
                image_pixels[y,x] = 0
            elif r3 < 23:
                image_pixels[y,x] = 0
            elif r4 < 33:
                image_pixels[y,x] = 0
            elif r5 < 27:
                image_pixels[y,x] = 0
            elif r6 < 52:
                image_pixels[y,x] = 0
            elif r7 < 33:
                image_pixels[y,x] = 0
            elif r8 < 62:
                image_pixels[y,x] = 0
            elif r9 < 48:
                image_pixels[y,x] = 0
            elif r10 < 47:
                image_pixels[y,x] = 0
            elif r11 < 34:
                image_pixels[y,x] = 0
            elif r12 < 32:
                image_pixels[y,x] = 0
            else:
                pass
    
    image_pixels[0:4611, 1425:1455] = 0
    image_pixels[115:150, 1290:1535] = 0 
    image_pixels[218:267, 1391:1474] = 0 
    image_pixels[312:360, 1015:1711] = 0 
    image_pixels[425:480, 1100:1653] = 0 
    image_pixels[422:454,1025:1045] = 0 
    image_pixels[3204:3420, 772:779] = 0
    image_pixels[2703:2837, 969:977] = 0
    image_pixels[2223:2356, 901:908] = 0 #masking effect
    
    image_pixels[-220:,-410:] = 0
    image_pixels[-410:,-210:] = 0 #removing corner

    image_pixels = np.delete(image_pixels, np.s_[0:100], 0)
    image_pixels = np.delete(image_pixels, np.s_[-100:], 0)
    image_pixels = np.delete(image_pixels, np.s_[0:100], 1)
    image_pixels = np.delete(image_pixels, np.s_[-100:], 1) #removing noisy edges
    
    hdulist[0].data = image_pixels
    plt.figure()   
    hdulist.writeto('masked_image.fits')   # new image with masking applied
    plt.imshow(image_pixels)
    return image_pixels

def histogram(data):
    mu = np.round(np.mean(data),1)
    std = np.round(np.std(data),1)
    def gaussian(mu,std,data):
        return (np.exp(-((data-mu)**2)/(2*std**2)))/(np.sqrt(2*np.pi)*std)
        
    plt.figure()
    x = np.linspace(np.min(data), np.max(data), 100)
    plt.plot(x, gaussian(mu,std,x), 'k', label = 'Gaussian: $\mu$ = %s, $\sigma$ = %s' % (mu, std))
    plt.hist(data,bins=50, color='b', density = True, label = 'Histogram of pixel values in masked image')
    plt.title('Distribution of Background Noise')
    plt.legend()
    plt.xlabel('Pixel Value')
    plt.ylabel('Number')
    plt.xlim(3360,3465)
    plt.show()

image_pixels = mask(image_pixels)
flat_pixels = image_pixels.flatten()
flat_pixels = [d for d in flat_pixels if 3360 < d < 3480] # shorten Gaussian
histogram(flat_pixels)