## written by Sarah Betti 2019
## updated 5 November 2019

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.wcs import WCS

import glob
import os
import time
import re

from getrot import *
from coreSelection import *
from getBlobs import *
from noiseblobs import *
plt.ion()

############## CHANGEABLE PARAMETERS ##################

#PATH TO SIGNAL/WEIGHT/SN/PSF MAPS (for TolTEC data product, everything is included in one file) 

# Signal 
signal = 'toltec_simu_a1100_citlali.fits' 

# Weight, SN, PSF 
weight = signal
sn = signal
psf = signal 

# Noise realization 
noisepath = 'toltec_simu_a1100_noise_citlali.fits'

#OBSERVATION/TARGET INFORMATION
#wavelength = wavelength of observations in microns
#pcdist = distance to target in parsecs       
#sn_thresh = S/N threshold for image segmentation
#sim = define whether or not maps are simulations (True) or observations (False)

wavelength = 1100 # [um] 
pcdist = 860 # distance of MonR2 [pc] 
sn_thresh = 2.5 
sim = True 

#PATH TO TEMPERATURE/COLUMN DENSITY/YSO LISTS

#path to temperature maps, e.g.: temp = temperature.fits or temp= None 
temp = None 

#path to column density map 
coldens  = '0s_NH2_snapshot_050_860_smo36_padded.fits'

# path to YSO catalog or None 
ysoradec = None 

#LOCATION TO SAVE OUTPUT FILES
#outbase='/path/to/output/files' 
outbase = './outbase'

#PARTS OF REDUCTION TO RUN (can be True or False for each step)
run_getblobs = True 
run_getnoiseblobs = False #True
run_coreselection = True  

############### END CHANGEABLE PARAMETERS ########################

############### RUN IMSEG ########################
#make output file folder
if not os.path.exists(outbase):
    os.makedirs(outbase)

t1 = time.time()

print('starting core extraction')
print(f'path to signal files: {signal}')
print(f'path to noise files: {noisepath}')  
print(f'path to output files: {outbase} \n')
print('---------------o---------------\n')

mb = getBlobs('cores', pcdist, signal, weight, sn, psf, outbase, temp=temp, temp_hdr=temp, sim=sim)

if run_getblobs: 
    print('starting getBlobs\n')
    mb.nsigma = 5 
    mb.nmxsigma = 3.5 
    mb.areacut = -100
    mb.covlim = 0.4
    mb.contrast = 0.001
    # set number of connected pixels based on distance to target
    npix = int((206265*0.04)/pcdist)
    mb.npixels = npix 

    mb.getBlobCatalogue(sn_thresh,ysoradec, plot_blobs=True)
    print('---------------o---------------\n')
    
if run_getnoiseblobs:
    print('starting getNoiseBlobs\n')
    runNoiseBlobs(pcdist, noisepath, signal, weight, psf, outbase, sn_thresh,sim, cleanup=False)
    print('---------------o---------------\n')
    
fwhm_beam = mb.FWHM() 
#print(fwhm_beam) 
pickle_file = mb.bloblist_name()
#noise_pickle_file = noisebloblist_name() # --> Need to make this check if there is a noiseblobpkl file or not
noise_pickle_file= noisebloblist_name(outbase)

cs = coreSelection(outbase, wavelength, pcdist, coldens, pickle_file, noise_pickle_file=noise_pickle_file)
if run_coreselection:
    print('starting runCoreSelection\n')
    good_indices, goodmap = cs.runCoreSelection()
    cs.saveGoodCores(fwhm_beam, good_indices=good_indices)
    cs.plot_ratio('sim', sim=sim, good_indices=good_indices, goodmap=goodmap)
    print('---------------o---------------\n')
    
t2 = time.time()
print('finished total time: ', (t2-t1)/60., ' min')
