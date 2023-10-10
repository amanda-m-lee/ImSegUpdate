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

#PATH TO SIGNAL/WEIGHT/SN/PSF MAPS
#signal  = 'coadded_signal.fits'
#weight  = 'coadded_weight.fits'
#sn      = 'coadded_sn.fits'
#psf     = 'coadded_psf.fits'

signal  = '/Users/amandalee/research/c2c/simImSeg/simdata/reduced_sim_100190_lpT_neig3_approx_Gchanges_freqlow0.0/redu00/coadded/filtered/toltec_simu_a1100_filtered.fits'
weight  = '/Users/amandalee/research/c2c/simImSeg/simdata/reduced_sim_100190_lpT_neig3_approx_Gchanges_freqlow0.0/redu00/coadded/filtered/toltec_simu_a1100_filtered.fits'
sn      = '/Users/amandalee/research/c2c/simImSeg/simdata/reduced_sim_100190_lpT_neig3_approx_Gchanges_freqlow0.0/redu00/coadded/filtered/toltec_simu_a1100_filtered.fits'
psf     = '/Users/amandalee/research/c2c/simImSeg/simdata/reduced_sim_100190_lpT_neig3_approx_Gchanges_freqlow0.0/redu00/coadded/filtered/toltec_simu_a1100_filtered.fits'

#PATH TO NOISE REALIZATION FILES.  DO NOT INCLUDE NAMES OF FILES!
#noisepath = '/path/to/noise/realization/' OR None
noisepath = '/Users/amandalee/research/c2c/simImSeg/simdata/reduced_sim_100190_lpT_neig3_approx_Gchanges_freqlow0.0/redu00/coadded/filtered/toltec_simu_a1100_noise_filtered.fits'
 
#OBSERVATION/TARGET INFORMATION
#wavelength = wavelength of observations in microns
#pcdist = distance to target in parsecs       
#sn_thresh = S/N threshold for image segmenataion
#sim = define whether or not maps are simulations (True) or observations (False)
wavelength = 1100
pcdist = 140 # 0s_toltecsim_M2e3_mid_140.0.fits
sn_thresh = 2.5
sim = True #False AL: I checked and whether this is False/True does not influence the analysis (just labeling)

#PATH TO TEMPERATURE/COLUMN DENSITY/YSO LISTS
#path to temperature maps
#temp = 'temperature.fits'

## -- edited imageSegment.py to take in temp from the multi-ext fitsfile
## also, temp map does not have WCS, so assuming it has the same wcs as the coldens map 

#### USE MY INTERPOLATED TEMP MAP
temp = '/Users/amandalee/research/c2c/test/data/herschel/interpolated_monr2_temp.fits' 

#path to column density map

## -- edited coreSelection.py to take in coldens from the multi-ext fitsfile

### USE my interpolated version of column density map to see if now rerunning coreselection
### captures back that high SN core!
coldens = '/Users/amandalee/research/c2c/test/data/herschel/interpolated_monr2_colden.fits'

#load YSO list
#ysoradec = np.loadtxt('YSOs.txt')
# -- see the IPAC catalog format here:
# -- https://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html

# ysoradec = array([[ra1, dec1], [ra2, dec2], ... [ran, decn]])

#ysoradec = np.loadtxt('/Users/amandalee/research/c2c/test/data/Mon R2_catalog_ipac1.txt', skiprows=5, usecols=(4,5)) ## ra, dec

monr2_cat = np.loadtxt('/Users/amandalee/research/c2c/test/data/Mon R2_catalog_ipac1.txt', skiprows=5, usecols=(4,5,48)) ## ra, dec, class

good = np.where( (monr2_cat[:,2] >= 0) & # get only information for class 0, 1, 2, 3, YSOs
                 (monr2_cat[:,2] < 4) )

ysoradec = monr2_cat[:,[0,1]][good] # access only the (RA, DEC) columns 

#LOCATION TO SAVE OUTPUT FILES
#outbase='/path/to/output/files' 
outbase = '/Users/amandalee/research/c2c/simImSeg/outbase'

#PARTS OF REDUCTION TO RUN
run_getblobs = True #True
run_getnoiseblobs = False #True 
run_coreselection = False#True

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
    mb.nsigma = 5 #5
    mb.nmxsigma = 3.5 #3.5
    mb.areacut = -100
    mb.covlim = 0.4
    mb.contrast = 0.001
    #set number of connected pixels based on distance to target
    npix = int((206265*0.04)/pcdist)
    mb.npixels = npix 

    mb.getBlobCatalogue(sn_thresh,ysoradec, plot_blobs=True)
    print('---------------o---------------\n')
    
if run_getnoiseblobs:
    print('starting getNoiseBlobs\n')
    runNoiseBlobs(pcdist, noisepath, signal, psf, outbase, sn_thresh,sim, cleanup=False)
    print('---------------o---------------\n')
    
fwhm_beam = mb.FWHM() # -- gets the same beam calculated from the psf from science obs.
#print(fwhm_beam) # --> 12.100518486600036 
pickle_file = mb.bloblist_name()
noise_pickle_file = noisebloblist_name()
cs = coreSelection(outbase, wavelength, pcdist, coldens, pickle_file, noise_pickle_file=noise_pickle_file)
if run_coreselection:
    print('starting runCoreSelection\n')
    good_indices, goodmap = cs.runCoreSelection()
    cs.saveGoodCores(fwhm_beam, good_indices=good_indices)
    cs.plot_ratio('sim', sim=sim, good_indices=good_indices, goodmap=goodmap)
    ###cs.plot_CMF()
    print('---------------o---------------\n')
    
t2 = time.time()
print('finished total time: ', (t2-t1)/60., ' min')
