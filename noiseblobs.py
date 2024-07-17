## written by Sarah Betti 2019
## updated 5 November 2019

## --- AL NOTE: wherever has fil = f1.split('.nc')[0], etc.
# changed to fil = f1.split('_signal.fits')[0] as a hack, to not require keeping .nc files
# on local computer 

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.wcs import WCS

import time
import os
import glob

from getBlobs import *

import re 

def split_multiext_fits(noisepath):
    # AL: Split the TolTEC multi-extension noise file into individual noise files to follow
    #     old method from Sarah 
    
    # e.g., fitsfile: e.g. toltec_simu_a1100_noise_filtered.fits
    
    # split the multiextension noise signal files into individual ones
    # previous naming scheme: noiseXX_signal.fits
    #                         noiseXX_weight.fits
    #                         noiseXX_psf.fits
    #                         noiseXX_sn.fits
    
    multiext = fits.open(noisepath) # noisepath is defined in driver.py and is the fitsfile 
    multiext = multiext[1:] # skip the first row of HDUList which is just PRIMARY (not a noise extension)

    for ii in range(len(multiext)):
        extname = multiext[ii].header['EXTNAME'] # each EXTNAME is of form: signal_XX_I 


        # TOLTEC has mJy/beam, AzTEC was Jy/beam. Rescale it to Jy/beam
        
        data = multiext[ii].data  * 1E-3         # get the noise data (mJy/beam --> Jy/beam)
        header = multiext[ii].header             # get the header 
        header['UNIT'] = 'Jy/beam'              # rewrite brightness unit as Jy/beam
        # to indicate that a conversion from mJy/beam to Jy/beam was done: 
        header.comments['UNIT'] = 'Conversion from Jy/beam to mJy/beam was applied in noiseblobs.py'
        # get the number index XX                # in noiseblobs.py
        match = re.search(r'\d+', extname)
        if match:
            numindex = match.group() # returns numbering index as a string

        # use the same naming scheme as AzTEC, i.e., noiseXX_signal.fits 
        fitsname= 'noise' + numindex + '_signal.fits'

        if os.path.isfile(os.path.dirname(noisepath) + '/' + fitsname) == False:  
        
            fits.writeto(os.path.dirname(noisepath) + '/' + fitsname, data, header, overwrite=True)
        # --> AL: for this, I need the fitsfile to be included in the noisepath, unlike the version on Sarah's github
        
        
def setup_NoiseBlobs(ncfile, signalheader, psfdata, weightdata, weightheader, psfheader):
    #go through noise.nc files and create fits files
    f1 = ncfile ## now is fits 

    '''
    if (os.path.isfile(f1.split('.nc')[0] + '_weight.fits') == False) or (os.path.isfile(f1.split('.nc')[0] + '_signal.fits') == False):
        print(f'creating signal, weight, sn noise maps for noise map: {f1.split("/")[-1]}')
        os.system('fitswriter ' + f1)
    ''' 
    
    ## ========  write the weight, sn and psf files ========= 
    
    if os.path.isfile(f1.split('_signal.fits')[0] + '_psf.fits') == False:
        fil = f1.split('_signal.fits')[0]
        print('Writing PSF file: {}'.format(fil + '_psf.fits'))
        # == DID NOT rescale PSF from mJy/beam to Jy/beam since it is just for calculation of FWHM
        fits.writeto(fil + '_psf.fits', psfdata, psfheader, overwrite=True)

    # AL: write a weight file for each noise map (use weight maps from science maps for now)
    if os.path.isfile(f1.split('_signal.fits')[0] + '_weight.fits') == False:
        fil = f1.split('_signal.fits')[0]
        print('Writing weight file: {}'.format(fil + '_weight.fits'))

        # == weightdata has already been rescaled to 1/(Jy/beam)^2 at this point 
        fits.writeto(fil + '_weight.fits', weightdata, weightheader, overwrite=True)

    # AL: now create the s/n files for noise realizations with the weight and signal 
    if os.path.isfile(f1.split('_signal.fits')[0] + '_sn.fits') == False:  
        # fil = f1.split('.nc')[0]

        fil = f1.split('_signal.fits')[0]

        imwt = fits.getdata(fil + '_weight.fits') 
        imsig = fits.getdata(fil + '_signal.fits')

        if len(imsig.shape)  != 2:     imsig  = imsig[0,0,:,:]
        
        imsn = imsig*np.sqrt(imwt)
        
        #determine if any weights are less than 1.  if so, make the S/N at that location 0 and not negative. 
        snheader = signalheader
        snheader['UNIT'] = '' # no unit, since it is SN map
        
        w0 = np.where(imwt < 1.)
        n0 = len(w0[0])
        if n0 > 0:
            imsn[w0] = 0.
        print('Writing SN file: {}'.format(fil + '_sn.fits'))
        fits.writeto(fil + '_sn.fits', imsn, snheader, overwrite=True)
        
def noiseMaps(noisepath, psf, signal, weight):
    '''
    find all noise realization maps  
    
    Parameters
    ----------
    path : str
        location of coadd.nc
        
    Returns
    ---------
    globbed list of path location of all noise realization maps within /noise_maps/ folder.  
    
    '''

    # AL: split the multiextension noise file into individual noise files, following previous naming scheme and output them in the same directory as noisepath 
    split_multiext_fits(noisepath)
    
    f = glob.glob(os.path.dirname(noisepath) + '/noise*_signal.fits')    
    #print(f)

    #f = glob.glob(noisepath + '/noise*.nc') 

    ##### might parse out wrong number - maybe comment out next 3 lines 
    if len(f) > 0:
        ## reorder files so that they are in a natural/human ordering
        noise_number = [int(re.findall('\d+', i.split('/')[-1])[0]) for i in f] #pull out integer number of noise map 
        noise_num_ind = np.argsort(noise_number)  #reorder based on noisemap value
        f = np.array(f)[noise_num_ind] #sort initial list based on reordering
        
        # get coadded signal and psf

    if signal != psf: 
        psfdata = fits.getdata(psf) # -- THIS NEEDS TO BE CHANGED TO MULTIEXT FILE
        print(psf) # --> 8/14/23 is getting the signal path for some reason 
        signalheader = fits.getheader(signal) 
        print(signal) # --> 8/14/23 is getting the psf path for some reason 

    else: # AL: TolTEC sim support
        print(signal)
    # -- get the weight and psf (just take them from science map for now)
        hdu = fits.open(signal)
        signalheader = hdu['signal_I'].header
        psfdata = hdu['kernel_I'].data ### kernel is in mJy/beam but no rescaling since is for FWHM calculation. 
        psfheader = hdu['kernel_I'].header 
        # RESCALED WEIGHT to 1/(Jy/beam)^2 
        weightdata = hdu['weight_I'].data * 1E6 # 1/(mJy/beam)^2 --> 1/(Jy/beam)^2
        weightheader = hdu['weight_I'].header
        
        # to note that conversion from 1/(mJy/beam)^2 to 1/(Jy/beam)^2 was applied
        weightheader['UNIT'] = '1/(Jy/beam)^2'
        weightheader.comments['UNIT'] = 'Conversion from 1/(mJy/beam)^2 to 1/(Jy/beam)^2 was applied in noiseblobs.py' 

        if len(weightdata.shape)  != 2:     weightdata  = weightdata[0,0,:,:]
        if len(psfdata.shape)     != 2:     psfdata     = psfdata[0,0,:,:]

    for ncfile in f:
        # add in the psf header, since i do not change the units for it. 
        setup_NoiseBlobs(ncfile, signalheader, psfdata, weightdata, weightheader, psfheader)
        
    return f


def blobListSav(outbase):
    '''
    find all noise bloblist sav file
        
    Returns
    ---------
    globbed list of path location of all noise realization bloblist save files within /noise_maps/ folder.  
    '''
    f = glob.glob(outbase + '/bloblist_noise*.pkl')
    
    ## reorder files so that they are in a natural/human ordering
    noise_number = [int(re.findall('\d+', i.split('/')[-1])[0]) for i in f] #pull out integer number of noise map 
    noise_num_ind = np.argsort(noise_number)  #reorder based on noisemap value
    f = np.array(f)[noise_num_ind] #sort initial list based on reordering
    
    return f


def noise_maps(noisepath):
    signal = noisepath +  '_signal.fits'
    weight = noisepath + '_weight.fits'
    sn = noisepath + '_sn.fits'
    psf = noisepath +  '_psf.fits'

    #print(signal, weight, sn, psf)
    return signal, weight, sn, psf

        
#def runNoiseBlobs(pcdist, noisepath, signal, psf, outbase, sn_thresh, sim, cleanup=False):
def runNoiseBlobs(pcdist, noisepath, weight, signal, psf, outbase, sn_thresh, sim, cleanup=False):
    '''
    determine noise blobs in each noise map by running each noise .nc map through self.runMacanaBlobs() and save resulting properties in pickled pandas dataframe 
    
    Parameters
    ----------
    coaddpath : str
        location of coadd.nc
    noisepath : str
        location of noise*.nc files
    pcdist : float
        distance to source in pc
    cleanup : boolean 
        if True - delete all .fits files for noise*.nc 
        default is False
    
            
    Returns
    ---------
    saves pickled pandas dataframe of all noise core properties
    '''
    
    #get noise*.nc fies
    
    ######## f1 = noiseMaps(noisepath, signal, psf)
    # -->8/14/23  AL: THIS IS THE PROBLEM!!! The order of signal and psf are swapped here, compared
    # to the order that needs to be called for def noiseMaps()

    #f1 = noiseMaps(noisepath, psf, signal) # AL edit, this is now correct.
    f1 = noiseMaps(noisepath, psf, signal, weight) # AL edit, this is now correct. 
    
    #run macana_blobs on noise images
    for i in np.arange(len(f1)):
        ## AL CHANGED
        #fil = f1[i].split('.nc')[0]

        fil = f1[i].split('_signal.fits')[0]
        print(fil)
        
        #-- AL CHECK
        #print(fil) #  '/Users/amandalee/research/c2c/test/data/sokol19_field09/noise/noise44'
        num = fil.split('/noise')[-1]
        # -- AL CHECK
        #print(num, type(num)) 
        noisesignal, noiseweight, noisesn, noisepsf = noise_maps(fil)
        
       
        # -- AL: if bloblist_noiseXX.pkl file is false, then do noiseblob finding 
        if os.path.isfile(outbase + '/bloblist_noise' + str(num) + '.pkl') == False: 

            image_noise_sn = fits.getdata(noisesn)
            good = np.where(image_noise_sn >= sn_thresh)
            npix = int((206265*0.04)/pcdist)

            # AL: if there are "good" pixels, make sure there are enough for at least one blob
            if len(good[0]) >= npix: 
                gb = getBlobs('noise'+num,pcdist,noisesignal,noiseweight,noisesn,noisepsf,outbase,sim=sim)
                gb.npixels=npix

                gb.getBlobCatalogue(sn_thresh, ysoradec=None, plot_blobs=False)
            else:
                print('WARNING: there are not any pixels with sn >= {} in noise{}_sn.fits for noiseblob identification. No noise blob pickle file wias be made.'.format(sn_thresh, num))
                print()
                continue 
            '''    
            # -- create object of getBlobs class 
            gb = getBlobs('noise'+num, pcdist, noisesignal,noiseweight,noisesn,noisepsf,outbase, sim=sim)
            
            npix = int((206265*0.04)/pcdist)
            gb.npixels = npix

            # -- get blobs within noise realization maps, based on sn_thresh
            gb.getBlobCatalogue(sn_thresh, ysoradec = None, plot_blobs=False)
            print()   
            '''
        else:
            print('noise map ' + str(num) + ' have already been segmented: output files in: ' + outbase + '/bloblist_noise' + str(num) + '.pkl')

        if cleanup == True:
            os.system('rm ' + fil + '_*.fits')

    #open pickle files from macana blobs and combine into one noise blobs file
    f = blobListSav(outbase)
    print()
    print(f'number of bloblist.pkl noise files: {len(np.array(f))}')

    # if there are noiseblob pickle files, then make a total noiseblob list ("bloblist_total_noise.pkl")
    if len(np.array(f)) != 0: 

        total_noise = []
        for i in f:
            c = pd.read_pickle(i)
            total_noise.append(c)
        total_noise = pd.concat(total_noise, axis=0)
    
        df_total_noise = pd.DataFrame(data = total_noise)
        df_total_noise = df_total_noise.reset_index()

        # -- numnoise then is just the number of bloblist_noise*.pkl files included in the
        # total bloblist noise file 
        df_total_noise['numnoise'] = np.ones_like(df_total_noise['ra'])*len(f)
    
        df_total_noise.to_pickle(outbase + '/bloblist_total_noise.pkl')
    
        print(f'saved final noisebloblist to: {outbase}/bloblist_total_noise.pkl')

    else:
        print('There were not any bloblist.pkl noise files. No bloblist_total_noise.pkl file was made.')
    '''    
    total_noise = []
    for i in f:
        c = pd.read_pickle(i)
        total_noise.append(c)
    total_noise = pd.concat(total_noise, axis=0)
    
    df_total_noise = pd.DataFrame(data = total_noise)
    df_total_noise = df_total_noise.reset_index()

    # -- numnoise then is just the number of bloblist_noise*.pkl files included in the
    # total bloblist noise file 
    df_total_noise['numnoise'] = np.ones_like(df_total_noise['ra'])*len(f)
    
    df_total_noise.to_pickle(outbase + '/bloblist_total_noise.pkl')
    
    print(f'saved final noisebloblist to: {outbase}/bloblist_total_noise.pkl')
    ''' 
    
#def noisebloblist_name():
    #return 'bloblist_total_noise.pkl'

def noisebloblist_name(outbase):
    if os.path.isfile(outbase + '/bloblist_total_noise' + '.pkl') == False:
        return None
    else:
        return 'bloblist_total_noise.pkl'     
        
    
    
