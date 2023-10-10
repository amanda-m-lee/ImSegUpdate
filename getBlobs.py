## written by Sarah Betti 2019
## updated 5 November 2019

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.wcs import WCS

import time
import os
import glob
import re

from imageSegment import *
from getrot import *


class getBlobs:
    
    nsigma = 5      # sigma cut for cores
    nmxsigma = 3.5  # max sigma cut for peak of cores
    areacut = -100  # area cut
    covlim = 0.4  # coverage limit
    npixels = 10    # number of pixels touching to sbe considered a blob
    contrast=0.001  # fraction of the total source flux to that a local peak to be considered as a separate blob
    
    def __init__(self, name, pcdist, signal, weight, sn, psf, outbase, temp=None, temp_hdr=None, sim=False):
        
        '''
        Purpose 
        ----------
        Finds potential cores within molecular cloud 
        
        Parameters
        ----------
        name : string
            type of .fits file (cores, noise)
        pcdist : string
            distance to cloud
        signal : string
            path to signal map (.fits file)
        weight : string
            path to weight map (.fits file)
        sn : string
            path to S/N map (.fits file)
        psf : string
            path to psf map (.fits file)
        outbase : string
            path to folder to save output files
        sim : boolean (optional)
            whether the maps are from a simulation (True) or observed data (False)
        '''

        self.name       =   name 
        self.pcdist     =   pcdist
        
        self.signal     =   signal
        self.weight     =   weight
        self.sn         =   sn
        self.psf        =   psf
        
        self.outbase    =   outbase
        self.sim        =   sim
        
        if temp is not None:
            self.temp = temp
        if temp_hdr is not None:
            self.temp_hdr = temp_hdr
        if (temp is not None) and (temp_hdr is None):
            self.temp_hdr = temp
            
    
    def signal(self):
        '''signal path'''
        return self.signal
    
    
    def weight(self):
        '''weight path'''
        return self.weight
    
    
    def sn(self):
        '''sn path'''
        return self.sn
    
    
    def psf(self):
        '''psf path'''
        return self.psf

    
    def FWHM(self):
        '''
        calculates the fwhm from the psf 
        
        Parameters
        ----------
        normkernel : boolean
            normkernel = True to calculate the fwhm from the psf
            normkernel = False use a fwhm = 8.5 (for AzTEC on the LMT)
            
        Returns
        ---------
        if self.normkernel = True:
            fwhm calculated from psf 
        if self.normkernel = False     
            fwhm = 8.5  
        '''
        weightPath  = self.weight
        psfPath     = self.psf
            
        if weightPath != psfPath:
            psf         = fits.getdata(psfPath)
            cdelt2  = fits.getheader(weightPath)['CDELT2']       
        else: # AL: TolTEC sim support
            hdu     = fits.open(psfPath)
            psf     = hdu['kernel_I'].data # -- Yes it is Q: is kernel_I considered PSF of observation?
            if len(psf.shape)   != 2:      psf  = psf[0,0,:,:]

            cdelt2  = hdu['weight_I'].header['CDELT2']

        w       = np.where(psf > 0.5 * np.max(psf))
        nw      = len(w[0])
        fwhm    = np.sqrt(nw / np.pi) * 2*cdelt2*3600. #arcsec
        print('FWHM calculated from psf in def FWHM() is:', fwhm)
        return fwhm
    
    
    def setup_ImSeg(self):
        '''
        creates a clump find instance of the signal/psf/weight/sn fits map 
        
        Returns
        ----------
        clump : instance of the clumpFind object for the data
        
        '''
        weightPath  = self.weight
        snPath      = self.sn
        signalPath  = self.signal
           
        # open all fits files
        if signalPath != snPath:
            insig   = fits.getdata(signalPath)
            hdr     = fits.getheader(snPath)
            insn    = fits.getdata(snPath)
            wt      = fits.getdata(weightPath)
            wcs_info= WCS(hdr)
        else: # --AL: for TolTEC sim support 
            hdu     = fits.open(signalPath)
            insig   = hdu['signal_I'].data
            hdr     = hdu['signal_I'].header
            insn    = hdu['sig2noise_I'].data
            wt      = hdu['weight_I'].data
            wcs_info = WCS(hdr).sub(2)
            hdu.close()

        # --AL: for TolTEC sim support
        if len(insn.shape)  != 2:     insn  = insn[0,0,:,:]
        if len(insig.shape) != 2:    insig  = insig[0,0,:,:]
        if len(wt.shape)    != 2:       wt  = wt[0,0,:,:]
            
        #initialize image segmentation 
        clumps = clumpFind(insig, insn, hdr, wt=wt, w=wcs_info)
        return clumps
    
    def run_ImSeg(self, sn_thresh):
        '''
        Finds clumpFind.deblendSources attribute for the data
        
        Returns
        ----------
        clumps : instance of the clumpFind object
        fin_obj : object of the deblended sources 
        fin : deblended sources data (same size as the signal map.  (integer values of each blob, 0 for background)
        '''
        clumps = self.setup_ImSeg()
        #get an object class and array of blobs
        fin_obj, fin = clumps.deblendSources(sn_thresh,covlim=self.covlim, npixels=self.npixels, contrast=self.contrast)
        return clumps, fin_obj, fin
    
    def plot_ImSeg(self, sn_thresh):
        '''
        plot 2D segmentation image found from clumpFind.plotClumps
        
        '''
        clumps, fin_obj, fin = self.run_ImSeg(sn_thresh)
        clumps.plotClumps(fin, self.name, self.outbase)
        
    def properties_ImSeg(self, sn_thresh):
        '''
        Find properties of the blobs using clumpFind.sourceProperties
        
        Returns
        ----------
        cores : pandas dataframe of properties of the blobs
        '''
        clumps, fin_obj, fin = self.run_ImSeg(sn_thresh, self.covlim, self.npixels, self.contrast)
        fwhm = self.FWHM()
        
        if hasattr(self, 'temp'):
            cores = clumps.sourceProperties(fin_obj, fwhm, self.covlim, temp = self.temp, temp_hdr = self.temp_hdr)
        else:
            cores = clumps.sourceProperties(fin_obj, fwhm, self.covlim)
        return cores
        
    def getBlobCatalogue(self, sn_thresh, ysoradec, plot_blobs=True):
        '''
        runs through ImSeg:
            - saves 2D segmentation image to fits file
            - determine properties of all sources
            - determines which sources pass threshold cuts
            - finds YSO matches 
            - saves updated dataframe of core properties to pickle (.pkl)
        This is the main workhorse of the getBlobs class
        '''
        
        print(f'starting segmentation on: {self.name}')
        fwhm = self.FWHM()
        print('FWHM from def FWHM() but called in def getBlobCatalogue() is', fwhm)
        clumps = self.setup_ImSeg()

        t1 = time.time()
        #get an object class and array of blobs
        fin_obj, fin = clumps.deblendSources(sn_thresh,self.covlim, self.npixels, self.contrast)
        t2 = time.time()
        print(f'finished segmenting in : {t2-t1} sec')
        
        #open up header and WCS information to save blob ID image to fits
        hdr = fits.getheader(self.weight)
        
        #plot the blobs and write blob ID image to fits file
        if plot_blobs:
            outline_data = clumps.plotClumps(fin, self.name, self.outbase)
            fits.writeto(self.outbase + '/' + self.name + '_blobprints.fits', outline_data, hdr, overwrite=True)  
            
            print(f'saved bloblist fits to: {self.outbase}/{self.name}_blobprints.fits')

        
        if hasattr(self, 'temp'):
            if self.temp_hdr is None:
                self.temp_hdr = self.temp
            
            cores = clumps.sourceProperties(fin_obj, fwhm, self.covlim, temp = self.temp, temp_hdr = self.temp_hdr)
        else:
            cores = clumps.sourceProperties(fin_obj, fwhm, self.covlim)
       
        # determine where the sn is above threshold, if the blob is not on the edges, where the max sn is above threshold and the area is within some area cut. -- only looking for large significant cores.  
        cores['wgd'] = np.where((cores['sn'] > self.nsigma) & (cores['maskbit'] == True) & 
                                (cores['max_sn'] > self.nmxsigma) & (cores['area'] > self.areacut), 1, 0)
        
        
        print(f'number cores found: {len(cores["id"])}')
    
        if ysoradec is not None:
            print('Matching to YSO table')
            idnum, ysocounts = self.ysoBlobs(ysoradec, fin)
            
            cores['ycnt'] = ysocounts
            print(f"number cores matched to YSOs: {len(cores['ycnt'][cores['ycnt'] > 0])}")
        else:
            cores['ycnt'] = np.zeros_like(cores['id'])
            
        cores.to_pickle(self.outbase + '/bloblist_' + self.name + '.pkl')
        print(f'saved bloblist to: {self.outbase}/bloblist_{self.name}.pkl')
        
    def bloblist_name(self):
        '''
        name of bloblist .pkl file
        '''
        return 'bloblist_' + self.name + '.pkl'
                   
    
    def ysoBlobs(self, ysoradec, fin=None):
        ''' 
        determine number of unique blobs with YSOs and how many YSOs per blob
        must have 2d array of yso RA (column 1) & Dec (column 2)
        
        Returns
        --------
        ncores : integer
            number of unique cores
        ycnt : list
            number of ysos found in each blobprint.  List is ordered by core ID number.
        '''
        
        if fin is None:
            clumps, fin_obj, fin = self.run_ImSeg(sn_thresh)
            
        wcs = WCS(self.weight)

    
        ## AL check 
        ##print(np.unique(fin)) # -- just assigns an index for each blob (each has val as id)
                              # -- e.g., [0, 1, 2, ..., 561] for field09
        ## AL check 
        ##print(fin[0], fin[1])



        
        ncores = np.unique(fin)[1::]  # takes elements from index 1 to end, so [1... 561]


        ## AL check 
        ##print('ncores is:', ncores)
        ##print('length of ncores:', len(ncores))


        sz = np.shape(fin)
        
        if len(ncores) > 0:
            ycnt =  []
    
            if len(ysoradec) > 0:
                #find x and y pixel coords of yso ra and dec
                x1, y1 = wcs.wcs_world2pix(ysoradec[:,0], ysoradec[:,1], 0)

#                print(x1) 
#                print(y1)

#                print(y1[y1>0]) # why are these all negative? -was reading in SESNA catalog wrong. 
                
                #determine if they are within the bounds of our image
                w1 = np.where((x1 >= 0) & (x1 <= sz[1] - 1)&(y1 >= 0) &(y1 <= sz[0]-1))
                n1 = len(w1[0])
                
                #if they are within the bounds of image determine their location on flattened fin image 

                #-- AL personal check
                '''
                print('fin is:', fin) # this would be deblended map, 0 for background
                plt.title('Deblended image')
                plt.imshow(fin, origin = 'lower') # looks like values are blobid's 
                plt.plot(x1[w1],y1[w1], "r.", markersize=0.4, label='YSOs')
                plt.legend(loc='upper left')
                plt.show()
                ''' 

                if n1 > 0:                   
                    yy1 = np.int_(y1[w1]+ 0.5)
                    xx1 = np.int_(x1[w1] + 0.5)
                    X   = np.stack((yy1,xx1), axis=-1)
                    idx = np.ravel_multi_index(X.T,np.array([sz[0], sz[1]]))
                    fin2= fin.flatten()
                    fin1= fin2[idx]

                #go through list of cores and determine which blobs the YSOs belong to and how many YSOs there are per blob and append.  
                for z in ncores:
                    if n1 > 0:
                        w  = np.where(fin1 ==z)
                        nw = len(w[0])
                  
                    else:
                        nw = 0
                    ycnt.append(nw)
            else:
                for z in ncores:
                    ycnt.append(0)

            return ncores, ycnt
    
   
        
