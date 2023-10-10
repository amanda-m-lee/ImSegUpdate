## written by Sarah Betti 2019
## updated 5 November 2019

import numpy as np
import numpy.ma as ma
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mpl

from astropy.io import fits
from astropy.stats import SigmaClip
from astropy.visualization import ZScaleInterval
from astropy.wcs import WCS

from photutils import Background2D, SExtractorBackground
from photutils import detect_sources
from photutils import deblend_sources
from photutils import source_properties
from photutils import SegmentationImage

from getrot import *


class clumpFind:
    def __init__(self, data, sn, hdr, wt=None, w=None):
        '''
        Purpose 
        ----------
        clump finding using photutils image segmentation and deblending.  
        
        Parameters
        ----------
        data : 2d numpy.ndarray
            signal data
        sn : 2d numpy.ndarray
            signal to noise data
        hdr : fits file header
        wt : 2d numpy.ndarray (optional)
            weight data
            if wt = None, then the weight data will be calculated from (sn / data)
        w : wcs (optional)
            WCS information from hdr
            if w = None, then the wcs information will be found from the hdr
        '''
        self.data   =   np.nan_to_num(data)
        self.sn     =   np.nan_to_num(sn)
        self.hdr    =   hdr
        
        if wt is None:
            self.wt     =  (np.divide(self.sn, self.data, out=np.zeros_like(self.sn), where=self.data !=0)**2. )
        else:
            self.wt = np.nan_to_num(wt)
            
        if w is None:
            self.w   = WCS(self.hdr) ## WCS information
        else:
            self.w = w

                        
    def background(self):
        '''
        calculates the background of the data signal image where signal > 0.
        
        Returns
        ---------
        noise_map : 2d numpy.ndarray
            background noise map calculated from (data / sn)
            
        '''

        noise_map = np.divide(self.data, self.sn,out=np.zeros_like(self.data), where=self.sn !=0)
        return noise_map
    
    
    
    def detectSources(self, sn_thresh, covlim = 0.4, npixels=10):
        '''
        uses photutils detect_sources to determine blobs.  
        
        Parameters
        ----------
        sn_thresh : float
            minimum S/N to be considered part of a blob 
        covlim : float      
            coverage limit to consider when detecting blobs
        npixels : float
            number of connected pixesl for source matching
        
        Returns
        ----------
        segm : object class
            object of the segmentation of blobs.  can be accessed with segm.data.  should be used in deblendSources()
        
        
        '''
        #mask data before source detection
        #   - only choose souces with: weights above coverage limit*median 
        #   - sn above sn threshold
        #   - where data != 0
        mask0   = self.wt > 0
        mask00  = self.wt < covlim*np.median(self.wt[mask0])
        mask1   = self.sn < sn_thresh
        mask2   = self.data == 0

        #make mask
        # -- AL note: 
        # -- bitwise OR | --> only returns False is BOTH X and Y e.g., in X|Y are false
        # -- otherwise, it will be True 
    
        mask    =  (mask00 | mask1) | mask2  ## mask1 | mask2

        #calculate background
        threshold = self.background() # -- AL note: this is just the noise map, data/sn 

        #find sources that have npixels connected pixels that are each greater than the corresponding pixel-wise threshold level defined above by the noise map.  A value of zero is always reserved for the background.

        # -- AL note: anything in mask for photutils IS NOT included for detecting sources
        segm = detect_sources(self.data, threshold, npixels=npixels, mask=mask)
        return segm
        
    def deblendSources(self, sn_thresh, covlim=0.4, npixels=10, contrast=0.001):
        '''
        deblend sources which overlap.  
        uses photutils deblend_sources() that uses a combination of multi-thresholding and watershed segmentation.  Note that in order to deblend sources, they must be separated enough such that there is a saddle between them. 
        
        Parameters
        ----------
        sn_thresh : float
            minimum S/N to be considered part of a blob  
        covlim : float      
            coverage limit to consider when detecting blobs
        npixels : float
            number of connected pixesl.  must be the same as in detectSources()
        contrast : float 
            contrast is the fraction of the total source flux that a local peak must have to be considered as a separate object.
        
        Returns
        ---------
        segm : object
            deblended sources object
        segm.data : numpy.ndarray
            deblended source data array
        
        NOTE
        --------
        had to update photutils deblend.py with an update that is not implemented yet.  see https://github.com/saimn/photutils/commit/76a65ecc7f24c653f788fefb34931798296baa68 and follow green updates or contact sbetti@umass.edu for updates. 
        
        Further updates are necessary to change nlevels from input variable to variable based on sigma of each blob.  
        '''
        
        print('***If deblending does not occur: check to make sure all corrections to photutils deblend.py have been applied***')
        #runs detect sources
        segm    = self.detectSources(sn_thresh, covlim, npixels)
  
        # create the same masks to run deblending
        mask0   = self.wt > 0
        mask00  = self.wt < covlim*np.median(self.wt[mask0])
        mask1   = self.sn < sn_thresh
        mask2   = self.data == 0

        mask    = (mask00 | mask1) | mask2 ## mask1 | mask2 
        
        #copy data in order to mask without overwriting
        datt        = np.copy(self.data)
        datt[mask]  = 0
        
        
        snn = np.copy(self.sn)
        snn[mask]   =   0

        # -- AL note: how come did not also just use background() here like before for noise?
        noise       =   np.divide(datt, snn, out=np.zeros_like(datt), where=snn !=0)
        noise[noise == 0] = np.nan
        mednoise    =   np.nanmedian(noise)
        
        byte_im     =   np.divide(datt, mednoise, out=np.zeros_like(datt), where=mednoise !=0)
        byte_im[byte_im < 0]    =   0
        byte_data = 255 - np.int_(byte_im)
        
        #try deblending. 
        segm_deblend = deblend_sources(datt, byte_data, segm, npixels=npixels, contrast=contrast)  ### manually updated version of deblend_sources
            
        try:
            print('deblending complete')
            return segm_deblend, segm_deblend.data
        except:
            print('no deblending occurred')
            return segm, segm.data
                           
    def edge(self, edgecut):
        '''
        finds edges of image 
        '''
        #get the rotation and cdelt of the image
        rot, cdelt = getRot(self.hdr)
        if edgecut == True:
            near = 70./3600. / cdelt[1] # -- AL Q: what is 70 for? 
        else:
            near = 8./3600./cdelt[1] # AL Q: to reject cands. within 8 arcsec of coverage?   
        return near
    

    def sourceProperties(self, segm_deblend, fwhm, covlim =0.4, edgecut=True, **kwargs):
        '''
        Finds the properties of the blobs:
            - ID number
            - ra
            - dec
            - max sn
            - max signal
            - total sn
            - total signal
            - area
            - maskbit
            - half power area
            - half peak power ra
            - half peak power dec 
            - temperature (dependent if temp map is given)
            - fwhm major axis 
            - fwhm minor axis
            
        Parameters
        --------
        segm_deblend : object 
            deblended source array
        fwhm : float    
            fwhm of data
        covlim : float
            coverage limit 
        edgecut : boolean 
            if True : cut the edges at 70":
            if False: cut the edges at 8"
        kwargs : 
            temp : 2d numpy.ndarray
                temperature map 
            temp_hdr : .fits header
                header of temperature map
        
        Returns
        --------
        df : pandas dataframe
            dataframe of all the properties for each blobs
            
        '''
            
        #get rotation
        rot, cdelt = getRot(self.hdr)
        #extract core properties and sn properties
        cat     = source_properties(self.data, segm_deblend, wcs=self.w)
        cat2    = source_properties(self.sn, segm_deblend, wcs=self.w) 

        #get array of background subtracted pixel values within segment
        meann = [np.mean(i) for i in cat.values]
        # sum of the noise (signal/sn) in each segment

        # -- AL Q: how come squared the sum of noise in each segment?
        
        noisesum = [np.sum((cat.values[i] / cat2.values[i])**2.) for i in np.arange(len(cat.values))]

        #x, y, ra, dec values for each segment 
        xc  = cat.xcentroid.value
        yc  = cat.ycentroid.value
        tra = cat.sky_centroid.ra.degree
        tdec= cat.sky_centroid.dec.degree 
        
        #maximum S/N 
        msns    = cat2.max_value  
        # signal at maximum signal to noise
        msigs   = self.data[np.int_(cat2.maxval_ypos.value), np.int_(cat2.maxval_xpos.value )] 

        # area of segment in degrees^2
        areas = cat.area.value*(cdelt[1]**2. )
        
        # average signal in segment
        sigs_i = np.array([np.mean(i) for i in cat.values]) * areas/(1.1331*(fwhm/3600.)**2. )
        
        # average noise in segment
        noise = np.sqrt(noisesum)/np.sqrt(1.1331*(fwhm/3600./cdelt[1])**2)
        #average S/N
        sns = sigs_i/noise  

        #half power signal, noise, sn 
        #finding all the values greater than half the maximum S/N
        whp =[cat[ind].values > 0.5*self.data[np.int_(cat2[ind].maxval_ypos.value), np.int_(cat2[ind].maxval_xpos.value )] for ind in np.arange(len(cat.values)) ]


        # ------- Debugging with Grant ---------- 
        #breakpoint()
        
        #print([i[1] for i in cat.coords])

        #print(cat.to_table())
        
        
        # finding x and y coordinates of each segment 
        #xhp0 = np.array([i[1] for i in cat.coords])  -- removed np.array() 
        #yhp0 = np.array([i[0] for i in cat.coords]) 

        # ----------------------------------------

        xhp0 = [i[1] for i in cat.coords]
        yhp0 = [i[0] for i in cat.coords]

        # find x and y coordinates of half power for each segment 
        xhp = [xhp0[ind][whp[ind]] for ind in np.arange(len(whp))]
        yhp = [yhp0[ind][whp[ind]] for ind in np.arange(len(whp))]

        # find average center x/y coordinate of half power

        # -- AL note: looks like flux weighted hp coords. 
        xchp = [np.sum(self.data[yhp[ind], xhp[ind]] * xhp[ind]) / np.sum(self.data[yhp[ind], xhp[ind]]) for ind in np.arange(len(xhp)) ]
        ychp = [np.sum(self.data[yhp[ind], xhp[ind]] * yhp[ind]) / np.sum(self.data[yhp[ind], xhp[ind]]) for ind in np.arange(len(xhp)) ]

        # ra, dec at center half power
        hpra, hpdec = self.w.wcs_pix2world(xchp, ychp, 0)
        hpareas = [(len(xhp[ind])*cdelt[1]**2) for ind in np.arange(len(xhp))]

        #find edges and create mask where values are close to the edge
        near = self.edge(edgecut)
        mask = np.copy(self.wt)
        masked = np.where(self.wt >0)
        mask[self.wt > covlim*np.median(self.wt[masked])] =  1.
        mask[self.wt <= covlim*np.median(self.wt[masked])] = 0
        
        mask2 = [ mask[int(np.min(yhp0[ind])-near):int(np.max(yhp0[ind])+near+1), int(np.min(xhp0[ind])-near):int(np.max(xhp0[ind])+near+1)] for ind in np.arange(len(xhp0))]
        
        #give each segment 1/0 if it is near the edge or not.
        # -- AL note: this comment should be reversed. 1 is good, 0 is bad=near edge
        
        maskbit = []
        for i in np.arange(len(mask2)):
            summ = np.sum(mask2[i]) # -- AL: mask2 only has values of 1's or 0's
            if summ == len(mask2[i].reshape(-1)): # AL: -- len of FLATTENED(-1) mask2[i] 
                maskbit.append(1)                 # would be same as summ if all are included
            else:                                 # (hence within coverage limit)
                maskbit.append(0)
        maskbit = np.array(maskbit)
        
        ## get FWHM of semimajor and semiminor axes of sources
        const = 2.*np.sqrt(2.*np.log(2.))
        fwhm_major = const*cat.semimajor_axis_sigma.value #pixels
        fwhm_minor = const*cat.semiminor_axis_sigma.value
        
        #calculate temperatures of each segment using temperature map 
        if 'temp' in kwargs:
            temp = kwargs['temp']
            temp_hdr = kwargs['temp_hdr']

            ## use my interpolated temperature map 
            temp = fits.getdata(temp)

            # this is just for WCS information. I set the header to be the same as coldens
            # (for my interpolated temp map)
            temp_hdr = fits.getheader(temp_hdr)

            '''
            ## AL edit since temperature is part of multi-ext file 
            hdu = fits.open(temp)
            temp = hdu['TEMP'].data
            temp_hdr = hdu['COLDEN'].header ## somehow temp map does not have WCS info
                                            ## i'm going to assume WCS info is same as coldens
                                            ## map... 
            '''
            
            temp_wcs = WCS(temp_hdr)
            avg_temp = []
            for ind in np.arange(len(xhp0)):
                ra_in = []
                dec_in = []
                for ind2 in np.arange(len(xhp0[ind])):
                    ra2, dec2 = self.w.wcs_pix2world(xhp0[ind][ind2], yhp0[ind][ind2], 0)
                    ra_in.append(ra2)
                    dec_in.append(dec2)
                xxx, yyy = temp_wcs.wcs_world2pix(ra_in, dec_in, 0)
                avg_temp.append(np.median(temp[np.int_(yyy),np.int_(xxx)]) )
                
            d = {'id':cat.id,'ra': tra, 'dec':tdec, 'max_sn':msns, 'max_sig':msigs, 'sn':sns, 'sig':sigs_i, 'area':areas, 'maskbit': maskbit, 'hparea':hpareas, 'hpra':hpra, 'hpdec':hpdec, 'temp':np.array(avg_temp), 'fwhm_maj': fwhm_major, 'fwhm_min': fwhm_minor}
            df = pd.DataFrame(data=d)
#            
            return df
        
        else:
            d = {'id':cat.id,'ra': tra, 'dec':tdec, 'max_sn':msns, 'max_sig':msigs, 'sn':sns, 'sig':sigs_i, 'area':areas, 'maskbit': maskbit, 'hparea':hpareas, 'hpra':hpra, 'hpdec':hpdec, 'fwhm_maj': fwhm_major, 'fwhm_min': fwhm_minor}
            df = pd.DataFrame(data=d)
            
            return df
       
             
    def saveSourcePropertiesDF(self, segm_deblend, outbase, **kwargs):
        '''
        save all properties as csv file
        
        Parameters
        ---------
        segm_deblend : object
            deblended source data
        outbase : string
            location to save files
        **kwargs : 
            temp : 2d numpy.ndarray
                temperature data 
            temp_hdr : fits header
                header of temperature data
        '''
        
        if 'temp' in kwargs:
            d = self.sourceProperties(segm_deblend, kwargs['temp'], kwargs['temp_hdr'])
        else:
            d = self.sourceProperties(segm_deblend)
        
        df.to_csv(outbase+'/source_properties.csv', index=False)
    
    def plotClumps(self, segm_deblend, source, outbase):
        '''
        plot clumps and original image
        
        Parameters
        -------
        segm_deblend : object
            deblended source data
        source : string
            name of objects being plotted (eg: cores)
        outbase : string
            location to save files
        
        Returns
        --------
        im2 : 2d numpy.ndarray
            outlines segmenation image
        
        '''

        interval = ZScaleInterval()
        vmin, vmax = interval.get_limits(self.data)
      
        plt.close()
        fig = plt.figure(figsize=(10,6))
        ax1 = plt.subplot(121,projection=self.w)
        ax1.imshow(self.data, cmap='Greys_r', vmin=vmin, vmax=vmax)
        ax1.set_title('Data:' + source)

        
        ## silly method to make an outline around each segment to show divisons.  divisions are now a pixel width with value = 0
        im = SegmentationImage(segm_deblend).outline_segments() #makes pixel width outline with the same label as the segment.  everything but outline have a value = 0
        im[im > 0 ] = -10000.  # makes all the outlines a negative number.  Since no labels are negative, you will never have to worry about accidently overlapping with a core label
        im2 = im + segm_deblend  #add the original segmentation image to this new outline image.  All outlines have value -10000, everything else is 0, so you preserve the original labels except at the outlines.
        im2[im2 < 0] = 0.  # make all the outlines 0, the same as the background
        outline_image = im2  #rename for simplity
        im3 = np.copy(im2)
        im3[im3 > 0]=1
        ax2 = plt.subplot(122, projection=self.w)
        cmap = mpl.colors.ListedColormap(np.random.rand(256,3))
        ax2.imshow(im3, cmap='binary')#segm.cmap(random_state=12345))
        ax2.set_title('Segmentation Image:' + source)

        plt.savefig(outbase + '/clumpfind_cores.png')
        plt.savefig(outbase + '/clumpfind_cores.pdf')
        plt.show()
        plt.pause(5)
        plt.close()
        
        
        
        return outline_image
        
    
    
    
    
    
    
    
    
