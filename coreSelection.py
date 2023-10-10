import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS
from astropy.modeling.models import Gaussian2D
import astropy.units as u

from scipy import signal
from scipy import interpolate
import scipy.stats as st

from getrot import *
from get_flux_mass_conv import *

np.set_printoptions(suppress=True)

class coreSelection:

    '''
    determine cores that pass threshold tests and calculate mass and size
    
    mincut : float, optional
        minimum log(coldensratio) threshold
        (default: -1.25)
    minsn : float, optional
        minimum total s/n threshold
        (default: 3)
    gth : float, optionalco 
        "good" score minimum threshold
        (default: 0.75)
    '''
    
    mincut=-1.25  # -- could you explain how this max dilution threshold was chosen again? 
    minsn=3  
    gth=0.75
      
    def __init__(self, outbase, wavelength, pcdist, coldens, pickle_file, noise_pickle_file=None):
        '''
        Parameters
        ----------
        outbase : string
            location of files
        wavelength : float
            wavelength of observations in mm
        pcdist : float
            distance to object in pc
        coldens : string
            path to column density map
        pickle_file : string
            name of blob .pkl file made from getBlobs
        noise_pickle_file : string or None
            name of noise blob .pkl file made from noiseblobs 
        '''
        self.outbase       =   outbase
        self.pickle_file=   pickle_file
        if noise_pickle_file is not None:
            self.noise_pickle_file = noise_pickle_file
        else:
            self.noise_pickle_file = None
        self.coldens    =   coldens
        self.signal     =   signal
        self.wavelength =   wavelength  # mm
    
        self.pcdist     =   pcdist      # distance to object in pc
        

    def restorePickleFile(self):
        '''
        read pickle_file into pandas dataframe
        
        Returns
        ----------
        blob_list : pandas dataframe
            pandas dataframe of core properties
        '''
        blob_list = pd.read_pickle(self.outbase + '/' + self.pickle_file)
        return blob_list

    def restoreNoisePickleFile(self):
        '''
        read noise_pickle_file into pandas dataframe
        
        Returns
        ----------
        noise_blob_list : pandas dataframe
            pandas dataframe of noise core properties
        '''
        noise_blob_list = pd.read_pickle(self.outbase + '/' + self.noise_pickle_file)
        return noise_blob_list


    def getNH2(self):
        '''
        opens column density image
        
        Returns
        ---------
        im : 2D numpy.ndarray
            coldens data
        hdr : astropy.io.fits.Header 
            header of coldens fits files
        
        '''
        '''
        ## AL edit (since column density / temp maps are in multi-ext fitfiles now)
        print(self.coldens)
        hdu = fits.open(self.coldens)

        im = hdu['COLDEN'].data
        hdr = hdu['COLDEN'].header
        '''

        ## now using my interpolated column density map, so can use old code: 
    
        im  =   fits.getdata(self.coldens)
        hdr =   fits.getheader(self.coldens)
        
        return im, hdr

    def mass_conversion(self,temp):
        '''
        finds conversion factor to convert flux mass 
        
        Parameters
        -------
        temp : list of float
            list or float of temperatures each for blob

        Returns
        --------
        masscon : list or float
            conversion factors to multiply to flux to get mass for each blob
        '''
        
        masscon  =   get_flux_mass_conv(self.pcdist,temp,self.wavelength)
        
        return masscon
    
    
    def findNH2ratio(self, blob_info): 
        '''
        determine ratio of column density of blob in both the coldens image and the blobprints data image
        
        Parameters
        ---------
        blob_info : pandas dataframe
            dataframe of blob properties, can be access with self.restorePickleFile() or self.restoreNoisePickleFile()
            
        Returns
        --------
        tmp : 2d numpy.ndarray
            2D histogram of the S/N and number density ratio
        xp : list
            list of log(s/n) values for each core binned in 0.1 size bins
        yp : list
            ratio of colden NH2 to data NH2 value binned in 0.1 size bins
        '''
        
        #interpolates  refname (fake core image) with column density image; if refname=False, returns column density map (self.coldens) and header
        im, hdr = self.getNH2()
        im = np.nan_to_num(im) 
        wcs     = WCS(hdr)
        sz      = np.shape(im) #determine shape of column density map

        dn_to_nh2 = 1. #number density to H2 number density conversion
        
        #converts ra and dec to pixel coords on the column density map in order to match column density to mm maps
        x1, y1 = wcs.wcs_world2pix(blob_info['ra'], blob_info['dec'], 1)
        x1 = np.array(x1) #only select x pixels
        
        #determine indices on column density map for core ra/dec if map was flattened
        yy1 = np.array(np.int_(y1+ 0.5))
        xx1 = np.array(np.int_(x1 + 0.5))
    
        
        #if (any(q > sz[0] for q in xx1)) or (any(q < 0 for q in xx1)) or (any(q < 0 for q in yy1)) or (any(q > sz[1] for q in yy1)):
            #for i in np.arange(len(xx1)):
                #if (xx1[i] > sz[0]) or (xx1[i] < 0) or (yy1[i] < 0) or (yy1[i] > sz[1]):
                    #xx1[i] = 0
                    #yy1[i] = 0 
        X = np.stack((yy1,xx1), axis=-1)
        ind1 = np.ravel_multi_index(X.T,np.array([sz[0], sz[1]]))
        #determines were there is a good core in both core and noise maps
        w   = blob_info['maskbit']   == 1
        #w[ind1 == 0] = False

        
        #find the mass conversion for good cores in cores and noise
        if 'temp' in blob_info.columns:
                temp   = blob_info['temp']
        else:
            temp   = np.ones_like(blob_info['ra'])*12.
            
        massconvg   = self.mass_conversion(temp)
        
        #flattens column density map
        im2 = im.reshape(-1)
        
        #making a 2D histogram of the S/N and number density ratio calcuated above between (0, 1.9) and (-2, 0.4)
        tmp, x,y = np.histogram2d(np.log10(blob_info['sn'][w]), np.log10(im2[np.int_(ind1[w])]*dn_to_nh2/(1e21)*15)-np.log10(blob_info['sig'][w]*massconvg[w]/(blob_info['area'][w]*(np.pi/180.*self.pcdist)**2)), range=[[0, 1.9],[-2, 0.4]],density=False, bins=[20,25])
        
        xp = np.log10(blob_info['sn']) #log10(S/N of good cores) / 0.1 - determining S/N of good cores
        #calculating column density of Herschel / number density of mm at good cores
        #going into the flattened image and finding the values of the indices of good cores; taking the flux of the good cores and multiplying the mass conversion / area * distance
        yp = ((np.log10(im2[np.int_(ind1)]*dn_to_nh2/1e21*15)) - (np.log10(blob_info['sig']*massconvg/(blob_info['area']*(np.pi/180.*self.pcdist)**2))))
        
        
        return tmp, xp,yp
    
        
    def NH2(self, blob_info, fwhm_beam):
        
        im, hdr = self.getNH2()
        wcs     = WCS(hdr)
        sz      = np.shape(im) #determine shape of column density map

        dn_to_nh2 = 1. #number density to H2 number density conversion
        
        #converts ra and dec to pixel coords on the column density map in order to match column density to mm maps
        x1, y1 = wcs.wcs_world2pix(blob_info['ra'], blob_info['dec'], 1)
        x1 = np.array(x1) #only select x pixels
        
        #determine indices on column density map for core ra/dec if map was flattened
        yy1 = np.int_(y1+ 0.5)
        xx1 = np.int_(x1 + 0.5)
        
        #if (any(q > sz[0] for q in xx1)) or (any(q < 0 for q in xx1)) or (any(q < 0 for q in yy1)) or (any(q > sz[1] for q in yy1)):
            #for i in np.arange(len(xx1)):
                #if (xx1[i] > sz[0]) or (xx1[i] < 0) or (yy1[i] < 0) or (yy1[i] > sz[1]):
                    #xx1[i] = 0
                    #yy1[i] = 0 
        
        X = np.stack((yy1,xx1), axis=-1)
        ind1 = np.ravel_multi_index(X.T,np.array([sz[0], sz[1]]))
        
        #find the mass conversion for good cores in cores and noise
        if 'temp' in blob_info.columns:
                temp   = blob_info['temp']
        else:
            temp   = np.ones_like(blob_info['ra'])*12.
            
        massconvg   = self.mass_conversion(temp)
        
        #flattens column density map
        im2 = im.reshape(-1)
        
        NH2_colden = []
        for i in ind1:
            if i == 0:
                NH2_colden.append(0)
            else:
                NH2_colden.append(im2[np.int_(i)])
        
        mass = (blob_info['sig']*massconvg)*1.989e33
        dist = (self.pcdist*u.pc).to(u.cm)
        m_hydrogen = 1.67e-24*u.g
        mu = 2.3  # -- AL note: Betti+21 says took mu=2.8 (same as reference) but here it is 2.3
        beam_SA = blob_info['area'] * (np.pi/180.)**2. 
        
        NH2_cores = ((mass / (dist**2. * beam_SA * mu* m_hydrogen))) # cm^-2
        
        return np.array(NH2_colden), NH2_cores
        

    def mass(self, flux,temp):
        # calculate mass of core
        massconvg   = self.mass_conversion(temp)
        massdiv = 1
        mass_cores = flux*(massconvg/massdiv)
        return mass_cores
    
    def size(self, hparea, fwhm_beam):    
        # calculate fwhm deconvolved corrected size
        #hparea = np.array(hparea)
        #-- convert hparea from deg2 to arcsec2
        hparea = np.array(hparea) * (3600**2) 
        ##print(hparea)
        fwhm = (2.* np.sqrt((hparea)/np.pi)) #* (self.pcdist/206265.)   #pc 
        ##print(fwhm)
        hpbw = fwhm_beam#* (self.pcdist/206265.)  #pc
        ###print(hpbw)
        size_cores = np.sqrt(fwhm**2. - hpbw**2.)*(self.pcdist/206265.) # --in parsecs
        
        ##print(fwhm**2 - hpbw**2) # --> these were all negatives, hparea is in deg2 not arcsec2
        ##print(size_cores) # --> all nans; after fix, still some nans (but those are not final cores)
                          # --> I guess those cores are < size of the aztec beam 
        return size_cores

        
    def runCoreSelection(self):
        '''
        Determines good cores based on thresholds 
        
        Returns
        -------
        good : list
            list of indices of blobs which pass the thresholds to be considered a core
        goodmap : 2d surface map of the good S/N cores and number densities
        '''
        
        #open pickle pandas dataframes of blobs
        blob_info = self.restorePickleFile()
        
        # find the column density and histogram of blobs
        
        good_blobs   = blob_info['maskbit'][blob_info['maskbit']   == 1]
        print('number of blobs found: ', len(good_blobs))
        
        tmp1, xp, yp = self.findNH2ratio(blob_info)

        # -- again, doing stuff to the bin edges
        xp /= 0.1
        yp = (yp+2)/0.1
        
        # creating a good 2d surface map of the good S/N cores and number densities
        tempgood = np.copy(tmp1)
        #saying for any value in 2D hist of good cores that are less than 1---make them 1
        # -- why would any value in 2D hist be < 1 if it's a histogram (counting?) 
        tempgood[tempgood < 1] = 1

        #open the noise pickle file
        if self.noise_pickle_file is not None:
            noise_blob_info = self.restoreNoisePickleFile()
            good_noise_blobs   =  noise_blob_info['maskbit'][noise_blob_info['maskbit']   == 1]
            print('number of noise blobs found: ', len(good_noise_blobs))
        
            #find the column density and histogram of noise blobs
            tmp2, xnp, ynp = self.findNH2ratio(noise_blob_info)
            xnp /= 0.1
            ynp = (ynp+2)/0.1
            # create goodmap

            # -- so... np.unique(noise_blob_info['numnoise'])[0] would be just 99,  
            # if 99 bloblist_noise*.pkl files are included in the final total bloblist_noise file.
            # then, tmp2/tot_number_noise_included_from_realizations is just
            # average histogram of sn vs. column density ratio across those 99 noise files
            # and then taking tmp2/tempgood is rlly just the ratio between
            # (2dhist of SN vs. NH2ratio from noise) / (2dhist of SN vs. NH2 ratio from data)
            # then we take 1 - (ratio) to get a goodmap
            goodmap = (1-((tmp2)/np.unique(noise_blob_info['numnoise'])[0])/tempgood)

            ##print(np.unique(noise_blob_info['numnoise'])[0])
            
            #interpolate over the surface map just created and find the value in the surface at the location of the good (S/N, number density) found above
            goodY = np.arange(20)
            goodX = np.arange(25)
            f = interpolate.interp2d(goodX,goodY, goodmap, fill_value=1)

            interpolate_noise = [f(ynp.values[i], xnp.values[i])[0] for i in np.arange(len(xnp))]
            noise_blob_info['interpolate_noise'] = interpolate_noise 
            
            #these are the indices of the "official" noise cores that are good and not totally noise
            # use a threshold to only select cores that pass:
            #   - gth: "good" score minimum threshold = interpolated values must be above this threshold
            #   - minsn: minimum s/n threshold = sn of good cores must be greater than this threshold
            #   - mincut: minimum column density ratio threshold = the ratio of herschel or mm of good cores must be above this threshold
            #these are the indices of the "official" cores that are said to be cores and not noise


            # -- so, if the interpolated noise is basically just a regridded form of
            # 1 - (ratio between noise/data 2d histograms of sn vs. coldens_ratio)
            # which means if noise is only 0.25x of flux 2dhist (from data), 
            # it is considered a good core??
            # and this is all for the false core detection, then? 
            noise_blob_info['wng'] = np.where((noise_blob_info['maskbit'] == 1) &
                            (noise_blob_info['interpolate_noise']>= self.gth) & 
                            (noise_blob_info['sn'] >= self.minsn) & 
                            ((ynp * 0.1) - 2. >= self.mincut),1, 0)
            # boolean mask of "official" cores where 1 = official cores
            goodnoise = noise_blob_info['wng'] == 1
            nh = len(noise_blob_info['wng'][goodnoise])
        
            interpolate_flux = [f(yp.values[i], xp.values[i])[0] for i in np.arange(len(xp))]
            blob_info['interpolate_flux'] = interpolate_flux
        
            # use a threshold to only select cores that pass:
            #   - gth: "good" score minimum threshold = interpolated values must be above this threshold
            #   - minsn: minimum s/n threshold = sn of good cores must be greater than this threshold
            #   - mincut: minimum column density ratio threshold = the ratio of herschel or mm of good cores must be above this threshold
            #these are the indices of the "official" cores that are said to be cores and not noise
            blob_info['wg'] = np.where((blob_info['maskbit'] == 1) & 
                                (blob_info['interpolate_flux'] >= self.gth)& 
                                (blob_info['sn'] >= self.minsn) & 
                                ((yp * 0.1) - 2. >= self.mincut)
                                , 1,0)
            # boolean mask of "official" cores where 1 = official cores
            good = blob_info['wg'] == 1
            ng = len(blob_info['wg'][good])
        
        
        else:
            goodmap = tempgood
            blob_info['wg'] = np.where((blob_info['maskbit'] == 1) & 
                                (blob_info['sn'] >= self.minsn) & 
                                ((yp * 0.1) - 2. >= self.mincut)
                                , 1,0)
            # boolean mask of "official" cores where 1 = official cores
            good = blob_info['wg'] == 1
            ng = len(blob_info['wg'][good])
        
        print('number of candidate cores: ', ng)
        
        # save all parameters to numpy save file (.npz) 
        if self.noise_pickle_file is not None:
            wb = noise_blob_info['maskbit'] == 1
            w =  blob_info['maskbit'] == 1
            #save core mass selection criteria as numpy save file due to uneven size data files.  
            np.savez(self.outbase + '/core_mass_select', goodcores = blob_info['id'].values[good], goodnoisecores=noise_blob_info['id'].values[goodnoise],good_fluxvalues=blob_info['interpolate_flux'][good], wg=good[good], wh=goodnoise[goodnoise], w=w[w], wb=wb[wb], tmp1=tmp1, tmp2=tmp2, goodmap=goodmap, numnoise=noise_blob_info['numnoise'][goodnoise])
        else:
            w =  blob_info['maskbit'] == 1
            #save core mass selection criteria as numpy save file due to uneven size data files.  
            np.savez(self.outbase + '/core_mass_select', goodcores = blob_info['id'].values[good], wg=good[good], w=w[w], tmp1=tmp1,  goodmap=goodmap)
        
        return good, goodmap
    
    
    def saveGoodCores(self,fwhm_beam, arcsecond_per_pixel = 1, good_indices=None):
        '''
        save all good core parameters to .csv file 
        
        Parameters
        ---------
        fwhm_beam : float
            fwhm of data
        arcsecond_per_pixel : float
            number of arcseconds per 1 pixel
        good_indices : list or None
            list of indices of true/good cores
            if None, will find good_indices from runCoreSelection() attribute
        '''
        blob_info = self.restorePickleFile()
        if 'temp' in blob_info.columns:
            temp   = blob_info['temp']
        else:
            temp   = np.ones_like(blob_info['ra'])*12.
            
        if good_indices is None:
            good_indices, goodmap = self.runCoreSelection() 
            
        #calculate NH2
        NH2_colden, NH2_cores = self.NH2(blob_info, fwhm_beam)
        
        #calculate mass
        mass_cores = self.mass(blob_info['sig'], temp)
    
        # calculate size
        size_cores = self.size(blob_info['hparea'], fwhm_beam)

        
        # calculate geometric size
        geo_size_cores = np.sqrt(arcsecond_per_pixel**2. * blob_info['fwhm_maj']*(self.pcdist/206265.)**2. *blob_info['fwhm_min'])
        
            
        #save good core property information as pandas dataframe csv file.
        good_cores = {'coreid': blob_info['id'][good_indices], 
                    'ra': blob_info['ra'][good_indices],  #deg
                    'dec':blob_info['dec'][good_indices], #deg
                    'peak_flux':blob_info['max_sig'][good_indices], # Jy / beam
                    'peak_sn':blob_info['max_sn'][good_indices],   
                    'tot_flux':blob_info['sig'][good_indices],       # Jy / beam
                    'tot_sn':blob_info['sn'][good_indices],
                    'area':blob_info['area'][good_indices],     # deg 
                    'yso_count':blob_info['ycnt'][good_indices],
                    'hparea':blob_info['hparea'][good_indices],  # deg
                    'hpra':blob_info['hpra'][good_indices],   # deg
                    'hpdec':blob_info['hpdec'][good_indices], # deg
                    'temp':temp[good_indices], # K 
                    'fwhm_maj': blob_info['fwhm_maj'][good_indices]*arcsecond_per_pixel * (self.pcdist/206265.), # arcseconds
                    'fwhm_min': blob_info['fwhm_min'][good_indices]*arcsecond_per_pixel * (self.pcdist/206265.), # arcseconds
                    'mass':mass_cores[good_indices], # solar masses
                    'size':size_cores[good_indices], # pc
                    'geo_size':geo_size_cores[good_indices], # pc
                    'NH2_colden': NH2_colden[good_indices],  # cm^-2
                    'NH2_cores': NH2_cores[good_indices]}    # cm^-2
        
        good_cores_final = pd.DataFrame(data=good_cores)
        print(f'saving final core selection to: {self.outbase}/final_core_selection.csv')
        good_cores_final.to_csv(self.outbase+'/final_core_selection.csv', index=False)
        
        return


    def plot_ratio(self, instrument, sim=False, good_indices=None, goodmap=None, logSN=None, ratio=None):
        '''
        plot the ratio of column densities to the S/N
        
        Parameters
        ---------
        sim : boolean
            if True: data is from a simulation
            if False: data if from observations
        good_indices : list or None
            list of indices of true/good cores
            if None, will find good_indices from runCoreSelection() attribute
        goodmap : 2d numpy.ndarray or None
            2d histogram of cores
            if None, will find goodmap from runCoreSelection() attribute
            if None, will find good_indices from runCoreSelection() attribute
        logSN : list or None
            list of logged S/N values
            if None, will find logS/N from findNH2Blobs() attribute
        ratio : list or None
            list of NH2 ratio from cores
            if None, will find ratio from findNH2Blobs() attribute  
        
        '''
        
        # create figure
        plt.figure()
        xx,yy = np.mgrid[0:2:20j, -2:0.5:25j] 
        x = np.linspace(0,2,20)
        y = np.linspace(-2,0.5, 25)
        print('plotting comparison to Column density and CMF ')

        blob_info = self.restorePickleFile()
        
        if (good_indices is None) or (goodmap is None):
            good_indices, goodmap = self.runCoreSelection() 
        if (logSN is None) or (ratio is None):
            tmp1, logSN, ratio = self.findNH2ratio(blob_info)
        
        w = blob_info['maskbit'] == 1
        
        # plot S/N vs column density ratio
        plt.plot(logSN[w],ratio[w], 'rD', markersize=4)

        # plot contour of good map -- everything outside is good
        plt.contour(x,y,goodmap.T, colors='k', levels=[0.5, 0.75, 0.9])
        plt.plot([np.log10(self.minsn), np.log10(self.minsn),2],[10, self.mincut, self.mincut], 'k--')
        plt.xlim(.001,2)
        plt.ylim(-1.7,0.5)
        if sim == True:
            plt.ylabel('log {N(H$_2$ ; Sim) / N(H$_2$ ; Sim+AzTEC)}', size=13)
        else:
            plt.ylabel(f'log [N(H$_2$ ; Herschel) / N(H$_2$ ; {instrument}]', size=13)
        plt.xlabel('log S/N', size=13)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)

        #if the number of official cores is more than 0, plot the S/N vs column density of them. --these are good cores
        if len(good_indices) > 0:
            plt.plot(logSN[good_indices],ratio[good_indices], 'D', color='darkred',markersize=4)
            #plt.plot(logSN[good_indices],ratio[good_indices], 'D', color='darkred', markerfacecolor='None' , markersize=4)

        plt.minorticks_on()
        plt.tick_params(which='minor', size=2.5, direction='in', top=True, left=True, right=True, bottom=True)
        plt.tick_params(which='major', size=5, direction='in', top=True, left=True, right=True, bottom=True)

        plt.savefig(self.outbase + '/SNcolumndensity.pdf')
        
        plt.show()
        plt.pause(2)
        plt.close()
        
        # save all parameters to npz files for easy plotting
        np.savez(self.outbase + '/SNvsNH2ratio', sn = logSN[w], NH2ratio = ratio[w])
        np.savez(self.outbase + '/goodmap', x = x, y=y, goodmap = goodmap.T)
        np.savez(self.outbase + '/SNvsNH2ratio_good', sn = logSN[good_indices] ,NH2ratio = ratio[good_indices])


    def plot_CMF(self, mass_cores=None, corrected=False):
        '''
        Plot Core mass function
        
        Parameters
        --------
        mass_cores : list of None
            list of core masses
            if None, will calculate the masses
        corrected : boolean
            if True : assumes masses have been corrected for completeness
            if False: masses have not been corrected for completeness
        
        '''
        massdiv = 1
        
        if mass_cores is None:
            blob_info = self.restorePickleFile()
            if 'temp' in blob_info.columns:
                temp   = blob_info['temp']
            else:
                temp   = np.ones_like(blob_info['ra'])*12.
            #calculate mass
            mass_cores = self.mass(blob_info['sig'], temp)
        
        p2=np.log10(mass_cores)
        binwidth = (3.5*np.std(p2))/(len(p2)**(1./3.))
        
        bins = np.arange(np.min(p2), np.max(p2)+binwidth, binwidth)
        n1, bins1 = np.histogram(p2,bins=bins)#, facecolor='none', linewidth=1.5, edgecolor='k')
        #print(n1) # number in each bin
        #print(bins1)
        
        # plot CMF            
        xhist=np.concatenate((np.array([bins1[0]-.2*massdiv]),
                          bins1[:-1],np.array([bins1[-2]+.2*massdiv])))
        #print(xhist)
        yhist=np.concatenate( (np.array([.001]),n1,np.array([.001])))
        #print(yhist)
        
        plt.figure()
        plt.step(xhist, (yhist/(xhist[1]-xhist[0])),where='mid', color='k')
        plt.ylim(1,1000)
        plt.yscale('log')
        plt.ylabel('dN/dlog M/Msun')
        plt.xlabel('log M/Msun')
        plt.ylim(1,1000)
        plt.minorticks_on()
        plt.tight_layout()
        
        if corrected==False:
            plt.title('Core Mass Function')
            plt.savefig(self.outbase + '/histogram_cmf.pdf')
            plt.savefig(self.outbase + '/histogram_cmf.png')
        else:
            plt.title('Core Mass Function Corrected')
            plt.savefig(self.outbase + '/histogram_cmf_corr.pdf')
            plt.savefig(self.outbase + '/histogram_cmf_corr.png')
        plt.show()
        plt.pause(2)
        plt.close()
    

        
        
       
