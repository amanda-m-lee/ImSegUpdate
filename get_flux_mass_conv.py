## written by Sarah Betti 2019
## updated 5 November 2019

import numpy as np
from astropy import units as u
from astropy import constants as const

def get_flux_mass_conv(pcdist,T, wavelength=1100):
    
    '''
    PURPOSE
    --------
    Return mass conversion factor for total flux in Jy. 

    Parameters
    --------
    pcdist: float in pc
        distance to object in pc
    T: float in Kelvin 
        temperature of object
    wavelength: float in microns
        wavelength at which to find flux.  Can only be [350, 850, 1100] microns
        default is 1100 microns (1.1 mm)

    Returns
    --------
    m_XXum: float
        conversion factor to convert total flux to XX Jy -- AL: is this supposed to be solar mass?
        
    Raises
    -------
    ValueError if wavelength is not 350, 850, 1100
        
    Examples
    --------
        m_350um = get_aztec_mass_conv(800, 18, 350)

    Comments
    --------
        written in IDL by Rob Gutermith 
        converted to python by Sarah Betti 2019
    '''

    h   =   const.h.cgs.value     # g cm^2 / s OR  erg s
    k   =   const.k_B.cgs.value   # erg / K
    c   =   const.c.cgs.value     # cm /s
    d   =   pcdist *3.0857e18     # cm

    if wavelength == 350.:
        ### Ossenkopf & Henning 1994 #5 opacities; "OH5"
        kappa_350um     =   0.101                   # cm^2 / g  (includes gas-to-dust)

        ### Dunham et al. 2009 in prep.
        snu_350um       =   151.2*1.0*10**-23       # erg/cm^2/s/Hz
        snu_350um_alt1  =   81.3*1.0*10**-23        # erg/cm^2/s/Hz
        snu_350um_alt2  =   34.6*1.0*10**-23        # erg/cm^2/s/Hz

        nu_350um        =   (2.99792458e14)/350.    # Hz
        planck_350um    =   (2.*h*(nu_350um**3.)/(c**2.))/((np.exp((h*nu_350um)/(k*T)))-1.)

        m_350um         =   (((d**2.)*snu_350um)/(planck_350um*kappa_350um))/1.99e33
        m_350um_alt1    =   (((d**2.)*snu_350um_alt1)/(planck_350um*kappa_350um))/1.99e33
        m_350um_alt2    =   (((d**2.)*snu_350um_alt2)/(planck_350um*kappa_350um))/1.99e33

        return m_350um


    elif wavelength == 450.:
        ### Ossenkopf & Henning 1994 #5 opacities; "OH5"
        kappa_450um     =   0.0674               # cm^2 / g  (includes gas-to-dust)

        
        snu_450um       =   1.0*10**-23      # erg/cm^2/s/Hz

        nu_450um        =   (2.99792458e14)/450.    # Hz
        planck_450um    =   (2.*h*(nu_450um**3.)/(c**2.))/((np.exp((h*nu_450um)/(k*T)))-1.)

        m_450um         =   (((d**2.)*snu_450um)/(planck_450um*kappa_450um))/1.99e33

        return m_450um

    elif wavelength == 850.:
        ### Ossenkopf & Henning 1994 #5 opacities; "OH5"
        kappa_850um     =   0.0114                  # cm^2 / g  (includes gas-to-dust) (OH*4*, 845um)

        
        snu_850um       =   1.0*10**-23             # erg/cm^2/s/Hz

        nu_850um        =   (2.99792458e14)/850.    # Hz
        planck_850um    =   (2.*h*(nu_850um**3.)/(c**2.))/((np.exp((h*nu_850um)/(k*T)))-1.)

        m_850um         =   (((d**2.)*snu_850um)/(planck_850um*kappa_850um))/1.99e33

        return m_850um


    elif wavelength == 1100.:
        ### Ossenkopf & Henning 1994 #5 opacities; "OH5"
        #kappa_1100um   =   0.0114                  # cm^2 / g  (includes gas-to-dust) (OH*4*, 845um)
        kappa_1100um    =   0.0121                  # cm^2 / g  (includes gas-to-dust) (close enough?)
        #kappa_1100um   =   0.0069                  # cm^2 / g  (includes gas-to-dust) (OH*4*, 1.1mm)

       
        snu_1100um      =   1.0*10**-23             # erg/cm^2/s/Hz

        nu_1100um       =   (2.99792458e14)/1100.0  # Hz
        # -- e.g., c = 3x10^8 m/s
        # -- wavelength is in micron, so multiply c*10^6 micron/m
        # hence, gives freq. in hertz. 
        planck_1100um   =   (2.*h*(nu_1100um**3.)/(c**2.))/((np.exp((h*nu_1100um)/(k*T)))-1.)
        # -- in units: ergs/cm2/s/steradian/Hz
        
        m_1100um        =   (((d**2.)*snu_1100um)/(planck_1100um*kappa_1100um))/1.99e33
        # -- becomes units: Msun * sr??? 
        return m_1100um

    else:
        raise ValueError('wavelength does not exist.  Use either 350, 345, 850, 1100 microns instead.')
