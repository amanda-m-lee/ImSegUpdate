## written by Sarah Betti 2019
## updated 5 November 2019

import numpy as np
from astropy.io import fits

def getRot(hdr):
    '''
    NAME:
        getRot

    PURPOSE:
        Return the rotation and plate scale of an image from its FITS header. 

    EXPLANATION:
        Derive the counterclockwise rotation angle, and the X and Y scale factors of an image, from a FITS image header.

    INPUTS:
        hdr = header of a fits file

    OUTPUTS:
        ROT = Scalar giving the counterclockwise rotation of NORTH in DEGREES from the +Y axis of the image.
        CDELT = 2 element numpy array giving the scale factors in DEGREES/PIXEL in the X and Y directions.   CDELT[1] is always positive, whereas CDELT[0] is negative for a normal left-handed coordinate system, and positive for a right-handed system. 

    EXAMPLE:
        header = fits.getheader('<filename>.fits')
        rot, cdelt = getRot(header)

    COMMENTS:
        Rewritten into python from astroIDL library function getrot.pro by Sarah Betti 2019
    '''

    radeg = 180. / np.pi
    cd = np.array([[hdr['CDELT1'], 0], [0, hdr['CDELT2']]])

    if ((cd[1,0] == 0) & (cd[0,1] == 0)): #Unrotated coordinates?
        
        # Amanda note: cd by construction is a 2x2 matrix:
        #              labeling indices, cd = [ [(0,0) (0,1)],
        #                                      [(1,0) (1,1)]]
        
        rot = 0.
        rot2 = 0.

        cdelt1 = hdr['CDELT1']
        cdelt2 = hdr['CDELT2']

    else:
        det = cd[0,0]*cd[1,1] - cd[0,1]*cd[1,0]

        if det < 0:
            sgn = -1
        else:
            sgn = 1

        rot  = np.arctan(sgn*cd[0,1], sgn*cd[0,0])
        rot2 = np.arctan(-cd[1,0], cd[1,1])


        if rot != rot2:          #Handle unequal rotations

            if abs(rot - rot2)*radeg < 2.:
                rot = (rot + rot2)/2.
            elif abs(rot - rot2- tpi)*radeg < 2.:
                rot = (rot +rot2 -tpi)/2.
            elif abs(rot - rot2 +tpi)*radeg < 2.:
                rot = (rot +rot2 +tpi)/2.
            else:
                raise ValueError('WARNING: X and Y axis rotations differ by more than 2 degrees')


        cdelt1 =   sgn*sqrt(cd[0,0]**2 + cd[0,1]**2)
        cdelt2 =   sqrt(cd[1,1]**2 + cd[1,0]**2)

    rot = rot*radeg
    return rot, np.array([cdelt1, cdelt2])
