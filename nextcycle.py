import numpy as np
from tqdm import tqdm
from scipy.ndimage import gaussian_filter, median_filter

"""
This module contains functions to smooth the inversion results
before using it in the next cycle.
"""

# ========================= PHI CORRECTION
def corrphi(azimuthmap):
    # Fix the azimuth values so that they are in the range [0,180]
    sin_az = np.sin(np.deg2rad(azimuthmap)*2.0)
    cos_az = np.cos(np.deg2rad(azimuthmap)*2.0)    
    azimuthmap = np.rad2deg(np.arctan2(sin_az, cos_az))/2.0
    azimuthmap[azimuthmap<0] = azimuthmap[azimuthmap<0]+180
    return azimuthmap


# ====================================================================
def smooth(fileinput, fwhm_gaussian=0.0, size_median=0, suffix='_smoothed', skip=1):
    """
    Smoothes the inversion results using a gaussian filter and/or a median filter.
    """
    
    # Read the inversion results:
    inversion_model = np.load(fileinput, allow_pickle=True)
    # The inversion model has dimensions:  [ny, nx, ntau, npar]

    # Before smoothing, we need to make sure that the parameters are within the limits:
    par = 6 # The inclination angle
    inversion_model[:, :, :, par] = np.clip(inversion_model[:, :, :, par], 0.0, 180.0)
    par = 7 # The azimuth angle
    inversion_model[:, :, :, par] = corrphi(inversion_model[:, :, :, par])
    
    # Run the smoothing for all optical depths for all parameters:
    for par in tqdm(range(1,inversion_model.shape[3])):
        for tau in tqdm(range(inversion_model.shape[2]), leave=False):
            
            if np.sum(np.isnan(inversion_model[:, :, tau, par])) > 0:
                # Now we fix the Nans with the "replace" astropy function:
                from astropy.convolution import interpolate_replace_nans
                from astropy.convolution import Gaussian2DKernel
                inversion_model[:, :, tau, par] = interpolate_replace_nans(inversion_model[:, :, tau, par], Gaussian2DKernel(1.0),boundary='wrap')
                # And if there are still Nans, we replace them with the median:
                inversion_model[:, :, tau, par] = np.nan_to_num(inversion_model[:, :, tau, par], nan=np.nanmean(inversion_model[:, :, tau, par]))
            
            # Now we can clip them to the percentiles 0.1 and 99.9 (to avoid extreme values):
            inversion_model[:, :, tau, par] = np.clip(inversion_model[:, :, tau, par], np.percentile(inversion_model[:, :, tau, par], 0.1), np.percentile(inversion_model[:, :, tau, par], 99.9))


            # The azimuth is an angle, so we need to smooth it in a different way:
            if par == 7:
                # The azimuth is in degrees, so we convert it to radians:
                inversion_model[:, :, tau, par] = np.deg2rad(inversion_model[:, :, tau, par])*2.0 
                # We multiply by 2 to avoid problems with the 0-360 deg transition
                
                # We smooth the sin and cos of the azimuth:
                sin_az = np.sin(inversion_model[:, :, tau, par])
                cos_az = np.cos(inversion_model[:, :, tau, par])
                if fwhm_gaussian > 0.0:
                    sin_az = gaussian_filter(sin_az, fwhm_gaussian)
                    cos_az = gaussian_filter(cos_az, fwhm_gaussian)
                if size_median > 0:
                    sin_az = median_filter(sin_az, size_median)
                    cos_az = median_filter(cos_az, size_median)
                # We reconstruct the azimuth:
                inversion_model[:, :, tau, par] = np.arctan2(sin_az, cos_az)
                # And we convert it back to degrees:
                inversion_model[:, :, tau, par] = np.rad2deg(inversion_model[:, :, tau, par])/2.0 # We divide by 2 to recover the original values
                # Add 180 to the negative values:
                inversion_model[:, :, tau, par][inversion_model[:, :, tau, par]<0] = inversion_model[:, :, tau, par][inversion_model[:, :, tau, par]<0]+180
            
            else:
                if fwhm_gaussian > 0.0:
                    inversion_model[:, :, tau, par] = gaussian_filter(inversion_model[:, :, tau, par], fwhm_gaussian)
                if size_median > 0:
                    inversion_model[:, :, tau, par] = median_filter(inversion_model[:, :, tau, par], size_median)
            
    # Finally we can interpolate the model to a higher resolution proporcional to the skip parameter
    # in the x and y directions:
    if skip > 1:
        inversion_model = np.repeat(inversion_model, skip, axis=0)
        inversion_model = np.repeat(inversion_model, skip, axis=1)
        print('The model has been interpolated to a higher resolution of {}x{} pixels'.format(inversion_model.shape[0], inversion_model.shape[1]))
            

    # Save the smoothed model adding the suffix '_smoothed':
    np.save(fileinput[:-4]+suffix, inversion_model.astype(np.float32))



if __name__ == "__main__":
    fileinput = 'inv_original_newgrid/finalSIR_cycle1_model.npy'
    fwhm_gaussian = 0.5
    size_median = 0
    smooth(fileinput, fwhm_gaussian, size_median, suffix='_smoothed', skip=1)
    