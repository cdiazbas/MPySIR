import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

"""
This module contains functions to smooth the inversion results
before using it in the next cycle.
"""

# ====================================================================
def smooth(fileinput, fwhm_gaussian=0.0, size_median=0):
    """
    It smooths the inversion results using a gaussian filter and/or a median filter.
    """
    
    # Read the inversion results:
    inversion_model = np.load(fileinput, allow_pickle=True)

    # Define the guassian and median filter from ndimage:
    from scipy.ndimage import gaussian_filter, median_filter

    # The inversion model has dimensions:  [ny, nx, ntau, npar]

    # Run the smoothing for all optical depths for all parameters:
    for tau in tqdm(range(inversion_model.shape[2])):
        for par in tqdm(range(inversion_model.shape[3]), leave=False):
            
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
            
            else:
                if fwhm_gaussian > 0.0:
                    inversion_model[:, :, tau, par] = gaussian_filter(inversion_model[:, :, tau, par], fwhm_gaussian)
                if size_median > 0:
                    inversion_model[:, :, tau, par] = median_filter(inversion_model[:, :, tau, par], size_median)
                                
                
    # Save the smoothed model adding the suffix '_smoothed':
    np.save(fileinput[:-4]+'_smoothed.npy', inversion_model)



if __name__ == "__main__":
    fileinput = 'finalSIR_model.npy'
    fwhm_gaussian = 2.0
    size_median = 0
    smooth(fileinput, fwhm_gaussian, size_median)
    