import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from astropy.convolution import interpolate_replace_nans, Gaussian2DKernel

"""
This module contains functions to fix inversion results. If mode is 'interpolate', the missing values are interpolated. If mode is 'replace', the missing values are replaced by the values of the best chi2 model.
"""


# ========================= PHI CORRECTION
def corrphi(azimuthmap):
    # Fix the azimuth values so that they are in the range [0,180]
    sin_az = np.sin(np.deg2rad(azimuthmap)*2.0)
    cos_az = np.cos(np.deg2rad(azimuthmap)*2.0)    
    azimuthmap = np.rad2deg(np.arctan2(sin_az, cos_az))/2.0
    azimuthmap[azimuthmap<0] = azimuthmap[azimuthmap<0]+180
    return azimuthmap


# ========================= MERGE
inversion_results = 'finalSIR_cycle1_model.npy' 
outputname = '_fixed.npy'
mode = 'interpolate' # 'interpolate' or 'replace' (replace not yet implemented)
chi2_limit = 100.0 # Limit for the chi2 to be considered a good inversion

# Load the inversion model
inversion_model = np.load(inversion_results)

# Same for the Stokes profiles:
stokes_input = np.load(inversion_results[:-9]+'profiles.npy')

# Azimuth
inversion_model[:,:,:,7] = corrphi(inversion_model[:,:,:,7])

# Run the code: 
if mode == 'interpolate':
    
    # Calculate the chi2 map:
    chi2 = inversion_model[:,:,0,11]
    
    print('Interpolating',np.sum(chi2>chi2_limit),'pixels')
   
    # Create mask (1 if chi2<chi2_limit, 0 otherwise):
    mask = np.zeros(chi2.shape)
    mask[chi2<chi2_limit] = 1.0
    
    # We now go through all parameters and interpolate the missing values
    # using the astropy function that takes into account the NaN values:
    for ip in tqdm(range(inversion_model.shape[3])):
        for it in tqdm(range(inversion_model.shape[2]),leave=False):
            model_tau = inversion_model[:,:,it,ip]*1.0
            model_tau[mask==0.0] = np.nan
                    
            # We interpolate the missing values:
            kernel = Gaussian2DKernel(x_stddev=1.0)
            inversion_model[:,:,it,ip] = interpolate_replace_nans(model_tau, kernel,boundary='wrap')


"""
# Not yet implemented:
elif mode == 'replace':
    
    # Calculate the chi2 map:
    chi2 = inversion_model[:,:,0,11]
    
    print('Replacing',np.sum(chi2>chi2_limit),'pixels')
    
    # For each pixel where chi2>chi2_limit, we replace the values 
    # with the values of the best chi2 model:
    for ix in tqdm(range(inversion_model.shape[0]),leave=True):
        for iy in tqdm(range(inversion_model.shape[1]),leave=False):
            if chi2[ix,iy] > chi2_limit:
                # We find the index of the best chi2 model:
                index_min_chi2 = np.argmin(chi2[ix,iy])
                # We convert the index to 2D coordinates:
                index_min_chi2 = np.unravel_index(index_min_chi2, chi2.shape)
                
                # We replace the model and stokes values:
                inversion_model[ix,iy,:,:] = inversion_model[index_min_chi2[0],index_min_chi2[1],:,:]
                stokes_input[ix,iy,:,:] = stokes_input[index_min_chi2[0],index_min_chi2[1],:,:]
"""

# Save the new model:
np.save(inversion_results[:-4]+outputname, inversion_model)


# Save the merged stokes:
np.save(inversion_results[:-9]+'profiles'+outputname, stokes_input)