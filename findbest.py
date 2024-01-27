import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from astropy.io import fits
"""
This module contains functions to fix bad pixels in the inversion results by finding a better model from
another pixel.
"""

# ========================= FINDBEST =========================
# The first one is the baseline
inversion_results = 'inv_original_newgrid/finalSIR_cycle3_model.npy'
outputname = '_fbest.npy'


# Observed profiles
directory = "/mn/stornext/d20/RoCS/carlosjd/projects/wSPRESOL/data"
observed_stokes = np.load(directory+"/sunspot_jmb_sir_synth_newgrid_profiles.npy")
observed_stokes = observed_stokes.transpose(0,1,2,3) # (x,y,lambda,stokes)
print('observed_stokes.shape = ',observed_stokes.shape)

# Number of pixels to fix:
npix = 4*4000
# npix = observed_stokes.shape[0]*observed_stokes.shape[1]
print('npix = ',npix)

# Load the inversion model as baseline
inversion_model = np.load(inversion_results, allow_pickle=True)
# Same for the Stokes profiles:
stokes = np.load(inversion_results[:-9]+'profiles.npy')


# Calculate the chi2 maps:
print("Calculate the initial chi2 map...")
chi2map = np.sum((observed_stokes[:,:,:,:]-stokes[:,:,:,:])**2.0,axis=(2,3))/observed_stokes.shape[3]

# Start fixing the worst pixels:
print("Fixing the worst pixels...")
chi2map_ordered = np.sort(chi2map.flatten())
for i in tqdm(range(npix)):
    # Find the pixel with the highest chi2:
    index = np.where(chi2map == chi2map_ordered[-(i+1)])
    # Take only one pixel:
    index = (index[0][0],index[1][0])
    
    # Calculate the new chi2map for this pixel:
    ichi2map = np.sum((observed_stokes[index[0],index[1],:,:]-stokes[:,:,:,:])**2.0,axis=(2,3))/observed_stokes.shape[3]
    # Find the pixel with the lowest chi2:
    index_min_chi2 = np.where(ichi2map == np.min(ichi2map))
    # Take only one pixel:
    index_min_chi2 = (index_min_chi2[0][0],index_min_chi2[1][0])
    
    print("Fixing pixel: ",index, "(chi2: {0:1.1e})".format(chi2map[index[0],index[1]]))
    print("with model for pixel: ",index_min_chi2, "(chi2: {0:1.1e})".format(ichi2map[index_min_chi2[0],index_min_chi2[1]]))

    # Set the inversion results 
    inversion_model[index[0],index[1],:,:] = inversion_model[index_min_chi2[0],index_min_chi2[1],:,:]
    stokes[index[0],index[1],:,:] = stokes[index_min_chi2[0],index_min_chi2[1],:,:]


# Save the merged model:
np.save(inversion_results[:-4]+outputname, inversion_model.astype(np.float32))

# Save the merged stokes:
np.save(inversion_results[:-9]+'profiles'+outputname, stokes.astype(np.float32))