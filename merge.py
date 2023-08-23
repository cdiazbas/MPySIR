import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from astropy.io import fits
"""
This module contains functions to merge inversion results.
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
chi2location = 11
# The first one is the baseline
inversion_results = ['finalSIR_cycle3_model.npy','finalSIR_cycle2_model.npy'] 
outputname = '_merged.npy'

# Observed profiles
directory = "/mn/stornext/d20/RoCS/carlosjd/projects/wSPRESOL/data"
observed_stokes =  fits.open(directory+"/test_field_1_630_0.fits")[0].data
observed_stokes = observed_stokes.transpose(2,3,1,0) # (x,y,lambda,stokes)
observed_stokes.shape


# Load any inversion model as baseline
inversion_model_list = []
for i in range(len(inversion_results)):
    inversion_model_list.append(np.load(inversion_results[i], allow_pickle=True))

# Same for the Stokes profiles:
stokes_list = []
for i in range(len(inversion_results)):
    stokes_list.append(np.load(inversion_results[i][:-9]+'profiles.npy'))


# Load all the chi2 maps:
chi2maps = []
for i in tqdm(range(len(inversion_results))):
    ichi2map = np.sum((observed_stokes[:,:,:,1:]-stokes_list[i][:,:,:,1:])**2.0,axis=(2,3))/observed_stokes.shape[3]
    chi2maps.append(ichi2map)
chi2maps = np.array(chi2maps)

# Compute the index of the minimum chi2 of each pixel across all the maps:
index_min_chi2 = np.argmin(chi2maps, axis=0)

# Set first model as baseline:
inversion_model = np.copy(inversion_model_list[0])

# Fill the baseline model with the values of all the variables:
for ip in tqdm(range(inversion_model.shape[3])):
    for ix in tqdm(range(inversion_model.shape[0]),leave=False):
        for iy in tqdm(range(inversion_model.shape[1]),leave=False):
            inversion_model[ix,iy,:,ip] = inversion_model_list[index_min_chi2[ix,iy]][ix,iy,:,ip]

# Azimuth
inversion_model[:,:,:,7] = corrphi(inversion_model[:,:,:,7])
            
# Save the merged model:
np.save(inversion_results[0][:-4]+outputname, inversion_model)


# Set first file as baseline:
stokes = np.copy(stokes_list[0])

# Fill the baseline stokes with the values of all the variables:
for ix in tqdm(range(stokes.shape[0])):
    for iy in tqdm(range(stokes.shape[1]),leave=False):
        stokes[ix,iy,:,:] = stokes_list[index_min_chi2[ix,iy]][ix,iy,:,:]

# Save the merged stokes:
np.save(inversion_results[0][:-9]+'profiles'+outputname, stokes)