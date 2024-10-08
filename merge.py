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
# The first one is the baseline
inversion_results = ['inv_100k_5E-3_5mA_1line/100k_5E-3_5mA_1line_conf2_model_fbest.npy','inv_100k_5E-3_5mA_1line/100k_5E-3_5mA_1line_conf1_model.npy'] 
outputname = '_merged.npy'

# Observed profiles
directory = "/mn/stornext/d20/RoCS/carlosjd/projects/wSPRESOL/data"
observed_stokes = np.load(directory+"/sunspot_jmb_sir_synth_profiles_R_100k_5E-3.npy")
observed_stokes = observed_stokes.transpose(0,1,2,3) # (x,y,lambda,stokes)
print('observed_stokes.shape = ',observed_stokes.shape)

observed_stokes = observed_stokes[:,:,range(200,401),:] # Range of wavelength points to be used in the inversion

# Load any inversion model as baseline
inversion_model_list = []
for i in range(len(inversion_results)):
    inversion_model_list.append(np.load(inversion_results[i]))

# Same for the Stokes profiles:
stokes_list = []
for i in range(len(inversion_results)):
    stokes_list.append(np.load(inversion_results[i].replace('model','profiles')))


# Calculate the chi2 maps:
print("Calculating chi2 maps...")
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
print("Merging models...")
for ip in tqdm(range(inversion_model.shape[3])):
    for ix in tqdm(range(inversion_model.shape[0]),leave=False):
        for iy in tqdm(range(inversion_model.shape[1]),leave=False):
            inversion_model[ix,iy,:,ip] = inversion_model_list[index_min_chi2[ix,iy]][ix,iy,:,ip]

# Azimuth
inversion_model[:,:,:,7] = corrphi(inversion_model[:,:,:,7])
            
# Save the merged model:
np.save(inversion_results[0].replace('.npy',outputname), inversion_model.astype(np.float32))


# Set first file as baseline:
stokes = np.copy(stokes_list[0])

# Fill the baseline stokes with the values of all the variables:
print("Merging Stokes profiles...")
for ix in tqdm(range(stokes.shape[0])):
    for iy in tqdm(range(stokes.shape[1]),leave=False):
        stokes[ix,iy,:,:] = stokes_list[index_min_chi2[ix,iy]][ix,iy,:,:]

# Save the merged stokes:
np.save(inversion_results[0].replace('model','profiles').replace('.npy',outputname), stokes.astype(np.float32))