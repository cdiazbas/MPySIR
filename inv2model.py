import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

"""
This module contains functions to map the inversion results to
2 files: one with the maps of the parameters and another with the
synthetic profiles.
"""


# ====================================================================
def readSIRMap(outputSir, parameter, tau):
    """
    It returns the map of a given parameter from the inversion at a given optical depth
    """
    heightMap = outputSir.shape[0]
    widthMap = outputSir.shape[1]
    parmap = np.zeros((widthMap, heightMap))
    for pix_y in range(0, heightMap):
        for pix_x in range(0, widthMap):
            # For vmac, fill, stray and chi2 we need to take the first value
            if parameter == 8 or parameter == 9 or parameter == 10 or parameter == 11:
                parmap[pix_x, pix_y] = outputSir[pix_y,pix_x][1][0][parameter]
            else:
                parmap[pix_x, pix_y] = outputSir[pix_y,pix_x][1][0][parameter][tau]
    return parmap.T


# ====================================================================
def create_modelmap(inversion_file, npar = 12):
    """
    It creates a file with the model parameters [ny, nx, ntau, npar] from the inversion 
    """
    # Read the inversion file:
    inversion = np.load(inversion_file,allow_pickle=True)
    
    # It should have the shape: [ny, nx, ntau, npar]
    logtau = inversion[0,0][1][0][0]
    ntau = len(logtau)
    
    # The height and width of the map:
    ny = inversion.shape[0]
    nx = inversion.shape[1]
    
    # Create the file:
    modelmap = np.zeros((ny, nx, ntau, npar))
    for tau in tqdm(range(ntau)):
        for par in tqdm(range(npar), leave=False):
            modelmap[:, :, tau, par] = readSIRMap(inversion, par, tau)
            
    # Before smoothing, we need to make sure that the parameters are within the limits:
    par = 6 # The inclination angle
    modelmap[:, :, :, par] = np.clip(modelmap[:, :, :, par], 0.0, 180.0)
        
    # Save the file:
    np.save(inversion_file[:-4]+'_model.npy', modelmap)
    

# ====================================================================
def readSIRProfileMap(outputSir, Nstoke):
    """
    It returns the synthetic profiles from the inversion for a given Stokes parameter
    """
    cont = 0
    ny = outputSir.shape[0]
    nx = outputSir.shape[1]
    nwav = len(outputSir[0,0][2][0])
    syn_map = np.zeros((ny, nx, nwav))
    for ypix in range(0, ny):
        for xpix in range(0, nx):
            syn_map[ypix, xpix, :] = outputSir[ypix,xpix][2][1][Nstoke]
    return syn_map


# ====================================================================
def create_profilemap(inversion_file):
    """
    It creates a file with the synthetic profiles from the inversion
    """
    # Read the inversion file:
    inversion = np.load(inversion_file,encoding='latin1',allow_pickle=True)

    # It should have the shape: [ny, nx, nwav, nstokes]
    ny = inversion.shape[0]
    nx = inversion.shape[1]
    nwav = inversion[0,0][2][0].shape[0]
    
    # Create the file:
    profilemap = np.zeros((ny, nx, nwav, 4))
    for stoke in tqdm(range(4)):
        profilemap[:, :, :, stoke] = readSIRProfileMap(inversion, stoke)
    
    # Save the file:
    np.save(inversion_file[:-4]+'_profiles.npy', profilemap)



if __name__ == "__main__":

    filename = 'finalSIR_cycle1.npy'
    create_modelmap(filename)
    create_profilemap(filename)
