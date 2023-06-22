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
    cont = 0
    heightMap = outputSir[-1][-1][0][0]+1
    widthMap = outputSir[-1][-1][0][1]+1
    mapa = np.zeros((widthMap, heightMap))
    for fila in range(0, heightMap):
        for columna in range(0, widthMap):
            punto = cont % outputSir.shape[1]
            veces = int(cont/outputSir.shape[1])
            if parameter == 8 or parameter == 9 or parameter == 10 or parameter == 11:
                mapa[columna,fila] = outputSir[veces][punto][1][0][parameter]
            else:
                mapa[columna,fila] = outputSir[veces][punto][1][0][parameter][tau]
            cont += 1
    return mapa.T


# ====================================================================
def readSIRProfileMap(outputSir, Nstoke):
    """
    It returns the synthetic profiles from the inversion for a given Stokes parameter
    """
    cont = 0
    [height, width, rangoLambda] = shapeSIRMap(outputSir)
    synMapa = np.zeros((width, height, rangoLambda))
    for fila in range(0, height):
        for columna in range(0, width):
            punto = cont % outputSir.shape[1]
            veces = int(cont/outputSir.shape[1])
            synMapa[columna, fila, :] = outputSir[veces][punto][2][1][Nstoke][:]
            cont += 1
    return synMapa

# ====================================================================
def shapeSIRMap(resultadoSir):
    """
    It returns the shape of the SIR map
    """
    height = resultadoSir[-1][-1][0][0]+1
    width  = resultadoSir[-1][-1][0][1]+1
    nlambda = len(resultadoSir[0][0][2][1][0])
    return [height, width, nlambda]


# ====================================================================
def create_profilemap(inversion_file):
    """
    It creates a file with the synthetic profiles from the inversion
    """
    # Read the inversion file:
    inversion = np.load(inversion_file,encoding='latin1',allow_pickle=True)

    # It should have the shape: [ny, nx, nwav, nstokes]
    [ny, nx, nwav] = shapeSIRMap(inversion)
    
    # Create the file:
    profilemap = np.zeros((nx, ny, nwav, 4))
    for stoke in tqdm(range(4)):
        profilemap[:, :, :, stoke] = readSIRProfileMap(inversion, stoke)
    
    # Save the file:
    np.save(inversion_file[:-4]+'_profiles.npy', profilemap)


# ====================================================================
def create_modelmap(inversion_file, npar = 12):
    """
    It creates a file with the model parameters from the inversion
    """
    # Read the inversion file:
    inversion = np.load(inversion_file,encoding='latin1',allow_pickle=True)
    
    # It should have the shape: [ny, nx, ntau, npar]
    logtau = inversion[0][0][1][0][0]
    ntau = len(logtau)
    
    # The height and width of the map:
    ny = inversion.shape[0]*(inversion[0][-1][0][0]+1)
    nx = (inversion[0][-1][0][1]+1)
    
    # Create the file:
    modelmap = np.zeros((nx, ny, ntau, npar))
    for tau in tqdm(range(ntau)):
        for par in tqdm(range(npar), leave=False):
            modelmap[:, :, tau, par] = readSIRMap(inversion, par, tau)
    
    # Save the file:
    np.save(inversion_file[:-4]+'_model.npy', modelmap)
    

if __name__ == "__main__":
    filename = 'finalSIR.npy'
    create_modelmap(filename)
    create_profilemap(filename)
