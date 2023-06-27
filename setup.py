###############################################################
#  MPySIR: MPI python script for SIR
#
#  CALL:    mpirun -n 10 python setup.py
##############################################################

"""
This script runs the SIR inversion in parallel using MPI.

SIRMODE:
====================
1.- 'perPixel': Each pixel is inverted independently
2.- 'continue': The inversion is continued from a previous one

# Not tested yet:
3.- 'gammaV': The inclination is inferred form the Stokes V profile
4.- 'addFullProfile': The synthetic profiles are in another grid

"""

# ================================================= TIME
import time; start_time = time.time()

# ================================================= LIBRARIES
import numpy as np
from mpi4py import MPI
import os
import sirutils
from sirutils import pprint
import sirtools
import sys
import time
import datetime


# ================================================= MPI INIT - CLEAN
comm = MPI.COMM_WORLD
widthT = 1
if comm.rank == 0:
	(widthT, heightT) = sirutils.getTerminalSize()
	print('-'*widthT)
	print('Running on %d cores' % comm.size)
	print('-'*widthT)
	try:
		pncore()
	except:
		pass
	from clean import clean; clean()
comm.Barrier()


# ================================================= INPUT
imagefits = '../../data/sunspot_jmb_sir_synth.fits'
original_axis = 'ns nw ny nx'
fov = None #'4,4' # We extract a 20x20 pixels for the inversion
skip = 8 # We skip pixels in the 2D grid to reduce the number of pixels

wavefile = '../../data/wav.npy'
dictLines = {'atom':'200,201'}  # Line Number in LINEAS file
rango = range(0,401) # Range to invert
modeloFin = 'hsraB_3.mod'
sirmode = 'perPixel' # 'continue'
sirmode = 'continue' # 'continue'
continuemodel = 'finalSIR_cycle1_model_smoothed.npy'

chi2map = True # By default, we compute the chi2 map
lambdaRef = 6301.5080 # This is extracted from LINEAS file
test1pixel = False # Test the inversion of one pixel
outputfile = 'finalSIR_cycle1.npy'
verbose = True


# ================================================= INVERSION PARAMETERS
Nodes_temperature = '2,3,5'
Nodes_magneticfield = '1,2,3'
Nodes_LOSvelocity = '1,2,3'
Nodes_gamma = '1,1,2'                   # Magnetic field inclination
Nodes_phi = '1,1,2'                     # Magnetic field azimuth
Invert_macroturbulence = '0'
Initial_vmacro = 0.0 # km/s


# Cycle 2:
# Nodes_temperature = '5,7,9'
# Nodes_magneticfield = '3,5,5'
# Nodes_LOSvelocity = '3,5,5'
# Nodes_gamma = '2,3,3'                   # Magnetic field inclination
# Nodes_phi = '2,3,3'                     # Magnetic field azimuth
# Invert_macroturbulence = '0'
# Initial_vmacro = 0.0 # km/s


# ================================================= LOAD DATA


# Load wavelength (this is the same for all nodes):
xlambda = np.load(wavefile)
x = (xlambda[rango] -lambdaRef)*1e3  # Wavelength in mA


# Updates malla.grid and sir.trol if master:
if comm.rank == 0:
    
    # Modify the "malla.grid" file to change the wavelength range.
    sirutils.modify_malla(dictLines, x)
    
    # Modify the "sir.trol" file to change the inversion parameters.
    sirutils.modify_sirtrol(Nodes_temperature, Nodes_magneticfield, Nodes_LOSvelocity, Nodes_gamma, Nodes_phi, Invert_macroturbulence)

    # Modify the initial model with the initial macro velocity:
    sirutils.modify_vmacro(Initial_vmacro)


# Now only the master node reads the data and broadcasts it to the rest of the nodes:
if comm.rank == 0:

    # ================================================= LOAD INPUT DATA
    # Check if image is fits or npy:
    if imagefits[-3:] == 'npy':
        isfits = False
    elif imagefits[-4:] == 'fits':
        isfits = True
    else:
        print('ERROR: Image format not recognized.')
        sys.exit()


    # Load image:
    FileName = imagefits
    if isfits:
        from astropy.io import fits
        image = fits.open(imagefits)[0].data
    else:
        image = np.load(imagefits)
        

    # We now swap the axes to [ny, nx, ns, nw]:
    pprint('[INFO] Before - Image shape: '+str(image.shape))
    from einops import rearrange
    image = rearrange(image, original_axis+'-> ny nx ns nw')
    pprint('[INFO] After - Image shape: '+str(image.shape))


    # If fov is not None, we extract a portion of the image:
    if fov is not None:
        image = image[0:int(fov.split(',')[0]),0:int(fov.split(',')[1]),:,:]
    
    # If skip is not 1, we skip pixels:
    if skip != 1:
        image = image[::skip,::skip,:,:]

    # Data dimensions:
    height, width, nStokes, nLambdas = image.shape

    # Now we divide the image in portions and send them to the nodes. For that,
    # we can flatten the X&Y dimensions and divide them with array_split:
    totalpixels = height*width
    print('[INFO] Total pixels: '+str(totalpixels))
    listofpixels = np.arange(totalpixels)
    listofparts = np.array_split(listofpixels, comm.size)
    
    # We need to flatten the image to send it to the nodes, by
    # moving the X&Y dimensions to the first axis and then flattening:
    image = image.reshape((totalpixels, nStokes, nLambdas))    
    
    # Print the list of parts for each node:
    for nodei in range(comm.size):
        if verbose:
            print('[INFO] Node '+str(nodei)+' will receive: from pixel '+str(listofparts[nodei][0])+' to '+str(listofparts[nodei][-1]))

    # Divide the image in small portions and broadcast them to the nodes:
    for nodei in range(1, comm.size):
        myrange = listofparts[nodei]
        comm.send(image[myrange,:,:], dest = nodei, tag = 100)

    # The master node keeps the first portion:
    myrange = listofparts[0]
    myPart = np.copy(image[myrange,:,:])
    del image # We delete the original image to save memory.
    if verbose: print('Node 0 received data -> '+str(myPart.shape))

    
    # ================================================= LOAD CONTINUE MODEL
    if sirmode == 'continue' and continuemodel is not None:
        init_model = np.load(continuemodel)
        # If it was created with inv2model routine it should have the axes: [ny, nx, ntau, npar]
        
        # We now also flatten the model to send it to the nodes, by
        # flattening, as model is already in the correct axes:
        init_model = init_model.reshape((init_model.shape[0]*init_model.shape[1], init_model.shape[2], init_model.shape[3]))
        
        # We now broadcast the model to the rest of the nodes like the data:
        for nodei in range(1, comm.size):
            myrange = listofparts[nodei]
            comm.send(init_model[myrange,:,:], dest = nodei, tag = 200)
            
        # The master node keeps the first portion:
        myrange = listofparts[0]
        myInit_model = np.copy(init_model[myrange,:,:])
        if verbose: print('Node 0 received new init model -> '+str(myInit_model.shape))
        del init_model # We delete the original image to save memory.
        
    
if comm.rank != 0:
    
    # ================================================= LOAD INPUT DATA
    # The rest of the nodes receive their portion:
    myPart = comm.recv(source = 0, tag = 100)
    if verbose:
        print('Node '+str(comm.rank)+' received data -> '+str(myPart.shape),flush=True)

    # ================================================= LOAD CONTINUE MODEL
    if sirmode == 'continue' and continuemodel is not None:
        myInit_model = comm.recv(source = 0, tag = 200)
        if verbose:
            print('Node '+str(comm.rank)+' received new init model -> '+str(myInit_model.shape),flush=True)


comm.Barrier()
pprint('==> Data loaded ..... {0:2.3f} s'.format(time.time() - start_time))
pprint('-'*widthT)


# ================================================= COPY FILES
# We prepare the directory with all necesary
os.system('cp -R invDefault node'+str(comm.rank))
# if verbose: print('Node '+str(comm.rank)+' prepared.',flush=True)
comm.Barrier()
time.sleep(0.1) # We wait a bit to avoid problems with the prints


# ================================================= INVERSION
curr = os.getcwd()  # Current path
os.chdir(curr+'/node'+str(comm.rank))  # enter to each node_folder

comm.Barrier()
pprint('==> Ready to start! ..... {0:2.3f} s'.format(time.time() - start_time))



# We execute SIR according to the OS:
import platform
_platform = platform.system() 
if _platform == "Linux": # Linux OS
	sirfile = './sir.x'
	pprint('Platform: Linux OS')
elif _platform == "Darwin": # MAC OS X
	sirfile = './sir_mac.x'
	pprint('Platform: Mac OS X')


resultadoSir = []
totalPixel = myPart.shape[0]
if comm.rank == 0: print(f'\r... {0:4.2f} % ...'.format(0.0), end='', flush=True)


# We invert only one pixel to test the code:
if test1pixel:
    totalPixel = 1
    comm.Barrier()
    pprint('==> Testing 1 pixel ..... {0:2.3f} s'.format(time.time() - start_time))


# Ensure that the first column has a single line when writting the "data.per" file
if len(dictLines['atom'].split(',')) > 1:
    wperfilLine = dictLines['atom'].split(',')[0]
else:
    wperfilLine = dictLines['atom']


# If we are continuing a inversion from a previous model, we modify the initial model:
if sirmode == 'continue' and continuemodel is not None:
    # Use the hsraB.mod as the baseline model for changing the initial model:
    [tau_init, model_init] = sirtools.lmodel12('hsraB.mod')
    

# ================================================= INVERSION
# Start the inversion:
for currentPixel in range(0,totalPixel):
    mapa = myPart[currentPixel,:,:]
    stokes = [mapa[0,rango],mapa[1,rango],mapa[2,rango],mapa[3,rango]]
    sirtools.wperfil('data.per',wperfilLine,x,stokes)
    
    if sirmode == 'continue' and continuemodel is not None:
        # We write the initial model as hsraB.mod which is the default name for the initial model in SIR [ny, nx, ntau, npar]
        init_pixel = myInit_model[currentPixel,:,:]
        sirutils.write_continue_model(tau_init, model_init, init_pixel, final_filename='hsraB.mod')

    
    # Run SIR:
    sirutils.sirexe(0,0,0,comm.rank,sirfile, modeloFin, resultadoSir, sirmode, chi2map)

    if test1pixel:
        sirutils.plotper()  # Plots the profiles if we are testing 1 pixel
        sirutils.plotmfit() # Plots the model if we are testing 1 pixel
    
    # Print the percentage of the inversion:
    porcentage = float(currentPixel)/float(totalPixel)*100.
    if comm.rank == 0: print('\r... {0:4.2f} % ...'.format(porcentage), end='', flush=True)


comm.Barrier()
pprint('==> Inversion finished. Now gathering the results ..... {0:2.3f} s'.format(time.time() - start_time))
pprint('-'*widthT)


# ================================================= SAVING RESULTS
# The results are transformed to a numpy array to be able to save them using MPI:
resultadoSir = np.array(resultadoSir, dtype=object)

# We could save the results in a file if needed by uncommenting the following lines:
# np.save('modelos.npy', resultadoSir)


# Gather the results in the master node:
if comm.rank != 0:
    # We send the results to the master node:
	comm.send(resultadoSir, dest = 0 , tag = 0)
else:
    # We receive the results from the rest of the nodes:
    finalSir = []
    finalSir.append(resultadoSir)
    for nodei in range(1, comm.size):
        vari = comm.recv(source = nodei, tag = 0)
        finalSir.append(vari)
    os.chdir(curr)

    # Now that we have all the results for all the pixels, we can concatenate them
    # and reshape them to the original shape of the input data:
    finalSir = np.concatenate(finalSir, axis=0)
    finalSir = finalSir.reshape((height, int(finalSir.shape[0]/height),finalSir.shape[1],finalSir.shape[2]))

    np.save(outputfile, finalSir)

comm.Barrier()
pprint('==> MPySIR <==')

# Print the total time in the format HH:MM:SS (hours, minutes, seconds):
total_time = time.time() - start_time
if comm.rank == 0:
    print('Total time: '+str(datetime.timedelta(seconds=total_time)))
    print(' Output file: '+outputfile)
pprint('-'*widthT)


# We clean the directory with all the temporary files, except when we are testing the code:
if comm.rank == 0 and not test1pixel:
    clean()
    
