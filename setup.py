###############################################################
#  MPySIR: MPI python script for SIR
#
#  CALL:    mpirun -n 10 python setup.py
##############################################################

"""
SCRIPT: setup.py

SIRMODE:
====================
1.- 'perPixel'
2.- 'continue'

3.- 'beforePixel'
4.- 'gammaV'
5.- 'medianFilter'
6.- 'addFullProfile'
7.- 'gammVaddFullProfile'

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
imagefits = '../../data/test_field_1_630_0.fits'
original_axis = 'ns nw ny nx'
fov = None #'4,4' # We extract a 20x20 pixels for the inversion

wavefile = '../../data/wav.npy'
dictLines = {'atom':'200,201'}  # Line Number in LINEAS file
rango = range(0,401) # Range to invert
modeloInit = 'hsraB.mod'
modeloFin = 'hsraB_3.mod'
# sirmode = 'gammVaddFullProfile'
sirmode = 'perPixel' # 'continue'
continuemodel = 'finalSIR_model.npy'

chi2map = True
lambdaRef = 6301.5080
test1pixel = False
outputfile = 'finalSIR_cycle1.npy'


# ================================================= INVERSION PARAMETERS
Nodes_temperature = '2,3,5'
Nodes_magneticfield = '1,2,3'
Nodes_LOSvelocity = '1,2,3'
Nodes_gamma = '1,1,2'
Nodes_phi = '1,1,2'
Invert_macroturbulence = '0'
Initial_vmacro = 0.0 # km/s




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
        

    # We now swap the axes to [ny, ns, nx, nw] from original_axis:
    pprint('[INFO] Before - Image shape: '+str(image.shape))
    from einops import rearrange
    image = rearrange(image, original_axis+'-> ny ns nx nw')
    pprint('[INFO] After - Image shape: '+str(image.shape))


    # If fov is not None, we extract a portion of the image:
    if fov is not None:
        image = image[0:int(fov.split(',')[0]),:,0:int(fov.split(',')[1]),:]

    # Data dimensions:
    height, nStokes, width, nLambdas = image.shape

    # Now we divide the image in portions and send them to the nodes:
    listofheights = np.arange(height)   
    listofparts = np.array_split(listofheights, comm.size)
    
    # Print the list of parts for each node:
    for nodei in range(comm.size):
        print('[INFO] Node '+str(nodei)+' will receive: from column '+str(listofparts[nodei][0])+' to '+str(listofparts[nodei][-1]))

    # Divide the image in small portions and broadcast them to the nodes:
    for nodei in range(1, comm.size):
        myHeight = listofparts[nodei]
        comm.send(image[myHeight,:,:,:], dest = nodei, tag = 100)

    # The master node keeps the first portion:
    myHeight = listofparts[0]
    myPart = np.copy(image[myHeight,:,:,:])
    del image # We delete the original image to save memory.
    
    
    # ================================================= LOAD CONTINUE MODEL
    if sirmode == 'continue' and continuemodel is not None:
        init_model = np.load(continuemodel)
        # If it was created with inv2model routine it should have the axes: [ny, nx, ntau, npar]
        
        # We now broadcast the model to the rest of the nodes like the data:
        for nodei in range(1, comm.size):
            myHeight = listofparts[nodei]
            comm.send(init_model[myHeight,:,:,:], dest = nodei, tag = 200)
            
        # The master node keeps the first portion:
        myHeight = listofparts[0]
        myInit_model = np.copy(init_model[myHeight,:,:,:])
        del init_model # We delete the original image to save memory.
        
    
if comm.rank != 0:
    
    # ================================================= LOAD INPUT DATA
    # The rest of the nodes receive their portion:
    myPart = comm.recv(source = 0, tag = 100)
    print('Node '+str(comm.rank)+' received data -> '+str(myPart.shape),flush=True)

    # ================================================= LOAD CONTINUE MODEL
    if sirmode == 'continue' and continuemodel is not None:
        myInit_model = comm.recv(source = 0, tag = 200)
        print('Node '+str(comm.rank)+' received new init model -> '+str(myInit_model.shape),flush=True)


comm.Barrier()
pprint('==> Data loaded ..... {0:2.3f} s'.format(time.time() - start_time))
pprint('-'*widthT)


# ================================================= COPY FILES
# We prepare the directory with all necesary
os.system('cp -R invDefault node'+str(comm.rank))
print('Node '+str(comm.rank)+' prepared.',flush=True)
comm.Barrier()
time.sleep(0.1) # We wait a bit to avoid problems with the prints


# ================================================= INVERSION
curr = os.getcwd()  # Current path
os.chdir(curr+'/node'+str(comm.rank))  # enter to each node_folder

comm.Barrier()
pprint('==> Inside! ..... {0:2.3f} s'.format(time.time() - start_time))



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
myHeight = myPart.shape[0]
width = myPart.shape[2]
totalPixel = myHeight*width
if comm.rank == 0: print('... {0:4.2f} % ...'.format(0./totalPixel*100))

if test1pixel:
    myHeight = 1
    width = 1
    comm.Barrier()
    pprint('==> Test 1 pixel ..... {0:2.3f} s'.format(time.time() - start_time))


# Ensure that the first column has a single line when writting the "data.per" file
if len(dictLines['atom'].split(',')) > 1:
    wperfilLine = dictLines['atom'].split(',')[0]
else:
    wperfilLine = dictLines['atom']

if sirmode == 'continue' and continuemodel is not None:
    # Use the hsraB.mod as the baseline model for changing the initial model:
    [tau_init, model_init] = sirtools.lmodel12('hsraB.mod')
    


# ================================================= INVERSION
# Start the inversion:
for fila in range(0,myHeight):
    for columna in range(0,width):
        mapa = myPart[fila,:,columna,:]
        stokes = [mapa[0,rango],mapa[1,rango],mapa[2,rango],mapa[3,rango]]
        sirtools.wperfil('data.per',wperfilLine,x,stokes)
        
        if sirmode == 'continue' and continuemodel is not None:
            # We write the initial model as hsraB.mod which is the default name for the initial model in SIR [ny, nx, ntau, npar]
            init_pixel = init_model[fila,columna,:,:]
    
            model_init[0] = init_pixel[:,1]
            model_init[1] = init_pixel[:,2]
            model_init[2] = init_pixel[:,3]
            model_init[3] = init_pixel[:,4]
            model_init[4] = init_pixel[:,5]
            model_init[5] = init_pixel[:,6]
            model_init[6] = init_pixel[:,7]
            model_init[7] = vmac[0,8]
            model_init[8] = fill[0,9]
            model_init[9] = stray[0,10]
            # new_model_init[10] = zz
            # new_model_init[10] = pgas
            # new_model_init[11] = rho
            
            wmodel12([tau_init, model_init], filename, verbose=False)
            [tau, todoPlot]
        
        
        # Run SIR:
        sirutils.sirexe(fila,columna,myHeight,comm.rank,sirfile, modeloFin, resultadoSir, sirmode, chi2map)

        if test1pixel:
            sirutils.plotper()  # Realiza los plots de los perfiles de Stokes
            sirutils.plotmfit() # Realiza los plots de los modelos ajustados
            
    actualPixel = (fila+1)*(width)
    if comm.rank == 0: print('... {0:4.2f} % ...'.format(float(actualPixel)/float(totalPixel)*100.))

comm.Barrier()
pprint('==> Done ..... {0:2.3f} s'.format(time.time() - start_time))
pprint('-'*widthT)


# ================================================= SAVING RESULTS
resultadoSir = np.array(resultadoSir, dtype=object)
np.save('modelos.npy', resultadoSir)


if comm.rank != 0:
	comm.send(resultadoSir, dest = 0 , tag = 0)
else:
	finalSir = []
	finalSir.append(resultadoSir)
	for nodei in range(1, comm.size):
		vari = comm.recv(source = nodei, tag = 0)
		finalSir.append(vari)
	os.chdir(curr)
	np.save(outputfile, finalSir)

comm.Barrier()
pprint('==> MPySIR <== ..... {0:2.3f} s'.format(time.time() - start_time))
pprint('-'*widthT)

