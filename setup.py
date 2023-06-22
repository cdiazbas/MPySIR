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
2.- 'beforePixel'
3.- 'gammaV'
4.- 'medianFilter'
5.- 'addFullProfile'
6.- 'gammVaddFullProfile'
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
wavefile = '../../data/wav.npy'
dictLines = {'atom':'200,201'}  # Line Number in LINEAS file
rango = range(0,401) # Range to invert
modeloInit = 'hsraB.mod'
modeloFin = 'hsraB_3.mod'
# sirmode = 'gammVaddFullProfile'
sirmode = 'perpixel'
chi2map = True
lambdaRef = 6301.5080
test1pixel = False
continueInversion = None # If not None, it should be a string with the name of the file to continue the inversion

# ================================================= LOAD DATA
xlambda = np.load(wavefile)
x = (xlambda[rango] -lambdaRef)*1e3  # Wavelength in mA


# Modify the "malla.grid" file to change the wavelength range:
# Only if master:
if comm.rank == 0:
    sirutils.modify_malla(dictLines, x)

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
    pprint('[Before] Image shape: '+str(image.shape))
    
    # Swap axes from [ns, nw, ny, nx] to [ny, ns, nx, nw]
    image = image.transpose(2,0,3,1)
    pprint('[After] Image shape: '+str(image.shape))

    
else:
    image = np.load(imagefits)

# Data dimensions:
height, nStokes, width, nLambdas = image.shape
pprint((height, nStokes, width, nLambdas))

divImage = height % comm.size
if divImage == 0:
	pprint('Height is exactly divisible')
else :
	pprint('Not divisible')
	import sys
	sys.exit()


myHeight = int(height / comm.size)
pprint('==> Dividida en '+str(comm.size)+' trozos de altura '+str(myHeight)+'.')
myPart = np.copy(image[myHeight*comm.rank:myHeight*comm.rank+myHeight,:,:,:])
print(myPart.shape)


comm.Barrier()
pprint('==> Data loaded ..... {0:2.3f} s'.format(time.time() - start_time))
pprint('-'*widthT)


# ================================================= COPY FILES
# We prepare the directory with all necesary
os.system('cp -R invDefault node'+str(comm.rank))
print('Node '+str(comm.rank)+' prepared.')
comm.Barrier()



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


#from cjd_pysir.sirtools import lmodel8
#tau, magnitudes = lmodel8(modeloInit,verbose=False)
#comm.Barrier()


resultadoSir = []
totalPixel = myHeight*width
if comm.rank == 0: print('... {0:4.2f} % ...'.format(0./totalPixel*100))

if test1pixel:
    myHeight = 1
    width = 1
    comm.Barrier()
    pprint('==> Test 1 pixel ..... {0:2.3f} s'.format(time.time() - start_time))


# Ensure that the line is a single line when writting the "data.per" file
if len(dictLines['atom'].split(',')) > 1:
    wperfilLine = dictLines['atom'].split(',')[0]
else:
    wperfilLine = dictLines['atom']


# Start the inversion:
for fila in range(0,myHeight):
    for columna in range(0,width):
        mapa = myPart[fila,:,columna,:]
        stokes = [mapa[0,rango],mapa[1,rango],mapa[2,rango],mapa[3,rango]]
        sirtools.wperfil('data.per',wperfilLine,x,stokes)

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
	np.save('finalSIR.npy', finalSir)

comm.Barrier()
pprint('==> MPySIR <== ..... {0:2.3f} s'.format(time.time() - start_time))
pprint('-'*widthT)

