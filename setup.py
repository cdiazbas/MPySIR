###############################################################
#  MPySIR: MPI python script for SIR
#
#  CALL:    mpirun -n 10 python setup.py
##############################################################

"""
# Author: cdiazbas@iac.es
# Date: 09.06.2015
# Version: 1.7
"""

"""
Archivos necesarios:
====================
imagen.npy / imagen.fits
xLambda.npy


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
from sirutils import *
import pySir.sirtools as sirtools
import pyLib.imtools as imtools

# ================================================= MPI INIT - CLEAN
comm = MPI.COMM_WORLD
widthT = 1
if comm.rank == 0:
	(widthT, heightT) = getTerminalSize()
	print('-'*widthT)
	print('Running on %d cores' % comm.size)
	print('-'*widthT)
	try:
		pncore()
	except:
		pass
	from clean import clean; clean()
comm.Barrier()


# We prepare the directory with all necesary
os.system('cp -R invDefault node'+str(comm.rank))
print('Node '+str(comm.rank)+' prepared.')
comm.Barrier()


# ================================================= INPUT
imagefits = 'mapa34B2.npy'
dictLines = {'SiI':300}  # Line Number in LINEAS file
rango = range(0,120) # Where the spectral line is
modeloInit = 'hsraB.mod'
modeloFin = 'hsraB_3.mod'
sirmode = 'gammVaddFullProfile'
chi2map = True
lambdaRef = 10827.110

# ================================================= LOAD DATA
xlambda = np.load('xLambdaBin.npy')
x = (xlambda[rango] -lambdaRef)*1e3      # Longitud de onda en mA

FileName = imagefits
imagen = np.load(imagefits)
#imagen = np.copy(imagen1[:,:,48:50,:])
height, nStokes, width, nLambdas = imagen.shape
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
myPart = np.copy(imagen[myHeight*comm.rank:myHeight*comm.rank+myHeight,:,:,:])
print(myPart.shape)


comm.Barrier()
pprint('==> Data loaded ..... {0:2.3f} s'.format(time.time() - start_time))
pprint('-'*widthT)


# ================================================= INVERSION
curr = os.getcwd()  # Current path
os.chdir(curr+'/node'+str(comm.rank))  # enter to each node_folder
# # os.system('mkdir node_'+str(comm.rank))

comm.Barrier()
pprint('==> Inside! ..... {0:2.3f} s'.format(time.time() - start_time))


import platform; _platform = platform.system() # We execute SIR according to the OS:
if _platform == "Linux": # Linux OS
	sirfile = 'sir.x'
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

for fila in range(0,myHeight):
	for columna in range(0,width):
		mapa = myPart[fila,:,columna,:]
		stokes = [mapa[0,rango],mapa[1,rango],mapa[2,rango],mapa[3,rango]]
		sirtools.wperfil('data.per',dictLines['SiI'],x,stokes)

		sirexe(fila,columna,myHeight,comm.rank,sirfile, modeloFin, resultadoSir, sirmode, chi2map)

		# plotper()  # Realiza los plots de los perfiles de Stokes
		# plotmfit() # Realiza los plots de los modelos ajustados
		
	actualPixel = (fila+1)*(width)
	if comm.rank == 0: print('... {0:4.2f} % ...'.format(float(actualPixel)/float(totalPixel)*100.))

#print('Soy {0} a tiempo {1:2.3f} s'.format(comm.rank,time.time() - start_time))
comm.Barrier()
pprint('==> Inversion terminada ..... {0:2.3f} s'.format(time.time() - start_time))
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
	np.save('finalSir.npy', finalSir)

comm.Barrier()
pprint('==> MPySIR <== ..... {0:2.3f} s'.format(time.time() - start_time))
pprint('-'*widthT)

