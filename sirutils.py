###############################################################
#  MPySIR: MPI python script for SIR
#
#  SCRIPT: sirutils.py
###############################################################


from mpi4py import MPI
from sirtools import lmodel8, wmodel8
from sirtools import lperfil
import matplotlib.pyplot as plt
import numpy as np
import os


#=============================================================================
def sirexe(fila, columna, myHeight, rank, sirfile, modeloFin, resultadoSir, sirmode, chi2map = True):

    finalProfile = 'hsraB_3.per'

    if sirmode == 'gammaV' or sirmode == 'gammVaddFullProfile':
        gammaV()

    if sirmode == 'medianFilter':
        medianFilter(fila, columna)


    # We run SIR
    import os
    os.system('echo sir.trol | '+sirfile+' > pylog.txt')

    # We now read the model parameters after the fitting:
    tau, magnitudes = lmodel8(modeloFin,verbose=False)
    ERRtau, ERRmagnitudes = lmodel8(modeloFin[0:-4]+'.err',verbose=False)
    magnitudes.insert(0,tau)
    ERRmagnitudes.insert(0,ERRtau)

    # We add the chi2 value to the list of parameters:
    if chi2map:
        chifile = open('sir.chi','r')
        for line in chifile:
            pass
        chi2 = float(line.split()[1])
        lenmag = len(magnitudes)
        magnitudes.insert(lenmag,chi2)
        ERRmagnitudes.insert(lenmag,chi2)

    if sirmode == 'addFullProfile' or sirmode == 'gammVaddFullProfile':
        addFullProfile(sirfile)
        finalProfile = 'dataFull.per'

    # We read the final synthetic profiles after the fitting:
    xFull, stokesFull, [nL,posi,nN] = lperfil(finalProfile,verbose=False)
    perfiles = [xFull,stokesFull]
    modelos = [magnitudes, ERRmagnitudes]
    punto = [fila + myHeight*rank, columna]
    resultadoSir.append([punto,modelos,perfiles])

    if sirmode == 'beforePixel':
        os.system('rm hsraB.mod'); os.system('cp hsraB_3.mod hsraB.mod')



#=============================================================================
def modify_malla(dictLines, x):
    """
    Modifies the "malla.grid" file to change the wavelength range.
    """
    # Read the file:
    f = open('invDefault/malla_.grid','r')
    lines = f.readlines()
    f.close()

    # Only modify the last line 
    step = x[1]-x[0]
    space = 10*' '
    lines[-1] = '{0}     :     {1:6.4f},     {2:6.4f},     {3:6.4f}'.format(dictLines['atom'],x[0],step,x[-1])+'\n'
    print('[INFO] malla.grid updated: ',lines[-1][:-1])

    # Write the file:
    f = open('invDefault/malla.grid','w')
    f.writelines(lines)
    f.close()


#=============================================================================
def modify_sirtrol(Nodes_temperature, Nodes_magneticfield, Nodes_LOSvelocity, Nodes_gamma, Nodes_phi, Invert_macroturbulence):
    """
    Modifies the "sir.trol" file to change the number of nodes.
    """
    # Read the file:
    f = open('invDefault/sir_.trol','r')
    lines = f.readlines()
    f.close()

    """The lines appear in the following order:
    Nodes for temperature 1      :2,3,5
    Nodes for electr. press. 1   :0
    Nodes for microturb. 1       :0
    Nodes for magnetic field 1   :1,2,2
    Nodes for LOS velocity 1     :1,2,2
    Nodes for gamma 1            :1,1,2
    Nodes for phi 1              :1,1,2
    Invert macroturbulence 1?    :0
    """

    # Modify the lines:
    lines[14] = 'Nodes for temperature 1      :'+str(Nodes_temperature)+'\n'
    lines[17] = 'Nodes for magnetic field 1   :'+str(Nodes_magneticfield)+'\n'
    lines[18] = 'Nodes for LOS velocity 1     :'+str(Nodes_LOSvelocity)+'\n'
    lines[19] = 'Nodes for gamma 1            :'+str(Nodes_gamma)+'\n'
    lines[20] = 'Nodes for phi 1              :'+str(Nodes_phi)+'\n'
    lines[21] = 'Invert macroturbulence 1?    :'+str(Invert_macroturbulence)+'\n'

    # Write the file:
    f = open('invDefault/sir.trol','w')
    f.writelines(lines)
    f.close()
    print('[INFO] sir.trol updated')



#=============================================================================
def modify_vmacro(initial_vmacro):
    """
    Modifies the guess model with the initial macroturbulence.
    """
    # Read the file:
    f = open('invDefault/hsraB_.mod','r')
    lines = f.readlines()
    f.close()
    
    # The vmacro is in the index 0 (of 3 elements) of the 1st line:
    firstLine = lines[0].split()
    newLine = '  '+str(initial_vmacro)+'  '+firstLine[1]+'  '+firstLine[2]+'\n'
    lines[0] = newLine
    
    # Write the file:
    f = open('invDefault/hsraB.mod','w')
    f.writelines(lines)
    f.close()
    print('[INFO] Initial model updated with vmacro = ',initial_vmacro,' km/s')
    
    


#=============================================================================
def addFullProfile(sirfile):
	import os
	os.system('echo sirFull.trol | '+sirfile+' > pylogFull.txt')




#=============================================================================
def pncore():
		import platform; _platform = platform.system() # We execute SIR according to the OS:
		from subprocess import PIPE, Popen
		if _platform == "Linux": # Linux OS
				proceso = Popen(['nproc'], stdout=PIPE, stderr=PIPE)
				ncores = proceso.stdout.read().split('\n')[0]
				print('Available cores = '+ncores)
		elif _platform == "Darwin": # MAC OS X
				proceso = Popen(['sysctl','hw.ncpu'], stdout=PIPE, stderr=PIPE)
				ncores = proceso.stdout.read().split('\n')[0]
				print('Available cores = '+ncores.split(':')[-1])


#=============================================================================
def pprint(ini='', end='', comm=MPI.COMM_WORLD):
		"""Print for MPI parallel programs: Only rank 0 prints *str*."""
		if comm.rank == 0:
				print(str(ini)+end)

#=============================================================================
def getTerminalSize():
	import os
	env = os.environ
	def ioctl_GWINSZ(fd):
		try:
			import fcntl, termios, struct, os
			cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,
			'1234'))
		except:
			return
		return cr
	cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
	if not cr:
		try:
			fd = os.open(os.ctermid(), os.O_RDONLY)
			cr = ioctl_GWINSZ(fd)
			os.close(fd)
		except:
			pass
	if not cr:
		try:
			cr = (env.get('LINES', 25), env.get('COLUMNS', 80))
			# import os
			# cr[0],cr[1] = os.popen('stty size', 'r').read().split()
		except:
			pass
			# print("Default")
			# cr = (env.get('LINES', 25), env.get('COLUMNS', 80))
	return int(cr[1]), int(cr[0])


#=============================================================================
def plotper():
    """
    Plots the observed and synthetic profiles.
    """

    filename='hsraB_3.per'
    LineName = ['FeI 6301.5']#,'FeI 6302.5','SiI 10827.1','FeI 15648.5','FeI 15652.9']
    NumeroLineas = len(LineName)
    MainFile = 'data.per'
    Color2 = 'm'
    Color1='k'

    YRANGEMAX = np.array([1.1,3.,3.,3.,3.])
    YRANGEMIN = -YRANGEMAX; YRANGEMIN[0] = 0


    # Abrimos los ficheros:
    x0, stokes0, [nL,posi,nN] = lperfil(MainFile)

    if nL == 1:

        x0 = np.array(x0)
        StokeI0=stokes0[0]
        StokeQ0=stokes0[1]
        StokeU0=stokes0[2]
        StokeV0=stokes0[3]

        x, stokes, [nL,posi,nN] = lperfil(filename)

        x = np.array(x)
        StokeI=stokes[0]
        StokeQ=stokes[1]
        StokeU=stokes[2]
        StokeV=stokes[3]

    lennN = len(nN)
    NumeroLineas = nL
    PosiNn0T = list(posi); PosiNn0T.append(lennN-1)
    PosiNn1T = list(posi); PosiNn1T.append(lennN-1)
    x0A=x0/1000.
    xA=x/1000.


    sParam = 4
    plt.figure(figsize=(15,5))
    for Index in range(0,NumeroLineas):
        
        plt.subplot(NumeroLineas,sParam,Index+1)
        plt.plot(x0A[PosiNn0T[Index]:PosiNn0T[Index+1]-1],StokeI0[PosiNn0T[Index]:PosiNn0T[Index+1]-1],Color1,lw=1.0)
        plt.tick_params(axis='y', direction='in')
        plt.tick_params(axis='x', direction='in')
        plt.xlabel(r'$\Delta\lambda$ [$\AA$]', fontsize=15)
        plt.ylabel(r'$I/I_c$', fontsize=15)
        plt.xlim(x0A[PosiNn0T[Index]],x0A[PosiNn0T[Index+1]-1])
        plt.ylim(YRANGEMIN[0],YRANGEMAX[0])
        plt.grid(alpha=0.2,linestyle='-')
        plt.locator_params(axis = 'x', nbins = 4)
        plt.locator_params(axis = 'y', nbins = 6)


        plt.subplot(NumeroLineas,sParam,Index+1+NumeroLineas)
        plt.plot(x0A[PosiNn0T[Index]:PosiNn0T[Index+1]-1],100*StokeQ0[PosiNn0T[Index]:PosiNn0T[Index+1]-1],Color1,lw=1.0)
        plt.tick_params(axis='y', direction='in')
        plt.tick_params(axis='x', direction='in')
        plt.xlabel(r'$\Delta\lambda$ [$\AA$]', fontsize=15)
        plt.ylabel(r'$Q/I_c$ $[\%]$', fontsize=15)
        plt.ylim(YRANGEMIN[1],YRANGEMAX[1])
        plt.xlim(x0A[PosiNn0T[Index]],x0A[PosiNn0T[Index+1]-1])
        plt.grid(alpha=0.2,linestyle='-')
        plt.locator_params(axis = 'x', nbins = 4)
        plt.locator_params(axis = 'y', nbins = 6)



        plt.subplot(NumeroLineas,sParam,Index+2+NumeroLineas)
        plt.plot(x0A[PosiNn0T[Index]:PosiNn0T[Index+1]-1],100*StokeU0[PosiNn0T[Index]:PosiNn0T[Index+1]-1],Color1,lw=1.0)
        plt.tick_params(axis='y', direction='in')
        plt.tick_params(axis='x', direction='in')
        plt.xlabel(r'$\Delta\lambda$ [$\AA$]', fontsize=15)
        plt.ylabel(r'$U/I_c$ $[\%]$', fontsize=15)
        plt.grid(alpha=0.2,linestyle='-')
        plt.locator_params(axis = 'x', nbins = 4)
        plt.xlim(x0A[PosiNn0T[Index]],x0A[PosiNn0T[Index+1]-1])
        plt.ylim(YRANGEMIN[2],YRANGEMAX[2])
        plt.locator_params(axis = 'y', nbins = 6)

        


        plt.subplot(NumeroLineas,sParam,Index+3+NumeroLineas)
        plt.plot(x0A[PosiNn0T[Index]:PosiNn0T[Index+1]-1],100*StokeV0[PosiNn0T[Index]:PosiNn0T[Index+1]-1],Color1,lw=1.0)
        plt.tick_params(axis='y', direction='in')
        plt.tick_params(axis='x', direction='in')
        plt.xlabel(r'$\Delta\lambda$ [$\AA$]', fontsize=15)
        plt.ylabel(r'$V/I_c$ $[\%]$', fontsize=15)
        plt.ylim(YRANGEMIN[3],YRANGEMAX[3])
        plt.xlim(x0A[PosiNn0T[Index]],x0A[PosiNn0T[Index+1]-1])
        plt.grid(alpha=0.2,linestyle='-')
        plt.locator_params(axis = 'x', nbins = 4)
        plt.locator_params(axis = 'y', nbins = 6)


    # ========================================================================================================

        plt.subplot(NumeroLineas,sParam,Index+1)
        plt.plot(xA[PosiNn1T[Index]:PosiNn1T[Index+1]-1],StokeI[PosiNn1T[Index]:PosiNn1T[Index+1]-1],Color2,lw=1.0)

        plt.subplot(NumeroLineas,sParam,Index+1+NumeroLineas)
        plt.plot(xA[PosiNn1T[Index]:PosiNn1T[Index+1]-1],100*StokeQ[PosiNn1T[Index]:PosiNn1T[Index+1]-1],Color2,lw=1.0)


        plt.subplot(NumeroLineas,sParam,Index+2+NumeroLineas)
        plt.plot(xA[PosiNn1T[Index]:PosiNn1T[Index+1]-1],100*StokeU[PosiNn1T[Index]:PosiNn1T[Index+1]-1],Color2,lw=1.0)

        plt.subplot(NumeroLineas,sParam,Index+3+NumeroLineas)
        plt.plot(xA[PosiNn1T[Index]:PosiNn1T[Index+1]-1],100*StokeV[PosiNn1T[Index]:PosiNn1T[Index+1]-1],Color2,lw=1.0)

    # ========================================================================================================

    plt.tight_layout()

    plt.savefig('P'+filename[0:-4]+'.pdf', bbox_inches='tight')#, pad_inches=0)
    print('P'+filename[0:-4]+'.pdf'+':: SAVED')

    return

#=============================================================================
def plotmfit():
# def plotmfit(initModel,outModel,errorModel=False,whichMag,verbose=False):
	import matplotlib.pyplot as plt
	import os
	import numpy as np

	# import seaborn as sns
	# sns.set_style("ticks")

	filename =  'hsraB_3.mod'
	errorModel = True
	#====================================================================

	#; 0:= temp , 1:= pres, 2:= vmic, 3:= B, 4:= vlos 5:=gamma
	PosiPlot = [0,3,4,5]

	LabelPlot= ['$T$ $[K]$',r'$P_e$'+' [dyn cm^-3]',r'$v_{mic}$'+' [cm/s]','$B$ $[G]$',r'$v_{LOS}$'+' $[m/s]$',r'$\gamma$ $[deg]$']
	TEXTAU = r'$\tau$'

	NumPlots = len(PosiPlot)

	MainFile = 'hsraB.mod'

	Color1='k'
	Color2 = 'm'


	tau, TodoPlot = lmodel8(MainFile,verbose=False)
	tau2, TodoPlot2 = lmodel8(filename,verbose=False)


	if errorModel:
		ERRtau2, TodoPlot3 = lmodel8(filename[0:-4]+'.err')



	MMargen=[0.2,0.3,0.3,0.3]

	plt.figure(figsize=(4*NumPlots,5))
	for i in range(0,len(PosiPlot)):
		
		plt.subplot(1,len(PosiPlot),i+1)
		Cantidad0=TodoPlot[PosiPlot[i]]
		Cantidad1=TodoPlot2[PosiPlot[i]]
		plt.plot(tau,Cantidad0,Color1,lw=1.0)
		plt.plot(tau2,Cantidad1,Color2,lw=1.0)
		
		if errorModel:
			ERRCantidad1=TodoPlot3[PosiPlot[i]]
			plt.fill_between(tau2, Cantidad1-ERRCantidad1, Cantidad1+ERRCantidad1,facecolor='m', alpha=0.2)
		#ylim(YRANGEMIN[Index],YRANGEMAX[Index])
		plt.xlim(min(tau2),max(tau2))
		
		Min0 = min(Cantidad0)
		Min1 = min(Cantidad1)
		Max0 = max(Cantidad0)
		Max1 = max(Cantidad1)

		if Min0 == Min1 : MinAbs = Min0
		if Min0 != Min1 : MinAbs = min([Min0,Min1])
		if Max0 == Max1 : MaxAbs = Max0
		if Max0 != Max1 : MaxAbs = max([Max0,Max1])
		NuevoMargen= abs(MinAbs-MaxAbs)*MMargen[i]
		plt.ylim(MinAbs-NuevoMargen,MaxAbs+NuevoMargen)
		#ylim(0.,2.5)

		plt.tick_params(axis='y', direction='in')#,labelsize=20)
		plt.tick_params(axis='x', direction='in')#,labelsize=20)
		plt.xlabel(r'$log(\tau)$', fontsize=20)
		plt.ylabel(LabelPlot[PosiPlot[i]], fontsize=20)
		plt.locator_params(axis = 'y', nbins = 6)
		plt.locator_params(axis = 'x', nbins = 4)
		plt.grid(alpha=0.2,linestyle='-')



	plt.tight_layout()
	# plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.4, hspace=0.3)

	plt.savefig('M'+filename[0:-4]+'.pdf', bbox_inches='tight')#, pad_inches=0)
	print('M'+filename[0:-4]+'.pdf'+':: GUARDADO')
	return


#=============================================================================
def cerca(number, array):
	from numpy import argmin, abs
	indice = argmin(abs(number-array))
	return indice


#=============================================================================
def gammaV():
	from scipy import integrate
	from scipy.interpolate import interp1d

	MainFile = 'data.per'

	x0, stokes0, [nL,posi,nN] = lperfil(MainFile)

	indice = cerca(0.,x0)
	distmin = min([abs(indice),abs(len(x0)-indice)])
	centro = indice
	x = x0
	y = stokes0[3]
	xlobuloazul = x[centro+1-distmin:centro+1]  
	ylobuloazul = y[centro+1-distmin:centro+1] 
	xlobulorojo = x[centro:centro+distmin] 
	ylobulorojo = y[centro:centro+distmin]
	int_roja = integrate.simps(ylobulorojo,xlobulorojo)
	int_azul = integrate.simps(ylobuloazul,xlobuloazul)

	if int_azul > int_roja: gamma = 45.0
	if int_azul < int_roja: gamma = 135.0

	from numpy import ones
	tau, magnitudes = lmodel8('hsraB.mod',verbose=False)
	modelo = [tau, magnitudes]
	magnitudes[5] = gamma*ones(len(magnitudes[5]))
	wmodel8(modelo,'hsraB.mod',verbose=False)


