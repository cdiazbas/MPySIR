from mpi4py import MPI
from sirtools import lmodel8, wmodel8, wmodel12
from sirtools import lperfil
import matplotlib.pyplot as plt
import numpy as np
import os

"""
This module contains functions to run SIR and modify the SIR files.
"""


#=============================================================================
def sirexe(fila, columna, myHeight, rank, sirfile, modeloFin, resultadoSir, sirmode, chi2map = True):
    """
    Runs SIR for a given pixel.
    """
    
    # By default, we use 3 cycles so the final name of the model is hsraB_3.mod
    finalProfile = 'hsraB_3.per'

    if sirmode == 'gammaV' or sirmode == 'gammVaddFullProfile':
        gammaV()

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
def write_continue_model(tau_init, model_init, continue_model, final_filename='hsraB.mod'):
    """
    Writes an input model (from a previous inversion) to be used as a starting model
    """
    # Modify the model:
    model_init[0] = continue_model[:,1] # temperature
    model_init[1] = continue_model[:,2] # electron pressure
    model_init[2] = continue_model[:,3] # microturbulence
    model_init[3] = continue_model[:,4] # magnetic field
    model_init[4] = continue_model[:,5] # LOS velocity
    model_init[5] = continue_model[:,6] # inclination
    model_init[6] = continue_model[:,7] # azimuth
    model_init[7] = continue_model[0,8] # macro velocity
    model_init[8] = continue_model[0,9] # filling factor
    model_init[9] = continue_model[0,10] # stray light
    wmodel12([tau_init, model_init], 'hsraB.mod', verbose=False)




#=============================================================================
def addFullProfile(sirfile):
    """
    Having the option of synthetic profiles with other properties
    """
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


# ========================================================================================================
def plotper(main_file='data.per',
            synth_file='hsraB_3.per',
            color1='k',
            color2='m',
            y_range_max=np.array([1.1, 3., 3., 3., 3.])):
    """
    Plot the observed and synthetic profiles.
    """

    y_range_min = -y_range_max
    y_range_min[0] = 0

    # Load data
    x0, stokes0, [num_lines, pos, num_points] = lperfil(main_file)
    if num_lines == 1:
        x, stokes, [_, _, _] = lperfil(synth_file)

    # Prepare positions for slicing
    pos_new = list(pos) + [len(num_points) - 1]
    x0A, xA = np.array(x0) / 1000., np.array(x) / 1000.

    # Initialize the figure
    plt.figure(figsize=(15, 5 * num_lines))

    # Iterate through lines and Stokes parameters
    for line_idx in range(num_lines):
        for sParam in range(4):
            plt_idx = line_idx * 4 + sParam + 1
            plt.subplot(num_lines, 4, plt_idx)

            # Plot the main data
            x_range = slice(pos_new[line_idx], pos_new[line_idx + 1] - 1)
            data = stokes0[sParam][x_range]
            if sParam > 0:  # Apply scaling factor of 100 to Q, U, and V parameters
                data = data * 100
            plt.plot(x0A[x_range], data, color1, lw=1.0)

            # Customize the plot
            plt.xlabel(r'$\Delta\lambda$ [$\AA$]', fontsize=15)
            plt.ylabel(['I/Ic', 'Q/Ic [%]', 'U/Ic [%]', 'V/Ic [%]'][sParam], fontsize=15)
            plt.xlim(x0A[pos_new[line_idx]], x0A[pos_new[line_idx + 1] - 1])
            plt.ylim(y_range_min[sParam], y_range_max[sParam])
            plt.grid(alpha=0.2, linestyle='-')
            plt.locator_params(axis='both', nbins=4)

            # Plot synthetic profiles
            synth_data = stokes[sParam][x_range]
            if sParam > 0:  # Apply scaling factor of 100 to Q, U, and V parameters
                synth_data = synth_data * 100
            plt.plot(xA[x_range], synth_data, color2, lw=1.0)

    # Save the figure
    plt.tight_layout()
    output_file = 'P' + synth_file[:-4] + '.pdf'
    plt.savefig(output_file, bbox_inches='tight')
    print(output_file + ':: SAVED')




#=============================================================================
def plotmfit(main_file='hsraB.mod',
             synth_file='hsraB_3.mod',
             error_model=True,
             indices_to_plot=[0, 3, 4, 5],
             labels=['$T$ $[K]$', r'$P_e$' + ' [dyn cm^-3]', r'$v_{mic}$' + ' [cm/s]', '$B$ $[G]$', r'$v_{LOS}$' + ' $[m/s]$', r'$\gamma$ $[deg]$'],
             color1='k',
             color2='m',
             margin=[0.2, 0.3, 0.3, 0.3]):
    """
    Plot the initial and final model parameters.
    """

    num_plots = len(indices_to_plot)

    # Load data
    tau, data = lmodel8(main_file, verbose=False)
    tau2, data2 = lmodel8(synth_file, verbose=False)
    if error_model:
        _, error_data = lmodel8(synth_file[:-4] + '.err')

    plt.figure(figsize=(4 * num_plots, 5))

    for i, index in enumerate(indices_to_plot):
        # Plot the data
        plt.subplot(1, num_plots, i + 1)
        quantity0, quantity1 = data[index], data2[index]
        plt.plot(tau, quantity0, color1, lw=1.0)
        plt.plot(tau2, quantity1, color2, lw=1.0)

        # Plot the error model data
        if error_model:
            error_quantity = error_data[index]
            plt.fill_between(tau2, quantity1 - error_quantity, quantity1 + error_quantity, facecolor='m', alpha=0.2)

        # Set the limits for x and y axes with margin adjustment
        min_val, max_val = min(min(quantity0), min(quantity1)), max(max(quantity0), max(quantity1))
        margin_adjustment = abs(min_val - max_val) * margin[i]
        plt.gca().update(dict(xlim=(min(tau2), max(tau2)), ylim=(min_val - margin_adjustment, max_val + margin_adjustment)))

        # Customize the plot
        plt.tick_params(axis='both', direction='in')
        plt.xlabel(r'$log(\tau)$', fontsize=20)
        plt.ylabel(labels[index], fontsize=20)
        plt.locator_params(axis='both', nbins=4)
        plt.grid(alpha=0.2, linestyle='-')

    # Save the figure
    plt.tight_layout()
    output_file = 'M' + synth_file[:-4] + '.pdf'
    plt.savefig(output_file, bbox_inches='tight')
    print(output_file + ':: SAVED')


#=============================================================================
def cerca(number, array):
	from numpy import argmin, abs
	indice = argmin(abs(number-array))
	return indice


#=============================================================================
def gammaV():
    """
    Compute the inclination from the Stokes V profile to have a good guess
    """
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


