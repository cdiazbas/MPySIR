from mpi4py import MPI
from sirtools import lmodel8, wmodel8, wmodel12
from sirtools import lperfil
import matplotlib.pyplot as plt
import numpy as np
import os
from tqdm import tqdm
import glob
import sys

"""
This module contains functions to run SIR and modify the SIR files.
"""


#=============================================================================
def sirexe(rank, sirfile, resultadoSir, sirmode, chi2map = True, x=None):
    """
    Runs SIR for a given pixel.
    """

    # We run SIR
    os.system('echo sir.trol | '+sirfile+' > pylog.txt')

    # ====================================================
    # INVERSION MODE
    if sirmode != 'synthesis':
        # Last files written by SIR:
        file_list = glob.glob("*.mod")
        file_list.sort(key=os.path.getmtime, reverse=True)
        finalModel = os.path.basename(file_list[0])
        finalProfile = finalModel[0:-4]+'.per'
 
        tau, magnitudes = lmodel8(finalModel,verbose=False)
        ERRtau, ERRmagnitudes = lmodel8(finalModel[0:-4]+'.err',verbose=False)
        magnitudes.insert(0,tau)
        ERRmagnitudes.insert(0,ERRtau)    

        # We add the chi2 value to the list of parameters:
        if chi2map:
            try:
                with open('sir.chi', 'r') as file:
                    last_line = file.readlines()[-1]
                chi2 = float(last_line.split()[1])
            except:
                # If SIR does not write the chi2 file
                chi2 = 1000.0
            # Add the chi2 value to the list of parameters:
            lenmag = len(magnitudes)
            magnitudes.insert(lenmag,chi2)
            ERRmagnitudes.insert(lenmag,chi2)


        # We read the final synthetic profiles after the fitting:
        xFull, stokesFull, [nL,posi,nN] = lperfil(finalProfile,verbose=False)
        perfiles = [xFull,stokesFull]
        modelos = [magnitudes, ERRmagnitudes]
        punto = [0,0]
        resultadoSir.append([punto,modelos,perfiles])

        if sirmode == 'beforePixel':
            os.system('rm hsraB.mod'); os.system('cp hsraB_3.mod hsraB.mod')

        # Clean some files:
        try:
            os.remove('sir.chi')
        except:
            print('[INFO] sir.chi was not created by SIR.')



    # ====================================================
    # SYNTHESIS MODE
    else:
        # Here we are in synthesis mode, so we only read the synthetic profiles:
        finalProfile = 'data.per'
        try:
            xFull, stokesFull, [nL,posi,nN] = lperfil(finalProfile,verbose=False)
            perfiles = [xFull,stokesFull]
        except:
            print('[INFO] SIR failed to create the synthetic profiles.')
            # Create fake perfiles:
            xFull = x.copy()
            stokesFull = [np.zeros_like(x) for i in range(4)]
            perfiles = [xFull,stokesFull]

        punto = [0,0]
        resultadoSir.append([punto,punto,perfiles])






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
def modify_sirtrol_synthesis(Linesfile, Abundancefile, mu_obs):
    """
    Modifies the "sir.trol" file to change the number of nodes.
    """
    # Read the file:
    f = open('invDefault/sir_.trol','r')
    lines = f.readlines()
    f.close()

    # Modify the lines:
    lines[0] = 'Number of cycles             :'+str(0)+'\n'
    lines[5] = 'Atomic parameters file       :'+str(Linesfile)+'\n'
    lines[6] = 'Abundances file              :'+str(Abundancefile)+'\n'
    lines[32] = 'mu=cos (theta)               :'+str(mu_obs)+'\n'

    # Write the file:
    f = open('invDefault/sir.trol','w')
    f.writelines(lines)
    f.close()
    print('[INFO] sir.trol updated with Linesfile and Abundancefile.')


#=============================================================================
def modify_sirtrol(Nodes_temperature, Nodes_magneticfield, Nodes_LOSvelocity, Nodes_gamma, Nodes_phi, Invert_macroturbulence, Linesfile, Abundancefile,mu_obs,Nodes_microturbulence,weightStokes):
    """
    Modifies the "sir.trol" file to change the number of nodes.
    """
    # Read the file:
    f = open('invDefault/sir_.trol','r')
    lines = f.readlines()
    f.close()
    
    # Number of cycles:
    ncycles = np.max([len(Nodes_temperature.split(',')),len(Nodes_magneticfield.split(',')),len(Nodes_LOSvelocity.split(',')),
                      len(Nodes_gamma.split(',')),len(Nodes_phi.split(',')),len(Nodes_microturbulence.split(','))])
    print('[INFO] Number of cycles:',ncycles,'with weights:',weightStokes)

    # Modify the lines:
    lines[0] = 'Number of cycles             :'+str(ncycles)+'\n'
    lines[5] = 'Atomic parameters file       :'+str(Linesfile)+'\n'
    lines[6] = 'Abundances file              :'+str(Abundancefile)+'\n'
    lines[9] = 'Weight for Stokes I          :'+str(weightStokes.split(',')[0])+'\n'
    lines[10] = 'Weight for Stokes Q          :'+str(weightStokes.split(',')[1])+'\n'
    lines[11] = 'Weight for Stokes U          :'+str(weightStokes.split(',')[2])+'\n'
    lines[12] = 'Weight for Stokes V          :'+str(weightStokes.split(',')[3])+'\n'
    lines[14] = 'Nodes for temperature 1      :'+str(Nodes_temperature)+'\n'
    lines[16] = 'Nodes for microturb. 1       :'+str(Nodes_microturbulence)+'\n'
    lines[17] = 'Nodes for magnetic field 1   :'+str(Nodes_magneticfield)+'\n'
    lines[18] = 'Nodes for LOS velocity 1     :'+str(Nodes_LOSvelocity)+'\n'
    lines[19] = 'Nodes for gamma 1            :'+str(Nodes_gamma)+'\n'
    lines[20] = 'Nodes for phi 1              :'+str(Nodes_phi)+'\n'
    lines[21] = 'Invert macroturbulence 1?    :'+str(Invert_macroturbulence)+'\n'
    lines[32] = 'mu=cos (theta)               :'+str(mu_obs)+'\n'

    # Write the file:
    f = open('invDefault/sir.trol','w')
    f.writelines(lines)
    f.close()
    print('[INFO] sir.trol updated with nodes, Linesfile and Abundancefile.')


#=============================================================================
def modify_vmacro(initial_vmacro, filename_base = 'invDefault/hsraB_.mod', filename_final='invDefault/hsraB.mod', verbose=True):
    """
    Modifies the guess model with the initial macroturbulence.
    """
    # Read the file:
    f = open(filename_base,'r')
    lines = f.readlines()
    f.close()

    # The vmacro is in the index 0 (of 3 elements) of the 1st line:
    firstLine = lines[0].split()
    newLine = '  '+str(initial_vmacro)+'  '+firstLine[1]+'  '+firstLine[2]+'\n'
    lines[0] = newLine

    # Write the file:
    f = open(filename_final,'w')
    f.writelines(lines)
    f.close()
    
    if verbose:
        print('[INFO] Initial model updated with vmacro = ',initial_vmacro,' km/s')


#=============================================================================
def get_ntau():
    """
    Get the number of points in the wavelength axis
    """
    # Read the file:
    f = open('invDefault/hsraB.mod','r')
    lines = f.readlines()
    f.close()
    
    # The number of points is the lenght of the column:
    ntau = len(lines)-1
    return ntau

#=============================================================================
def calculate_divisors(n):
    # Calculate the divisors of n:
    divisors = []
    for i in range(1, n + 1):
        if n % i == 0:
            divisors.append(i)
    return divisors

#=============================================================================
def calculate_nodes():
    """
    The number of nodes is calculated as the minimum number of divisors of ntau-1
    """
    ntau = get_ntau()
    divisors = calculate_divisors(ntau-1)
    
    # Transform to array + 1:
    return np.array(divisors)+1

#=============================================================================
def modify_vmicro(initial_vmicro, filename_base = 'invDefault/hsraB.mod', filename_final='invDefault/hsraB.mod', verbose=True):
    """
    Modifies the guess model with the initial microturbulence.
    """
    # Read the file:
    f = open(filename_base,'r')
    lines = f.readlines()
    f.close()

    # The vmicro is the 4th column starting from the 2nd line:
    for i in range(1,len(lines)):
        line = lines[i].split()
        line[3] = '{0:1.2e}'.format(initial_vmicro)
        lines[i] = '  '.join(line)+'\n'

    # Write the file:
    f = open(filename_final,'w')
    f.writelines(lines)
    f.close()
    if verbose:
        print('[INFO] Initial model updated with vmicro = ',initial_vmicro/1e5,' km/s')    


#=============================================================================
def write_continue_model(tau_init, model_init, continue_model, final_filename='hsraB.mod', apply_constraints=True):
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
    
    if apply_constraints:
        # Temperature cannot be larger than 12000 K:
        loc = np.where(model_init[0] > 12000.0)[0]
        model_init[0][loc] = 12000.0

        # All the values are the same before -3 or after 0.5 (except for the temperature):    
        loc = np.where(tau_init < -3.0)[0]
        for j in range(3,7):
            model_init[j][loc] = model_init[j][loc[0]]
        
        loc = np.where(tau_init > 0.0)[0]
        for j in range(3,7):
            model_init[j][loc] = model_init[j][loc[-1]]
        
        # Smooth in 1D from astropy:
        from astropy.convolution import convolve, Gaussian1DKernel
        for j in range(3,7):
            model_init[j] = convolve(model_init[j], Gaussian1DKernel(2),boundary='extend')
        # Also for the temperature:
        model_init[0] = convolve(model_init[0], Gaussian1DKernel(2), boundary='extend')

    wmodel12([tau_init, model_init], final_filename, verbose=False)




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
            # For vmac, filling factor, stray-light and chi2 we take the single values:
            if parameter == 8 or parameter == 9 or parameter == 10 or parameter == 11:
                parmap[pix_x, pix_y] = outputSir[pix_y,pix_x][1][0][parameter]
            else:
                parmap[pix_x, pix_y] = outputSir[pix_y,pix_x][1][0][parameter][tau]
    return parmap.T


# ====================================================================
def create_modelmap(inversion, inversion_file, npar = 12):
    """
    It creates a file with the model parameters [ny, nx, ntau, npar] from the inversion 
    """

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
            
    # Make sure that the parameters are within the limits:
    par = 6 # The inclination angle
    modelmap[:, :, :, par] = np.clip(modelmap[:, :, :, par], 0.0, 180.0)
        
    # Save the file:
    if inversion_file[-4:] == '.npy':
        np.save(inversion_file[:-4]+'_model.npy', modelmap.astype(np.float32))
    else:
        np.save(inversion_file+'_model.npy', modelmap.astype(np.float32))
    

# ====================================================================
def readSIRProfileMap(outputSir, Nstoke):
    """
    It returns the synthetic profiles from the inversion for a given Stokes parameter
    """
    ny = outputSir.shape[0]
    nx = outputSir.shape[1]
    nwav = len(outputSir[0,0][2][0])
    syn_map = np.zeros((ny, nx, nwav))
    for ypix in range(0, ny):
        for xpix in range(0, nx):
            syn_map[ypix, xpix, :] = outputSir[ypix,xpix][2][1][Nstoke]
    return syn_map


# ====================================================================
def create_profilemap(inversion, inversion_file):
    """
    It creates a file with the synthetic profiles from the inversion
    """
    
    # It should have the shape: [ny, nx, nwav, nstokes]
    ny = inversion.shape[0]
    nx = inversion.shape[1]
    nwav = inversion[0,0][2][0].shape[0]
    
    # Create the file:
    profilemap = np.zeros((ny, nx, nwav, 4))
    for stoke in tqdm(range(4)):
        profilemap[:, :, :, stoke] = readSIRProfileMap(inversion, stoke)
    
    # Save the file:
    if inversion_file[-4:] == '.npy':
        np.save(inversion_file[:-4]+'_profiles.npy', profilemap.astype(np.float32))
    else:
        np.save(inversion_file+'_profiles.npy', profilemap.astype(np.float32))



#=============================================================================
def addFullProfile(sirfile):
    """
    TODO: Having the option of generating synthetic profiles with other properties
    """
    import os
    os.system('echo sirFull.trol | '+sirfile+' > pylogFull.txt')


#=============================================================================
def total_cores():
    """
    Get the total number of cores in the machine
    """
    try:
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        print('[INFO] Available cores = '+str(num_cores))
    except:
        pass

#=============================================================================
def pprint(ini='', end='', comm=MPI.COMM_WORLD):
    """Print for MPI parallel programs: Only rank 0 prints *str*."""
    if comm.rank == 0:
        print(str(ini)+end)


#=============================================================================
def getTerminalSize():
    """
    Get the terminal size in characters
    """
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
        except:
            pass
    return int(cr[1]), int(cr[0])


#=============================================================================
def plotper(main_file='data.per',
            synth_file=None,
            color1='k',
            color2='m',
            y_range_max=np.array([1.1, 3., 3., 3., 3.])):
    """
    Plot the observed and synthetic profiles.
    """
    if synth_file is None:
        # Read the last model written by SIR:
        file_list = glob.glob("*.per")
        file_list.sort(key=os.path.getmtime, reverse=True)
        synth_file = os.path.basename(file_list[0])

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
            plt.xlabel(r'$\Delta\lambda$ [$\AA$]', fontsize=12)
            plt.ylabel(['I/Ic', 'Q/Ic [%]', 'U/Ic [%]', 'V/Ic [%]'][sParam], fontsize=12)
            plt.xlim(x0A[pos_new[line_idx]], x0A[pos_new[line_idx + 1] - 1])
            plt.grid(alpha=0.2, linestyle='-')
            plt.locator_params(axis='both', nbins=4)
            plt.minorticks_on()

            # Plot synthetic profiles
            synth_data = stokes[sParam][x_range]
            if sParam > 0:  # Apply scaling factor of 100 to Q, U, and V parameters
                synth_data = synth_data * 100
            plt.plot(xA[x_range], synth_data, color2, lw=1.0)
            
            if sParam == 0:
                plt.ylim(y_range_min[sParam], np.max(data) * 1.05)
            else:
                plt.ylim(y_range_min[sParam], y_range_max[sParam])


    # Save the figure
    plt.tight_layout()
    output_file = 'P' + synth_file[:-4] + '.pdf'
    plt.savefig(output_file, bbox_inches='tight')
    print(output_file + ':: SAVED')


#=============================================================================
def plotmfit(main_file='hsraB.mod',
             synth_file=None,
             error_model=True,
             index_to_plot=[0, 3, 4, 5],
             labels=['$T$ [K]', r'$P_e$' + ' [dyn cm^-3]', r'$v_{mic}$' + ' [cm/s]', '$B$ [G]', r'$v_{LOS}$' + ' [m/s]', r'$\gamma$ [deg]'],
             color1='k',
             color2='m',
             margin=[0.2, 0.3, 0.3, 0.3]):
    """
    Plot the initial and final model parameters.
    """
    if synth_file is None:
        # Read the last model written by SIR:
        file_list = glob.glob("*.mod")
        file_list.sort(key=os.path.getmtime, reverse=True)
        synth_file = os.path.basename(file_list[0])

    num_plots = len(index_to_plot)

    # Load data
    tau, data = lmodel8(main_file, verbose=False)
    tau2, data2 = lmodel8(synth_file, verbose=False)
    if error_model:
        _, error_data = lmodel8(synth_file[:-4] + '.err')

    plt.figure(figsize=(4 * num_plots, 5))

    for i, index in enumerate(index_to_plot):
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
        plt.xlabel(r'$\log(\tau_{500})$', fontsize=12)
        plt.ylabel(labels[index], fontsize=12)
        plt.locator_params(axis='both', nbins=4)
        plt.grid(alpha=0.2, linestyle='-')
        plt.minorticks_on()

    # Save the figure
    plt.tight_layout()
    output_file = 'M' + synth_file[:-4] + '.pdf'
    plt.savefig(output_file, bbox_inches='tight')
    print(output_file + ':: SAVED')


#=============================================================================
def close_index(number, array):
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

    indice = close_index(0.,x0)
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


#=============================================================================
def checkParamsfile(Paramsfile):
    """
    Check that the a file exists inside the inversions folder:
    """
    import copy
    import shutil

    OParamsfile = copy.copy(Paramsfile)
    # Extract the name of the file:
    Paramsfile = Paramsfile.split('/')[-1]
    
    # Check if the file exists in the folder invDefault, if not copy it:
    if not os.path.exists('invDefault/'+Paramsfile):
        # Try to copy the file:
        try:
            shutil.copy(Paramsfile, 'invDefault/'+Paramsfile)
            print('[INFO] Copied to invDefault/'+Paramsfile)
        except FileNotFoundError:
            print("[INFO] File not found! Exitting...")
            sys.exit(1)
                
    else:
        print('[INFO] Already exists in invDefault/'+Paramsfile)
    return Paramsfile


#=============================================================================
def getLambdaRef(dictLines,Linesfile):
    """
    Open the Linesfile goes to the location of the line number
    and extract the reference wavelength:
    """
    with open('invDefault/'+Linesfile) as f:
        for i, line in enumerate(f):
            atomindex = line.split('=')[0].strip()
            if atomindex == dictLines['atom'].split(',')[0]:
                lambdaRef = float(line.split()[2])
                break
        print('[INFO] lambdaRef = %.4f' % lambdaRef)
    return lambdaRef


#=============================================================================
def loadanyfile(inpufile):
    """
    Detect and load any file with the correct format
    """
    # Check if image is fits or npy:
    isfits = inpufile.split('.')[-1] == 'fits'
    isfits = isfits or inpufile.split('.')[-1] == 'FITS'
    
    # Load image:
    if isfits:
        from astropy.io import fits
        inputdata = fits.open(inpufile)[0].data
    else:
        inputdata = np.load(inpufile)
    return inputdata

#=============================================================================
def varname(p):
    """
    Returns the name of a variable as a string
    """
    import inspect, re
    for line in inspect.getframeinfo(inspect.currentframe().f_back)[3]:
        m = re.search(r'\bvarname\s*\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*\)', line)
        if m:
            return m.group(1)

#=============================================================================
def check_nodes(nodes_allowed, node_variable_list, node_names):
    """
    Check that the nodes are allowed, as SIR will change them internally later
    """
    for j, node_variable_i in enumerate(node_variable_list):
    
        # Check if any of the elements of node_variable_i is different from the allowed:
        if any([int(i) not in nodes_allowed for i in node_variable_i.split(',')]):
            for i in range(len(node_variable_i.split(','))):
                if int(node_variable_i.split(',')[i]) not in nodes_allowed:
                    node_variable_i = node_variable_i.split(',')
                    closest_lower_node = np.max(np.array(nodes_allowed)[np.array(nodes_allowed)<int(node_variable_i[i])])
                    node_variable_i[i] = str(closest_lower_node)
                    node_variable_i = ','.join(node_variable_i)
            print('[INFO] Nodes in '+node_names[j]+' not allowed. Changing to '+str(node_variable_i))