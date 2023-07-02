import matplotlib.pyplot as plt
import numpy as np

"""
This module contains functions to plot the inversion results.
"""


# ========================= PHI colormap
import matplotlib.colors as mcolors
def make_colormap(seq):
	seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
	cdict = {'red': [], 'green': [], 'blue': []}
	for i, item in enumerate(seq):
		if isinstance(item, float):
			r1, g1, b1 = seq[i - 1]
			r2, g2, b2 = seq[i + 1]
			cdict['red'].append([item, r1, r2])
			cdict['green'].append([item, g1, g2])
			cdict['blue'].append([item, b1, b2])
	return mcolors.LinearSegmentedColormap('CustomMap', cdict)
c = mcolors.ColorConverter().to_rgb
phimap = make_colormap([c('white'), c('tomato'), 0.33, c('tomato'), c('deepskyblue'), 0.66, c('deepskyblue'),c('white')])


# ========================= PHI CORRECTION
def corrphi(azimuthmap):
    # Fix the azimuth values so that they are in the range [0,180]
    azimuthmap[azimuthmap<0] = (azimuthmap[azimuthmap<0]+360) % 360
    azimuthmap[azimuthmap>180] = (azimuthmap[azimuthmap>180]-180)
    return


# ========================= PLOT 1 MAP
def plot1map(indexlogTau, parameter, inversion_model = 'finalSIR_model.npy', extra=''):

    cmapArray = ['gray','gray','gray','bone','bone','seismic','Spectral_r',phimap,'bone','gray','gray','cubehelix']
    magTitle = [r'${\rm log(\tau)=}$',r'${\rm T\ [K]}$','p',r'${\rm v\ [cm/s]}$',r'${\rm B\ [G]}$',r'${\rm v\ [cm/s]}$',r'${\rm \gamma\ [d]}$',r'${\rm \phi\ [d]}$','vmacro','filling factor','stray-light (alpha)',r'${\rm \chi^2}$']
    magFile = ['_LOGTAU','_TEMP','_PGAS','_VMICRO','_B','_VLOS','_INCLINATION','_AZIMUTH','_VMACRO','_FILLING','_ALPHA','_CHI2']

    # ========================= MAP
    # Load the inversion model
    inversion_model = np.load(inversion_model, allow_pickle=True)
    param2plot = inversion_model[:,:,indexlogTau,parameter]
    
    # Print the logtau value at that index:
    logTau = inversion_model[0,0,indexlogTau,0]
    print('logTau: {0:2.2f}'.format(logTau))
    
    # Fix the azimuth values so that they are in the range [0,180]
    if parameter == 7: corrphi(param2plot)

    # Vmin and Vmax for the plot:
    vmini = np.mean(param2plot)-3*np.std(param2plot)
    if np.min(param2plot) >= 0.0 and parameter != 1: vmini = 0.
    vmaxi = np.mean(param2plot)+3*np.std(param2plot)
    if parameter == 1 or parameter == 4: vmini = np.min(param2plot); vmaxi = np.max(param2plot)
    if parameter == 6: vmaxi = 180.
    if parameter == 7: vmaxi = 180.;vmini = 0.
    if parameter == 11: vmaxi = np.mean(param2plot)+6*np.std(param2plot); vmini = 0.
    if parameter == 5: vmini = np.mean(param2plot)-4*np.std(param2plot); vmaxi = -vmini
    if parameter == 11: vmini = 0.0; vmaxi = 30.0
    if parameter == 4: vmini = 0.0; vmaxi = np.percentile(param2plot,99.0)

    # Plot the map associated to the parameter:
    ratio = param2plot.shape[1]/param2plot.shape[0]
    plt.figure(figsize=(ratio*4,4))
    plt.imshow(param2plot,cmap=cmapArray[parameter],origin='lower',interpolation='None',vmin=vmini,vmax=vmaxi)
    plt.xlabel('X Axis [pix]')
    plt.ylabel('Y Axis [pix]')
    cb = plt.colorbar()
    loglabel = r'${\rm log(\tau)=}$'
    cb.set_label(r""+magTitle[parameter]+r", "+loglabel+"{0}".format(logTau), labelpad=8., y=0.5, fontsize=12.)

    plt.savefig(magFile[parameter]+'_log{0:02.2f}{1}.pdf'.format(logTau,extra), bbox_inches='tight')
    
    print(magFile[parameter]+'_log{0:02.2f}{1}.pdf SAVE'.format(logTau,extra))
    print('-----------------------'+str(parameter))
    plt.clf()




inversion_model = 'finalSIR_cycle2_model.npy'
extra = '_cycle2'
index = 14 # 14 corresponds to logtau = 0.0, 24 to logtau=-1.0

rangeparams = [1,2,4,5,6,7,11]
for parameter in rangeparams:
    plot1map(index,parameter, inversion_model = inversion_model, extra=extra)
