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
    sin_az = np.sin(np.deg2rad(azimuthmap)*2.0)
    cos_az = np.cos(np.deg2rad(azimuthmap)*2.0)    
    azimuthmap = np.rad2deg(np.arctan2(sin_az, cos_az))/2.0
    azimuthmap[azimuthmap<0] = azimuthmap[azimuthmap<0]+180
    return azimuthmap

# ==================================================================== COLORBAR
def add_colorbar(im, aspect=20, pad_fraction=0.5, nbins=5, **kwargs):
    """Add a vertical color bar to an image plot."""
    from mpl_toolkits import axes_grid1
    divider = axes_grid1.make_axes_locatable(im.axes)
    width = axes_grid1.axes_size.AxesY(im.axes, aspect=1./aspect)
    pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
    current_ax = plt.gca()
    cax = divider.append_axes("right", size=width, pad=pad)
    plt.sca(current_ax)
    from matplotlib import ticker
    cb = im.axes.figure.colorbar(im, cax=cax, **kwargs)
    if nbins is not None:
        tick_locator = ticker.MaxNLocator(nbins=nbins)
        cb.locator = tick_locator
        cb.update_ticks()
    return cb


# ========================= PLOT 1 MAP
def plot1map(indexlogTau, parameter, inversion_model = 'finalSIR_model.npy', extra=''):

    cmapArray = ['gray','gray','gray','bone','bone','seismic','Spectral_r',phimap,'bone','gray','gray','cubehelix']
    magTitle = [r'${\rm log(\tau)=}$',r'${\rm T\ [K]}$',r'$P_e$ [dyn/cm$^2$]',r'${\rm v_{micro}\ [cm/s]}$',r'${\rm B\ [G]}$',r'${\rm v_{los}\ [km/s]}$',r'${\rm \Theta_B\ [d]}$',r'${\rm \Phi_B\ [d]}$',r'${\rm v_{macro}\ [km/s]}$','filling factor','stray-light (alpha)',r'${\rm \chi^2}$']
    magFile = ['_LOGTAU','_TEMP','_PGAS','_VMICRO','_B','_VLOS','_INCLINATION','_AZIMUTH','_VMACRO','_FILLING','_ALPHA','_CHI2']

    # ========================= MAP
    # Load the inversion model
    inversion_model = np.load(inversion_model, allow_pickle=True)
    # Check that the parameter is in range:
    if parameter in range(inversion_model.shape[3]):
        pass
    else:
        print('Error: parameter out of range')
        return 0
    
    # Search for the index of the logtau value:
    indexlogTau = np.argmin(np.abs(inversion_model[0,0,:,0]-indexlogTau))
    
    # Extract the parameter to plot:
    param2plot = inversion_model[:,:,indexlogTau,parameter]
    
    # Print the logtau value at that index:
    logTau = inversion_model[0,0,indexlogTau,0]
    print('logTau: {0:2.2f}'.format(logTau))
    
    # Fix the NaN values before plotting:
    param2plot = np.nan_to_num(param2plot, nan=np.nanmean(param2plot))

    # Fix the azimuth values so that they are in the range [0,180]
    if parameter == 7: param2plot = corrphi(param2plot)
    
    # LOS velocity in km/s
    if parameter == 5: param2plot = param2plot/1e5

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
    if parameter == 8: vmini = np.percentile(param2plot,1.0)
    
    print('Average: {0:2.2f}'.format(np.mean(param2plot)))

    # Plot the map associated to the parameter:
    ratio = param2plot.shape[1]/param2plot.shape[0]
    plt.figure(figsize=(ratio*4,4))
    im = plt.imshow(param2plot,cmap=cmapArray[parameter],origin='lower',interpolation='nearest',vmin=vmini,vmax=vmaxi)
    plt.xlabel('X axis [pix]')
    plt.ylabel('Y axis [pix]')
    cb = add_colorbar(im)
    loglabel = r'${\rm log(\tau)=}$'
    cb.set_label(r""+magTitle[parameter]+r", "+loglabel+"{0}".format(logTau), labelpad=8., y=0.5, fontsize=12.)

    plt.savefig(magFile[parameter]+'_log{0:02.2f}{1}.pdf'.format(logTau,extra), bbox_inches='tight')
    
    print(magFile[parameter]+'_log{0:02.2f}{1}.pdf SAVE'.format(logTau,extra))
    print('-----------------------'+str(parameter))
    plt.clf()




inversion_model = 'finalSIR_cycle1_model.npy'
extra = '_cycle1'
logtau = 0.0

rangeparams = [1,2,4,5,6,7,8,11]
for parameter in rangeparams:
    plot1map(logtau,parameter, inversion_model = inversion_model, extra=extra)
