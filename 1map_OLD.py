# Author: cdiazbas@iac.es

import matplotlib.pyplot as plt
import pyLib.imtools as imtools
import numpy as np


# # ========================= CREANDO DICCIONARIO

# cdict1={'red':	((0.0, 0.0, 0.0),
# 				(0.5, 0.0, 0.1),
# 				(1.0, 1.0, 1.0)),
# 		'green':((0.0, 0.0, 0.0),
# 				(1.0, 0.0, 0.0)),
# 		'blue': ((0.0, 0.0, 0.0),
# 				(0.5, 0.0, 0.1),
# 				(1.0, 1.0, 1.0))
# 				}
				
import matplotlib.colors as mcolors
# #blue_red1 = mcolors.LinearSegmentedColormap('BlueRed1', cdict1)


def make_colormap(seq):
		"""Return a LinearSegmentedColormap
		seq: a sequence of floats and RGB-tuples. The floats should be increasing
		and in the interval (0,1).
		"""
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
# #phimap = make_colormap(
# 		#[c('blue'), c('white'), 0.33, c('white'), 0.66, c('white'),c('blue')])
# #phimap = make_colormap(
# 		#[c('grey'), c('white'),0.5, c('white'), c('grey')])
phimap = make_colormap([c('white'), c('tomato'), 0.33, c('tomato'), c('deepskyblue'), 0.66, c('deepskyblue'),c('white')])
# phimap = make_colormap([c('white'), c('tomato'), 0.33, c('tomato'), c('steelblue'), 0.66, c('steelblue'),c('white')])

# phimap = make_colormap([c('red'), 0.33, c('red'), c('blue'), 0.66, c('blue')])

# phimap = make_colormap([c('tomato'), c('gold'), 0.25, c('gold'), c('deepskyblue'), 0.50, c('deepskyblue'),c('hotpink'), 0.75, c('hotpink'),c('tomato')])
# phimap = make_colormap([c('tomato'), 0.33, c('gold'), 0.66, c('deepskyblue')])


import numpy as np
from matplotlib.colors import LinearSegmentedColormap as lsc


def cmap_map(function, cmap, name='colormap_mod', N=None, gamma=None):
    """
    Modify a colormap using `function` which must operate on 3-element
    arrays of [r, g, b] values.

    You may specify the number of colors, `N`, and the opacity, `gamma`,
    value of the returned colormap. These values default to the ones in
    the input `cmap`.

    You may also specify a `name` for the colormap, so that it can be
    loaded using plt.get_cmap(name).
    """
    if N is None:
        N = cmap.N
    if gamma is None:
        gamma = cmap._gamma
    cdict = cmap._segmentdata
    # Cast the steps into lists:
    step_dict = {key: map(lambda x: x[0], cdict[key]) for key in cdict}
    # Now get the unique steps (first column of the arrays):
    step_list = np.unique(sum(step_dict.values(), []))
    # 'y0', 'y1' are as defined in LinearSegmentedColormap docstring:
    y0 = cmap(step_list)[:, :3]
    y1 = y0.copy()[:, :3]
    # Go back to catch the discontinuities, and place them into y0, y1
    for iclr, key in enumerate(['red', 'green', 'blue']):
        for istp, step in enumerate(step_list):
            try:
                ind = step_dict[key].index(step)
            except ValueError:
                # This step is not in this color
                continue
            y0[istp, iclr] = cdict[key][ind][1]
            y1[istp, iclr] = cdict[key][ind][2]
    # Map the colors to their new values:
    y0 = np.array(map(function, y0))
    y1 = np.array(map(function, y1))
    # Build the new colormap (overwriting step_dict):
    for iclr, clr in enumerate(['red', 'green', 'blue']):
        step_dict[clr] = np.vstack((step_list, y0[:, iclr], y1[:, iclr])).T
    return lsc(name, step_dict, N=N, gamma=gamma)

def cmap_map(function,cmap):
    """ Applies function (which should operate on vectors of shape 3:
    [r, g, b], on colormap cmap. This routine will break any discontinuous     points in a colormap.
    """
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ('red','green','blue'):         step_dict[key] = map(lambda x: x[0], cdict[key])
    step_list = sum(step_dict.values(), [])
    step_list = array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : array(cmap(step)[0:3])
    old_LUT = array(map( reduced_cmap, step_list))
    new_LUT = array(map( function, old_LUT))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i,key in enumerate(('red','green','blue')):
        this_cdict = {}
        for j,step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j,i]
            elif new_LUT[j,i]!=old_LUT[j,i]:
                this_cdict[step] = new_LUT[j,i]
        colorvector=  map(lambda x: x + (x[1], ), this_cdict.items())
        colorvector.sort()
        cdict[key] = colorvector

    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)



def dimMap(resultadoSir):
	height = resultadoSir.shape[0]*(resultadoSir1[0][-1][0][0]+1)
	width = (resultadoSir1[0][-1][0][1]+1)
	return [height,width]


def readmapa(resultadoSir, mapa, magnitud):
	cont = 0
	for fila in range(0, height):
		for columna in range(0, width):
			punto = cont % resultadoSir.shape[1]
			veces = int(cont/resultadoSir.shape[1])
			if magnitud == 8 or magnitud == 9 or magnitud == 10 or magnitud == 11:
				mapa[columna,fila] = resultadoSir[veces][punto][1][0][magnitud]
			else:
				mapa[columna,fila] = resultadoSir[veces][punto][1][0][magnitud][index]
			cont += 1
	return mapa


def corrphi(mapa):
	mapa[mapa<0] = (mapa[mapa<0]+360) % 360; mapa[mapa>180] = (mapa[mapa>180]-180)


# ==============================================================================================
global index
global magnitud

import matplotlib
#hsv
# from numpy import array
# phimap = cmap_map(lambda x: x/2+0.5, matplotlib.cm.jet)

# ========================= INPUT
invSir1 = 'finalSir.npy'
logTau = 0.0
magnitud = 7
cmapArray = ['','gray','gray','bone','bone','seismic','Spectral_r',phimap,'bone','gray','gray','cubehelix']
magTitle = ['TAU','$T$ $[kK]$','p','$v$ $[km/s]$','$B$ $[kG]$','$v$ $[km/s]$','$\gamma$ $[d]$','$\phi$ $[d]$','vmacro','fillingf','difusa','$\chi^2$']
magFile = ['TAU','TEMP','PRESION','VMICRO','CAMPO','VLOS','GAMMA','PHI','VMACRO','FILLING','DIFUSA','CHI2']

# ========================= MAP
resultadoSir1 = np.load(invSir1)
height, width = dimMap(resultadoSir1)
print('height:',height,'width:',width)

mapa = np.zeros((height, width))
index = np.where(resultadoSir1[0][0][1][0][0] == logTau)[0][0]
print('logTau: '+str(logTau)+' -> index: '+str(index))
readmapa(resultadoSir1, mapa.T ,magnitud)

# Limites en la escala de color
if magnitud == 7: corrphi(mapa)
print('3sigma_map: {0:2.2f}'.format(3*np.std(mapa)))
print('Mean_map: {0:2.2f}'.format(np.mean(mapa)))
print('Min_map: {0:2.2f}'.format(np.min(mapa)))
print('Max_map: {0:2.2f}'.format(np.max(mapa)))

vmini = np.mean(mapa)-3*np.std(mapa)
if np.min(mapa) >= 0.0 and magnitud != 1: vmini = 0.
vmaxi = np.mean(mapa)+3*np.std(mapa)
if magnitud == 1: vmini = np.min(mapa); vmaxi = np.max(mapa)
if magnitud == 6: vmaxi = 180.
if magnitud == 7: vmaxi = 180.;vmini = 0.
if magnitud == 5: vmaxi = np.mean(mapa)+3*np.std(mapa); vmini = -vmaxi


from matplotlib.colors import LogNorm
plt.imshow(mapa,cmap=cmapArray[magnitud],origin='lower',interpolation='None',vmin=vmini,vmax=vmaxi)#norm=LogNorm()
plt.title('Map 17jun14.006 (1)')
plt.xlabel('Slit Axis [pix]')
plt.ylabel('Time Axis [pix]')
cb = plt.colorbar(shrink=.46)#, ticks=[0.6, 0.8, 1., 1.2])
#cb = plt.colorbar(shrink=.46, ticks=[0.3, 0.6, 0.9, 1.2, 1.5])
# cb.set_label(r'Intensity HeI ({0:4.1f}) /$I_{{qs}}$({1:4.1f})'.format(xLambda[341],xLambda[posicontinuo]), labelpad=5., y=0.5, fontsize=12.)
cb.set_label(r""+magTitle[magnitud]+", $log(\\tau)$={0}".format(logTau), labelpad=8., y=0.5, fontsize=12.)

# plt.show()
plt.savefig(magFile[magnitud]+'_log{0:02d}.pdf'.format(int(logTau)), bbox_inches='tight')
print(magFile[magnitud]+'_log{0:02d}.pdf SAVE'.format(int(logTau)))



