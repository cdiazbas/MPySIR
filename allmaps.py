# Author: cdiazbas@iac.es

import matplotlib.pyplot as plt
import pyLib.imtools as imtools
import numpy as np


# ========================= CREANDO PHIMAP

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


def dimMap(resultadoSir):
    height = resultadoSir.shape[0]*(resultadoSir[0][-1][0][0]+1)
    width = (resultadoSir[0][-1][0][1]+1)
    return [height, width]


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


def do1map(logTau, magnitud):
    # ==============================================================================================
    # global index
    # global magnitud

    # ========================= INPUT
    invSir1 = 'finalSir.npy'
    # logTau = 0.0
    # magnitud = 2
    # hsv
    cmapArray = ['gray','gray','gray','bone','bone','seismic','Spectral_r',phimap,'bone','gray','gray','cubehelix']
    magTitle = ['TAU',r'${\rm T\ [kK]}$','p',r'${\rm v\ [km/s]}$',r'${\rm B\ [kG]}$',r'${\rm v\ [km/s]}$',r'${\rm \gamma\ [d]}$',r'${\rm \phi\ [d]}$','vmacro','fillingf','difusa',r'${\rm \chi^2}$']
    magFile = ['TAU','TEMP','PRESION','VMICRO','CAMPO','VLOS','GAMMA','PHI','VMACRO','FILLING','DIFUSA','CHI2']

    # ========================= MAP
    resultadoSir1 = np.load(invSir1)
    # height, width = dimMap(resultadoSir1)
    # print('height:',height,'width:',width)

    # mapa = np.zeros((height, width))
    index = np.where(resultadoSir1[0][0][1][0][0] == logTau)[0][0]
    print('logTau: '+str(logTau)+' -> index: '+str(index))
    # readmapa(resultadoSir1, mapa.T ,magnitud)

    from pySir import sirtools as st
    mapa = st.readSIRMap(resultadoSir1, magnitud, index)

    # Limites en la escala de color
    if magnitud == 7: corrphi(mapa)
    print('3sigma_map: {0:2.2f}'.format(3*np.std(mapa)))
    print('Mean_map: {0:2.2f}'.format(np.mean(mapa)))
    print('Min_map: {0:2.2f}'.format(np.min(mapa)))
    print('Max_map: {0:2.2f}'.format(np.max(mapa)))

    vmini = np.mean(mapa)-3*np.std(mapa)
    if np.min(mapa) >= 0.0 and magnitud != 1: vmini = 0.
    vmaxi = np.mean(mapa)+3*np.std(mapa)
    if magnitud == 1 or magnitud == 4: vmini = np.min(mapa); vmaxi = np.max(mapa)
    if magnitud == 6: vmaxi = 180.
    if magnitud == 7: vmaxi = 180.;vmini = 0.
    if magnitud == 11: vmaxi = np.mean(mapa)+6*np.std(mapa); vmini = 0.
    if magnitud == 5: vmini = np.mean(mapa)-4*np.std(mapa); vmaxi = -vmini

    from matplotlib.colors import LogNorm
    plt.imshow(mapa,cmap=cmapArray[magnitud],origin='lower',interpolation='None',vmin=vmini,vmax=vmaxi)#norm=LogNorm()
    plt.title('Map 17jun14.006 (3-4)')
    plt.xlabel('Slit Axis [pix]')
    plt.ylabel('Time Axis [pix]')
    cb = plt.colorbar(shrink=.46)#, ticks=[0.6, 0.8, 1., 1.2])
    #cb = plt.colorbar(shrink=.46, ticks=[0.3, 0.6, 0.9, 1.2, 1.5])
    # cb.set_label(r'Intensity HeI ({0:4.1f}) /$I_{{qs}}$({1:4.1f})'.format(xLambda[341],xLambda[posicontinuo]), labelpad=5., y=0.5, fontsize=12.)
    loglabel = r'${\rm log(\tau)=}$'
    cb.set_label(r""+magTitle[magnitud]+r", "+loglabel+"{0}".format(logTau), labelpad=8., y=0.5, fontsize=12.)

    # plt.show()
    plt.savefig(magFile[magnitud]+'_log{0:02d}.pdf'.format(int(logTau)), bbox_inches='tight')
    print(magFile[magnitud]+'_log{0:02d}.pdf SAVE'.format(int(logTau)))
    print('-----------------------'+str(magnitud))
    plt.clf()


for magnitud in range(12):
    do1map(0.0,magnitud)
