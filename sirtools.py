# cdiazbas@iac.es
# Code: SIR-files tools

"""
This file include:

1.-  lambda_mA, stokesIQUV, [nL,posi,nN] = lperfil(filename)

2.-  wperfil(filename, numberLine, lambda_mA, stokes)

3.-  [tau, todoPlot] = lmodel8(filename, verbose=True)

4.-  wmodel8(modelo, filename, verbose=False)

5.-  mapa = readSIRMap(resultadoSir, magnitud)

6.-  [height, width, nlambda] = shapeSIRMap(resultadoSir)

7.-  mapa = readSIRProfileMap(resultadoSir, Nstoke)



MODEL ATMOSPHERE FILES TO BE USED WITH SIR
Each model file contains the macroturbulent velocity (km/s), the
filling factor (only to be used with two-component models, ranging
from 0 to 1), and the stray light factor (in percent) in the first line.
Then, eight columns follow:
 
Column 1: log tau_5 (logarithm of the continuum optical depth at 5000 A)
Column 2: Temperature (K)
Column 3: Electron pressures (dyn/cm^2)
Column 4: Microturbulent velocity (cm/s)
Column 5: Magnetic field strength (G)
Column 6: Line-of-sight velocity (cm/s)
Column 7: Inclination angle of the magnetic field vector in deg 
          from 0 (pointing to the observer) to 180 (pointing away from the
          observer)
Column 8: Azimuthal angle of the magnetic field vector in deg. 
Column 9: Geometrical scale (km)
Column 10: Gas presure (dyn/cm^2)
Column 11: Gas density (gr/cm^3)

"""


# ====================================================================
def circular_mean(alpha):
    import numpy as np
    return np.arctan2(np.sum(np.sin(alpha*np.pi/180.)), np.sum(np.cos(alpha*np.pi/180.)))*180./np.pi

# ====================================================================
def circular_map_smooth(mapa, cuanto=1):
    import numpy as np
    new_mapa = np.copy(mapa)
    for i in range(new_mapa.shape[0]):
        for j in range(new_mapa.shape[1]):
            if i > cuanto and j>cuanto:
                new_mapa[i,j] = (circular_mean(2*mapa[i-cuanto:i+cuanto,j-cuanto:j+cuanto])/2. +360. ) %180.
            else:
                new_mapa[i,j] = (circular_mean(2*mapa[i:i+cuanto,j:j+cuanto])/2. +360. ) %180.
    return new_mapa

# ====================================================================
def vectorMapa(phiMap, sep, color, suma, difu,xscale=1.,yscale=1.):
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.ndimage as sn
    plt.autoscale(False)
    newPhi = sn.filters.gaussian_filter(phiMap, difu)
    for j in range(0, phiMap.shape[0], sep):
        for i in range(0, phiMap.shape[1], sep):
            plt.plot(np.array([i, i+1.*np.cos((newPhi[j, i]+suma)/180.*np.pi)])*xscale, np.array([j, j+1.*np.sin((newPhi[j, i]+suma)/180.*np.pi)])*yscale, color=color, lw=0.5)

# ====================================================================
def corrphi(mapa):
    mapa[:] = (mapa[:]+360.) %180.
    pass


# ====================================================================
def lperfil(filename, verbose=False):
    """Read SIR Stokes profile

    Args:
        filename (string)

    Returns:
        lambda_mA, stokesIQUV, [nL,posi,nN] = lperfil(filename)
    """
    from numpy import array
    fo = open(filename, 'r')

    Nn=[]; NumeroLineas = 1; PosiNn0T= []
    x0=[]
    StokeI0=[];  StokeQ0=[];  StokeU0=[];  StokeV0=[]
    x=[]
    StokeI=[];  StokeQ=[];  StokeU=[];  StokeV=[]
    for ii in fo:
        linea_split = ii.split()
        Nn.append(linea_split[0]) # do not convert to float as it can be blends
        x0.append(float(linea_split[1]))
        StokeI0.append(float(linea_split[2]))
        StokeQ0.append(float(linea_split[3]))
        StokeU0.append(float(linea_split[4]))
        StokeV0.append(float(linea_split[5]))

    # Conversion a array:
    x0=array(x0)
    StokeI0 = array(StokeI0)
    StokeQ0 = array(StokeQ0)
    StokeU0 = array(StokeU0)
    StokeV0 = array(StokeV0)
    lenNn = len(Nn)

    # Posiciones de las distintas lineas del filename
    PosiNn0T.append(0)
    try:
        NnInit = Nn[0]
        PosiNn0 = 0
        for NextI in range(0,lenNn-1):
            if (Nn[NextI] != NnInit):
                PosiNn0T.append(NextI)
            NnInit = Nn[NextI]

    except:
        print('Only1Line') #REVISAR!!
    PosiNn0T.append(lenNn-1)
    NumeroLineas = len(PosiNn0T)-1

    # Almaceno las lineas dentro del array
    for Index in range(NumeroLineas):
        StokeI.append(StokeI0[PosiNn0T[Index]:PosiNn0T[Index+1]-1])
        StokeQ.append(StokeQ0[PosiNn0T[Index]:PosiNn0T[Index+1]-1])
        StokeU.append(StokeU0[PosiNn0T[Index]:PosiNn0T[Index+1]-1])
        StokeV.append(StokeV0[PosiNn0T[Index]:PosiNn0T[Index+1]-1])
        x.append(x0[PosiNn0T[Index]:PosiNn0T[Index+1]-1])

    PosiNn0T = PosiNn0T[:-1]

    # Si hay una linea sola
    if len(x) == 1:
        x = x0; StokeI = StokeI0; StokeQ = StokeQ0;
        StokeU = StokeU0; StokeV = StokeV0
    if verbose:
        print('NumeroLineas:'+str(NumeroLineas))
        print('Info: lambda in mA')
        print('lambda_mA, stokesIQUV, [nL,posi,nN]')
    fo.close()
        
    return [x, [StokeI, StokeQ, StokeU, StokeV], [NumeroLineas, PosiNn0T, Nn]]


# ====================================================================
def wperfil(filename, numberLine, lambda_mA, stokes):
    """Write SIR Stokes profile in a file
    """
    si = stokes[0]
    sq = stokes[1]
    su = stokes[2]
    sv = stokes[3]
    fo = open(filename, 'w')
    for i in range(len(lambda_mA)):
        fo.write('   {0}  {1:3.4f}  {2:2.6E}  {3:2.6e}  {4:2.6e}  {5:2.6e}\n'\
            .format(numberLine,lambda_mA[i],si[i],sq[i],su[i],sv[i]))
    fo.close()
    return


# ====================================================================
def lmodel8(modelo, verbose=False):

    from numpy import array
    fo = open(modelo, 'r')
    tau = []
    temp = []
    Pres = []
    vmic = []
    BMag = []
    vlos = []
    gamma = []
    phi = []
    c = 0
    for ii in fo:
        linea_split = ii.split()
        if c == 0:
            vmac = float(linea_split[0])
            fill = float(linea_split[1])
            stray = float(linea_split[2])
        if c != 0:
            tau.append(float(linea_split[0]))
            temp.append(float(linea_split[1]))
            Pres.append(float(linea_split[2]))
            vmic.append(float(linea_split[3]))
            BMag.append(float(linea_split[4]))
            vlos.append(float(linea_split[5]))
            gamma.append(float(linea_split[6]))
            phi.append(float(linea_split[7]))
        c += 1
    # Conversion a array:
    tau = array(tau)
    temp = array(temp)
    Pres = array(Pres)
    vmic = array(vmic)
    BMag = array(BMag)
    vlos = array(vlos)
    gamma = array(gamma)
    phi = array(phi)
    lenTau = len(tau)

    todoPlot = [temp,Pres,vmic,BMag,vlos,gamma,phi,vmac,fill,stray]
    fo.close()
    if verbose:
        print('temp[kK], Pres[dyn cm^-3], vmic[km/s], BMag[kG], vlos[km/s], gamma[deg], phi[deg], vmac[km/s], fill, stray')
        print('Out: {tau, magnitudes}')
    return [tau, todoPlot]




# ====================================================================
def wmodel8(modelo, filename, verbose=False):

    [tau, todoPlot] = modelo
    temp = todoPlot[0]
    Pres = todoPlot[1]
    vmic = todoPlot[2]
    Bmag = todoPlot[3]
    vlos = todoPlot[4]
    gamma = todoPlot[5]
    phi = todoPlot[6]
    vmac = todoPlot[7]
    fill = todoPlot[8]
    stray = todoPlot[9]

    fo = open(filename, 'w')
    for i in range(-1,len(temp)):
        if i == -1:
            fo.write('   {0:3.4f}  {1:3.4f}  {2:3.4f}\n'.format(vmac, fill, stray))
        if i != -1:
            fo.write('   {0:2.4f}  {1:2.6e}  {2:2.6e}  {3:2.6e}  {4:2.6e}  {5:2.6e}  {6:2.6e}  {7:2.6e}\n'
                     .format(tau[i], temp[i], Pres[i], vmic[i], Bmag[i], vlos[i], gamma[i], phi[i]))
    fo.close()
    return


# ====================================================================
def lmodel12(modelo, verbose=False):
    '''
    MODEL ATMOSPHERE FILES TO BE USED WITH SIR

    Each model file contains the macroturbulent velocity (km/s), the
    filling factor (only to be used with two-component models, ranging
    from 0 to 1), and the stray light factor (in percent) in the first line.
    Then, eight columns follow:
     
    Column 1: log tau_5 (logarithm of the continuum optical depth at 5000 A)
    Column 2: Temperature (K)
    Column 3: Electron pressures (dyn/cm^2)
    Column 4: Microturbulent velocity (cm/s)
    Column 5: Magnetic field strength (G)
    Column 6: Line-of-sight velocity (cm/s)
    Column 7: Inclination angle of the magnetic field vector in deg 
              from 0 (pointing to the observer) to 180 (pointing away from the
              observer)
    Column 8: Azimuthal angle of the magnetic field vector in deg. 
    Column 9: Geometrical scale (km)
    Column 10: Gas presure (dyn/cm^2)
    Column 11: Gas density (gr/cm^3)
    '''

    from numpy import array
    fo = open(modelo, 'r')
    tau = []
    temp = []
    Pres = []
    vmic = []
    BMag = []
    vlos = []
    gamma = []
    phi = []
    zz = []
    pgas = []
    rho = []
    c = 0
    for ii in fo:
        linea_split = ii.split()
        if c == 0:
            vmac = float(linea_split[0])
            fill = float(linea_split[1])
            stray = float(linea_split[2])
        if c != 0:
            tau.append(float(linea_split[0]))
            temp.append(float(linea_split[1]))
            Pres.append(float(linea_split[2]))
            vmic.append(float(linea_split[3]))
            BMag.append(float(linea_split[4]))
            vlos.append(float(linea_split[5]))
            gamma.append(float(linea_split[6]))
            phi.append(float(linea_split[7]))
            zz.append(float(linea_split[8]))
            pgas.append(float(linea_split[9]))
            rho.append(float(linea_split[10]))
        c += 1
    # Conversion a array:
    tau = array(tau)
    temp = array(temp)
    Pres = array(Pres)
    vmic = array(vmic)
    BMag = array(BMag)
    vlos = array(vlos)
    gamma = array(gamma)
    phi = array(phi)
    zz = array(zz)
    pgas = array(pgas)
    rho = array(rho)
    # lenTau = len(tau)

    todoPlot = [temp,Pres,vmic,BMag,vlos,gamma,phi,vmac,fill,stray,zz,pgas,rho]
    fo.close()
    if verbose:
        print('temp[K], Pres[dyn cm^-2], vmic[cm/s], BMag[G], vlos[cm/s], gamma[deg], phi[deg], vmac[km/s], fill, stray,zz,pgas,rho')
        print('Out: {tau, magnitudes}')
    return [tau, todoPlot]



# ====================================================================
def wmodel12(modelo, filename, verbose=False):
    """
    Write a model file for SIR with the format of lmodel12
    
    Each model file contains the macroturbulent velocity (km/s), the
    filling factor (only to be used with two-component models, ranging
    from 0 to 1), and the stray light factor (in percent) in the first line.
    Then, eight columns follow:
     
    Column 1: log tau_5 (logarithm of the continuum optical depth at 5000 A)
    Column 2: Temperature (K)
    Column 3: Electron pressures (dyn/cm^2)
    Column 4: Microturbulent velocity (cm/s)
    Column 5: Magnetic field strength (G)
    Column 6: Line-of-sight velocity (cm/s)
    Column 7: Inclination angle of the magnetic field vector in deg 
              from 0 (pointing to the observer) to 180 (pointing away from the
              observer)
    Column 8: Azimuthal angle of the magnetic field vector in deg. 
    Column 9: Geometrical scale (km)
    Column 10: Gas presure (dyn/cm^2)
    Column 11: Gas density (gr/cm^3)
    '''
    """

    [tau, todoPlot] = modelo
    temp = todoPlot[0]
    Pres = todoPlot[1]
    vmic = todoPlot[2]
    Bmag = todoPlot[3]
    vlos = todoPlot[4]
    gamma = todoPlot[5]
    phi = todoPlot[6]
    vmac = todoPlot[7]
    fill = todoPlot[8]
    stray = todoPlot[9]
    zz = todoPlot[10]
    pgas = todoPlot[11]
    rho = todoPlot[12] 

    fo = open(filename, 'w')
    for i in range(-1,len(temp)):
        if i == -1:
            fo.write('   {0:3.4f}  {1:3.4f}  {2:3.4f}\n'.format(vmac, fill, stray))
        if i != -1:
            fo.write('   {0:2.4f}  {1:2.1f}  {2:2.3e}  {3:2.3e}  {4:2.3e}  {5:2.3e}  {6:2.3e}  {7:2.3e}  {8:2.3e}  {9:2.3e}  {10:2.4e}\n'
                     .format(tau[i], temp[i], Pres[i], vmic[i], Bmag[i], vlos[i], gamma[i], phi[i], zz[i], pgas[i], rho[i]))
    fo.close()
    return



