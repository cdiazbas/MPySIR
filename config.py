"""
Configuration file for the inversion using MPySIR
"""


# ================================================= INPUT
# We define the input parameters
imagefits = '../../data/sunspot_jmb_sir_synth.fits' # Input image (can be a FITS file or a numpy NPY file).
original_axis = 'ns nw ny nx' # Axis of the input image (ns: spectral points, nw: wavelength points, ny: y-direction, nx: x-direction)
fov = None # Using '20,20' we extract a 20x20 pixels for the inversion
skip = 1 # We can skip pixels in the 2D grid to reduce the number of pixels

wavefile = '../../data/wav.npy' # Wavelength axis in Angstroms absolute reference (the code will use the reference in LINEAS file)
dictLines = {'atom':'200,201'}  # Line Number in LINEAS file
wavrange = None # Range of wavelength points to be used in the inversion. Default: None, Option: range(0,401)
sirmode = 'perPixel' # If 'perPixel', SIR uses a FALC model. If 'continue', SIR starts from a previous model
# sirmode = 'continue' # 'continue'
continuemodel = 'finalSIR_cycle1_model.npy'

apply_constraints = True # Apply smooth extrapolation outside [-3.5,0]
chi2map = True # By default, we include the chi2 map in the output model
test1pixel = False # Test the inversion of one pixel
outputfile = 'finalSIR_cycle2.npy'
verbose = True # Verbose mode


# TODO: move the lambdaRef to the LINEAS file
lambdaRef = 6301.5080 # This is extracted from LINEAS file
# TODO: move the modeloFin to a variable that will change depending on the number of cycles
modeloFin = 'hsraB_3.mod'


# ================================================= INVERSION PARAMETERS
# We define the nodes for the inversion
Nodes_temperature = '5,7,9'            # Temperature
Nodes_magneticfield = '3,5,5'          # Magnetic field strength
Nodes_LOSvelocity = '3,5,5'            # Line-of-sight velocity
Nodes_gamma = '2,3,3'                  # Magnetic field inclination
Nodes_phi = '2,3,3'                    # Magnetic field azimuth
Invert_macroturbulence = '0'           # Macroturbulence (cannot be inverted if vmacro=0)
Initial_vmacro = 0.0 # km/s            # Initial macroturbulence
