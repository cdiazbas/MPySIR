# MPySIR
SIR inversions: MPI implementation with python

## Requirements

- Open the file `requirements.txt` and install the required packages

## Setting up the inversion

- Open the file `config.py.example` file and save it as `config.py`
- Modify the parameters in the `config.py` file to suit your needs
- If you do not have the source code of SIR, you can download it from [here](https://github.com/cdiazbas/SIRcode) and copy it to the `invDefault` folder.
- Run the inversion in MPI mode with `mpirun -np <nproc> python setup.py`

## Input and output formats

- The input data is expected to be FITS or NUMPY files with no specific order of the dimensions. This information is added in the `config.py` file.
- The output data is now a NUMPY file under the name given in the `config.py` file, of size [ny,nx,ntau,nparam] with 12 parameters in the following order: logtau, temperature (K), Pe (cm^-3), v_micro (cm/s), B (G), vlos (cm/s), inclination (deg), azimuth (deg), v_mac (cm/s), filling factor, stray light (fraction), chi2.


## Repository Structure

The repository structure is organized as follows:

- `LICENSE`: Contains the license details for this project.
- `README.md`: This file, containing details about the project and instructions for setting it up.
- `config.py.example`: An example configuration file to guide users in setting up their own inversion.
- `sirutils.py`: Python script containing utilities for this MPI implementation.
- `sirtools.py`: Python script containing utilities for reading and writing SIR files.
- `clean.py`: Python script to clean up the output files if needed.
- `allmaps.py`: Python script to produce some quick plots of the inversion results.
- `invDefault`: Folder containing the default inversion configuration files.
- `nextcycle.py`: Python script to filter the inversion results and prepare for the next cycle.
- `requirements.txt`: Lists all Python dependencies required.
- `setup.py`: Main Python script to run the inversion.


## Commit messages

Use following tags for commit messages:

       [ADD] : Adding new feature
       [DEV] : Improve existing feature
       [DEL] : Removing files, routines
       [FIX] : Fixes existing features
       [OPT] : Optimisation
       [DOC] : Documentation only
