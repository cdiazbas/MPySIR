# MPySIR
a parallel LTE inversion code for spectropolarimetric observations based on the *Stokes Inversion based on Response functions* (SIR) code.

> The original SIR code (Ruiz Cobo & del Toro Iniesta, 1992) is not parallelized, and inverting large datasets can be computationally demanding. To address this, we have developed a parallel version of the SIR code, named MPySIR. This version is written in Python and utilizes the MPI library to parallelize the inversion process. As a result, the inversion tasks can be distributed across multiple processors. It's important to note that the core implementation remains unaltered, with all MPI calls being made from within Python.

> This new implementation integrates in a single configuration file the previous functionalities of the SIR code and new functionalities related to the parallelization. Through this file, users can control various aspects, such as the input/output files, abundances, the mode of synthesis or inversion, the number of nodes for each physical parameter, and more. Additional features include debugging tools, the option to perform inversions only within a specified region of the dataset, the ability to combine different inversion results, the option to use previous inversion results as inputs for subsequent cycles, along with numerous other possibilities. The only feature that we did not carry over is multi-component inversions, mostly because we are interested in very high-resolution observations where we deem that feature unnecessary. 

## Requirements

- Open the file `requirements.txt` and install the required packages

## Setting up the inversion

- Open the file `config.py.example` file and save it as `config.py`
- Modify the parameters in the `config.py` file to suit your needs
- If you do not have the source code of SIR, you can download it from [here](https://github.com/cdiazbas/SIRcode), compile it, and copy the `sir.x` executable inside the `invDefault` folder.
- Run the inversion in MPI mode with `mpirun -np <nproc> python setup.py`

## Input and output formats

- The input data is expected to be FITS or NUMPY files with no specific order of the dimensions. This information is added in the `config.py` file.
- The output data is now a NUMPY file under the name given in the `config.py` file, of size [ny,nx,ntau,nparam] with 12 parameters in the following order: $\log(\tau_{500})$ (logarithm of the continuum optical depth at 500nm), Temperature ($\rm K$), Electron pressure ($\rm dyn/cm^2$), Microturbulent velocity ($\rm cm/s$), Magnetic field strength ($\rm G$), Line-of-sight velocity ($\rm cm/s$), Inclination angle of the magnetic field vector (deg), Azimuthal angle of the magnetic field vector (deg), Macroturbulent velocity ($\rm km/s$), Filling factor, Stray light (fraction), $\chi^2$.

## Telegram Notifications

MPySIR includes a feature to send Telegram notifications when the inversion is finished. To enable this feature, you need to set the following environment variables:

- `TELEGRAM_TOKEN`: Your Telegram bot token.
- `TELEGRAM_CHATID`: Your Telegram chat ID.

The function `notify_telegram` in `sirutils.py` is used to send Telegram notifications. It checks for the environment variables `TELEGRAM_TOKEN` and `TELEGRAM_CHATID` to send the notification. If these variables are not set, the notification will not be sent.

## Repository Structure

The repository structure is organized as follows:

- `invDefault`: Folder containing the default inversion configuration files.
- `config.py.example`: An example configuration file to guide users in setting up their own inversion. Remove the `.example` extension to use it.
- `setup.py`: Main Python script to run the inversion.

The rest of the files are utilities and scripts to help with the inversion process:
- `sirutils.py`: Python script containing utilities for this MPI implementation.
- `sirtools.py`: Python script containing utilities for reading and writing SIR files.
- `clean.py`: Python script to clean up the output files if needed.
- `allmaps.py`: Python script to produce some quick plots of the inversion results.
- `merge.py`: Python script to combine different inversion results according to the quality of the fit.
- `findbest.py`: Python script to improve the inversion results by finding better solutions in other pixels.
- `nextcycle.py`: Python script to filter the inversion results and prepare for the next cycle.
- `requirements.txt`: Lists all Python dependencies required.
- `LICENSE`: Contains the license details for this project.
- `README.md`: This file contains details about the project and instructions for setting it up.


## Commit messages

Use the following tags for commit messages:

    [ADD] : Adding a new feature
    [DEV] : Improving an existing feature
    [DEL] : Removing files, routines
    [FIX] : Fixing existing features
    [OPT] : Optimisation
    [DOC] : Documentation only
