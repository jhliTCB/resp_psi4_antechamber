# resp_psi4_antechamber
use psi4  to calculate properties, python "resp" package to fit ESP, antechamber to fit RESP



## Requirements

### psi4: version 1.9.1
An open source quantum chemical program.  https://psicode.org/installs/v191/

### antechamber:
The program for preparing molecular dynamics for Amber. Written in C. see: ambermd.org

### Python packages:
resp:  https://github.com/cdsgroup/resp
fortran formatter: https://pypi.org/project/fortranformat/
 

## Files
- **driver.py**   it is modified from the resp package, add lines to the end for generation the esp file
- **psi4_esp.py**  setup the options for the esp calculations
- **gesp2pdb.py** write the the esp charges of the atoms and the fitting points into the beta-factor column in a PDB format. For visualization.
- **run.sh** A bash wrapper for the submitting the python scripts (*SLURM*)
- **respfit.sh** A bash wrapper for using antechamber
- **fit_cov_resp.py**  A wrapper to fit the resp charge from the ligand-covalent bonded residue with ACE and NME capping atoms. It will generates two sets of parameters for different purposes of MD simulations.

