# GW-GC Cross Correlation
Contains the pipeline for the analysis presented in Canas-Herrera, Contigiani, Vardanyan [arXiv:1910.08353]. Please cite this paper when making use of the repository. 


Main flowchart
----------

* /Setup : this creates some required files for further computations.
* /MCMC_run/Fiducial : the notebook here specifies the fiducial model and generates mock data.
* /MCMC_run/Chain : contains the likelihood.
* /Plotting : contains utilities for producing the plots in the paper.

Main Files
----------

* CrossCorr_Cells.f95 : Fortran integrator which outputs the C_ell's.
* bassels.py : some auxiliary functions for the Fortran computation.
* aux.py : auxiliary functions for the likelihood computation.
* chain.py : likelihood computation.
* check-fit.ipynb : used to check the status of the chains and maybe show some temporary results

Large Files
----------
* Download some massive files from here: https://www.dropbox.com/sh/h88esrgfyl25ewr/AADMIqX2QhU1tLx830bA5HYsa?dl=0
This includes 
  * a) the EFTCAMB outputs, 
  * b) the post-processed versions of the EFTCAMB outputs,
  * c) the used fiducials in all the 10 runs,
  * d) all the chains.
