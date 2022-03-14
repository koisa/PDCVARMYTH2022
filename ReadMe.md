# ReadMe

February 1, 2022

The PDCVARMYTH2022 package is MATLAB routines and functions comprising Supplemental Material in following the article that generates all its figures: 

* Luiz A. Baccalá and Koichi Sameshima (2022) ‘Partial Directed Coherence and the Vector Autoregressive Modelling Myth and
  a Caveat. ](https://www.frontiersin.org/articles/10.3389/fnetp.2022.845327) 

## Installation and usage

The PDCVARMYTH2022 package contains MATLAB mfiles and subfolders you may copy into your local preferred working directory to execute them. To begin with you should run the startup.m script in the MATLAB command line window to set path and check for the requirements.

`> startup`

In addition to adding the paths, it will also check for the existence of the required MATLAB toolboxes. All other routines, including data simulation and total PDC estimation are included. This is a standalone version.

To run all Examples, choose [1], and [0] otherwise when requested on `startup.m` script execution.

The figures it generates are similar to those in the paper. Differences are due to possibly different random seeds.

To play with the scripts,  you may run each Example individually

`> Example1 | Example2 | Example3 | Example4`

This material will be incorporated into future releases of the AsympPDC package:

* [Sameshima K, Baccal´a LA. Asymp PDC Package (2014).[Click here](https://www.lcs.poli.usp.br/~baccala/pdc/CRCBrainConnectivity/AsympPDC/index.html). [Accessed: 2022-01-29.] and [here in Github](https://github.com/koisa/asympPDC))

## License

These routines are distributed under GNU General Public License v3.0 under
authorship of Koichi Sameshima and Luiz A. Baccal - January 2022.
