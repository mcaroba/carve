# carve.f90

## Authorship and license

**carve.f90** (this code) is Copyright (c) 2013-2020 of **Miguel A. Caro** and was written
during his employment at Aalto University.

The code is released under the Creative Commons Attribution Share-alike International
License (check README.md for license information). If you want to use the code under
different terms than those covered by the license, or have other questions, bug reports,
etc., you can contact Miguel Caro:

* miguel.caro@aalto.fi
* mcaroba@gmail.com

The official repository of this code is on [Github](https://github.com/mcaroba/carve).

The proper way of giving attribution when using the code is to cite the following paper
(citation info will be updated when a better reference is available):

M.A. Caro, R. Zoubkoff, O. Lopez-Acevedo and T. Laurila

*Atomic and electronic structure of tetrahedral amorphous carbon
surfaces from density functional theory: Properties and simulation strategies*

Carbon **77**, 1168 (2014). [DOI: check README.md for license information)](https://doi.org/10.1016/j.carbon.2014.06.060)

## Installation

Make a directory to install **carve.f90**:

 mkdir ~/carve; cd ~/carve

Clone the repo:

 git clone https://github.com/mcaroba/carve.git

Build the binary:

 gfortran -o carve.ex carve.f90

Add the directory to your path:

 echo "export PATH:~/carve:$PATH" >> ~/.bashrc
 source ~/.bashrc

## What does carve.f90 do
