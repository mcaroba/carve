# carve.f90

## Authorship and license

**carve.f90** (this code) is Copyright (c) 2013-2020 of **Miguel A. Caro** and was written
during his employment at Aalto University.

The code is released under the Creative Commons Attribution Share-Alike International
License (check LICENSE.md for license information). If you want to use the code under
different terms than those covered by the license, or have other questions, bug reports,
etc., you can contact Miguel Caro:

* miguel.caro@aalto.fi
* mcaroba@gmail.com

The official repository of this code is on [Github](https://github.com/mcaroba/carve).

The proper way of giving attribution when using the code is to cite the following papers
(citation info will be updated when a better reference is available):

* M.A. Caro, R. Zoubkoff, O. Lopez-Acevedo and T. Laurila. *Atomic and electronic
structure of tetrahedral amorphous carbon surfaces from density functional theory:
Properties and simulation strategies*. Carbon **77**, 1168 (2014).
[DOI: 10.1016/j.carbon.2014.06.060](https://doi.org/10.1016/j.carbon.2014.06.060).

* M.A. Caro, G. Csányi, T. Laurila and V.L. Deringer. *Machine learning driven simulated
deposition of carbon films: From low-density to diamondlike amorphous carbon*.
Phys. Rev. B **102**, 174201 (2020).
[DOI: 10.1103/PhysRevB.102.174201](https://doi.org/10.1103/PhysRevB.102.174201).

## Installation

Clone the repo:

    git clone https://github.com/mcaroba/carve.git

Build the binary:

    cd carve
    gfortran -o carve.ex carve.f90

Add the directory to your path:

    echo "export PATH=$(pwd):\$PATH" >> ~/.bashrc
    source ~/.bashrc

## What does carve.f90 do

**carve.f90** cuts out a spherical chunk of an atomic structure, centered on
a specific atom, and it passivates broken bonds and follows some simple chemical
rules to avoid unphysical resulting structures. Therefore, the carved
structures may extend beyond the boundary of the cutoff sphere. At the moment,
the code can only handle systems containing C, H and O (CHO).

Why on earth would you want to do this? The most common scenario would be that
you want to compute some local properties of the original atomic structure,
but it is so large that your method of choice cannot handle it. For example,
you generate a structure using a machine learning force field containing thousands
of atoms, and then want to calculate some bonding properties using DFT or some
other (significantly more expensive) electronic structure method where running the
calculation on the whole structure is impractical, like in
[DOI: 10.1103/PhysRevB.102.174201](https://doi.org/10.1103/PhysRevB.102.174201).

The code is currently computationally inefficient, mostly due to inefficient
neighbor list builds. For the typical systems we look at, the code runs in
under a second in a single CPU. If you want to use the code on multimillion
atom structures, you may start experiencing some serious slow down. If you
want to fix it, I recommend you start by taking a look at the neighbor lists.
You can also wait until I have the need for a more efficient code and implement
changes, but that may never actually happen!

## Usage (aka codumentation)

Using **carve.f90** is very easy. It takes VASP structure files (POSCAR, CONTCAR) as
input and produces either an XYZ or a POSCAR file. You can use ASE to convert to
VASP format from a number of different formats, if your original file to carve is
not in VASP format. These are the input parameters that **carve.f90** takes:

* input atoms file, `atoms_file`
* cutoff radius, `rcut`
* central atom index , `central_atom`, *in Fortran convention* (first atom is 1 [also the default] *not* 0, as it would be in Python)
* `dHC`: distance between a C atom and the H atom added to passivate a broken C-C bond (in Angstrom, default 0.9 Angstrom)
* `dHH_min`: minimum distance between two passivating H atoms (in Angstrom, default 1.1 Ansgtrom; if put too
close to one another, they will form an H2 molecule!)
* `format`: it can be "xyz", "vasp", "fix_in" or "fix_out" (default is "xyz", see below for details)
* `rfix`: radius for fixing atoms (in Angstrom, no default)

An example command is:

    carve.ex atoms_file="CONTCAR" rcut=5. format="xyz" dHC=0.9 dHH_min=1.1 > carved.xyz

### Fixing atoms

If outputting with VASP format, you can also choose to free some of the atoms. In that case use
either `format="fix_in"`, if you can to fix the atoms within an inner cutoff sphere, or
`format="fix_out"`, if you want to fix the atoms within an outer cutoff sphere. The command
line argument `rfix` determines where the boundary of the sphere lies:

    carve.ex atoms_file="CONTCAR" rcut=5. format="fix_out" rfix=3. dHC=0.9 dHH_min=1.1 > carved.vasp

This option is useful if you want to do relaxation of the inner portion of the carved structure
but keep the outer atoms fixed, as we did in
[DOI: 10.1103/PhysRevB.102.174201](https://doi.org/10.1103/PhysRevB.102.174201).
