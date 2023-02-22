## Table of contents
* [General info](#general-info)
* [Technologies](#technologies)
* [Setup](#setup)
* [How to contribute](#how-to-contribute)


## General info

This is a fortran package related to the field of astrophysics.  It recover limb darkening of distant stars from the gravitational microlensing light curves.

The input data files have 3 columns: time   Magnitude  error of magnitude
Also, the primary parameters of the microlensing event is needed and is stored in events.f90 file.
The output files have 3 columns    : 
Column 1: r
Column 2: Surface brightness at r 
Column 3: Estimated error of surface brightness

To know more about the theory behind this package and its algorithm please see "Golchin, Rahvar, 2020, MNRAS, 494, p584-597. https://doi.org/10.1093/mnras/staa743 "   
 
## Technologies
Project is created with:
 Fortran 90

## Setup
To run this project, copy all the files and put them all in the same folder.
 Then, run the command below in terminal, in the same folder that you copied the package:

sh limb.sh
 
## Support
To get help with this project email me at laya.golchin@gmail.com


## How to contribute
You can contribute to the project by uploading more recent data files (2010+) to the Data folder and upgrade the events.f90 file accordingly.
