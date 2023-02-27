## Table of contents
* [General info](#general-info)
* [Technologies](#technologies)
* [Setup](#setup)
* [How to contribute](#how-to-contribute)


## General info

This is a fortran package related to the field of astrophysics.  It recovers limb darkening of distant stars from the gravitational microlensing light curves.

The input data files have 4 columns: 
Column1: time   Column2: Apparent magnitude Column3:Error of time  Column4: Error of magnitude

Also, the primary parameters of the microlensing event is needed and is stored in events.f90 file.

The output files have 3 columns    : 
Column 1: r, 
Column 2: Surface brightness at r,  
Column 3: Estimated error of surface brightness.

To know more about the theory behind this package and its algorithm please see "Golchin, Rahvar, 2020, MNRAS, 494, p584-597. https://doi.org/10.1093/mnras/staa743 "   
 
## Technologies
Project is created with: Fortran 95, compiled with gfortran compiler. 
 <br>
OS: Macosx big sur
 <br>
Requierments: gcc, Xcode
 
## Setup
To run this project, put all the files in the same folder.<br>
Change the DATA_PATH string in the Events.f95 file.
 Then, run the command below in terminal, in the same folder that you copied the package:

sh limb.sh
 
## Support
To get help with this project email me at laya.golchin@gmail.com


## How to contribute
You can contribute to the project by uploading more recent data files (2010+) to the Data folder and upgrade the events.f90 file accordingly.
It would be much appreciated.
