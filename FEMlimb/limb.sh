#!/usr/bin/env bash

gfortran -c nrtype.f95
gfortran -c nrutil.f95
gfortran -c Events.f95
gfortran -c elliptics.f95
gfortran -c lensFs2.f95
gfortran -c Mesh.f95
gfortran -c femFs2.f95
gfortran -c limb_fem2.f95

gfortran -o limb2.out  nrtype.o nrutil.o Events.o elliptics.o lensFs2.o Mesh.o femFs2.o limb_fem2.o

