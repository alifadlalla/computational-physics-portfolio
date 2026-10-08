# Heat diffusion in a copper bar

I wrote this Fortran code during M1 to solve the 1D heat equation with forward and backward Euler schemes. The implicit step uses the Thomas algorithm. [My report](report.pdf) compares the simulations with copper-bar measurements.

`functions.f90` contains the numerical routines. `QualityControl.f90` checks both methods against an analytical solution. A small array allocation error was fixed in this file.

```sh
gfortran -std=legacy -ffree-line-length-none -fcheck=all functions.f90 QualityControl.f90 -o quality-control
./quality-control
```

The first two cells of `plot.ipynb` plot the resulting `error.txt`. The supplied copper-bar measurements are in `Cu_isol.txt` and `Cu_vent.txt`. The radiation branch in `heat_1d.f90` still needs checking.
