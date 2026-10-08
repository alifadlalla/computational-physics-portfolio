# N-body simulation with OpenMP

In this coursework, I studied a gravitational N-body simulation, velocity-Verlet integration, energy conservation and OpenMP performance. [My report](report.pdf) contains the calculations and timings. The reported run with 24 threads was about 8 times faster than the sequential run.

The program started from Barnabé Deforet's teaching code. I have included the version saved with my report.

The Fortran program, Makefile, benchmark script and Python plotting tools are in `src/`. `input.yaml` contains the simulation settings. `earlier-work/` contains other versions I found.

On Linux, compile and run from `src/`:

```sh
gfortran -O3 -fopenmp settings.f95 main.f95 -o main
./main
```