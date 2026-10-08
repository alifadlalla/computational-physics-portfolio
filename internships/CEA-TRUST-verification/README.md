# TRUST verification, 2024

I did this internship at CEA Saclay, STMF, from April to August 2024, with Alexiane Plessier and Michael Ndjinga.

I compared the low-Mach compressible solver with analytical solutions for steady 1D flows: gravity, regular and singular friction, and two parallel channels. [My report](report.pdf) contains the derivations and results.

The [regular-friction notebook](regular-friction/verification.ipynb) is an example of the code I used during the internship. It prepares input files, runs TRUST, extracts the results and compares them with the analytical solution. Its input template is in `regular-friction/src/`.

To run it, open Jupyter from `regular-friction/` in a TRUST environment with `trustutils`.

[My presentation](presentation.pdf) and [an example generated report](regular-friction/example-report.pdf) are also included.
