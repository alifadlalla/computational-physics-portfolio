# Quantum dynamics

I completed these notebooks as a group project.

- [Numerical representations](numerical-representations.ipynb): finite-difference and DVR Hamiltonians for a particle in a box, a harmonic oscillator and a molecular potential.
- [Wave-packet propagation](wave-packet-propagation.ipynb): matrix-exponential propagation and kinetic/potential splitting, with barrier and time-dependent potential examples.

They use NumPy, SciPy and Matplotlib. Install the root `requirements.txt` and open them in Jupyter.

The supplied H₂ reference values and their comparisons are included in the numerical-representations notebook. In the final Fermi-accelerator sweep, `simulation(p)` fixes `p = 1`, so those calls repeat the same initial state.
