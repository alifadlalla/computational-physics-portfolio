# Hamiltonian chaos

I studied the Chirikov standard map and plotted phase portraits for different values of K. [My report](report.pdf) contains the discussion.

`main.py` uses `functions.py`. Both use NumPy and Matplotlib. The saved plots are in `figures/`, and other versions are in `earlier-work/`.

The Lyapunov calculation needs correcting: the saved routine does not advance the orbit properly. The Jacobian sign in the report also needs checking. I have kept this as coursework from the time.
