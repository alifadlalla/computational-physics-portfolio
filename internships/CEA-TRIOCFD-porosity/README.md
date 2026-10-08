# TrioCFD porosity validation, 2025

I did this internship at CEA Saclay from April to August 2025, supervised by Alexiane Plessier. I studied a porosity-turbulence model for tube-support-plate clogging in steam generators.

I automated runs to compare the porous model with explicit geometry and varied the Darcy–Forchheimer coefficients `c0` and `c1`. [My report](report.pdf) gives the results and discusses the asymmetric velocity field that appeared at higher inlet velocities.

These notebooks are examples of the code I used during the internship to test different coefficients and compare the results. Each folder has its own input templates.

They need the internship's TRUST/TrioCFD environment, `trustutils`, and the `Frottement_Poreux` and `K_eps_conv` source terms. Open Jupyter from the notebook's folder.

[My presentation](presentation.pdf) and [an example generated report](case_2/full_analysis/example-report.pdf) are also included.
