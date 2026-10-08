# U-Net dose prediction, 2023

I did this M1 internship at Laboratoire Chrono-environnement, Montbéliard, France, from April to June 2023, supervised by Pierre-Emmanuel Leni. I tested a TensorFlow/Keras U-Net for VMAT prostate dose prediction. I used images prepared from DICOM files and GATE simulations.

[My report](report.pdf) describes the earlier runs: black outputs and a validation MAE stuck at `0.0378` over 10 epochs. After sending the report, I got more promising results and added them to slides 19–20 of [my presentation](presentation.pdf). The reconstruction had visible structure, but still needed further testing.

I wrote `unet.py`, `data_utils.py`, `train.py` and `paths.py`. The preprocessing scripts `dicom.py` and `gate.py` were given to me and I adapted them; they are not included here.

The dataset and trained model are not included. The patient mapping in `data_utils.py` is empty, so you cannot run training directly.
