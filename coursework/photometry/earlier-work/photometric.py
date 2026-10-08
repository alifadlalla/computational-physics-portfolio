"""
This script performs photometric analysis of astronomical data for 
Kuiper Belt Object 29981 (1999 TD10). The analysis is based on data 
from October 8th and 9th, including magnitude calculation and animations.

Directory structure:
- Data_for_astro_data_processing/
    - 8_Oct/
    - 9_Oct/
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from photutils.centroids import centroid_quadratic
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
from astropy.visualization import simple_norm
from astropy.io import fits
from math import log10


def photometric(day, folder_path, ref_i, ref_f, obj_i, obj_f, create_animation=False):
    """
    Perform photometric analysis for a given day's dataset.

    Args:
        day (str): '8_OCT' or '9_OCT'.
        folder_path (str): Path to the folder containing FITS files.
        ref_i (list): Initial coordinates of the reference star.
        ref_f (list): Final coordinates of the reference star.
        obj_i (list): Initial coordinates of the target object.
        obj_f (list): Final coordinates of the target object.
        showplot (bool): Whether to display individual plots.
        create_animation (bool): Whether to create an animation.
        animation_path (str): Path to save the animation file.

    Returns:
        Tuple containing positions and magnitudes of the objects and stars.
    """
    if day not in ['8_OCT', '9_OCT']:
        raise ValueError(f"Invalid value for day. Use '8_OCT' or '9_OCT'.")

    # File list and data storage
    files = os.listdir(folder_path)
    Mag_obj, Mag_star, MJD_array = [], [], []
    star_x, star_y, obj_x, obj_y, old_obj_x, old_obj_y, image_data = [], [], [], [], [], [], []

    # Initial distances
    dist_x = [abs(ref_i[0] - obj_i[0]), abs(ref_f[0] - obj_f[0])]
    dist_y = [abs(ref_i[1] - obj_i[1]), abs(ref_f[1] - obj_f[1])]

    k=0
    for file in files:
        if file.startswith('td'):
            print(f"iteration={k}, Processing file: {file}" )
            # Open FITS file
            hdulist = fits.open(os.path.join(folder_path, file))
            header = hdulist[0].header
            scidata = hdulist[0].data
            mjd = header['MJD-OBS']
            MJD_array.append(mjd)
            exptime = header['EXPTIME']

            # Photometric parameters
            if day == '9_OCT':
                r, r_in, r_out = 8, 15, 25
            elif day == '8_OCT':
                r, r_in, r_out = 5, 10, 15

            # Object and star tracking
            dt = mjd - MJD_array[0]
            vx = (dist_x[0] - dist_x[1]) / MJD_array[0]
            vy = (dist_y[0] - dist_y[1]) / MJD_array[0]

            if k == 0:
                Xcent, Ycent = ref_i[0], ref_f[1]
            else:
                Xcent, Ycent = star_x[-1], star_y[-1]

            radius = 150
            mask = np.ones_like(scidata, dtype=bool)
            mask[int(Ycent) - radius:int(Ycent) + radius, int(Xcent) - radius:int(Xcent) + radius] = False
            xstar, ystar = centroid_quadratic(scidata, mask=mask)
            star_x.append(round(xstar))
            star_y.append(round(ystar))

            # Objet

            if k == 0:
                xobj, yobj = dist_x[0] + r_out + vx * dt, -dist_y[0] + vy * dt
            else:
                xobj, yobj = xobj - vx * dt, yobj + vy * dt

            data_center = scidata[int(ystar)-300-(r+2):int(ystar)+r_out+10, int(xstar)-r_out -(r+2):int(xstar)+600]
            
            radius = 10
            obj_mask = np.ones_like(data_center, dtype=bool)
            obj_mask[int(yobj) - radius:int(yobj) + radius, int(xobj) - radius:int(xobj) + radius] = False
            xobj, yobj = centroid_quadratic(data_center, mask=obj_mask)
            print(xobj, yobj)
            obj_x.append(xobj+int(xstar)-r_out-(r+2))
            obj_y.append(yobj+int(ystar)-300-(r+2))

            # Save data for animation
            image_data.append(data_center)
            print("data for image:",xobj, yobj)
            old_obj_x.append(xobj)
            old_obj_y.append(yobj)

            # Photometric calculations Object
            positions_obj = [(xobj, yobj)]
            aperture_obj = CircularAperture(positions_obj, r)
            annulus_obj = CircularAnnulus(positions_obj, r_in, r_out)

            # Object flux and magnitude
            phot_table = aperture_photometry(data_center, aperture_obj)
            Sum_target_raw = phot_table['aperture_sum']
            phot_table = aperture_photometry(data_center, annulus_obj)
            Sky_background = phot_table['aperture_sum']
            Sum_target = Sum_target_raw - Sky_background / (np.pi * (r_out**2 - r_in**2)) * np.pi * r**2
            object_Magnitude = -2.5 * log10(Sum_target[0])
            Mag_obj.append(object_Magnitude)

             # Photometric calculations Star
            positions_obj = [(r_out+(r+2), 300 +(r+2))]
            aperture_obj = CircularAperture(positions_obj, r)
            annulus_obj = CircularAnnulus(positions_obj, r_in, r_out)

            # Object flux and magnitude
            phot_table = aperture_photometry(data_center, aperture_obj)
            Sum_target_raw = phot_table['aperture_sum']
            phot_table = aperture_photometry(data_center, annulus_obj)
            Sky_background = phot_table['aperture_sum']
            Sum_target = Sum_target_raw - Sky_background / (np.pi * (r_out**2 - r_in**2)) * np.pi * r**2
            star_Magnitude = -2.5 * log10(Sum_target[0])
            Mag_star.append(star_Magnitude)

            k+=1
            hdulist.close()

    # Calibrate magnitudes
    r_mag_st = -16.374
    exposure_time = 960
    r = r_mag_st + 25 + 2.5 * log10(exposure_time)
    Mag = np.array(Mag_obj) - np.array(Mag_star) + r + 0.25

    # Animation generation
    if create_animation:
        create_animation_frames(image_data, old_obj_x, old_obj_y, r, r_in, r_out, day)

    return star_x, star_y, obj_x, obj_y, Mag, MJD_array


def create_animation_frames(data, obj_x, obj_y, r, r_in, r_out, title):
    """
    Create an animation of object movement across FITS frames.

    Args:
        data (list): List of image arrays.
        obj_x (list): X-coordinates of the object.
        obj_y (list): Y-coordinates of the object.
        r, r_in, r_out (int): Radii for aperture and annulus.
        title (str): Title for the animation.
    """
    def update_plot(i, ax):
        ax.clear()
        norm = simple_norm(data[i], 'linear', percent=99.5)
        ax.imshow(data[i], norm=norm, origin='lower')
        
        
        aperture = CircularAperture([(obj_x[i], obj_y[i])], r)
        annulus = CircularAnnulus([(obj_x[i], obj_y[i])], r_in, r_out)
        aperture.plot(color='red', lw=2)
        annulus.plot(color='red', lw=2)


        aperture = CircularAperture([(r_out+ r+2, 300 +r+2 )], r)
        annulus_aperture = CircularAnnulus([(r_out + r-2, 300 + r-2 )], r_in, r_out)
        aperture.plot(color='white', lw=2)
        annulus_aperture.plot(color='white', lw=2)
        
        
        ax.set_title(f"{title} Frame {i + 1}")
        ax.set_xlabel("X Pixels")
        ax.set_ylabel("Y Pixels")

    fig, ax = plt.subplots(figsize=(8, 8))
    ani = FuncAnimation(fig, update_plot, frames=len(data), fargs=(ax,), interval=500)
    plt.show()


for day in ['8_OCT', '9_OCT']:
    if day == '8_OCT':
        folder_path = 'Data_for_astro_data_processing/8_Oct'
        ref_i, ref_f = [1141, 1188], [1049, 1184]
        obj_i, obj_f = [1149, 1156], [1173, 1104]
    elif day == '9_OCT':
        folder_path = 'Data_for_astro_data_processing/9_Oct'
        ref_i, ref_f = [737, 1394], [748, 1268]
        obj_i, obj_f = [1169, 1168], [1309, 992]

    star_x, star_y, obj_x, obj_y, Mag, MJD_array = photometric(day, folder_path, ref_i, ref_f, obj_i, obj_f, create_animation=True)
