import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from photutils.centroids import centroid_quadratic
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
from astropy.visualization import simple_norm
from astropy.io import fits
from math import log10

def photometric(day, folder_path, ref_i, ref_f, obj_i, obj_f, showplot=False, create_animation=False, animation_path="animation.mp4"):
    if day not in ['8_OCT', '9_OCT']:
        raise ValueError(f"Invalid value for day.")

    files = os.listdir(folder_path)
    k = 0
    Mag_obj, Mag_star, MJD_array = [], [], []
    star_x, star_y, obj_x, obj_y = [], [], [], []
    old_obj_x, old_obj_y, image_data = [], [], []

    dist_x = [abs(ref_i[0] - obj_i[0]), abs(ref_f[0] - obj_f[0])]
    dist_y = [abs(ref_i[1] - obj_i[1]), abs(ref_f[1] - obj_f[1])]

    for file in files:
        if file.startswith('td'):
            print(file)
            k += 1
            hdulist = fits.open(os.path.join(folder_path, file))
            header = hdulist[0].header
            scidata = hdulist[0].data
            mjd = header['MJD-OBS']
            MJD_array.append(mjd)
            print('MJD = ', mjd)

            if day == '9_OCT':
                r, r_in, r_out = 8, 15, 25
            elif day == '8_OCT':
                r, r_in, r_out = 5, 10, 15

            dt = mjd - MJD_array[0]
            vx = (dist_x[0] - dist_x[1]) / MJD_array[0]
            vy = (dist_y[0] - dist_y[1]) / MJD_array[0]

            if k == 1:
                Xcent, Ycent = ref_i[0], ref_f[1]
            else:
                Xcent, Ycent = star_x[-1], star_y[-1]

            print("Current reference star position:", Xcent, Ycent)
            radius = 150
            X_pixels, Y_pixels = scidata.shape
            mask = np.array([[True] * Y_pixels] * X_pixels)
            mask[int(Ycent) - radius:int(Ycent) + radius, int(Xcent) - radius:int(Xcent) + radius] = False

            xstar, ystar = centroid_quadratic(scidata, mask=mask)
            print("New reference star Position (after Centroid):", xstar, ystar)
            data_center = scidata[int(ystar) - 300 - (r + 2):int(ystar) + r_out + 10, 
                                  int(xstar) - r_out - (r + 2):int(xstar) + 600]

            star_x.append(round(xstar))
            star_y.append(round(ystar))

            if k == 1:
                xobj, yobj = dist_x[0] + r_out + (r + 2) + vx * dt, -dist_y[0] + 300 + (r + 2) + vy * dt
            else:
                xobj, yobj = xobj - vx * dt, yobj + vy * dt

            print("Current object location:", xobj, yobj)
            radius = 10
            x_pixel, y_pixel = data_center.shape
            mask = np.array([[True] * y_pixel] * x_pixel)
            mask[int(yobj) - radius:int(yobj) + radius, int(xobj) - radius: int(xobj) + radius] = False
            xobj, yobj = centroid_quadratic(data_center, mask=mask)
            xobj_new, yobj_new = xobj + int(xstar) - r_out - (r + 2), yobj + int(ystar) - 300 - (r + 2)
            obj_x.append(xobj_new)
            obj_y.append(yobj_new)
            image_data.append(data_center)
            old_obj_x.append(xobj)
            old_obj_y.append(yobj)
            print("New object location (after Centroid):", xobj_new, yobj_new)

            # Photometry
            obj_magnitude = calculate_magnitude(data_center, [(xobj, yobj)], r, r_in, r_out)
            Mag_obj.append(obj_magnitude)

            star_magnitude = calculate_magnitude(data_center, [(r_out + r + 2, 300 + r + 2)], r, r_in, r_out)
            Mag_star.append(star_magnitude)

            hdulist.close()

    Mag = compute_absolute_magnitude(Mag_obj, Mag_star)
    if create_animation:
        create_animation_frames(image_data, old_obj_x, old_obj_y, r, r_in, r_out, day)

    return star_x, star_y, obj_x, obj_y, Mag, MJD_array

def calculate_magnitude(data, positions, r, r_in, r_out):
    aperture = CircularAperture(positions, r)
    annulus_aperture = CircularAnnulus(positions, r_in, r_out)

    phot_table = aperture_photometry(data, aperture)
    phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
    target_flux = (phot_table['aperture_sum'])
    sky_background = aperture_photometry(data, annulus_aperture)['aperture_sum'][0]
    sky_corrected_flux = target_flux - sky_background / (r_out**2 - r_in**2)*(r**2)

    return -2.5 * log10(sky_corrected_flux)

def compute_absolute_magnitude(Mag_obj, Mag_star):
    r_mag_st = -16.374
    exposure_time = 960
    r = r_mag_st + 25 + 2.5 * log10(exposure_time)
    return np.array(Mag_obj) - np.array(Mag_star) + r + 0.25

def create_animation_frames(data, obj_x, obj_y, r, r_in, r_out, title):
    def update_plot(i, ax):
        ax.clear()
        norm = simple_norm(data[i], 'linear', percent=99.5)
        ax.imshow(data[i], norm=norm, origin='lower')
        
        # Object
        positions = [(obj_x[i], obj_y[i])]
        aperture = CircularAperture(positions, r)
        annulus_aperture = CircularAnnulus(positions, r_in, r_out)
        aperture.plot(color='red', lw=2)
        annulus_aperture.plot(color='red', lw=2)

        # Star
        positions = [(r_out + r+2, 300 + r+2 )]
        aperture = CircularAperture(positions, r)
        annulus_aperture = CircularAnnulus(positions, r_in, r_out)
        aperture.plot(color='white', lw=2)
        annulus_aperture.plot(color='white', lw=2)


        ax.set_title(f"{title} Frame {i + 1}")
        ax.set_xlabel("X Pixels")
        ax.set_ylabel("Y Pixels")

    fig, ax = plt.subplots(figsize=(8, 8))
    ani = FuncAnimation(fig, update_plot, frames=len(data), fargs=(ax,), interval=500)
    plt.show()

def save_results(day, Mag, MJD_array, star_x, star_y, obj_x, obj_y):
    Mag_mjd = np.column_stack([Mag, MJD_array])
    with open(f"{day}_mag_mjd.txt", "w") as o:
        for line in Mag_mjd:
            o.write(f"{line[0]} {line[1]}\n")

    pos = np.column_stack([star_x, star_y, obj_x, obj_y])
    with open(f"{day}_pos_ref_obj.txt", "w") as o:
        for line in pos:
            o.write(f"{line[0]} {line[1]} {line[2]} {line[3]}\n")

def main():
    for day in ['8_OCT', '9_OCT']:
        if day == '8_OCT':
            folder_path = 'Data_for_astro_data_processing/8_Oct'
            ref_i, ref_f = [1141, 1188], [1049, 1184]
            obj_i, obj_f = [1149, 1156], [1173, 1104]
        elif day == '9_OCT':
            folder_path = 'Data_for_astro_data_processing/9_Oct'
            ref_i, ref_f = [737, 1394], [748, 1268]
            obj_i, obj_f = [1169, 1168], [1309, 992]

        star_x, star_y, obj_x, obj_y, Mag, MJD_array = photometric(
            day, folder_path, ref_i, ref_f, obj_i, obj_f, create_animation=True
        )
        save_results(day, Mag, MJD_array, star_x, star_y, obj_x, obj_y)

if __name__ == "__main__":
    main()
