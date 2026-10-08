"""
there should be a directory named Data_for_astro_data_processing containng the dolfer for 8 and 9  oct

"""

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from photutils.centroids import centroid_com, centroid_1dg, centroid_2dg, centroid_quadratic
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
#from photutils.aperture import CircularAperture
from astropy.visualization import simple_norm
from astropy.io import fits
from math import log10
import os

def photometric(day, folder_path, ref_i, ref_f, obj_i, obj_f, showplot=False, create_animation=False, animation_path="animation.mp4"):


    # Check if the variable is in the allowed values
    if day not in ['8_OCT', '9_OCT']:
        raise ValueError(f"Invalid value for day.")

    # The files
    files = os.listdir(folder_path)

    k = 0
    Mag_obj=[]
    Mag_star = []
    MJD_array=[]
    star_x=[]
    star_y=[]
    obj_x=[]
    obj_y=[]
    old_obj_x=[]
    old_obj_y=[]
    image_data = []

    dist_x=[abs(ref_i[0]-obj_i[0]), abs(ref_f[0]-obj_f[0])]
    dist_y=[abs(ref_i[1]-obj_i[1]), abs(ref_f[1]-obj_f[1])]

    #print("inital distance between object and star:", dist_x[0], dist_y[0])
    #print("final distance between object and star:", dist_x[1], dist_y[1])

    for file in files:
        if file.startswith('td'):
            print(file)
            k += 1
            # Open the FITS image
            hdulist = fits.open(folder_path+'/'+file)
            header = hdulist[0].header
            scidata = hdulist[0].data
            mjd = header['MJD-OBS']
            MJD_array.append(mjd)
            print('MJD = ',mjd)
            exptime = header['EXPTIME']

            # Define parameters for photometry
            if day=='9_OCT':
                r = 8  # Radius for measuring target flux
                r_in = 15  # Inner radius for the annulus used for measuring sky background
                r_out = 25  # Outer radius for the annulus used for measuring sky background
            elif day=='8_OCT':
                r = 5  # Radius for measuring target flux
                r_in = 10  # Inner radius for the annulus used for measuring sky background
                r_out = 15  # Outer radius for the annulus used for measuring sky background

            dt = mjd-MJD_array[0] # time
            vx =(dist_x[0]-dist_x[1])/MJD_array[0]
            vy =(dist_y[0]-dist_y[1])/MJD_array[0]

            #print("vx, vy:", vx, vy)
            #print("mjd", mjd)
            #print("dt",dt)
            # Approximate target coordinates (x, y)
            #Xcent = 1083
            #Ycent = 660
            if k == 1:
                Xcent = ref_i[0]
                Ycent = ref_f[1]
            else:
                Xcent = star_x[-1]
                Ycent = star_y[-1]

            print("Current reference star postion:", Xcent, Ycent)

            radius = 150 # radius of the centroid
            X_pixels, Y_pixels = scidata.shape
            mask = np.array([[True]*Y_pixels]*X_pixels)
            mask[int(Ycent) - radius:int(Ycent) + radius, int(Xcent) - radius:int(Xcent) + radius] = False

            xstar, ystar = centroid_quadratic(scidata, mask=mask) 
            #print("star, ystar:", xstar, ystar)
            print("New referencce star Position (after Centroid)", xstar, ystar)
            data_center = scidata[int(ystar)-300-(r+2):int(ystar)+r_out+10, int(xstar)-r_out -(r+2):int(xstar)+600]
            
            #print(data_center)
            star_x.append(round(xstar))
            star_y.append(round(ystar))


            # object detection
            if k ==1:
                xobj, yobj= dist_x[0]+r_out+(r+2)+vx*dt , -dist_y[0]+300+(r+2)+vy*dt 
            else:
                xobj, yobj= xobj-vx*dt , yobj+vy*dt

            
            print("current object location:", xobj, yobj)
            radius = 10
            print(k)
            x_pixel, y_pixel = data_center.shape
            mask =  np.array([[True]*y_pixel]*x_pixel)
            mask[int(yobj)-radius:int(yobj)+radius, int(xobj)-radius: int(xobj)+radius] = False
            xobj, yobj = centroid_quadratic(data_center, mask=mask)
            xobj_new , yobj_new= xobj+int(xstar) - r_out-(r+2), yobj + int(ystar)-300 - (r+2) 
            obj_x.append(xobj_new)
            obj_y.append(yobj_new)

            # Save data for animation
            image_data.append(data_center) 
            print("data for image:", xobj, yobj)
            old_obj_x.append(xobj)
            old_obj_y.append(yobj)     


            #print('Star position = ', xstar_new,ystar_new)
            print('New object location (after Centroid):', xobj_new, yobj_new)



            positions = [(xobj, yobj)]
            aperture = CircularAperture(positions, r)
            annulus_aperture = CircularAnnulus(positions, r_in, r_out)

            # magnitude of the object

            phot_table = aperture_photometry(data_center, aperture)
            phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
            #print(phot_table)
            Sum_target_raw=(phot_table['aperture_sum'])
            #print('\nSum_target_raw=%.3f\n'%(Sum_target_raw))

            phot_table = aperture_photometry(data_center, annulus_aperture)
            phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
            #print(phot_table)
            Sky_background=(phot_table['aperture_sum'])

            #print('\nSky_background=%.3f\n'%(Sky_background))
            Sum_target=Sum_target_raw-Sky_background/(r_out*r_out-r_in*r_in)*(r*r)
            #print('Target flux without sky=%.3f\n'%(Sum_target))
            object_Magnitude=-2.5*log10(Sum_target[0])
            Mag_obj.append(object_Magnitude)
            print('instrumental magnitude of the object=%.3f\n'%(object_Magnitude))
            

            positions = [(r_out+(r+2), 300 +(r+2))]
            #print (positions)
            aperture = CircularAperture(positions, r)
            annulus_aperture = CircularAnnulus(positions, r_in, r_out)

            
            # Magnitude of the star
            phot_table = aperture_photometry(data_center, aperture)
            phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
            #print(phot_table)
            Sum_target_raw=(phot_table['aperture_sum'])
            #print('\nSum_target_raw=%.3f\n'%(Sum_target_raw[0]))

            phot_table = aperture_photometry(data_center, annulus_aperture)
            phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
            #print(phot_table)
            Sky_background=(phot_table['aperture_sum'])

            #print('\nSky_background=%.3f\n'%(Sky_background[0]))
            Sum_target=Sum_target_raw-Sky_background/(r_out*r_out-r_in*r_in)*(r*r)
            #print('Target flux without sky=%.3f\n'%(Sum_target[0]))
            star_Magnitude=-2.5*log10(Sum_target[0])
            Mag_star.append(star_Magnitude)
            
            #print('instrumental magnitude of the star=%.3f\n'%(star_Magnitude))
            #hdulist.close()

            #print( 40*'-' )
            #print("\n \n \n \n")

            hdulist.close()
    
    ## From the analysis of Pansstar file
    if day=='8_OCT':
        r_mag_st = -16.374      # 8 oct
    elif day=='9_OCT':
        r_mag_st = -16.374      # 9 oct
    exposure_time = 960

    r = r_mag_st + 25 +2.5*log10(exposure_time)

    Mag = np.array(Mag_obj) - np.array(Mag_star) + r + 0.25

    # Generate animation if requested
    if create_animation:
        create_animation_frames(image_data,  old_obj_x, old_obj_y, r, r_in, r_out, day, animation_path)

    return star_x, star_y, obj_x, obj_y, Mag, MJD_array

def create_animation_frames(data, obj_x, obj_y, r, r_in, r_out, title, save_path):
    def update_plot(i, ax):
        ax.clear()
        norm = simple_norm(data[i], 'linear', percent=99.5)
        ax.imshow(data[i], norm=norm, origin='lower')
        positions = [(obj_x[i], obj_y[i])]
        aperture = CircularAperture(positions, r)
        annulus_aperture = CircularAnnulus(positions, r_in, r_out)
        aperture.plot(color='red', lw=2)
        annulus_aperture.plot(color='red', lw=2)
        #ax.plot(star_x[i], star_y[i], 'bo', label="Reference Star")
        #ax.plot(obj_x[i], obj_y[i], 'ro', label="Target Object")
        positions = [(r_out + r-2, 300 + r-2 )]
        #print (positions)
        aperture = CircularAperture(positions, r)
        annulus_aperture = CircularAnnulus(positions, r_in, r_out)
        aperture.plot(color='white', lw=2)
        annulus_aperture.plot(color='white', lw=2)
        #ax.legend()
        ax.set_title(f"{title} Frame {i + 1}")
        ax.set_xlabel("X Pixels")
        ax.set_ylabel("Y Pixels")

    fig, ax = plt.subplots(figsize=(8, 8))
    ani = FuncAnimation(fig, update_plot, frames=len(data), fargs=(ax,), interval=500)
    #ani.save(save_path, writer="ffmpeg")
    plt.show()


for day in ['8_OCT', '9_OCT']:
    if day=='8_OCT':
        # 8 oct
        folder_path = 'Data_for_astro_data_processing/8_Oct'
        # refernce star
        ref_i =[1141, 1188]
        ref_f =[1049, 1184]

        # target object
        obj_i=[1149, 1156]
        obj_f=[1173, 1104]
    if day=='9_OCT':
        # 9 OCT
        folder_path = 'Data_for_astro_data_processing/9_Oct'

        # refernce star
        ref_i =[737, 1394]
        ref_f =[748, 1268]

        # target object
        obj_i=[1169, 1168]
        obj_f=[1309, 992]

    star_x, star_y, obj_x, obj_y, Mag, MJD_array = photometric(day, folder_path, ref_i, ref_f, obj_i, obj_f,create_animation=True)

    # write the files
    for i in range( len(Mag) ) : 
        print( "Magnitude of the Object, MDJ = ", Mag[i]  , MJD_array[i]  ) 
    
    Mag_mjd = np.zeros( [len(Mag),3] ) 
    Mag_mjd[:,0] = Mag 
    Mag_mjd[:,1] = MJD_array 

    with open(day+"_oct_mag_mjd.txt", "w") as o:
        for line in Mag_mjd:
            print("{} {}".format(line[0] , 
                                    line[1]), 
                                    file=o) 

    pos = np.zeros( [len(Mag),4] ) 
    pos[:,0] = star_x
    pos[:,1] = star_y
    pos[:,2] = obj_x
    pos[:,3] = obj_y

    with open(day+"_oct_pos_ref_obj.txt", "w") as o:
        for line in pos :
            print("{} {} {} {}".format(line[0]  , 
                                    line[1]  , 
                                    line[2]  ,
                                    line[3] ), 
                                    file=o) 

