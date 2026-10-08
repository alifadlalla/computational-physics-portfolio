# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np

import matplotlib as plt

# For the up mirror

# the E components
def method(st="up"):

    def inp( file ):
        
     
        file_amp = file[:,0]
        file_phase = file[:,1]*3.14/180
        
        return file_amp* np.exp(1j*file_phase)
        
        
    orientation = st 
    Ex_up = np.loadtxt(f"fwtmp_"+orientation+"_f2_ex.dat",skiprows=2)
    #Ey_up = np.loadtxt(f"fwtmp_"+orientation+"_f2_ey.dat")
    Ez_up = np.loadtxt(f"fwtmp_"+orientation+"_f2_ez.dat",skiprows=2)
    
    
    
    #Hx_up = np.loadtxt(f"fwtmp_"+orientation+"_f2_hx.dat")
    Hy_up = np.loadtxt(f"fwtmp_"+orientation+"_f2_hy.dat",skiprows=2)
    #Hz_up = np.loadtxt(f"fwtmp_"+orientation+"_f2_hz.dat")
    
     
    Ex_up = inp(Ex_up)
    #Ey_up = inp(Ey_up)
    Ez_up = inp(Ez_up)
    
    #Hx_up = inp(Hx_up)
    Hy_up = inp(Hy_up)
    #Hz_up = inp(Hz_up)
    
    E0 = [Ex_up, np.zeros(len(Ex_up)), Ez_up]
    H0 = [np.zeros(len(Hy_up)), Hy_up, np.zeros(len(Hy_up))]
    
    E0 = np.transpose(E0)
    H0 = np.transpose(H0)
   
    return 0.5*np.real(np.cross(E0,np.conjugate(H0)))

rho_up=method(st="up")
rho_down=method(st="down")
rho_right=method(st="right")
rho_left=method(st="left")

#print(rho_up)
P=np.sum(rho_up[:,2]) + np.sum(rho_right[:,0]) - np.sum(rho_down[:,2]) - np.sum(rho_left[:,0])
print(P*1e-3*2/376.7)
#xprint(np.sum(P)*2/376700)
#print(np.average(P)*2/376700)