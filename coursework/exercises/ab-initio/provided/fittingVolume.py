### Data fitting ###

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def part(x,v0):
    return (v0/x)**(2/3)-1

def partbis(x,v0):
    return 6-4*(v0/x)**(2/3)

def birchmurn(x,e0,v0,bm0,bm0d):
    return e0+(9*v0*bm0/16)*((part(x,v0)**3)*bm0d+(part(x,v0)**2)*partbis(x,v0))

#e0=-15.809
#v0=22.0149
#bm0=6.8262
#bm0d=0.181201


data = np.loadtxt('VOLUME.dat')
xdft = data[:,1]
ydft = data[:,2]

popt,pcov = curve_fit(birchmurn,xdft,ydft,p0=[-15.8,22,6,1])

xint = np.linspace(15,30,100)
plt.plot(xdft,ydft,marker="o",label="DFT results")
plt.plot(xint,birchmurn(xint,*popt),label="fitting curve")
plt.xlabel('Volume (Ang**3)')
plt.ylabel('Potential energy (eV)')
plt.title('Relaxation volume of Fe bulk')
plt.legend(loc="upper right")
plt.show()

print(80*"-")
print("Fitted parameters are:",popt)
print("\nMinimum volume is:",popt[1],"Ang**3")
print("corresponding to a lattice parameter:",popt[1]**(1/3),"Ang")
print("\nThe bulk modulus is:",popt[2]*160.2,"GPa")

