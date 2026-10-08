### HF calculation for He atom ###

from math import *
import numpy as np
import matplotlib.pyplot as plt


#### Basis set definition ####
def slater(r,alpha):        
    return np.sqrt(alpha**3/np.pi)*np.exp(-alpha*r)


n=int(input("Number of Slater function for the basis set definition: "))
print(str(n) + " Slater functions should be defined")
x=np.linspace(-5,5,100)
r=abs(x)

alpha=[]
for i in range(n):
    alp=float(input("Effective charge of the nucleus alpha for function " + str(i) + ": "))
    alpha.append(alp)
    plt.xlabel("r")
    plt.ylabel("STO(r)")
    plt.plot(x,slater(r,alp),label='alpha='+str(alp))
    plt.legend()
plt.show()

#### Integral calculations for the Slater basis set: Sij, Iij, Jij ####

def Sint(palpha,salpha):
    return (2*np.sqrt(palpha)/salpha)**3

###      Other solution but Spq is defined locally in that case...
#def Sint2(alpha,n):
#    Spq=np.zeros((n,n))
#    for i in range(n):
#        for j in range(n):
#            Spq[i,j]=(2*np.sqrt(alpha[i]*alpha[j])/(alpha[i]+alpha[j]))**3
#    return Spq
#print(Sint2(alpha,n))
    
def Iint(palpha,salpha):
    return 4*(np.sqrt(palpha)/salpha)**3*(palpha-2*salpha)
    
def Jint(salpha,sbialpha,ptotalpha,stotalpha):
    return 32*(np.sqrt(ptotalpha))**3*(1/(salpha**3*sbialpha**2)-1/(salpha**3*stotalpha**2)-1/(salpha**2*stotalpha**3))

Sij=np.zeros((n,n))
Iij=np.zeros((n,n))
Jij=np.zeros((n,n,n,n))
for i in range(n):
    for j in range(n):
        palpha=alpha[i]*alpha[j]
        salpha=alpha[i]+alpha[j]
        Sij[i,j]=Sint(palpha,salpha)
        Iij[i,j]=Iint(palpha,salpha)
        for k in range(n):
            for l in range(n):
                sbialpha=alpha[k]+alpha[l]
                ptotalpha=palpha*alpha[k]*alpha[l]
                stotalpha=salpha+alpha[k]+alpha[l]
                Jij[i,j,k,l]=Jint(salpha,sbialpha,ptotalpha,stotalpha)

Svp=np.linalg.eig(Sij)
Sp1=Svp[1].dot(np.diag(Svp[0]**(1/2))).dot(np.linalg.inv(Svp[1]))    # S^(1/2)
Sp2=Svp[1].dot(np.diag(Svp[0]**(-1/2))).dot(np.linalg.inv(Svp[1]))   # S^(-1/2)
#print(Sij)
#print(Iij)
#print(Jij)

#### Initial Ci coefficients ####
Ci=np.zeros((n))
for i in range(n):
    Ci[i]=float(input("Value of the C" + str(i) + " coefficient: "))
    
print("\n")

#### Beginning of the iteration ####
maxit=20      # maximum number of iterations
iterat=0
energy=0.0
enerit=[]
crit=1e-8      # convergence criterium in eV
conv=27.2114   # conversion of Hartree in eV

while(iterat<maxit):
    iterat +=1
    print("\nIteration n°",iterat)
    print(15*"-")
    
#### Calculations of the Fock matrix Fij according to Ci coefficients ####
    Fij=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            Fij[i,j]=Iij[i,j]
            for k in range(n):
                for l in range(n):
                    Fij[i,j]=Fij[i,j]+Ci[k]*Ci[l]*Jij[i,j,k,l]
#    print(Fij)
#### Calculation of the electronic energy ####
    energyini=energy
    energy=0.0
    energyI=0.0
    energyJ=0.0
    for i in range(n):
        for j in range(n):
            energyI=energyI+Ci[i]*Ci[j]*Iij[i,j]
            energyJ=energyJ+Ci[i]*Ci[j]*Fij[i,j]-Ci[i]*Ci[j]*Iij[i,j]
            energy=(energy+Ci[i]*Ci[j]*(Iij[i,j]+Fij[i,j]))
    enerit.append(energy*conv)
    print("Total electronic energy:",energy*conv,"eV")
    # print("Monoelectronic energy:",energyI*conv,"eV")
    # print("Coulomb repulsion energy:",energyJ*conv,"eV")

#### Determination of new Ci coefficients ####
    Fp=Sp2.dot(Fij).dot(Sp2)
    Fpdiag=np.linalg.eig(Fp)
    vp=Fpdiag[0]
    Cp=Fpdiag[1]
    enew=np.min(vp)
    Ci=Sp2.dot(Cp[:,np.argmin(vp)])
    print("New Ci coefficients:",Ci)
    print("New orbital energy:",enew*conv,"eV")
    
    if(abs(energy*conv-energyini*conv)<crit):
        print("\n")
        print(40*"-")
        print("Convergence is reached after " + str(iterat) + " iterations")
        print("Final electronic energy:",energy*conv,"eV")
        print("Experimental value is: -79.0 eV")
        error=abs(energy*conv+79.0)/79.0*100
        print("Relative errror with the exp. value: {:.5f} %".format(error))
        print("Monoelectronic energy I:",energyI*conv,"eV")
        print("Coulomb repulsion energy J:",energyJ*conv,"eV")
        break

#### Figure of convergence of electronic energy ####
x=np.linspace(0,iterat,iterat)

plt.xlabel("Number of iterations")
plt.ylabel("Electronic energy (eV)")
plt.ylim(-80,-75)
plt.plot(x,enerit,marker="*",label="HF calc.")
plt.axhline(y=-77.879,color="green",ls="--",label="HF limit")
plt.axhline(y=-79.0,color="red",ls="--",label="Exp. value")
plt.legend()
#plt.savefig("Convergence-He.png")
plt.show()





       
