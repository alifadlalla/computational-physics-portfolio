#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 14:41:03 2023

@author: Quentin Ansel
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy import linalg as ln

#------------------------------------------------------------------------------
#
# Parameters of the physical system
#
#------------------------------------------------------------------------------

delta=5.5;                                                                     # Offset of the spin.
tf=2*np.pi/np.sqrt(1**2 + delta**2);                                           # Final time, given by the minimum time whith a control amplitude of 1.                            
Nstep = 200;                                                                   # Number of time points
dt = tf/Nstep;                                                                 # time-setp

psi0=np.array([1.0,0.0]);                                                      # Initial state
psitarget=np.array([0.0,1.0]);                                                 # Target state
H0=delta*0.5*np.array([[1,0],[0,-1]]);                                         # free Hamiltonian
H1=0.5*np.array([[0,1],[1,0]]);                                                # control hamiltonian


#------------------------------------------------------------------------------
#
# Functions
#
#------------------------------------------------------------------------------

# Evolution operator for a time step dt
def Propagator(un):
    return ln.expm(-1j*dt*(H0 + un*H1))  

# Compute the evolution operator at each time-step
def Uevol(u):
    Nt        = u.size;                                                         
    U         = np.zeros((Nt,2,2), dtype = 'complex_');
    for n in range(Nt):
        U[n,:,:] = Propagator(u[n]);
    return U

# Compute the quantum state at each time step using the initial state and the evolution operator
def SolveSchrodingerEq(psi0,U):
    #INITIALISATION
    Nt        = U.shape[0];                                                    # Number of time points
    C         = np.zeros((2,Nt+1), dtype = 'complex_');                        # Initialisation of the matrix
    C[:,0]    = psi0;                                                          # Initial conditions
    
    #INSTRUCTION
    for n in range(Nt):
        C[:,n+1] = U[n] @ C[:,n];                                              # Forward propagation at each time step                     
    return C

# Compute the adjoint state at each time step using the target state and the evolution operator
def SolveReverseSchrodingerEq(chif,U):
    #INITIALISATION
    Nt        = U.shape[0];                                                    # Number of time points
    C         = np.zeros((2,Nt+1), dtype = 'complex_');                        # Initialisation of the matrix
    C[:,Nt]   = np.conjugate(chif);                                            # final condition
    
    #INSTRUCTION
    for n in range(Nt,0,-1):
        C[:,n-1] = C[:,n] @ U[n-1];                                            # Backward propagation at each time step                     
    return C


# Compute the cost associated with a control field u
def Cost(u):
    Un      = Uevol(u);
    psit    = SolveSchrodingerEq(psi0,Un);                                     # Computation of the final state  
    dyncost = dt*sum(u*u) /2;                                                  # Computation of the dynamical cost 
    return -np.absolute(np.dot(np.conjugate(psitarget),psit[:,Nstep]))**2 + 0.1*dyncost/tf # cost with p0 = 0.1/tf.

# Compute the gradient
def GradientGRAPE(u):
    Nt   = u.size; 
    dF   = np.zeros(Nt, dtype = 'float');
    Un   = Uevol(u);
    psit = SolveSchrodingerEq(psi0,Un);
    chit = SolveReverseSchrodingerEq(psitarget,Un);
    c    = np.conjugate(chit[:,Nt]@ psit[:,Nt]);
    for n in range(Nt):
        dF[n] = -dt*2*np.imag(c*(chit[:,n+1]@ (H1 @ psit[:,n+1]))) + 0.1* u[n]*dt/tf
    return dF



#------------------------------------------------------------------------------
#
# Solve the optimal control problem
#
#------------------------------------------------------------------------------


u0=np.cos(np.linspace(0, np.pi, Nstep));                                       # Initial guess of the control field.

t = np.linspace(0, tf, Nstep+1)                                                # vector of time points, used for the plots

sol=minimize(Cost,u0,method="BFGS",jac=GradientGRAPE)                          # Minimization algorithm with a specified gradient

uopt=sol.x                                                                     # save the optimal control field into another variable

PSIT=SolveSchrodingerEq(psi0,Uevol(uopt));                                     # Compute the trajectory of the quantum state
z = np.abs(PSIT)**2                                                            # Compute the population at each time step

print(z[:,Nstep])                                                              # Print the population at the final time

plt.plot(t, z.T[::,0:2])                                                       # plot the population as a function of time

plt.xlabel('t')

plt.legend(['$|\psi_1|²$', '$|\psi_2|²$'], shadow=True)

plt.title('Probability to find the excited and the ground states')

plt.show()

plt.plot(t[0:Nstep], uopt.T)                                                   # Plot the optimal control field and the initial one.
plt.plot(t[0:Nstep], u0.T)

plt.xlabel('t')
plt.ylabel('u(t)')

plt.show()