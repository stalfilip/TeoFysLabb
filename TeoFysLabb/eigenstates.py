# -*- coding: utf-8 -*-
# Python simulation of an electron in a 1d infinite box potential
# Integrate time independent SE using the Verlet method
# Boundary conditions are found by shooting
# MW 220325

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq


a=1.0e-9           # well width a=1 nm
hbar=1.0545718e-34 # Plancks constant
m=9.10938356e-31   # electron mass
e=1.60217662e-19   # electron charge=-e
c=2.0*m/hbar**2    # constant in Schrödinger equation
EeV=(hbar*np.pi/a)**2/(2.0*m)/e

#print('Exact solution for infinite box potential')
print('E1=',(hbar*np.pi/a)**2/(2*m)/e,'eV')
print('E2=',(hbar*2*np.pi/a)**2/(2*m)/e,'eV')
print('E3=',(hbar*3*np.pi/a)**2/(2*m)/e,'eV')
print('E4=',(hbar*4*np.pi/a)**2/(2*m)/e,'eV')
print('E5=',(hbar*5*np.pi/a)**2/(2*m)/e,'eV')
estimated_energies = [0.376,1.50,3.38,6.02,9.40]


def CompareN(E):
    lnerrors = []
    iter = [10,10**2,10**3,10**4,10**5]
    theoretical_err = [f(N,10**(-9)) for N in iter]
    for N in iter:
        psi_tab,x_tab = Energy(N,E)
        lnerrors.append(np.abs(psi_tab[-1]))
    plt.loglog(iter, lnerrors,'o', label='Numerical Error')
    plt.loglog(iter, theoretical_err,'-', label='Theoretical Error')
    plt.loglog()
    plt.xlabel('N')
    plt.ylabel('Error')
    plt.legend()
    plt.show()

def f(N,const):
    return const/N**2
      
# potential energy function
def V(x):
    #y = 0.0
    y = x**2/2 # harmonic oscillator
    #y = x**2/2 + x**4 # anharmonic oscillator
    return y

# initial values and lists


def Energy1(N,EeV):
    E = EeV*e          # input energy in J
    dx=a/N             # step length
    dx2=dx**2          # step length squared
    psi = 0.0           # wave function at initial position
    dpsi = 1.0          # derivative of wave function at initial position
    psi_tab = []        # list to store wave function for plot
    psi_tab.append(psi)
    x = 0               # initial value of position x
    x_tab = []          # list to store positions for plot
    x_tab.append(x/a)
    for i in range(N) :
        d2psi = c*(V(x)-E)*psi
        d2psinew = c*(V(x+dx)-E)*psi
        psi += dpsi*dx + 0.5*d2psi*dx2
        dpsi += 0.5*(d2psi+d2psinew)*dx
        x += dx
        x_tab.append(x/a)
        psi_tab.append(psi)
    plot(psi_tab,x_tab)
    return psi_tab, x_tab

def HarmonicEnergy(N,E):
    dx=4/N             # step length
    dx2=dx**2          # step length squared
    psi = 0.0           # wave function at initial position
    dpsi = 1.0          # derivative of wave function at initial position
    psi_tab = []        # list to store wave function for plot
    psi_tab.append(psi)
    x = 0               # initial value of position x
    x_tab = []          # list to store positions for plot
    x_tab.append(x/a)
    for i in range(N) :
        d2psi = c*(V(x)-E)*psi
        d2psinew = c*(V(x+dx)-E)*psi
        psi += dpsi*dx + 0.5*d2psi*dx2
        dpsi += 0.5*(d2psi+d2psinew)*dx
        x += dx
        x_tab.append(x/a)
        psi_tab.append(psi)
    plot(psi_tab,x_tab)
    return psi_tab, x_tab

def plot(psi_tab,x_tab):
    plt.close()
    plt.figure(num=None, figsize=(8,8), dpi=80, facecolor='w', edgecolor='k')
    print('E=',EeV,'eV , psi(x=a)=',psi_tab[-1])    
    plt.plot(x_tab, psi_tab, linewidth=1, color='red',label='E(2)')
    plt.legend()
    #plt.xlim(0, 1)
    #limit=1.e-9
    #plt.ylim(0, limit)
    #plt.ylim(-limit, limit)
    #plt.autoscale(False)
    plt.xlabel('x/a')
    plt.ylabel('$\psi$')
    #plt.savefig('psi.pdf')
    plt.show()

#print(find_zero_energy(0.2,0.5))
def main():
    HarmonicEnergy(100,3/2)

if __name__ == "__main__":
    main()

