import numpy as np
import matplotlib.pyplot as plt

def V(x):
    return 0.5 * x**2 + x**4

a=4          # well width a=1 nm
hbar=1.0 # Plancks constant
m=1.0   # electron mass
e=1.0   # electron charge=-e
c=2.0*m/hbar**2    # constant in Schrödinger equation
N=10000               # number of mesh points
dx=a/N             # step length
dx2=dx**2          # step length squared

def solve_se_at_a(E_initial,Even = True):
    x = -a/2               # initial value of position x
    if Even:
        psi = 1.0           # wave function at initial position
        dpsi = 0.0          # derivative of wave function at initial position
    else:
        psi = 0.0           # wave function at initial position
        dpsi = 1.0          # derivative of wave function at initial position

    x_tab = []          # list to store positions for plot
    psi_tab = []        # list to store wave function for plot
    x_tab.append(x/a)
    psi_tab.append(psi)

    for i in range(N) :
        d2psi = c*(V(x)-E_initial)*psi
        d2psinew = c*(V(x+dx)-E_initial)*psi
        psi += dpsi*dx + 0.5*d2psi*dx2
        dpsi += 0.5*(d2psi+d2psinew)*dx
        x += dx
        x_tab.append(x/a)
        psi_tab.append(psi)
    
    return [x_tab, psi_tab]

def interval_halving(E_low, E_high, tol):
    """
    Finds the energy eigenvalue of a particle in a box using the interval halving method.
    """
    while True:
        
        psi_at_a_high = solve_se_at_a(E_high*e)
        psi_at_a_low = solve_se_at_a(E_low*e)

        if abs(psi_at_a_low[1][-1]) < abs(psi_at_a_high[1][-1]):
            E_high = (E_high + E_low)/2
        else:
            E_low = (E_high + E_low)/2
            
        if abs(E_high - E_low) < tol:
            break

    return (E_high + E_low)/2, psi_at_a_high

def plotEvenAndOdd():
    E = 30
    psi_at_a = solve_se_at_a(E, True)
    psi_at_a_odd = solve_se_at_a(E, False)

    plt.plot(np.array(psi_at_a[0]), psi_at_a[1], label='$\Psi(x)$')
    plt.plot(np.array(psi_at_a_odd[0]), psi_at_a_odd[1], label='$\Psi(x)$')
    plt.xlabel('x [nm]')
    plt.ylabel('$\Psi(x)$')
    plt.title('Numerical wavefunction for E = 30')
    plt.legend()
    plt.grid(True)

    plt.show()

def main():
    # plot for the wavefunction
    EeV_bounds = np.array([[0.0, 1.0], [2.0, 3.0], [4.0, 6.0], [7.0, 8.3], [10.0, 12.0]])
    EeV_bounds = np.array([[30.0, 33.0]])
    
    for n, EeV_bound in enumerate(EeV_bounds):
        E_low, E_high = EeV_bound
        E_mid, psi_tab = interval_halving(E_low, E_high, 1e-10)
        print('E_' + str(n+1) + ' =', E_mid)

        # Plot the wavefunction
        plt.plot(np.array(psi_tab[0]), psi_tab[1], label='$\Psi_{' + str(n+1) + '}(x)$')
        plt.xlabel('x [nm]')
        plt.ylabel('$\Psi(x)$')

    plt.title('Numerical wavefunction for E = 30')
    plt.legend()
    plt.grid(True)

    plt.show()

plotEvenAndOdd()
#main()