import numpy as np
import matplotlib.pyplot as plt

# Parameters
N = 1000
dx = 4 / N
dx2 = dx**2
x = np.linspace(-2, 2, N)

# Harmonic oscillator potential
def V(x):
    return 0.5 * x**2

def AnHo(x):
    return x**2/2 + x**4 # anharmonic oscillator

# Initial guesses for the energy
E_harmonic = 0.8
c = 2

# Verlet method
def verlet(E, V, dx, dx2, c, Even = True):
    x = 0
    if Even:
        psi = 0.0
        dpsi = 1.0
    else:
        psi = 1.0
        dpsi = 0.0
    x_tab = []
    psi_tab = []
    x_tab.append(x)
    psi_tab.append(psi)
    for i in range(N):
        d2psi = c * (V(x) - E) * psi
        d2psinew = c * (V(x + dx) - E) * psi
        psi += dpsi * dx + 0.5 * d2psi * dx2
        dpsi += 0.5 * (d2psi + d2psinew) * dx
        x += dx
        x_tab.append(x)
        psi_tab.append(psi)
    return x_tab, psi_tab

x_harmonic, psi_harmonic = verlet(E_harmonic, AnHo, dx, dx2, c)

# Normalize wavefunctions
def normalize(psi, dx):
    norm = np.sum(np.abs(psi)**2) * dx
    return psi / np.sqrt(norm)

psi_harmonic_normalized = normalize(psi_harmonic, dx)

# Plot wavefunctions and potentials
plt.figure(figsize=(10, 5))

plt.plot(x_harmonic, psi_harmonic_normalized, label="Harmonic oscillator ground state")
plt.plot(x, V(x), label="Harmonic potential", linestyle="--", color="gray")

plt.xlabel("x")
plt.ylabel("Ψ(x)")
plt.legend()
plt.grid(True)
plt.show()
