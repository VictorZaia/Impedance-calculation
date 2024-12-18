import numpy as np
import matplotlib.pyplot as plt

# Constants
c0 = 343  # Speed of sound in air (m/s)
rho0 = 1.21  # Air density (kg/m^3)
mu = 1.81e-5  # Dynamic viscosity Pa.s
nu = mu / rho0  # Kinematic viscosity m2/s

# Plate's dimensions
sigma = 0.05  # Porosity
e = 1.0e-3  # Thickness in meters
d = 1.0e-3  # Diameter in meters

# Cavity
L = 0.08  # Cavity length in meters

# Frequency range
frequencies = np.linspace(0.1, 5000, 100)

def Initiate_wave(frequencies):
    """
    Initiates wave properties: omega and k.
    """
    omega = 2 * np.pi * frequencies  # Angular frequency
    k = omega / c0  # Wave number
    return omega, k

def compute_impedance_plate(omega, k):
    """
    Computes the acoustic resistance and reactance for the plate.
    """
    # Acoustic resistance
    r_visc = np.sqrt(8 * nu * omega) / (c0 * sigma) * (1 + e / d)
    r_rad = 1 / (8 * sigma) * (k * d)**2
    r_tot = r_visc + r_rad

    # Acoustic reactance
    chi_mass = omega / (sigma * c0) * (e + (8 * d) / (3 * np.pi) * (1 - 0.7 * np.sqrt(sigma)))
    chi_visc = omega / (sigma * c0) * (np.sqrt(8 * nu / omega) * (1 + e / d))
    chi_tot = chi_mass + chi_visc
    
    return r_visc, r_rad, r_tot, chi_mass, chi_visc, chi_tot

def compute_impedance_cavity(L, k):
    """
    Computes the cavity impedance.
    """
    return -1 / np.tan(k * L)

# Compute wave properties
omega, k = Initiate_wave(frequencies)

# Compute impedance components for plate and cavity
r_visc, r_rad, r_tot, chi_mass, chi_visc, chi_tot = compute_impedance_plate(omega, k)
chi_cavity = compute_impedance_cavity(L, k)

# Plotting Resistance and Reactance
fig, axs = plt.subplots(2, 1, figsize=(10, 12))

# Plot acoustic resistance
axs[0].plot(frequencies, r_visc, color='blue', label='Viscous effects')
axs[0].plot(frequencies, r_rad, color='red', label='Radiation effects')
axs[0].plot(frequencies, r_tot, color='green', label='Total resistance')
axs[0].set_xlabel('Frequency (Hz)')
axs[0].set_ylabel('Acoustic Resistance (r)')
axs[0].legend()
axs[0].grid(True)

# Plot acoustic reactance
axs[1].plot(frequencies, chi_visc, color='blue', label='Viscous effects')
axs[1].plot(frequencies, chi_mass, color='red', label='Mass effects')
axs[1].plot(frequencies, chi_tot, color='green', label='Total reactance')
axs[1].set_xlabel('Frequency (Hz)')
axs[1].set_ylabel('Acoustic Reactance (χ)')
axs[1].legend()
axs[1].grid(True)

# Plot cavity impedance
plt.figure(figsize=(10, 6))
plt.plot(frequencies, chi_cavity, color='blue', label='Cavity impedance')
plt.plot(frequencies, chi_tot + chi_cavity, color='red', label='Cavity and plate impedance')
plt.xlabel('Frequency (Hz)')
# plt.xlim([0, 3.1])
plt.ylabel('Coefficient B')
plt.ylim([-20, 50])
plt.legend()
plt.grid(True)

# Non-linear contribution (Resistance)
def compute_resistance_guess_formule():
    """
    Returns the resistance from the guess formula.
    """
    B = (1 - sigma**2) / (sigma * c0)
    return B

def compute_resistance_Melling_formule(C):
    """
    Returns the resistance using Melling's formula.
    """
    B = 4 / (3 * np.pi) * (1 - sigma**2) / (sigma * c0 * C**2)
    return B

v0 = np.linspace(0, 5, 100)  # Acoustic velocity (m/s)
C_D = 0.7

plt.figure(figsize=(10, 6))
plt.plot(v0, compute_resistance_guess_formule() * v0, color='blue', label='Guess')
plt.plot(v0, compute_resistance_Melling_formule(C_D) * v0, color='red', label='Melling')
plt.xlabel('Acoustic velocity (m/s)')
plt.ylabel('Acoustic resistance')
plt.legend()
plt.grid(True)
plt.show()
