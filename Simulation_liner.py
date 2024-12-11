import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


""" Constants """

# Air's constants
c0 = 343  # Speed of sound in air (m/s)
rho0 = 1.21  # Air density (kg/m^3)
mu = 1.81 * 10**-5  # Dynamic viscosity Pa.s
nu = mu / rho0  # Kinematic viscosity m2/s

# Plate's dimensions
sigma = 0.0139  # Porosity
e = 1.02e-3  # Thickness in meters (converted from mm)
d = 0.68e-3  # Diameter in meters (converted from mm)

# Cavity
L = 10e-3  # Cavity length in meters

# Frequency range
frequencies = np.linspace(0.1, 8000, 100)


""" Function to compute resistance and reactance """

def Initiate_wave(frequencies):
    omega = 2 * np.pi * frequencies  # Angular frequency
    k = omega / c0  # Wave number
    return omega, k

def compute_impedance_plate(omega, k):
    # Acoustic resistance (r)
    r_visc = np.sqrt(8 * nu * omega) / (c0 * sigma) * (1 + e / d)
    r_rad = 1 / (8 * sigma) * (k * d)**2
    r_tot = r_visc + r_rad

    # Acoustic reactance (χ)
    chi_mass = omega / (sigma * c0) * (e + (8 * d) / (3 * np.pi) * (1 - 0.7 * np.sqrt(sigma)))
    chi_visc = omega / (sigma * c0) * (np.sqrt(8 * nu / omega) * (1 + e / d))
    chi_tot = chi_mass + chi_visc
    
    return r_visc, r_rad, r_tot, chi_mass, chi_visc, chi_tot

def compute_impedance_cavity(L, k):
    chi_cavity = - 1 / np.tan(k * L)
    return chi_cavity


""" Compute resistance and reactance """

omega, k = Initiate_wave(frequencies)
r_visc, r_rad, r_tot, chi_mass, chi_visc, chi_tot = compute_impedance_plate(omega, k)
chi_cavity = compute_impedance_cavity(L, k)

# Nonlinear coefficient
B = (1 - sigma**2) / (sigma * c0)  # Example value (adjust based on system)

# Define the equation to solve for r
def resistance_eq(r, omega, p_bar, r_tot, chi_tot, chi_cavity):
    a = r_tot  # Use the total linear resistance A(omega) already calculated
    impedance_magnitude = abs(r + (1j * (chi_tot + chi_cavity)))  # Correct imaginary handling
    return r - (a + B * p_bar / (rho0 * c0 * sigma * impedance_magnitude))

# Define pressure amplitude values (Pa)
p_bar_values = [2, 63, 200, 1125, 2000]  # Example values in Pascals

# Solve r for each frequency and pressure amplitude
r_nonlinear = {}

for p_bar in p_bar_values:
    r_values = []
    for i in range(len(frequencies)):
        omega_i = omega[i]
        chi_tot_i = chi_tot[i]
        chi_cavity_i = chi_cavity[i]
        r_tot_i = r_tot[i]
        
        # Initial guess for r
        r_initial = 0.5
        
        # Solve numerically for r using fsolve
        r_solution = fsolve(resistance_eq, r_initial, args=(omega_i, p_bar, r_tot_i, chi_tot_i, chi_cavity_i))[0]
        r_values.append(r_solution)
    
    r_nonlinear[p_bar] = r_values

# Plot the nonlinear resistance for different pressure amplitudes
plt.figure(figsize=(10, 6))
for p_bar, r_values in r_nonlinear.items():
    plt.plot(frequencies, r_values, label=f"$|\\bar{{p}}| = {p_bar}$ Pa")

plt.xlabel("Frequency (Hz)")
plt.ylabel("Nonlinear Resistance $r$ (Pa·s/m)")
plt.title("Nonlinear Resistance vs Frequency for Different $|\\bar{p}|$")
plt.legend()
plt.grid(True)
plt.show()
