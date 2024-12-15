import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

""" Constants """

# Air's constants
c0 = 343  # Speed of sound in air (m/s)
rho0 = 1.21  # Air density (kg/m^3)
mu = 1.81 * 10**-5  # Dynamic viscosity Pa.s
nu = mu / rho0  # Kinematic viscosity m2/s

""" Functions to compute resistance and reactance """

def Initiate_wave(frequencies):
    omega = 2 * np.pi * frequencies  # Angular frequency
    k = omega / c0  # Wave number
    return omega, k

def compute_impedance_plate(omega, k, sigma, d, e):
    r_visc = np.sqrt(8 * nu * omega) / (c0 * sigma) * (1 + e / d)
    r_rad = 1 / (8 * sigma) * (k * d)**2
    r_tot = r_visc + r_rad

    chi_mass = omega / (sigma * c0) * (e + (8 * d) / (3 * np.pi) * (1 - 0.7 * np.sqrt(sigma)))
    chi_visc = omega / (sigma * c0) * (np.sqrt(8 * nu / omega) * (1 + e / d))
    chi_tot = chi_mass + chi_visc
    
    return r_tot, chi_tot

def compute_impedance_cavity(L, k):
    chi_cavity = - 1 / np.tan(k * L)
    return chi_cavity

def resistance_eq(r, omega, p_bar, r_tot, chi_tot, chi_cavity, sigma):
    B = (1 - sigma**2) / (sigma * c0)  # Nonlinear coefficient
    impedance_magnitude = abs(r + (1j * (chi_tot + chi_cavity)))
    return r - (r_tot + B * p_bar / (rho0 * c0 * sigma * impedance_magnitude))

""" General Function to Vary Parameters """

def compute_impedance_varying_param(frequencies, L, d, sigma, p_bar, e):
    # Detect the varying parameter
    parameters = {'L': L, 'd': d, 'sigma': sigma, 'p_bar': p_bar, 'e': e}
    varying_param = [k for k, v in parameters.items() if isinstance(v, (list, np.ndarray))]

    varying_param_name = varying_param[0]
    varying_param_values = parameters[varying_param_name]
    print(f"Varying parameter: {varying_param_name}")

    # Set fixed values for other parameters
    L_fixed, d_fixed, sigma_fixed, p_bar_fixed, e_fixed = L, d, sigma, p_bar, e
    omega, k = Initiate_wave(frequencies)

    # Store results
    results = []

    # Iterate over the varying parameter
    for value in varying_param_values:
        # Update the parameter based on which one is varying
        if varying_param_name == 'L':
            L_fixed = value
        elif varying_param_name == 'd':
            d_fixed = value
        elif varying_param_name == 'sigma':
            sigma_fixed = value
        elif varying_param_name == 'p_bar':
            p_bar_fixed = value
        elif varying_param_name == 'e':
            e_fixed = value

        # Compute impedance
        r_tot, chi_tot = compute_impedance_plate(omega, k, sigma_fixed, d_fixed, e_fixed)
        chi_cavity = compute_impedance_cavity(L_fixed, k)

        # Solve for nonlinear resistance
        r_nonlinear = []
        for i in range(len(frequencies)):
            omega_i = omega[i]
            chi_tot_i = chi_tot[i]
            chi_cavity_i = chi_cavity[i]
            r_tot_i = r_tot[i]
            
            r_initial = 0.5  # Initial guess for resistance
            r_solution = fsolve(resistance_eq, r_initial, args=(omega_i, p_bar_fixed, r_tot_i, chi_tot_i, chi_cavity_i, sigma_fixed))[0]
            r_nonlinear.append(r_solution)

        # Save the results
        results.append((value, r_nonlinear))

    # Plot the results
    plt.figure(figsize=(10, 6))
    for value, r_values in results:
        plt.plot(frequencies, r_values, label=f"{varying_param_name} = {value:.3e}")

    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Nonlinear Resistance $r$ (Pa·s/m)")
    plt.title(f"Nonlinear Resistance vs Frequency for Varying {varying_param_name}")
    plt.legend()
    plt.grid(True)
    plt.show()

""" Example Usage """

# Liner's geometry and properties
L = 10e-3
d = 0.5e-3
sigma = 0.015
p_bar = 1000
e = 1.0e-3 

# Define frequencies
# frequencies = np.linspace(0.1, 8000, 100)

# L_values = np.linspace(5e-3, 20e-3, 5)  # Vary L between 5mm and 20mm
# compute_impedance_varying_param(frequencies, L_values, d, sigma, p_bar, e)

# d_values = np.linspace(0.5e-3, 1.0e-3, 5)  # Vary d between 0.5mm and 1.0mm
# compute_impedance_varying_param(frequencies, L, d_values, sigma, p_bar, e)

# sigma_values = np.linspace(0.01, 0.02, 5)  # Vary sigma between 0.01 and 0.02
# compute_impedance_varying_param(frequencies, L, d, sigma_values, p_bar, e)

# p_bar_values = np.linspace(20, 1000, 5)  # Vary p_bar between 100Pa and 1000Pa
# compute_impedance_varying_param(frequencies, L, d, sigma, p_bar_values, e)

e_values = np.linspace(0.5e-3, 2e-3, 5)  # Vary e between 0.5mm and 2mm
compute_impedance_varying_param(frequencies, L, d, sigma, p_bar, e_values)


