import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import matplotlib.cm as cm

# %%
"""Constants """

# Air's constants
c0 = 343  # Speed of sound in air (m/s)
rho0 = 1.21  # Air density (kg/m^3) 
mu = 1.81 * 10**-5  # Dynamic viscosity Pa.s
nu = mu / rho0  # Kinematic viscosity m2/s

# %%
""" Functions to compute resistance and reactance """

def Initiate_wave(frequencies):
    """
    Initiates wave properties: omega and k.
    """
    omega = 2 * np.pi * frequencies  # Angular frequency
    k = omega / c0  # Wave number
    return omega, k

def compute_impedance_plate(omega, k, sigma, d, e, M):
    """
    Computes the acoustic resistance and reactance for the plate.
    """
    r_visc = np.sqrt(8 * nu * omega) / (c0 * sigma) * (1 + e / d)
    r_rad = 1 / (8 * sigma) * (k * d)**2
    r_tot_plate = r_visc + r_rad

    eps = 1 / (1 + 305 * M**3) # Correction factor when considering airflow, M is the mach number
    chi_mass = omega / (sigma * c0) * (e + eps * (8 * d) / (3 * np.pi) * (1 - 0.71 * np.sqrt(sigma)))
    chi_visc = omega / (sigma * c0) * (np.sqrt(8 * nu / omega) * (1 + e / d))
    chi_tot_plate = chi_mass + chi_visc
    
    return r_tot_plate, chi_tot_plate

def compute_reactance_cavity(L, k):
    """
    Computes the cavity impedance.
    """
    chi_cavity = - 1 / np.tan(k * L)
    return chi_cavity

def compute_resistance_tangencial_airflow(sigma, M):
    """
    Computes the resistance term for a tangencial airflow
    M: mach number 
    """
    r_airflow = 0.3 * (1 - sigma**2) / sigma * M
    return r_airflow

def calculate_absorption_coefficient(resistance, reactance, theta=0):
    
    Z_surface = resistance + 1j * reactance

    R_surface = (Z_surface * np.cos(theta) - 1)/(Z_surface * np.cos(theta) + 1)
    
    return 1 - np.abs(R_surface)**2

# %%
""" Function of the resistance using the acoustic pressure rather than acoustic velocity """

def resistance_eq(r, omega, p_acous_pa, r_tot_plate, chi_tot_plate, chi_cavity, sigma, M):
    """ 
    The resistance can be given by a linear part, called A(w) and a non linear part which depends on the acoustic velocity |v|, multiplied by a coefficient B.
    r(w,|v|) = A(w) + B |v|/sigma + (1 -sigma²)/sigma kM
    |v| can be linked to |p| yielding the equation
    see Malmary if necessary
    """
    B = (1 - sigma**2) / (sigma * c0)
    impedance_magnitude = abs(r + (1j * (chi_tot_plate + chi_cavity)))
    r_airflow = compute_resistance_tangencial_airflow(sigma, M)
    return r - (r_tot_plate + B * p_acous_pa / (rho0 * c0 * sigma * impedance_magnitude) + r_airflow)

# %%
""" General Function to Vary Parameters """

def compute_impedance_varying_param(frequencies, L, d, sigma, p_acous_pa, e, M):
    """ 
    This function verifies which parameter is a list and solve the system for differnet values.
    """
    parameters = {'L': L, 'd': d, 'sigma': sigma, 'p_acous_pa': p_acous_pa, 'e': e, 'M':M}

    for k, v in parameters.items():
        if isinstance(v, np.ndarray):
            varying_param = k

    varying_param_values = parameters[varying_param]
    print(f"Varying parameter: {varying_param}")

    L_fixed, d_fixed, sigma_fixed, p_acous_pa_fixed, e_fixed, M_fixed = L, d, sigma, p_acous_pa, e, M
    omega, k = Initiate_wave(frequencies)

    results_resistance = {}
    results_alpha = {}

    for value in varying_param_values:
        # Check the name of the varying parameter and updates the value
        if varying_param == 'L':
            L_fixed = value
        elif varying_param == 'd':
            d_fixed = value
        elif varying_param == 'sigma':
            sigma_fixed = value
        elif varying_param == 'p_acous_pa':
            p_acous_pa_fixed = value
        elif varying_param == 'e':
            e_fixed = value
        elif varying_param == 'M':
            M_fixed = value

        r_tot_plate, chi_tot_plate = compute_impedance_plate(omega, k, sigma_fixed, d_fixed, e_fixed, M_fixed)
        chi_cavity = compute_reactance_cavity(L_fixed, k)

        r_nonlinear = []
        alpha = []

        for i in range(len(frequencies)):
            omega_i = omega[i]
            chi_tot_plate_i = chi_tot_plate[i]
            chi_cavity_i = chi_cavity[i]
            r_tot_plate_i = r_tot_plate[i]
            
            r_initial = 0.5
            r_solution = fsolve(resistance_eq, r_initial, args=(omega_i, p_acous_pa_fixed, r_tot_plate_i, chi_tot_plate_i, chi_cavity_i, sigma_fixed, M_fixed))[0]
            
            alpha.append(calculate_absorption_coefficient(r_solution, chi_tot_plate_i + chi_cavity_i))
            r_nonlinear.append(r_solution)

        # Save the results
        results_resistance[value] = r_nonlinear
        results_alpha[value] = alpha

    return results_resistance, results_alpha

# %%
""" Example Usage """

# Liner's geometry and properties
L = 15e-3
d = 1.5e-3
sigma = 0.15
e = 1.5e-3

p_acous_pa = 1000


# Define frequencies
frequencies = np.linspace(0.1, 5000, 100)

M_values = np.linspace(0, 0.5, 5)
L_values = np.linspace(10e-3, 20e-3, 5)

cmap = cm.get_cmap('viridis', M_values.size)

plt.figure(figsize=(10, 6))
for i in range(M_values.size):
    results_resistance, results_alpha = compute_impedance_varying_param(frequencies, L_values, d, sigma, p_acous_pa, e, M_values[i])

    color = cmap(i)

    value = list(results_alpha.keys())
    r_values = list(results_alpha.values())
    freq_max = []
    for j in range(len(value)):
        freq_max.append(frequencies[np.argmax(r_values[j])])
    plt.plot(value, freq_max, color=color, label=f"{M_values[i]}")
    plt.xlabel("cavity size")
    plt.ylabel("Frequency (Hz)")
    plt.legend()
    plt.grid(True)
plt.show()



















































"""
plt.figure(figsize=(10, 6))
for value, r_values in results_resistance:
    plt.plot(frequencies, r_values, label=f"{varying_param} = {20 * np.log10(value / 20e-6):.0f} dB")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Acoustic Resistance (r)")
plt.title(f"Nonlinear Resistance vs Frequency for Varying {varying_param}")
plt.legend()
plt.grid(True)

plt.figure(figsize=(10, 6))
for value, r_values in results_alpha:
    plt.plot(frequencies, chi_cavity + chi_tot_plate, label=f"{varying_param} = {20 * np.log10(value / 20e-6):.0f} dB")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Acoustic Reactance (χ)")
plt.ylim([-20, 50])
plt.legend()
plt.grid(True)

plt.figure(figsize=(10, 6))
for value, r_values in results_alpha:
    plt.plot(frequencies, r_values, label=f"{varying_param} = {20 * np.log10(value / 20e-6):.0f} dB")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Absorption coefficient")
plt.legend()
plt.grid(True)
plt.show()
"""