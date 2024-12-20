import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

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
    r_tot = r_visc + r_rad

    eps = 1 / (1 + 305 * M**3) # Correction factor when considering airflow, M is the mach number
    chi_mass = omega / (sigma * c0) * (e + eps * (8 * d) / (3 * np.pi) * (1 - 0.71 * np.sqrt(sigma)))
    chi_visc = omega / (sigma * c0) * (np.sqrt(8 * nu / omega) * (1 + e / d))
    chi_tot = chi_mass + chi_visc
    
    return r_tot, chi_tot

def compute_impedance_cavity(L, k):
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


# %%
""" Function of the resistance using the acoustic pressure rather than acoustic velocity """

def resistance_eq(r, omega, p_pa, r_tot, chi_tot, chi_cavity, sigma, M):
    """ 
    The resistance can be given by a linear part, called A(w) and a non linear part which depends on the acoustic velocity |v|, multiplied by a coefficient B.
    r(w,|v|) = A(w) + B |v|/sigma + (1 -sigma²)/sigma kM
    |v| can be linked to |p| yielding the equation
    see Malmary if necessary
    """
    B = (1 - sigma**2) / (sigma * c0)
    impedance_magnitude = abs(r + (1j * (chi_tot + chi_cavity)))
    r_airflow = compute_resistance_tangencial_airflow(sigma, M)
    return r - (r_tot + B * p_pa / (rho0 * c0 * sigma * impedance_magnitude) + r_airflow)


# %%
""" General Function to Vary Parameters """

def compute_impedance_varying_param(frequencies, L, d, sigma, p_pa, e, M):
    """ 
    This function verifies which parameter is a list and solve the system for differnet values.
    """
    parameters = {'L': L, 'd': d, 'sigma': sigma, 'p_pa': p_pa, 'e': e, 'M':M}

    for k, v in parameters.items():
        if isinstance(v, np.ndarray):
            varying_param = k

    varying_param_values = parameters[varying_param]
    print(f"Varying parameter: {varying_param}")

    L_fixed, d_fixed, sigma_fixed, p_pa_fixed, e_fixed, M_fixed = L, d, sigma, p_pa, e, M
    omega, k = Initiate_wave(frequencies)

    results = []

    for value in varying_param_values:
        # Check the name of the varying parameter and updates the value
        if varying_param == 'L':
            L_fixed = value
        elif varying_param == 'd':
            d_fixed = value
        elif varying_param == 'sigma':
            sigma_fixed = value
        elif varying_param == 'p_pa':
            p_pa_fixed = value
        elif varying_param == 'e':
            e_fixed = value
        elif varying_param == 'M':
            M_fixed = value

        r_tot, chi_tot = compute_impedance_plate(omega, k, sigma_fixed, d_fixed, e_fixed, M_fixed)
        chi_cavity = compute_impedance_cavity(L_fixed, k)

        r_nonlinear = []
        for i in range(len(frequencies)):
            omega_i = omega[i]
            chi_tot_i = chi_tot[i]
            chi_cavity_i = chi_cavity[i]
            r_tot_i = r_tot[i]
            
            r_initial = 0.5
            r_solution = fsolve(resistance_eq, r_initial, args=(omega_i, p_pa_fixed, r_tot_i, chi_tot_i, chi_cavity_i, sigma_fixed, M_fixed))[0]
            r_nonlinear.append(r_solution)

        # Save the results
        results.append((value, r_nonlinear))

    # Plot the results
    plt.figure(figsize=(10, 6))
    for value, r_values in results:
        plt.plot(frequencies, r_values, label=f"{varying_param} = {value:.3e}")

    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Nonlinear Resistance $r$ (Pa·s/m)")
    plt.title(f"Nonlinear Resistance vs Frequency for Varying {varying_param}")
    plt.legend()
    plt.grid(True)
    plt.show()


# %%
""" Example Usage """

# Liner's geometry and properties
L = 10e-3
d = 0.5e-3
sigma = 0.015
p_pa = 1000
e = 1.0e-3
M = 0.1

# Define frequencies
frequencies = np.linspace(0.1, 8000, 100)

L_values = np.linspace(5e-3, 20e-3, 5)  # Vary L between 5mm and 20mm
compute_impedance_varying_param(frequencies, L_values, d, sigma, p_pa, e, M)

# d_values = np.linspace(0.5e-3, 1.0e-3, 5)  # Vary d between 0.5mm and 1.0mm
# compute_impedance_varying_param(frequencies, L, d_values, sigma, p_pa, e, M)

# sigma_values = np.linspace(0.01, 0.02, 5)  # Vary sigma between 0.01 and 0.02
# compute_impedance_varying_param(frequencies, L, d, sigma_values, p_pa, e, M)

# p_pa_values = np.linspace(20, 1000, 5)  # Vary p_pa between 100Pa and 1000Pa
# compute_impedance_varying_param(frequencies, L, d, sigma, p_pa_values, e, M)

#e_values = np.linspace(0.5e-3, 2e-3, 5)  # Vary e between 0.5mm and 2mm
#compute_impedance_varying_param(frequencies, L, d, sigma, p_pa, e_values, M)

# M_values = np.linspace(0.01, 0.5, 4)  # Vary e between 0.5mm and 2mm
# compute_impedance_varying_param(frequencies, L, d, sigma, p_pa, e, M_values)
