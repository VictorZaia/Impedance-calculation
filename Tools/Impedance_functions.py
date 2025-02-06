import numpy as np

def compute_resistance_plate(omega, k, sigma, d, e, nu, speed_of_sound):
    """
    Computes the acoustic resistance for the plate.
    """
    r_visc = np.sqrt(8 * nu * omega) / (speed_of_sound * sigma) * (1 + e / d)
    r_rad = 1 / (8 * sigma) * (k * d)**2
    r_tot_plate = r_visc + r_rad
    
    return r_tot_plate

def compute_reactance_plate(omega, sigma, d, e, nu, speed_of_sound, M):
    """
    Computes the acoustic reactance for the plate.
    """
    eps = 1 / (1 + 305 * M**3) # Correction factor when considering airflow, M is the mach number
    chi_mass = omega / (sigma * speed_of_sound) * (e + eps * (8 * d) / (3 * np.pi) * (1 - 0.71 * np.sqrt(sigma)))
    chi_visc = omega / (sigma * speed_of_sound) * (np.sqrt(8 * nu / omega) * (1 + e / d))
    chi_tot_plate = chi_mass + chi_visc
    
    return chi_tot_plate

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

def resistance_eq(r, p_acous_pa, r_tot_plate, chi_tot_plate, chi_cavity, sigma, rho, speed_of_sound, M):
    """ 
    The resistance can be given by a linear part, called A(w) and a non linear part which depends on the acoustic velocity |v|, multiplied by a coefficient B.
    r(w,|v|) = A(w) + B |v|/sigma + (1 -sigma²)/sigma M
    |v| can be linked to |p| yielding the equation
    see Malmary if necessary
    """
    B = (1 - sigma**2) / (sigma * speed_of_sound)
    impedance_magnitude = abs(r + (1j * (chi_tot_plate + chi_cavity)))
    r_airflow = compute_resistance_tangencial_airflow(sigma, M)

    return r - (r_tot_plate + B * p_acous_pa / (rho *speed_of_sound * sigma * impedance_magnitude) + r_airflow)

def calculate_absorption_coefficient(resistance, reactance, theta=0):
    """
    Calculate the absorption coefficient.
    """
    Z_surface = resistance + 1j * reactance
    R_surface = (Z_surface * np.cos(theta) - 1) / (Z_surface * np.cos(theta) + 1)

    return 1 - np.abs(R_surface)**2