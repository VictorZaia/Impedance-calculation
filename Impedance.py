import numpy as np
from scipy.optimize import fsolve

class Impedance:
    """Class to calculate the impedance (resistance and reactance) for a given wave and environment."""
    def __init__(self, Environment, Wave, Liner, p_acous_pa, M):

        self._resistance = None
        self._reactance = None
        self._absorption_coefficient = None

        self._Environment = Environment
        self._Wave = Wave
        self._Liner = Liner
        self._p_acous_pa = p_acous_pa
        self._M = M

    def set_resistance(self, resistance):
        self._resistance = resistance

    def set_reactance(self, reactance):
        self._reactance = reactance
    
    def set_absorption_coefficient(self, absorption_coefficient):
        self._absorption_coefficient = absorption_coefficient

    
    def compute_impedance_plate(self, omega, k, sigma, d, e, M):
        """
        Computes the acoustic resistance and reactance for the plate.
        """
        r_visc = np.sqrt(8 * self._Environment.nu * omega) / (self._Environment.c * sigma) * (1 + e / d)
        r_rad = 1 / (8 * sigma) * (k * d)**2
        r_tot_plate = r_visc + r_rad

        eps = 1 / (1 + 305 * M**3) # Correction factor when considering airflow, M is the mach number
        chi_mass = omega / (sigma * self._Environment.c) * (e + eps * (8 * d) / (3 * np.pi) * (1 - 0.71 * np.sqrt(sigma)))
        chi_visc = omega / (sigma * self._Environment.c) * (np.sqrt(8 * self._Environment.nu / omega) * (1 + e / d))
        chi_tot_plate = chi_mass + chi_visc
        
        return r_tot_plate, chi_tot_plate
    
    @staticmethod
    def compute_reactance_cavity(L, k):
        """
        Computes the cavity impedance.
        """
        chi_cavity = - 1 / np.tan(k * L)

        return chi_cavity
    
    @staticmethod
    def compute_resistance_tangencial_airflow(sigma, M):
        """
        Computes the resistance term for a tangencial airflow
        M: mach number 
        """
        r_airflow = 0.3 * (1 - sigma**2) / sigma * M

        return r_airflow
    
    @staticmethod
    def calculate_absorption_coefficient(resistance, reactance, theta=0):
        """Calculate the absorption coefficient."""
        Z_surface = resistance + 1j * reactance
        R_surface = (Z_surface * np.cos(theta) - 1) / (Z_surface * np.cos(theta) + 1)
        return 1 - np.abs(R_surface)**2

    def resistance_eq(self, r, p_acous_pa, r_tot_plate, chi_tot_plate, chi_cavity, sigma, M):
        """ 
        The resistance can be given by a linear part, called A(w) and a non linear part which depends on the acoustic velocity |v|, multiplied by a coefficient B.
        r(w,|v|) = A(w) + B |v|/sigma + (1 -sigma²)/sigma kM
        |v| can be linked to |p| yielding the equation
        see Malmary if necessary
        """
        B = (1 - sigma**2) / (sigma * self._Environment.c)
        impedance_magnitude = abs(r + (1j * (chi_tot_plate + chi_cavity)))
        r_airflow = Impedance.compute_resistance_tangencial_airflow(sigma, M)

        return r - (r_tot_plate + B * p_acous_pa / (self._Environment.rho * self._Environment.c * sigma * impedance_magnitude) + r_airflow)

    def compute_impedance_for_frequencies(self):
        """Compute the impedance for varying frequencies."""

        L = self._Liner._L
        d = self._Liner._d
        sigma = self._Liner._sigma
        e = self._Liner._e

        omega = self._Wave.omega
        K = self._Wave.K

        r_tot_plate, chi_tot_plate = self.compute_impedance_plate(omega, K, sigma, d, e, self._M)
        chi_cavity = Impedance.compute_reactance_cavity(L, K)

        r_nonlinear = []
        alpha = []

        for i in range(len(omega)):
            chi_tot_plate_i = chi_tot_plate[i]
            chi_cavity_i = chi_cavity[i]
            r_tot_plate_i = r_tot_plate[i]

            r_initial = 0.5
            r_solution = fsolve(self.resistance_eq, r_initial, args=(self._p_acous_pa, r_tot_plate_i, chi_tot_plate_i, chi_cavity_i, sigma, self._M))[0]

            alpha.append(Impedance.calculate_absorption_coefficient(r_solution, chi_tot_plate_i + chi_cavity_i))
            r_nonlinear.append(r_solution)

        self.set_resistance(r_nonlinear)
        self.set_reactance(chi_cavity + chi_tot_plate)
        self.set_absorption_coefficient(alpha)
