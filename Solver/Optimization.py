from Properties.Wave import Wave
from Properties.Flying_condition import Flying_condition
from Solver.Processor import Processor

from scipy.optimize import fsolve, minimize
import numpy as np
import matplotlib.pyplot as plt


class Optimizer(Processor):
     
    def __init__(self, environment, wave):
          super().__init__(environment, wave, "To be optimized")
    
    @staticmethod
    def compute_resistance_plate(omega, k, sigma, d, e, nu, speed_of_sound):
        """
        Computes the acoustic resistance for the plate.
        """
        r_visc = np.sqrt(8 * nu * omega) / (speed_of_sound * sigma) * (1 + e / d)
        r_rad = 1 / (8 * sigma) * (k * d)**2
        r_tot_plate = r_visc + r_rad
        
        return r_tot_plate
    
    @staticmethod
    def compute_reactance_plate(omega, sigma, d, e, nu, speed_of_sound, M):
        """
        Computes the acoustic reactance for the plate.
        """
        eps = 1 / (1 + 305 * M**3) # Correction factor when considering airflow, M is the mach number
        chi_mass = omega / (sigma * speed_of_sound) * (e + eps * (8 * d) / (3 * np.pi) * (1 - 0.71 * np.sqrt(sigma)))
        chi_visc = omega / (sigma * speed_of_sound) * (np.sqrt(8 * nu / omega) * (1 + e / d))
        chi_tot_plate = chi_mass + chi_visc
        
        return chi_tot_plate
    
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
    def resistance_eq(r, p_acous_pa, r_tot_plate, chi_tot_plate, chi_cavity, sigma, rho, speed_of_sound, M):
        """ 
        The resistance can be given by a linear part, called A(w) and a non linear part which depends on the acoustic velocity |v|, multiplied by a coefficient B.
        r(w,|v|) = A(w) + B |v|/sigma + (1 -sigma²)/sigma kM
        |v| can be linked to |p| yielding the equation
        see Malmary if necessary
        """
        B = (1 - sigma**2) / (sigma * speed_of_sound)
        impedance_magnitude = abs(r + (1j * (chi_tot_plate + chi_cavity)))
        r_airflow = Processor.compute_resistance_tangencial_airflow(sigma, M)

        return r - (r_tot_plate + B * p_acous_pa / (rho *speed_of_sound * sigma * impedance_magnitude) + r_airflow)
    
    @staticmethod
    def calculate_absorption_coefficient(resistance, reactance, theta=0):
        """
        Calculate the absorption coefficient.
        """
        Z_surface = resistance + 1j * reactance
        R_surface = (Z_surface * np.cos(theta) - 1) / (Z_surface * np.cos(theta) + 1)

        return 1 - np.abs(R_surface)**2
    
    
    @staticmethod
    def resistance_eq_opt(L, d, sigma, e, p_acous_pa, omega, K, nu, rho, speed_of_sound, M):
        """ 
        Resistance equation used for computing the resistance based on L. This is a pure function
        that doesn't depend on instance variables and can be used with JAX.
        """
        # Compute the resistance and reactance at the fixed frequency
        r_tot_plate = Optimizer.compute_resistance_plate(omega, K, sigma, d, e, nu, speed_of_sound)
        chi_tot_plate = Optimizer.compute_reactance_plate(omega, sigma, d, e, nu, speed_of_sound, M)
        chi_cavity = Optimizer.compute_reactance_cavity(L, K)
        
        # Solve the resistance equation for the fixed frequency
        r_initial = np.array(0.5)
        r_solution = fsolve(Optimizer.resistance_eq, r_initial, args=(p_acous_pa, r_tot_plate, chi_tot_plate, chi_cavity, sigma, rho, speed_of_sound, M))[0]
        
        alpha = Optimizer.calculate_absorption_coefficient(r_solution, chi_tot_plate + chi_cavity)
        return -alpha


    
    @staticmethod
    def find_L_with_gradient(d, sigma, e, p_acous_pa, omega, K, nu, rho, speed_of_sound, M):
        """
        Find the cavity height (L) that minimizes the resistance using gradient descent.
        """
        def resistance_for_L(L):
            
            return Optimizer.resistance_eq_opt(L, d, sigma, e, p_acous_pa, omega, K, nu, rho, speed_of_sound, M)

        # Initialize with an initial guess for L
        L_initial = 0.1
        
        result = minimize(resistance_for_L, L_initial, bounds=[(0.01  , 0.5)])

        optimal_L = result.x[0]  
        optimal_absorption = -result.fun  
        
        return optimal_L, optimal_absorption
    
    """
    @staticmethod
    def solve_opt(d, sigma, e, altitude, frequencie, p_acous_pa):
        flying_condition = Flying_condition(altitude)
        flying_condition.initialize_flying_condition()

        wave = Wave(frequencie, flying_condition._environment.speed_of_sound)

        optimizer = Optimizer(flying_condition._environment, wave)

        optimal_L, optimal_absorption = optimizer.find_L_with_gradient(d, sigma, e, p_acous_pa, optimizer._wave.omega, optimizer._wave.k, optimizer._environment.nu, optimizer._environment.rho, optimizer._environment.speed_of_sound, flying_condition._mach)

        return optimal_L, optimal_absorption
    """

    @staticmethod
    def plot_loss_function(L, d, sigma, e, altitude, frequencie, p_acous_pa):

        flying_condition = Flying_condition(altitude)
        flying_condition.initialize_flying_condition()

        wave = Wave(frequencie, flying_condition._environment.speed_of_sound)
        
        alpha = np.zeros((len(sigma), len(L)))

        for i, value_sigma in enumerate(sigma):
            for j, value_L in enumerate(L):
                alpha[i, j] = - Optimizer.resistance_eq_opt(value_L, d, value_sigma, e, p_acous_pa, wave.omega, wave.k, flying_condition._environment.nu, flying_condition._environment.rho, flying_condition._environment.speed_of_sound, flying_condition._mach)

        plt.figure(figsize=(10, 6))
        plt.pcolormesh(sigma*100,L*10**3, alpha, cmap='rainbow')
        plt.colorbar(label="Absorption Coefficient")
        plt.xlabel("sigma [%]")
        plt.ylabel("L [mm]")
        plt.grid(True)
        plt.show()

    @staticmethod
    def plot_3d_loss_function(L_values, d, sigma_values, e, altitude_values, frequencie, p_acous_pa):
        optimal_L = []
        optimal_sigma = []
        alpha_max = []

        for altitude in altitude_values:

            flying_condition = Flying_condition(altitude)
            flying_condition.initialize_flying_condition()

            wave = Wave(frequencie, flying_condition._environment.speed_of_sound)
            
            alpha = np.zeros((len(sigma_values), len(L_values)))

            for i, value_sigma in enumerate(sigma_values):
                for j, value_L in enumerate(L_values):
                    alpha[i, j] = - Optimizer.resistance_eq_opt(value_L, d, value_sigma, e, p_acous_pa, wave.omega, wave.k, flying_condition._environment.nu, flying_condition._environment.rho, flying_condition._environment.speed_of_sound, flying_condition._mach)
            alpha_max.append(np.max(alpha))
            alpha_max_index = np.where(alpha == np.max(alpha))
            optimal_L.append(L_values[alpha_max_index[1][0]])
            optimal_sigma.append(sigma_values[alpha_max_index[0][0]])
        
        plt.figure(figsize=(10, 6))

        # Plot optimal L vs altitude
        plt.subplot(2, 1, 1)
        plt.plot(altitude_values, optimal_L, label='Optimal L', color='b')
        plt.xlabel("Altitude (m)")
        plt.ylabel("Optimal L (m)")
        plt.title("Optimal Cavity Height (L) vs Altitude")
        plt.grid(True)

        # Plot optimal sigma vs altitude
        plt.subplot(2, 1, 2)
        plt.plot(altitude_values, optimal_sigma, label='Optimal Sigma', color='r')
        plt.xlabel("Altitude (m)")
        plt.ylabel("Optimal Sigma")
        plt.title("Optimal Sigma vs Altitude")
        plt.grid(True)

        plt.tight_layout()
        plt.show()

