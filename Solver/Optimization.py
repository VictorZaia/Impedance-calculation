from Tools.Progress_bar import *
from Tools.Impedance_functions import *

from Properties.Wave import Wave
from Properties.Flying_condition import Flying_condition

from scipy.optimize import fsolve
import numpy as np

class Optimizer():
     
    def __init__(self, environment, wave):
          super().__init__(environment, wave, "To be optimized")
        
    @staticmethod
    def resistance_eq_opt(L, d, sigma, e, p_acous_pa, omega, K, nu, rho, speed_of_sound, M):
        """ 
        Resistance equation used for computing the resistance based on L.
        """
        # Compute the resistance and reactance at the fixed frequency
        r_tot_plate = compute_resistance_plate(omega, K, sigma, d, e, nu, speed_of_sound)
        chi_tot_plate = compute_reactance_plate(omega, sigma, d, e, nu, speed_of_sound, M)
        chi_cavity = compute_reactance_cavity(L, K)

        # Solve the resistance equation for the fixed frequency
        r_initial = np.array(0.5)
        r_solution = fsolve(resistance_eq, r_initial, args=(p_acous_pa, r_tot_plate, chi_tot_plate, chi_cavity, sigma, rho, speed_of_sound, M))[0]
        
        alpha = calculate_absorption_coefficient(r_solution, chi_tot_plate + chi_cavity)
        return alpha

    @staticmethod
    def loss_function_values(L, d, sigma, e, altitude, frequencie, p_acous_pa):

        flying_condition = Flying_condition(altitude)
        flying_condition.initialize_flying_condition()

        wave = Wave(frequencie, flying_condition._environment.speed_of_sound)
        
        alpha = np.zeros((len(sigma), len(L)))

        for i, value_sigma in enumerate(sigma):
            for j, value_L in enumerate(L):
                alpha[i, j] = Optimizer.resistance_eq_opt(value_L, d, value_sigma, e, p_acous_pa, wave.omega, wave.k, flying_condition._environment.nu, flying_condition._environment.rho, flying_condition._environment.speed_of_sound, flying_condition._mach)
        
        return alpha

    @staticmethod
    def optimum_geometry(L_values, d, sigma_values, e, altitude_values, frequencie, p_acous_pa):

        optimal_L = []
        optimal_sigma = []
        alpha_max = []

        total_iterations = len(altitude_values)

        for idx, altitude in enumerate(altitude_values):

            flying_condition = Flying_condition(altitude)
            flying_condition.initialize_flying_condition()

            wave = Wave(frequencie, flying_condition._environment.speed_of_sound)
            
            alpha = np.zeros((len(sigma_values), len(L_values)))

            for i, value_sigma in enumerate(sigma_values):
                for j, value_L in enumerate(L_values):
                    alpha[i, j] = Optimizer.resistance_eq_opt(value_L, d, value_sigma, e, p_acous_pa, wave.omega, wave.k, flying_condition._environment.nu, flying_condition._environment.rho, flying_condition._environment.speed_of_sound, flying_condition._mach)
            alpha_max.append(np.max(alpha))
            alpha_max_index = np.where(alpha == np.max(alpha))
            optimal_L.append(L_values[alpha_max_index[1][0]])
            optimal_sigma.append(sigma_values[alpha_max_index[0][0]])

            progress_bar(idx + 1, total_iterations)

        return alpha_max, optimal_L, optimal_sigma

        