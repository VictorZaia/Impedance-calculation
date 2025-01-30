from Properties.Environment import *
from Properties.Wave import *
from Properties.Flying_condition import *
from Liner.Liner import *
from Solver.Impedance import *

import numpy as np
from scipy.optimize import fsolve
import copy

class Processor:

    __impedance = Impedance()

    def __init__(self, environment, wave, liner):
        self._environment = environment
        self._wave = wave
        self._liner = liner
    
    def compute_resistance_plate(self, omega, k, sigma, d, e):
        """
        Computes the acoustic resistance for the plate.
        """
        r_visc = np.sqrt(8 * self._environment.nu * omega) / (self._environment.speed_of_sound * sigma) * (1 + e / d)
        r_rad = 1 / (8 * sigma) * (k * d)**2
        r_tot_plate = r_visc + r_rad
        
        return r_tot_plate
    
    def compute_reactance_plate(self, omega, sigma, d, e, M):
        """
        Computes the acoustic reactance for the plate.
        """
        eps = 1 / (1 + 305 * M**3) # Correction factor when considering airflow, M is the mach number
        chi_mass = omega / (sigma * self._environment.speed_of_sound) * (e + eps * (8 * d) / (3 * np.pi) * (1 - 0.71 * np.sqrt(sigma)))
        chi_visc = omega / (sigma * self._environment.speed_of_sound) * (np.sqrt(8 * self._environment.nu / omega) * (1 + e / d))
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
    
    def resistance_eq(self, r, p_acous_pa, r_tot_plate, chi_tot_plate, chi_cavity, sigma, M):
        """ 
        The resistance can be given by a linear part, called A(w) and a non linear part which depends on the acoustic velocity |v|, multiplied by a coefficient B.
        r(w,|v|) = A(w) + B |v|/sigma + (1 -sigma²)/sigma kM
        |v| can be linked to |p| yielding the equation
        see Malmary if necessary
        """
        B = (1 - sigma**2) / (sigma * self._environment.speed_of_sound)
        impedance_magnitude = abs(r + (1j * (chi_tot_plate + chi_cavity)))
        r_airflow = Processor.compute_resistance_tangencial_airflow(sigma, M)

        return r - (r_tot_plate + B * p_acous_pa / (self._environment.rho * self._environment.speed_of_sound * sigma * impedance_magnitude) + r_airflow)
    
    @staticmethod
    def calculate_absorption_coefficient(resistance, reactance, theta=0):
        """
        Calculate the absorption coefficient.
        """
        Z_surface = resistance + 1j * reactance
        R_surface = (Z_surface * np.cos(theta) - 1) / (Z_surface * np.cos(theta) + 1)

        return 1 - np.abs(R_surface)**2

    @staticmethod
    def compute_impedance(environment, wave, liner, p_acous_pa, M):
        """
        Compute the impedance for varying frequencies.
        """
        processor = Processor(environment, wave, liner)

        L = processor._liner._L
        d = processor._liner._d
        sigma = processor._liner._sigma
        e = processor._liner._e

        omega = processor._wave.omega
        K = processor._wave.k

        r_tot_plate = processor.compute_resistance_plate(omega, K, sigma, d, e)
        chi_tot_plate = processor.compute_reactance_plate(omega, sigma, d, e, M)

        chi_cavity = Processor.compute_reactance_cavity(L, K)

        r_nonlinear = []
        alpha = []

        for i in range(len(omega)):
            chi_tot_plate_i = chi_tot_plate[i]
            chi_cavity_i = chi_cavity[i]
            r_tot_plate_i = r_tot_plate[i]

            r_initial = 0.5
            r_solution = fsolve(processor.resistance_eq, r_initial, args=(p_acous_pa, r_tot_plate_i, chi_tot_plate_i, chi_cavity_i, sigma, M))[0]

            alpha.append(Processor.calculate_absorption_coefficient(r_solution, chi_tot_plate_i + chi_cavity_i))
            r_nonlinear.append(r_solution)

        processor.__impedance.set_resistance(r_nonlinear)
        processor.__impedance.set_reactance(chi_cavity + chi_tot_plate)
        processor.__impedance.set_absorption_coefficient(alpha)

        return processor.__impedance
    
    @staticmethod
    def solve(L, d, sigma, e, altitude, frequencies, p_acous_pa):
        flying_condition = Flying_condition(altitude)
        flying_condition.initialize_flying_condition()

        liner = Liner(L, d, sigma, e)
        wave = Wave(frequencies, flying_condition._environment.speed_of_sound)

        impedance = Processor.compute_impedance(flying_condition._environment, wave, liner, p_acous_pa, flying_condition._mach)

        return impedance

    @staticmethod
    def compute_impedance_varying_param(param, L, d, sigma, e, altitude, frequencies, p_acous_pa):
        impedances = {}
        if param == "L":
            impedances[param] = []
            for value in L:
                var = Processor.solve(value, d, sigma, e, altitude, frequencies, p_acous_pa)
                impedances[param].append({value: copy.deepcopy(var)})

        elif param == "d":
            impedances[param] = []
            for value in d:
                var = Processor.solve(L, value, sigma, e, altitude, frequencies, p_acous_pa)
                impedances[param].append({value: copy.deepcopy(var)})
    
        elif param == "sigma":
            impedances[param] = []
            for value in sigma:
                var = Processor.solve(L, d, value, e, altitude, frequencies, p_acous_pa)
                impedances[param].append({value: copy.deepcopy(var)})

        elif param == "e":
            impedances[param] = []
            for value in e:
                var = Processor.solve(L, d, sigma, value, altitude, frequencies, p_acous_pa)
                impedances[param].append({value: copy.deepcopy(var)})

        return impedances
    
    @staticmethod
    def gradient_descent(f, L_init, learning_rate = 0.01, num_iterations = 1000):

        x = L_init
        
        # calculate the gradient of f at (x, y)
        grad_x = grad(f)
        grad_y = grad(f,argnums=(1))
        
        for i in range(num_iterations):
            # update x and y using gradient descent
            x = x - learning_rate * grad_x(x,y)
            y = y - learning_rate * grad_y(x,y)

        # print the final value of the loss at (x, y)
        print("iteration {}: f(x, y) = {}".format(i, f(x, y)))
            
        return x
    
    @staticmethod
    def resistance_eq_opt(r, L, d, sigma, e, M, wave, c, nu, who, p_acous_pa):

        r_visc = np.sqrt(8 * nu * wave.omega) / (c* sigma) * (1 + e / d)
        r_rad = 1 / (8 * sigma) * (wave.k * d)**2
        r_tot_plate = r_visc + r_rad

        eps = 1 / (1 + 305 * M**3) # Correction factor when considering airflow, M is the mach number
        chi_mass = wave.omega / (sigma * c) * (e + eps * (8 * d) / (3 * np.pi) * (1 - 0.71 * np.sqrt(sigma)))
        chi_visc = wave.omega / (sigma * c) * (np.sqrt(8 * nu / wave._omega) * (1 + e / d))
        chi_tot_plate = chi_mass + chi_visc

        chi_cavity = - 1 / np.tan(k* L)

        r_airflow = 0.3 * (1 - sigma**2) / sigma * M

        B = (1 - sigma**2) / (sigma * c)
        impedance_magnitude = abs(r + (1j * (chi_tot_plate + chi_cavity)))
        r_airflow = r_airflow = 0.3 * (1 - sigma**2) / sigma * M

        return r - (r_tot_plate + B * p_acous_pa / (rho * self._environment.speed_of_sound * sigma * impedance_magnitude) + r_airflow)
    
