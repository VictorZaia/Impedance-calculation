from Tools.Impedance_functions import *

from Properties.Wave import Wave
from Properties.Flying_condition import Flying_condition
from Properties.Impedance import Impedance
from Liner.Liner import Liner

from scipy.optimize import fsolve
import copy

class Processor:

    __impedance = Impedance()

    def __init__(self, environment, wave, liner):
        self._environment = environment
        self._wave = wave
        self._liner = liner
    
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

        rho = processor._environment.rho
        nu = processor._environment.nu
        speed_of_sound = processor._environment.speed_of_sound

        r_tot_plate = compute_resistance_plate(omega, K, sigma, d, e, nu, speed_of_sound)
        chi_tot_plate = compute_reactance_plate(omega, sigma, d, e, nu, speed_of_sound, M)

        chi_cavity = compute_reactance_cavity(L, K)

        r_nonlinear = []
        alpha = []

        for i in range(len(omega)):
            chi_tot_plate_i = chi_tot_plate[i]
            chi_cavity_i = chi_cavity[i]
            r_tot_plate_i = r_tot_plate[i]

            r_initial = 0.5
            r_solution = fsolve(resistance_eq, r_initial, args=(p_acous_pa, r_tot_plate_i, chi_tot_plate_i, chi_cavity_i, sigma, rho, speed_of_sound, M))[0]

            alpha.append(calculate_absorption_coefficient(r_solution, chi_tot_plate_i + chi_cavity_i))
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