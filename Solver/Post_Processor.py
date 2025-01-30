from Solver.Impedance import *

import matplotlib.pyplot as plt

class Post_Processor:

    @staticmethod
    def plot_resistance(frequencies, impedance):
        plt.figure(figsize=(10, 6))
        plt.plot(frequencies, impedance.get_resistance())
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Resistance (r)")
        plt.grid(True)
        plt.show()

    @staticmethod
    def plot_reactance(frequencies, impedance):
        plt.figure(figsize=(10, 6))
        plt.plot(frequencies, impedance.get_reactance())
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Reactance (xhi)")
        plt.ylim(-10, 10)
        plt.grid(True)
        plt.show()

    @staticmethod
    def plot_absorption_coefficient(frequencies, impedance):
        plt.figure(figsize=(10, 6))
        plt.plot(frequencies, impedance.get_absorption_coefficient())
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Absorption coefficient")
        plt.grid(True)
        plt.show()
    
    @staticmethod
    def plot_resistances(frequencies, impedances):
        plt.figure(figsize=(10, 6))
        impedances_param = next(iter(impedances.keys()))
        impedances_values = impedances.get(impedances_param)
        
        for impedance_dict in impedances_values:
            impedance = list(impedance_dict.values())[0].get_resistance()
            plt.plot(frequencies, impedance, label=f'{impedances_param}={list(impedance_dict.keys())[0]:.3e}')
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Absorption coefficient")
        plt.legend()
        plt.grid(True)
        plt.show()

    @staticmethod
    def plot_reactances(frequencies, impedances):
        plt.figure(figsize=(10, 6))
        impedances_param = next(iter(impedances.keys()))
        impedances_values = impedances.get(impedances_param)
        
        for impedance_dict in impedances_values:
            impedance = list(impedance_dict.values())[0].get_reactance()
            plt.plot(frequencies, impedance, label=f'{impedances_param}={list(impedance_dict.keys())[0]:.3e}')
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Absorption coefficient")
        plt.legend()
        plt.grid(True)
        plt.show()

    @staticmethod
    def plot_absorption_coefficients(frequencies, impedances):
        plt.figure(figsize=(10, 6))
        impedances_param = next(iter(impedances.keys()))
        impedances_values = impedances.get(impedances_param)
        
        for impedance_dict in impedances_values:
            impedance = list(impedance_dict.values())[0].get_absorption_coefficient()
            plt.plot(frequencies, impedance, label=f'{impedances_param}={list(impedance_dict.keys())[0]:.3e}')
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Absorption coefficient")
        plt.legend()
        plt.grid(True)
        plt.show()
        
        





