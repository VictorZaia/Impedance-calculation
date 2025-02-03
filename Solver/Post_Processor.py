from Properties.Impedance import *

import matplotlib.pyplot as plt
class Post_Processor:

    @staticmethod
    def plot_resistance(frequencies, impedance):

        plt.figure(figsize=(10, 6))
        plt.plot(frequencies, impedance.get_resistance())
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Resistance (r)")
        plt.grid(True)
        plt.show()

    @staticmethod
    def plot_reactance(frequencies, impedance):

        plt.figure(figsize=(10, 6))
        plt.plot(frequencies, impedance.get_reactance())
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Reactance (xhi)")
        plt.ylim(-10, 10)
        plt.grid(True)
        plt.show()

    @staticmethod
    def plot_absorption_coefficient(frequencies, impedance):

        plt.figure(figsize=(10, 6))
        plt.plot(frequencies, impedance.get_absorption_coefficient())
        plt.xlabel("Frequency [Hz]")
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
        plt.xlabel("Frequency [Hz]")
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
        plt.xlabel("Frequency [Hz]")
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
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Absorption coefficient")
        plt.legend()
        plt.grid(True)
        plt.show()

    @staticmethod
    def plot_loss_function(sigma, L, alpha):
        
        plt.figure(figsize=(10, 6))
        plt.pcolormesh(sigma * 100, L * 10**3, alpha, cmap='rainbow')
        plt.colorbar(label="Absorption Coefficient")
        plt.xlabel("sigma [%]")
        plt.ylabel("L [mm]")
        plt.grid(True)
        plt.show()
    
    @staticmethod
    def plot_optimum_geometry(altitudes, alpha_max, optimal_L, optimal_sigma):

        plt.figure(figsize=(10, 6))

        plt.subplot(3, 1, 1)
        plt.plot(altitudes, alpha_max, label='Optimal L', color='b')
        plt.xlabel("Altitude (m)")
        plt.ylabel("Optimal L (m)")
        plt.title("Maximum absorption vs Altitude")
        plt.grid(True)

        plt.subplot(3, 1, 2)
        plt.plot(altitudes, optimal_L, label='Optimal L', color='r')
        plt.xlabel("Altitude (m)")
        plt.ylabel("Optimal L (m)")
        plt.title("Optimal Cavity Height (L) vs Altitude")
        plt.grid(True)

        plt.subplot(3, 1, 3)
        plt.plot(altitudes, optimal_sigma, label='Optimal Sigma', color='g')
        plt.xlabel("Altitude (m)")
        plt.ylabel("Optimal Sigma")
        plt.title("Optimal Sigma vs Altitude")
        plt.grid(True)

        plt.tight_layout()
        plt.show()
        
        





