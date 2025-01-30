from Properties.Environment import *
from Properties.Wave import *
from Liner.Liner import *
from Solver.Impedance import *
from Solver.Processor import *
from Solver.Post_Processor import Post_Processor as post
from Properties.Flying_condition import *

"""Frequencies"""

frequencies = np.linspace(0.1, 8000, 100)

"""Noise level"""

p_acous_pa = 1000

"""Altitude"""

altitudes = 5000#np.linspace(0, 12000, 4)

"""Initialize the liner geometry"""

L = np.linspace(10e-3, 20e-3, 5)
d = 1.5e-3
sigma = 0.15
e = 1.5e-3

impedances = Processor.compute_impedance_varying_param("L", L, d, sigma, e, altitudes, frequencies, p_acous_pa)
print(impedances)
post.plot_absorption_coefficients(frequencies, impedances)



L = 15e-3
d  = np.linspace(1e-3, 2.e-3, 5)
sigma = 0.15
e = 1.5e-3

impedances = Processor.compute_impedance_varying_param("d", L, d, sigma, e, altitudes, frequencies, p_acous_pa)
print(impedances)
post.plot_absorption_coefficients(frequencies, impedances)



L = 15e-3
d  = 1.5e-3
sigma = np.linspace(0.1, 0.2, 5)
e = 1.5e-3

impedances = Processor.compute_impedance_varying_param("sigma", L, d, sigma, e, altitudes, frequencies, p_acous_pa)
print(impedances)
post.plot_absorption_coefficients(frequencies, impedances)



L = 15e-3
d  = 1.5e-3
sigma = 0.15
e = np.linspace(1e-3, 2.e-3, 5)

impedances = Processor.compute_impedance_varying_param("e", L, d, sigma, e, altitudes, frequencies, p_acous_pa)
print(impedances)
post.plot_absorption_coefficients(frequencies, impedances)































"""
plt.figure(figsize=(10, 6))
for Ls in L:
    liner_var = Liner(Ls, d, sigma, e)


    flight_conditions = Flying_condition(5000)
    flight_conditions.initialize_flying_condition()

    wave = Wave(frequencies, flight_conditions._environment.speed_of_sound)

    impedance = Processor.compute_impedance(flight_conditions._environment, wave, liner_var, p_acous_pa, flight_conditions._mach)

    plt.plot(frequencies, impedance.get_absorption_coefficient(), label=f"L={Ls}")
plt.legend()
plt.xlabel("Frequency (Hz)")
plt.ylabel("Absorption")
plt.ylim(0,1)
plt.grid(True)
plt.show()











L = 15e-3
d_values  = np.linspace(1e-3, 2.e-3, 5)
sigma = 0.15
e = 1.5e-3

plt.figure(figsize=(10, 6))
for ds in d_values:
    liner_var = Liner(L, ds, sigma, e)


    flight_conditions = Flying_condition(5000)
    flight_conditions.initialize_flying_condition()

    wave = Wave(frequencies, flight_conditions._environment.speed_of_sound)

    impedance = Processor.compute_impedance(flight_conditions._environment, wave, liner_var, p_acous_pa, flight_conditions._mach)

    plt.plot(frequencies, impedance.get_absorption_coefficient(), label=f"d={ds}")
plt.legend()
plt.xlabel("Frequency (Hz)")
plt.ylabel("Absorption")
plt.ylim(0,1)
plt.grid(True)
plt.show()


















L = 15e-3
d  = 1.5e-3
sigma_values = np.linspace(0.1, 0.2, 5)
e = 1.5e-3

plt.figure(figsize=(10, 6))
for sigmas in sigma_values:
    liner_var = Liner(L, d, sigmas, e)


    flight_conditions = Flying_condition(5000)
    flight_conditions.initialize_flying_condition()

    wave = Wave(frequencies, flight_conditions._environment.speed_of_sound)

    impedance = Processor.compute_impedance(flight_conditions._environment, wave, liner_var, p_acous_pa, flight_conditions._mach)

    plt.plot(frequencies, impedance.get_absorption_coefficient(), label=f"sigma={sigmas}")
plt.legend()
plt.xlabel("Frequency (Hz)")
plt.ylabel("Absorption")
plt.ylim(0,1)
plt.grid(True)
plt.show()














L = 15e-3
d  = 1.5e-3
sigma = 0.15
e_values = np.linspace(1e-3, 2.e-3, 5)

plt.figure(figsize=(10, 6))
for es in e_values:
    liner_var = Liner(L, d, sigma, es)


    flight_conditions = Flying_condition(5000)
    flight_conditions.initialize_flying_condition()

    wave = Wave(frequencies, flight_conditions._environment.speed_of_sound)

    impedance = Processor.compute_impedance(flight_conditions._environment, wave, liner_var, p_acous_pa, flight_conditions._mach)

    plt.plot(frequencies, impedance.get_absorption_coefficient(), label=f"e={es}")
plt.legend()
plt.xlabel("Frequency (Hz)")
plt.ylabel("Absorption")
plt.ylim(0,1)
plt.grid(True)
plt.show()"""