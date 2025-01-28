from Environment import *
from Wave import *
from Liner import *
from Impedance import *

import matplotlib.pyplot as plt
import matplotlib.cm as cm

"""Initialize the air properties"""

temperature = 300
pressure = 101325

env = Environment(temperature, pressure)

"""Initialize the wave"""

frequencies = np.linspace(0.1, 8000, 100)
wave = Wave(frequencies, env.c)


"""Initialize the liner geometry"""

L = 15e-3
d = 1.5e-3
sigma = 0.15
e = 1.5e-3

liner = Liner(L, d, sigma, e)

"""Flying condition"""
p_acous_pa = 1000
M_values = np.linspace(0, 0.22, 5)

cmap = cm.get_cmap('viridis', M_values.size)

plt.figure(figsize=(10, 6))
for i in range(M_values.size):
    im = Impedance(env, wave, liner, p_acous_pa, M_values[i])
    im.compute_impedance_for_frequencies()
    plt.plot(frequencies, im._resistance)
plt.xlabel("cavity size")
plt.ylabel("Frequency (Hz)")
plt.legend()
plt.grid(True)
plt.show()