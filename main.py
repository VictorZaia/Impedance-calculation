from Properties.Environment import *
from Properties.Wave import *
from Properties.Impedance import *
from Properties.Flying_condition import *
from Liner.Liner import *
from Solver.Processor import *
from Solver.Post_Processor import Post_Processor as post
from Solver.Optimization import Optimizer as opt

"""Frequencies"""

frequencies = np.linspace(0.1, 5000, 100)
frequencie = 2000

"""Noise level"""

p_acous_pa = 1000

"""Altitude"""

altitudes = np.linspace(0, 12000, 100)
altitude = 2000

"""Initialize the liner geometry"""

L = np.linspace(10e-3, 20e-3, 5)
d = 1.5e-3
sigma = 0.15#np.linspace(0.05, 0.20, 15) 0.15
e = 1.5e-3

impedances = Processor.compute_impedance_varying_param("L", L, d, sigma, e, altitude, frequencies, p_acous_pa)
post.plot_absorption_coefficients(frequencies, impedances)

L = 15e-3
d  = np.linspace(1e-3, 2.e-3, 5)
sigma = 0.15
e = 1.5e-3

impedances = Processor.compute_impedance_varying_param("d", L, d, sigma, e, altitude, frequencies, p_acous_pa)
post.plot_absorption_coefficients(frequencies, impedances)



L = 15e-3
d  = 1.5e-3
sigma = np.linspace(0.1, 0.2, 5)
e = 1.5e-3

impedances = Processor.compute_impedance_varying_param("sigma", L, d, sigma, e, altitude, frequencies, p_acous_pa)
post.plot_absorption_coefficients(frequencies, impedances)



L = 15e-3
d  = 1.5e-3
sigma = 0.15
e = np.linspace(1e-3, 2.e-3, 5)

impedances = Processor.compute_impedance_varying_param("e", L, d, sigma, e, altitude, frequencies, p_acous_pa)
post.plot_absorption_coefficients(frequencies, impedances)


L = np.linspace(10e-3, 50e-3, 15)
d = 1.5e-3
sigma = np.linspace(0.05, 0.30, 15)
e = 1.5e-3


opt.plot_loss_function(L, d, sigma, e, altitude, frequencie, p_acous_pa)

opt.plot_3d_loss_function(L, d, sigma, e, altitudes, frequencie, p_acous_pa)




















































































