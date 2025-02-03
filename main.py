# %%
from Solver.Processor import *
from Solver.Post_Processor import Post_Processor as post
from Solver.Optimization import *
# %% [markdown]
"""
# **Computing and ploting the absorption coefficient for varying parameters**
The impedance calculation will be given by: 

$$ r(w,|v|) = A(w) + B \frac{|v|}{\sigma} + 0.3 * \frac{(1 - \sigma^2)}{\sigma} * M $$

the resistance can be given by a linear part, called A(w) and a non linear part which depends on the acoustic velocity |v|, multiplied by a coefficient B.

$$ r = A(w) + B * \frac{|p|}{\rho  c  \sigma \sqrt{r^2 + (\chi - j \cot{kL})^2}} + 0.3 * \frac{(1 - \sigma^2)}{\sigma} * M $$

"""
# %% [markdown]
"""
## Frequencies
Used to compute the wave's properties.
"""
frequencies = np.linspace(0.1, 5000, 200)

# %% [markdown]
"""
## Noise level
Acoustic pressure that represents the sound level in Pa.
"""
p_acous_pa = 1000 # Pa

# %% [markdown]
"""
## Altitude
Altitude where the air's properties will be computed.
"""
altitude = 2000

# %% [markdown]
"""
## Liner's geometry
Geometry of the liner: size of the cavity, diameter of the hole, porosity and thickness of the plate.
"""
L = 15e-3
d = 1.5e-3
sigma = 0.15
e = 1.5e-3

# %%
L_var = np.linspace(10e-3, 20e-3, 5)

impedances = Processor.compute_impedance_varying_param("L", L_var, d, sigma, e, altitude, frequencies, p_acous_pa)
post.plot_absorption_coefficients(frequencies, impedances)

d_var  = np.linspace(1e-3, 2.e-3, 5)

impedances = Processor.compute_impedance_varying_param("d", L, d_var, sigma, e, altitude, frequencies, p_acous_pa)
post.plot_absorption_coefficients(frequencies, impedances)

sigma_var = np.linspace(0.1, 0.2, 5)

impedances = Processor.compute_impedance_varying_param("sigma", L, d, sigma_var, e, altitude, frequencies, p_acous_pa)
post.plot_absorption_coefficients(frequencies, impedances)

e_var = np.linspace(1e-3, 2.e-3, 5)

impedances = Processor.compute_impedance_varying_param("e", L, d, sigma, e_var, altitude, frequencies, p_acous_pa)
post.plot_absorption_coefficients(frequencies, impedances)

# %% [markdown]
"""
# **Studying alpha behaviour for varying heights and porosity**
## Frequency of interest
frequency which we want to optimize the abosorption coefficient.
"""
frequencie = 1000

# %% [markdown]
"""
## Altitudes during flight
"""
altitudes = np.linspace(0, 12000, 100)

#%% [markdown]
"""
## Varying parameters
Values of height and porosity that will be considered.
"""
L_var = np.linspace(10e-3, 20e-3, 50)
sigma_var = np.linspace(0.1, 0.2, 50)

# %%
alpha = Optimizer.loss_function_values(L_var, d, sigma_var, e, altitude, frequencie, p_acous_pa)
post.plot_loss_function(sigma_var, L_var, alpha)

#%%
alpha_max, optimum_L, optimum_sigma = Optimizer.optimum_geometry(L_var, d, sigma_var, e, altitudes, frequencie, p_acous_pa)
post.plot_optimum_geometry(altitudes, alpha_max, optimum_L, optimum_sigma)

# %%
