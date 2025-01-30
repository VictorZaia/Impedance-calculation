import numpy as np
import matplotlib.pyplot as plt

# Constants
T0 = 288.15  # Sea level temperature in K
P0 = 101325  # Sea level pressure in Pa
L = -0.0065  # Temperature lapse rate in K/m
g = 9.80665  # Gravitational acceleration in m/s^2
R = 287.05   # Specific gas constant for air in J/(kg·K)
gamma = 1.4  # Adiabatic index for air

def temperature_troposphere(h):
    """Temperature in the troposphere (0–11 km)"""
    return T0 + L * h

def pressure_troposphere(h):
    """Pressure in the troposphere (0–11 km)"""
    return P0 * (temperature_troposphere(h) / T0) ** (-g / (L * R))

def pressure_stratosphere(h):
    """Pressure in the stratosphere (11–20 km)"""
    P11 = pressure_troposphere(11000)
    return P11 * np.exp(-g * (h - 11000) / (R * 216.65))

def speed_of_sound(temperature):
    return np.sqrt(gamma * R * temperature)

def mach_number(airspeed, temperature):
    c = speed_of_sound(temperature)
    return airspeed / c

def airspeed_during_flight(altitude):
    if altitude < 1000:
        return 70 + (altitude / 1000) * 100  # From 70 m/s to 170 m/s during takeoff
    elif altitude < 10000:
        return 170 + (altitude - 1000) / 9000 * 80  # Increase airspeed during climb
    else: 
        return 250  # Constant cruise speed in m/s (approx Mach 0.85)

altitudes = np.linspace(0, 12000, 100) 

# Calculate temperature and pressure for each altitude
temperatures = np.piecewise(altitudes, [altitudes <= 11000, altitudes > 11000],
                            [temperature_troposphere, 216.65])
pressures = np.piecewise(altitudes, [altitudes <= 11000, altitudes > 11000],
                         [pressure_troposphere, pressure_stratosphere])

mach_numbers = []
airspeeds = []

for altitude in altitudes:
    airspeed = airspeed_during_flight(altitude) 
    temperature = temperatures[np.where(altitudes == altitude)[0]] 

    # Calculate Mach number and impedance
    mach = mach_number(airspeed, temperature)

    
    airspeeds.append(airspeed)
    mach_numbers.append(mach)

fig, ax1 = plt.subplots(figsize=(10, 6))

ax1.set_xlabel('Altitude (m)')
ax1.set_ylabel('Temperature (K)', color='tab:red')
ax1.plot(altitudes, temperatures, color='tab:red', label='Temperature')
ax1.tick_params(axis='y', labelcolor='tab:red')

ax2 = ax1.twinx() 
ax2.set_ylabel('Pressure (Pa)', color='tab:blue')
ax2.plot(altitudes, pressures, color='tab:blue', label='Pressure')
ax2.tick_params(axis='y', labelcolor='tab:blue')

ax3 = ax1.twinx() 
ax3.spines['right'].set_position(('outward', 70))
ax3.set_ylabel('Mach Number', color='tab:green')
ax3.plot(altitudes, mach_numbers, color='tab:green', label='Mach Number')
ax3.tick_params(axis='y', labelcolor='tab:green')

fig.tight_layout()
plt.show()


