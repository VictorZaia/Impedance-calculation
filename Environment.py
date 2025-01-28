import numpy as np

class Environment:
    """Attributes"""
    # Constants
    __R = 287.05  # Specific gas constant for air (J/kg·K)
    __gamma = 1.4  # Adiabatic index for air
    __mu_0 = 1.716e-5  # Reference viscosity at T0 (Pa·s)
    __T0 = 273.15  # Reference temperature (K)
    __C = 110.4  # Sutherland's constant (K)

    def __init__(self, temperature, pressure):

        self._temperature = temperature
        self._pressure = pressure

    @property
    def rho(self):

        return self._pressure / (self.__R * self._temperature)

    @property
    def c(self):

        return np.sqrt(self.__gamma * self.__R * self._temperature)

    @property
    def mu(self):

        return self.__mu_0 * (self.__T0 + self.__C) / (self._temperature + self.__C) * (self._temperature / self.__T0)**1.5

    @property
    def nu(self):

        return self.mu / self.rho

    def update_conditions(self, temperature, pressure):

        self._temperature = temperature
        self._pressure = pressure

    def __repr__(self):

        return (f"Environment(temperature={self._temperature}, pressure={self._pressure}, rho={self.rho}, c={self.c}, mu={self.mu}, nu={self.nu})")

