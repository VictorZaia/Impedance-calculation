import numpy as np

class Environment:
    """Attributes"""
    # Constants
    __R = 287.05  # Specific gas constant for air (J/kg·K)
    __gamma = 1.4  # Adiabatic index for air
    __mu_0 = 1.716e-5  # Reference viscosity at T0 (Pa·s)
    __T0 = 288.15   # Reference temperature (K)
    __P0 = 101325 # Reference pressure (Pa)
    __C = 110.4  # Sutherland's constant (K)
    __g = 9.81
    __L = -0.0065 # K/m

    @classmethod
    def get_R(cls):
        return cls.__R
    
    @classmethod
    def get_gamma(cls):
        return cls.__gamma
    
    @classmethod
    def get_mu0(cls):
        return cls.__mu_0

    @classmethod
    def get_T0(cls):
        return cls.__T0
    
    @classmethod
    def get_P0(cls):
        return cls.__P0
    
    @classmethod
    def get_C(cls):
        return cls.__C
    
    @classmethod
    def get_g(cls):
        return cls.__g
    
    @classmethod
    def get_L(cls):
        return cls.__L
    
    """Constructor"""

    def __init__(self, temperature, pressure):
        self._temperature = temperature
        self._pressure = pressure

    @property
    def rho(self):
        return self._pressure / (self.__R * self._temperature)

    @property
    def speed_of_sound(self):
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
        return (f"Environment(temperature={self._temperature}, pressure={self._pressure}, rho={self.rho}, c={self.speed_of_sound}, mu={self.mu}, nu={self.nu})")