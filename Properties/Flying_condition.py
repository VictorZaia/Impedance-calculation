from Properties.Environment import *

class Flying_condition:

    def __init__(self, _altitude):
        self._environment = None
        self._mach = None
        self._altitude = _altitude

    def temperature_troposphere(self):
        """
        Calculate temperature in the troposphere (0–11 km)
        """
        return Environment.get_T0() + Environment.get_L() * self._altitude

    def pressure_troposphere(self):
        """
        Calculate pressure in the troposphere (0–11 km)
        """
        temperature = self.temperature_troposphere()
        return Environment.get_P0() * (temperature / Environment.get_T0()) ** (-Environment.get_g() / (Environment.get_L() * Environment.get_R()))

    def pressure_stratosphere(self):
        """
        Calculate pressure in the stratosphere (11–20 km)
        """
        P11 = self.pressure_troposphere()
        return P11 * np.exp(-Environment.get_g() * (self._altitude - 11000) / (Environment.get_R() * 216.65))

    def mach_number(self, airspeed, c):
        """
        Calculate Mach number
        """
        return airspeed / c

    def airspeed_during_flight(self):
        """
        Calculate airspeed during flight based on altitude
        """
        if self._altitude < 1000:
            return 70 + (self._altitude / 1000) * 100
        elif self._altitude < 10000:
            return 170 + (self._altitude - 1000) / 9000 * 80
        else: 
            return 250

    def initialize_flying_condition(self):
        """
        Calculate and return the environment based on current altitude
        """
        if self._altitude < 11000:
            temperature = self.temperature_troposphere()
            pressure = self.pressure_troposphere()
        else:
            temperature = 216.65
            pressure = self.pressure_stratosphere()
        
        self._environment = Environment(temperature, pressure)

        airspeed = self.airspeed_during_flight()
        self._mach = self.mach_number(airspeed, self._environment.speed_of_sound)

    def __repr__(self):
        return (f"Flying Condition(Altitude={self._altitude}, Environment={self._environment}, Mach number={self._mach}")

    
    