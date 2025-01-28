import numpy as np
from scipy.optimize import fsolve


class Impedance:
    """Class to calculate the impedance (resistance and reactance) for a given wave and environment."""
    def __init__(self):

        self._resistance = None
        self._reactance = None
        self._absorption_coefficient = None

    def set_resistance(self, resistance):
        self._resistance = resistance

    def set_reactance(self, reactance):
        self._reactance = reactance
    
    def set_absorption_coefficient(self, absorption_coefficient):
        self._absorption_coefficient = absorption_coefficient
