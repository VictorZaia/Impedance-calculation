class Impedance:
    
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

    def get_resistance(self):
        return self._resistance

    def get_reactance(self):
        return self._reactance
    
    def get_absorption_coefficient(self):
        return self._absorption_coefficient
