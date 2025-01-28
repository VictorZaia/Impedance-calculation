import numpy as np

class Wave:

    def __init__(self, frequencies, c):

        self._frequencies = frequencies
        self._c = c
    
    @property
    def omega(self):

        return 2 * np.pi * self._frequencies
    
    @property
    def K(self):

        return self.omega / self._c
    
    def __repr__(self):

        return (f"Wave(Frequencies={self._frequencies}, speed of sound={self._c}, omega={self.omega}, wave number={self.K})")
        