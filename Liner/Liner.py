class Liner:

    def __init__(self, L, d, sigma, e):
        self._L = L
        self._d = d
        self._sigma = sigma
        self._e = e

    def __repr__(self):
        return (f"Liner(L={self._L}, d={self._d}, sigma={self._sigma}, e={self._e}")
