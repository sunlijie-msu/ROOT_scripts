from uncertainties import ufloat
from uncertainties.umath import log

hbar = 6.582119569e-16 # eV*s
half_life = ufloat(67.87e-9, 0.07e-9) # s

width = hbar / (half_life / log(2))

print(f"Decay Width (width): {width} eV")