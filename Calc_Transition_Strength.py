from uncertainties import ufloat
from uncertainties.umath import * # for functions such as sin(), etc.
import math


# E2 Transition strength calculator:
Branching_Ratio = ufloat(1.00, 0.01)  # Branching ratio (unitless)
E_gamma = ufloat(0.7449, 0.001)  # Gamma energy in MeV
Tau = ufloat(2.41, 0.30)  # Lifetime in picoseconds
BE2 = 0.0816 * Branching_Ratio / (E_gamma**5 * Tau)  # B(E2) in e^2*b^2
print(f"B(E2) = {BE2:.2f} e^2*b^2")

# E1 Transition strength calculator:
Branching_Ratio = ufloat(1.00, 0.01)  # Branching ratio (unitless)
E_gamma = ufloat(0.4903, 0.001)  # Gamma energy in MeV
Tau = ufloat(0.91, 0.11)  # Lifetime in picoseconds
BE1 = 1/(1.59e15 * E_gamma**3 * Tau)  # B(E1) in e^2*fm^2
print(f"B(E1) = {BE1:.2e} e^2*fm^2")