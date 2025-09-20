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

# E1 Transition strength calculator:
Branching_Ratio = ufloat(1.00, 0.01)  # Branching ratio (unitless)
E_gamma = ufloat(0.0581, 0.0001)  # Gamma energy in MeV
Tau = ufloat(135.6E3, 17.3E3)  # Lifetime in picoseconds converted from Lifetime in ns=135.6+/-17.3 ns
BE1 = 1/(1.59e15 * E_gamma**3 * Tau)  # B(E1) in e^2*fm^2
print(f"B(E1) = {BE1:.2e} e^2*fm^2")

# E1 Transition strength in Weisskopf units (W.u.)
# Formula: B(E1)(W.u.) = [6.764e-6 * BR] / [E_gamma^3 * A^(2/3) * T_half * (1 + alpha)]
# E_gamma in keV, T_half in seconds, BR unitless, alpha = conversion coefficient
# Branching_Ratio = ufloat(1.000, 0.001)  # Branching ratio (unitless)
Branching_Ratio = 1.0
E_gamma = ufloat(58.1, 0.1)  # keV
A = 105
T_half = ufloat(94e-9, 12e-9)     # seconds
alpha = ufloat(0.577, 0.009)
# alpha = 0
# Direct calculation for B(E1) (W.u.) using the formula:
# B(E1)(W.u.) = [6.764e-6 * BR] / [E_gamma**3 * A**(2/3) * T_half * (1 + alpha)]
BE1_Wu = (6.764e-6 * Branching_Ratio) / (E_gamma**3 * A**(2/3) * T_half * (1 + alpha))
print(f"B(E1) (W.u.) = {BE1_Wu:.3e}")
