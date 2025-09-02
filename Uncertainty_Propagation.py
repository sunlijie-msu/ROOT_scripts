from uncertainties import ufloat
from uncertainties.umath import * # for functions such as sin(), etc.
import math

# In Python, ufloat is used to represent a number with its uncertainty.
# It is part of the 'uncertainties' package, which allows for easy propagation of uncertainties through calculations.

# Define a variable with uncertainty
# example: c = ufloat(4.0, 0.3)  # 4.0 ± 0.3
g = ufloat(-1.2753, 0.0013)  # 2.0 ± 0.1
K = ufloat(6144.48, 3.7)
logft = ufloat(4.1144, 0.0237)


# Perform calculations
ft = 10**logft  # Convert log(ft) to ft
BGT = K / g**2 / ft

# Print results with uncertainties
print(f"BGT = {BGT:.4f}")

# Unrelated example:
# Given uncertainties
stat_uncertainty = 0.3  # ms
syst_uncertainty = 0.6  # ms

# Calculate total uncertainty
total_uncertainty = math.sqrt(stat_uncertainty**2 + syst_uncertainty**2)

print(f"Total uncertainty: ±{total_uncertainty:.2f} ms")

Gamma = ufloat(0.005, 0.0008)  # Width in eV.
hbar = 6.582119569e-16  # eV*s
tau = hbar / Gamma * 1e15  # Lifetime in femtoseconds
print(f"Lifetime (tau) = {tau:.4f} fs")