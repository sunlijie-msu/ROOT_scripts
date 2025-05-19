from uncertainties import ufloat
from uncertainties.umath import * # for functions such as sin(), etc.
import math

# In Python, ufloat is used to represent a number with its uncertainty.
# It is part of the 'uncertainties' package, which allows for easy propagation of uncertainties through calculations.

# Define a variable with uncertainty
a = ufloat(2.0, 0.1)  # 2.0 ± 0.1
b = ufloat(3.0, 0.2)  # 3.0 ± 0.2
c = ufloat(4.0, 0.3)  # 4.0 ± 0.3

# Perform calculations
x = a / (a + b)
y = b / (a + b)

# Print results with uncertainties
print(f"x = {x:.2f}")
print(f"y = {y:.2f}")

# Given uncertainties
stat_uncertainty = 0.3  # ms
syst_uncertainty = 0.6  # ms

# Calculate total uncertainty
total_uncertainty = math.sqrt(stat_uncertainty**2 + syst_uncertainty**2)

print(f"Total uncertainty: ±{total_uncertainty:.2f} ms")