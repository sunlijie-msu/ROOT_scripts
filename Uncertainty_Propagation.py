from uncertainties import ufloat
from uncertainties.umath import * # for functions such as sin(), etc.
import math

# In Python, ufloat is used to represent a number with its uncertainty.
# It is part of the 'uncertainties' package, which allows for easy propagation of uncertainties through calculations.

# Define a variable with uncertainty
# example: c = ufloat(4.0, 0.3)  # 4.0 ± 0.3

Sp = ufloat(2749.0, 0.4)
Ef28Si = ufloat(1779.030, 0.011)
p23 = ufloat(6656, 8)
p18 = ufloat(4870, 4)
E1 = p18 + Ef28Si + Sp
E2 = p23 + Sp
w1 = 1 / E1.s**2
w2 = 1 / E2.s**2
Eweighted = (w1 * E1 + w2 * E2) / (w1 + w2)
Eunweighted = (E1 + E2) / 2

print(f"E1: {E1:.3f} keV")
print(f"E2: {E2:.3f} keV")
print(f"Weighted mean: {Eweighted:.3f} keV")


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


# Transition strength calculator:
Branching_Ratio = ufloat(0.82, 0.00)  # Branching ratio (unitless)
E_gamma = ufloat(0.610, 0.001)  # Gamma energy in MeV
Tau = ufloat(1.34, 0.20)  # Lifetime in picoseconds
BE2 = 0.0816 * Branching_Ratio / (E_gamma**5 * Tau)  # B(E2) in e^2*b^2
print(f"B(E2) = {BE2:.4f} e^2*b^2")


# Calculate: 30973761.9977(8)/34976365.2518(8)
result = ufloat(30973761.9977, 0.0008) / ufloat(34976365.2518, 0.0008)
print(f"Result = {result:.9f}")


I_gamma_to_34S_ground = ufloat(129.3, 1.2) # From 1980Wi13
I_gamma_to_34Cl_ground = ufloat(94.4, 1.4) # From 1980Wi13 146gamma
I_gamma_to_34S_ground = ufloat(125.9, 1.2) # From 1975Va02
I_gamma_to_34Cl_ground = ufloat(101.0, 4.0) # From 1975Va02 146gamma
CC = ufloat(0.1656, 0.0024)  # Conversion coefficient
I_ce_to_34Cl_ground = I_gamma_to_34Cl_ground * CC
I_isomeric_transition_to_34Cl_ground = I_gamma_to_34Cl_ground * (1 + CC)
I_total = I_gamma_to_34S_ground + I_gamma_to_34Cl_ground + I_ce_to_34Cl_ground
BR_gamma_to_34S_ground = I_gamma_to_34S_ground / I_total
BR_isomeric_transition_to_34Cl_ground = I_isomeric_transition_to_34Cl_ground / I_total

print(f"I(CE to 34Cl ground) = {I_ce_to_34Cl_ground:.2f}")
print(f"I(isomeric transition to 34Cl ground) = {I_isomeric_transition_to_34Cl_ground:.2f}")
print(f"I(total) = {I_total:.2f}")
print(f"BR(gamma to 34S ground) = {BR_gamma_to_34S_ground:.4f}")
print(f"BR(isomeric transition to 34Cl ground) = {BR_isomeric_transition_to_34Cl_ground:.4f}")
