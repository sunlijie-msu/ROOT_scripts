from uncertainties import ufloat

# --- Inputs ---
# Spins
Jr = 3.5  # Resonance spin (7/2)
Jp = 0.5  # Proton spin (1/2)
JT = 0.0  # Target spin (0)

# Partial Widths [eV] with uncertainties
# Gp = ufloat(3.8, 1.3)  # Proton partial width 1975KE11
# Gg = ufloat(2.0, 1.3)  # Gamma partial width 1975KE11
Gp = ufloat(7.2, 1.5)  # Proton partial width 1974BI16
Gg = ufloat(3.8, 1.5)  # Gamma partial width 1974BI16
# Gp = ufloat(5.8, 0.7)  # Proton partial width 1971BRXT
# Gg = ufloat(4.7, 1.7)  # Gamma partial width 1971BRXT

# --- Calculation ---
# Statistical factor: w = (2Jr + 1) / [(2Jp + 1)(2JT + 1)]
stat_factor = (2 * Jr + 1) / ((2 * Jp + 1) * (2 * JT + 1))

# Resonance Strength: wy = w * (Gp * Gg) / (Gp + Gg)
wy = stat_factor * (Gp * Gg) / (Gp + Gg)

# --- Output ---
print(f"Statistical Factor (omega): {stat_factor}")
print(f"Proton Width (Gp):          {Gp} eV")
print(f"Gamma Width (Gg):           {Gg} eV")
print(f"Total Width (Gamma):        {Gp + Gg} eV")
print("-" * 40)
print(f"Resonance Strength (wy):    {wy:.4f} eV")
print(f"Resonance Strength (S):     {wy * 2:.4f} eV")
