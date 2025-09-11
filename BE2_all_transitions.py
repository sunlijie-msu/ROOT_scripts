from uncertainties import ufloat
from uncertainties.umath import * # for functions such as sin(), etc.
import math

print("B(E2) Transition Strength Calculator for 127I")
print("=" * 50)
print("Formula: B(E2) = 0.0816 * BR / (E_gamma^5 * Tau)")
print("Units: B(E2) in e²b², E_gamma in MeV, Tau in ps")
print("=" * 50)

# Band A transitions
print("\nBAND A:")
print("-------")

# 27/2− → 23/2− (982.1 keV)
Branching_Ratio = ufloat(0.99, 0.01)  # 99% branching ratio
E_gamma = ufloat(0.9821, 0.001)  # 982.1 keV in MeV
Tau = ufloat(1.32, 0.125)  # τAve=1.32 ps {I+12-13} → use average uncertainty (12+13)/2=12.5
BE2 = 0.0816 * Branching_Ratio / (E_gamma**5 * Tau)
print(f"27/2− (982.1 keV): B(E2) = {BE2:.2f} e²b²")

# 23/2− → 19/2− (431.2 keV)
Branching_Ratio = ufloat(0.98, 0.01)  # 98% branching ratio
E_gamma = ufloat(0.4312, 0.001)  # 431.2 keV in MeV
Tau = ufloat(2.00, 0.29)  # τAve=2.00 ps {I29} - CORRECTED: symmetric uncertainty
BE2 = 0.0816 * Branching_Ratio / (E_gamma**5 * Tau)
print(f"23/2− (431.2 keV): B(E2) = {BE2:.2f} e²b²")

# 19/2− → 15/2− (651.5 keV)
Branching_Ratio = ufloat(0.97, 0.01)  # 97% branching ratio
E_gamma = ufloat(0.6515, 0.001)  # 651.5 keV in MeV
Tau = ufloat(1.73, 0.21)  # τAve=1.73 ps {I+22-20} → use average uncertainty (22+20)/2=21
BE2 = 0.0816 * Branching_Ratio / (E_gamma**5 * Tau)
print(f"19/2− (651.5 keV): B(E2) = {BE2:.2f} e²b²")

# 15/2− → 11/2− (658.7 keV)
Branching_Ratio = ufloat(0.92, 0.01)  # 92% branching ratio
E_gamma = ufloat(0.6587, 0.001)  # 658.7 keV in MeV
Tau = ufloat(1.01, 0.13)  # τAve=1.01 ps {I+12-14} → use average uncertainty (12+14)/2=13
BE2 = 0.0816 * Branching_Ratio / (E_gamma**5 * Tau)
print(f"15/2− (658.7 keV): B(E2) = {BE2:.2f} e²b²")

print("\nBAND B:")
print("-------")

# 19/2+ → 15/2+ (877.0 keV)
Branching_Ratio = ufloat(0.94, 0.01)  # 94% branching ratio
E_gamma = ufloat(0.8770, 0.001)  # 877.0 keV in MeV
Tau = ufloat(1.02, 0.235)  # τAve=1.02 ps {I+24-23} → use average uncertainty (24+23)/2=23.5
BE2 = 0.0816 * Branching_Ratio / (E_gamma**5 * Tau)
print(f"19/2+ (877.0 keV): B(E2) = {BE2:.2f} e²b²")

# 17/2+ → 11/2+ (610.0 keV)
Branching_Ratio = ufloat(0.82, 0.01)  # 82% branching ratio
E_gamma = ufloat(0.6100, 0.001)  # 610.0 keV in MeV
Tau = ufloat(1.34, 0.185)  # τAve=1.34 ps {I+17-20} → use average uncertainty (17+20)/2=18.5
BE2 = 0.0816 * Branching_Ratio / (E_gamma**5 * Tau)
print(f"17/2+ (610.0 keV): B(E2) = {BE2:.2f} e²b²")

# 15/2+ → 11/2+ (763.3 keV)
Branching_Ratio = ufloat(0.94, 0.01)  # 94% branching ratio
E_gamma = ufloat(0.7633, 0.001)  # 763.3 keV in MeV
Tau = ufloat(0.79, 0.075)  # τAve=0.79 ps {I+6-9} → use average uncertainty (6+9)/2=7.5
BE2 = 0.0816 * Branching_Ratio / (E_gamma**5 * Tau)
print(f"15/2+ (763.3 keV): B(E2) = {BE2:.2f} e²b²")

# 11/2+ → 7/2+ (659.0 keV)
Branching_Ratio = ufloat(1.00, 0.01)  # 100% branching ratio
E_gamma = ufloat(0.6590, 0.001)  # 659.0 keV in MeV
Tau = ufloat(1.42, 0.105)  # τAve=1.42 ps {I+10-11} → use average uncertainty (10+11)/2=10.5
BE2 = 0.0816 * Branching_Ratio / (E_gamma**5 * Tau)
print(f"11/2+ (659.0 keV): B(E2) = {BE2:.2f} e²b²")

print("\n" + "=" * 50)
print("Note: Uncertainties include systematic errors from:")
print("- Branching ratios (estimated ±1-2%)")
print("- Gamma energies (assumed ±1 keV)")
print("- Lifetimes (from ENSDF asymmetric uncertainties)")
print("=" * 50)
