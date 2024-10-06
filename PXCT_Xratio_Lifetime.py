import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogFormatter

# Data as a multiline string
data = """
Ep	Ratio_Cu_Zn_min	Ratio_Cu_Zn_max	Lifetime_min	Lifetime_max
825	1.0528E-17	1.4634E-01	2.7743E+00	3.8564E+16
1300	6.0027E-02	2.8241E-01	1.4376E+00	6.7636E+00
1700	3.0179E-01	8.4569E-01	4.8008E-01	1.3453E+00
2100	1.6290E+00	2.5605E+00	1.5856E-01	2.4924E-01
2500	1.0823E+00	1.8735E+00	2.1671E-01	3.7514E-01
2900	2.5871E+00	6.4902E+00	6.2556E-02	1.5693E-01
3300	3.5106E+00	1.1390E+01	3.5645E-02	1.1565E-01
3700	4.4398E+00	8.6554E+00	4.6907E-02	9.1445E-02
4100	5.6108E+00	1.6507E+01	2.4596E-02	7.2361E-02
4500	4.7293E+00	1.2883E+01	3.1514E-02	8.5848E-02
4900	6.6937E+00	2.0159E+01	2.0140E-02	6.0654E-02
5300	5.1136E+00	1.1188E+01	3.6288E-02	7.9396E-02
5700	3.8970E+00	7.9412E+00	5.1126E-02	1.0418E-01
6100	8.2830E+00	3.9564E+01	1.0262E-02	4.9016E-02
6600	6.0998E+00	3.6548E+01	1.1109E-02	6.6560E-02
7200	4.8533E+00	1.7217E+16	2.3582E-17	8.3655E-02
"""

# Global settings for plots
plt.rcParams['axes.linewidth'] = 4.0
plt.rcParams['font.size'] = 60
plt.rcParams['font.family'] = ['Times New Roman', 'Georgia', 'Cambria', 'Courier New', 'serif']
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Times New Roman'
plt.rcParams['mathtext.it'] = 'Times New Roman:italic'
plt.rcParams['legend.frameon'] = False

# Parse the data
Ep = []
Ratio_Cu_Zn_min = []
Ratio_Cu_Zn_max = []
Lifetime_min = []
Lifetime_max = []

# Split the data into lines
lines = data.strip().split('\n')
header = lines[0].split()  # Extract header (optional)
for line in lines[1:]:
    tokens = line.split()
    if len(tokens) == 5:
        Ep_value = float(tokens[0])
        Ratio_min = float(tokens[1])
        Ratio_max = float(tokens[2])
        Lifetime_min_value = float(tokens[3])
        Lifetime_max_value = float(tokens[4])
        
        Ep.append(Ep_value)
        Ratio_Cu_Zn_min.append(Ratio_min)
        Ratio_Cu_Zn_max.append(Ratio_max)
        Lifetime_min.append(Lifetime_min_value)
        Lifetime_max.append(Lifetime_max_value)

# Convert lists to numpy arrays for calculations
Ep = np.array(Ep)
Ratio_Cu_Zn_min = np.array(Ratio_Cu_Zn_min)
Ratio_Cu_Zn_max = np.array(Ratio_Cu_Zn_max)
Lifetime_min = np.array(Lifetime_min)
Lifetime_max = np.array(Lifetime_max)

# Calculate central values and uncertainties
# For Ratio_Cu_Zn
Ratio_Cu_Zn = (Ratio_Cu_Zn_min + Ratio_Cu_Zn_max) / 2
Ratio_Cu_Zn_unc_lower = Ratio_Cu_Zn - Ratio_Cu_Zn_min
Ratio_Cu_Zn_unc_upper = Ratio_Cu_Zn_max - Ratio_Cu_Zn

# For Lifetime
Lifetime = (Lifetime_min + Lifetime_max) / 2
Lifetime_unc_lower = Lifetime - Lifetime_min
Lifetime_unc_upper = Lifetime_max - Lifetime

# Plotting
fig, ax1 = plt.subplots(figsize=(22, 11))
plt.subplots_adjust(top=0.95, bottom=0.19, left=0.13, right=0.87)

color = 'black'
ax1.set_xlabel('Proton Energy (keV)', fontsize=60, color=color, labelpad=26)
ax1.set_ylabel(r'$R_{\mathrm{Cu/Zn}}$', color=color, fontsize=60, labelpad=26)

# Plot Ratio_Cu_Zn with asymmetric error bars
ax1.errorbar(
    Ep,
    Ratio_Cu_Zn,
    yerr=[Ratio_Cu_Zn_unc_lower, Ratio_Cu_Zn_unc_upper],
    fmt='none',         # Do not plot markers
    linestyle='None',
    ecolor='black',
    capsize=11,        # how wide are the caps
    capthick=5,         # Thickness of the caps
    elinewidth=5,       # Thickness of the error bar lines
    label=r'$R_{\mathrm{Cu/Zn}}$'
)

ax1.set_yscale('log')  # Set Y-axis to logarithmic scale

# Adjust the tick parameters
tick_width = 4  # Thickness of the ticks
tick_length = 15  # Length of the ticks
ax1.tick_params(axis='both', which='major', direction='in', length=tick_length, width=tick_width, labelsize=60, colors=color, pad=10)
ax1.set_xlim(left=0, right=max(Ep) + 800)
ax1.set_ylim(bottom=2e-3, top=5e2)

# Create a secondary Y-axis that shares the same X-axis
def ratio_to_lifetime(ratio):
    return 0.406 / ratio

def lifetime_to_ratio(lifetime):
    return 0.406 / lifetime


# Desired lifetime tick values from bottom to top
lifetime_ticks = [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2]  # Increasing order
# Corresponding ratio tick values
ratio_ticks = [lifetime_to_ratio(lt) for lt in lifetime_ticks]

# No need to reverse the ticks
# ratio_ticks = ratio_ticks[::-1]
# lifetime_ticks = lifetime_ticks[::-1]

# Create the secondary Y-axis
ax2 = ax1.twinx()
ax2.set_yscale('log')
ax2.set_ylim(ax1.get_ylim())
ax2.set_yticks(ratio_ticks)
lifetime_labels = [f"$10^{{{int(np.log10(lt))}}}$" for lt in lifetime_ticks]
ax2.set_yticklabels(lifetime_labels)
# The secondary Y-axis (ax2) uses ratio_ticks for positions but labels them with lifetime_ticks.

ax2.tick_params(axis='both', which='major', direction='in', length=tick_length, width=tick_width, labelsize=60, colors=color, pad=10)

ax2.set_ylabel('Lifetime (fs)', color=color, fontsize=60, labelpad=30)

# Enable minor ticks
# ax1.minorticks_on()
# ax2.minorticks_on()
# Adjust minor ticks
# ax1.tick_params(axis='both', which='minor', direction='in', length=tick_length / 2, width=tick_width, labelsize=60, colors=color, pad=10)
# ax2.tick_params(axis='both', which='minor', direction='in', length=tick_length / 2, width=tick_width, labelsize=60, colors=color, pad=10)
# Add grid lines
ax1.grid(which='major', linestyle='--', linewidth=1)
# ax1.grid(which='minor', linestyle='--', linewidth=1)

plt.show()
# Save the figure
# plt.savefig(r'F:\e21010\pxct\Fig_PXCT_60Ga_Xratio.png')
# plt.savefig(r'F:\e21010\pxct\Fig_PXCT_60Ga_Xratio.eps')
