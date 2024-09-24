import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogFormatter

# Data as a multiline string
data = """
Ep      600     Cu_counts       1.55    8.2985  Zn_counts       28.55   10.7641 Real_Cu_counts  1.73689     9.29909 Real_Zn_counts  28.55   10.7641 Ratio_Cu_Zn     0.0608368       0.326519        Lifetime_proton_emitting_state      6.67359 35.818
Ep      1300    Cu_counts       12.9    12.0438 Zn_counts       75.9    15.782  Real_Cu_counts  14.4554     13.496  Real_Zn_counts  75.9    15.782  Ratio_Cu_Zn     0.190453        0.182169        Lifetime_proton_emitting_state      2.13175 2.03903
Ep      1500    Cu_counts       8.65    10.0978 Zn_counts       51.65   13.1646 Real_Cu_counts  9.69297     11.3154 Real_Zn_counts  51.65   13.1646 Ratio_Cu_Zn     0.187666        0.224239        Lifetime_proton_emitting_state      2.16341 2.58502
Ep      1700    Cu_counts       14.325  9.86102 Zn_counts       32.325  11.2752 Real_Cu_counts  16.0522     11.05   Real_Zn_counts  32.325  11.2752 Ratio_Cu_Zn     0.496589        0.383221        Lifetime_proton_emitting_state      0.817578        0.630931
Ep      1900    Cu_counts       38.7    12.1867 Zn_counts       22.7    11.063  Real_Cu_counts  43.3662 13.6561 Real_Zn_counts  22.7    11.063  Ratio_Cu_Zn     1.91041 1.1085  Lifetime_proton_emitting_state      0.21252 0.123313
Ep      2100    Cu_counts       255.1   30.0993 Zn_counts       141.1   26.8744 Real_Cu_counts  285.859 33.7285 Real_Zn_counts  141.1   26.8744 Ratio_Cu_Zn     2.02593 0.453908        Lifetime_proton_emitting_state      0.200402        0.0448999
Ep      2500    Cu_counts       124.85  21.9535 Zn_counts       97.85   20.9395 Real_Cu_counts  139.904 24.6006 Real_Zn_counts  97.85   20.9395 Ratio_Cu_Zn     1.42978 0.396009        Lifetime_proton_emitting_state      0.28396 0.0786492
Ep      2900    Cu_counts       113.175 18.0668 Zn_counts       34.175  14.2079 Real_Cu_counts  126.821 20.2452 Real_Zn_counts  34.175  14.2079 Ratio_Cu_Zn     3.71093 1.65261 Lifetime_proton_emitting_state      0.109407        0.0487225
Ep      3300    Cu_counts       97.625  16.1279 Zn_counts       18.625  11.6706 Real_Cu_counts  109.396 18.0725 Real_Zn_counts  18.625  11.6706 Ratio_Cu_Zn     5.87362 3.80623 Lifetime_proton_emitting_state      0.0691227       0.044793
Ep      3700    Cu_counts       238.6   23.3807 Zn_counts       43.6    15.9152 Real_Cu_counts  267.369 26.1998 Real_Zn_counts  43.6    15.9152 Ratio_Cu_Zn     6.13232 2.31773 Lifetime_proton_emitting_state      0.0662066       0.025023
Ep      4100    Cu_counts       405.3   29.6376 Zn_counts       61.3    19.1616 Real_Cu_counts  454.169     33.2112 Real_Zn_counts  61.3    19.1616 Ratio_Cu_Zn     7.40895 2.37847 Lifetime_proton_emitting_state      0.0547986       0.0175918
Ep      5100    Cu_counts       329.725 24.4759 Zn_counts       47.725  14.334  Real_Cu_counts  369.481     27.4271 Real_Zn_counts  47.725  14.334  Ratio_Cu_Zn     7.74189 2.39521 Lifetime_proton_emitting_state      0.052442        0.0162247
Ep      6100    Cu_counts       182.475 16.5272 Zn_counts       26.475  8.57413 Real_Cu_counts  204.477     18.52   Real_Zn_counts  26.475  8.57413 Ratio_Cu_Zn     7.72339 2.59725 Lifetime_proton_emitting_state      0.0525676       0.0176776
Ep      7100    Cu_counts       30.65   6.15937 Zn_counts       2.65    2.32366 Real_Cu_counts  34.3456     6.90204 Real_Zn_counts  2.65    2.32366 Ratio_Cu_Zn     12.9606 11.6592 Lifetime_proton_emitting_state      0.0313257       0.0281801
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
Ratio_Cu_Zn = []
Ratio_Cu_Zn_unc = []
Lifetime = []
Lifetime_unc = []

# Split the data into lines
lines = data.strip().split('\n')
for line in lines:
    tokens = line.split()
    # Find the indices of 'Ep', 'Ratio_Cu_Zn', and 'Lifetime_proton_emitting_state'
    try:
        # Ep value is right after 'Ep'
        Ep_index = tokens.index('Ep') + 1
        Ep_value = float(tokens[Ep_index])
        Ep.append(Ep_value)
        
        # Ratio_Cu_Zn value and uncertainty are right after 'Ratio_Cu_Zn'
        Ratio_index = tokens.index('Ratio_Cu_Zn') + 1
        Ratio_value = float(tokens[Ratio_index])
        Ratio_unc = float(tokens[Ratio_index + 1])
        Ratio_Cu_Zn.append(Ratio_value)
        Ratio_Cu_Zn_unc.append(Ratio_unc)
        
        # Lifetime value and uncertainty are right after 'Lifetime_proton_emitting_state'
        Lifetime_index = tokens.index('Lifetime_proton_emitting_state') + 1
        Lifetime_value = float(tokens[Lifetime_index])
        Lifetime_unc_value = float(tokens[Lifetime_index + 1])
        Lifetime.append(Lifetime_value)
        Lifetime_unc.append(Lifetime_unc_value)
    except (ValueError, IndexError):
        continue

# Convert lists to numpy arrays for calculations
Ep = np.array(Ep)
Ratio_Cu_Zn = np.array(Ratio_Cu_Zn)
Ratio_Cu_Zn_unc = np.array(Ratio_Cu_Zn_unc)
Lifetime = np.array(Lifetime)
Lifetime_unc = np.array(Lifetime_unc)

# Plotting
fig, ax1 = plt.subplots(figsize=(22, 11))
plt.subplots_adjust(top=0.96, bottom=0.19, left=0.13, right=0.87)

color = 'black'
ax1.set_xlabel('Proton Energy (keV)', fontsize=60, color=color, labelpad=26)
ax1.set_ylabel(r'$R_{\mathrm{Cu/Zn}}$', color=color, fontsize=60, labelpad=26)
ax1.errorbar(Ep, Ratio_Cu_Zn, yerr=Ratio_Cu_Zn_unc, fmt='s', linestyle='None', color=color, ecolor='black', capsize=14, capthick=4, label=r'$R_{\mathrm{Cu/Zn}}$', markersize=16, markeredgewidth=4, markeredgecolor='black', linewidth=4)
ax1.set_yscale('log')  # Set Y-axis to logarithmic scale
# Adjust the tick parameters
tick_width = 4  # Thickness of the ticks
tick_length = 15  # Length of the ticks
ax1.tick_params(axis='both', which='major', direction='in', length=tick_length, width=tick_width, labelsize=60, colors=color, pad=10)
# ax1.grid(True, which='both', linestyle='--', linewidth=2.0)
ax1.set_xlim(left=0, right=max(Ep) + 500)
ax1.set_ylim(bottom=min(Ratio_Cu_Zn) / 2, top=max(Ratio_Cu_Zn) * 2)

# Create a secondary Y-axis that shares the same X-axis
def ratio_to_lifetime(ratio):
    return 0.406 / ratio

def lifetime_to_ratio(lifetime):
    return 0.406 / lifetime

# Desired lifetime tick values from bottom to top
lifetime_ticks = [10, 5, 2, 1, 0.5, 0.2, 0.1, 0.02]
# Corresponding ratio tick values
ratio_ticks = [lifetime_to_ratio(lt) for lt in lifetime_ticks]

# Since ax1 (primary Y-axis) goes from bottom to top, we need ratio_ticks in increasing order
ratio_ticks = ratio_ticks[::-1]  # Reverse to get increasing order
lifetime_ticks = lifetime_ticks[::-1]  # Reverse to match the order

# Set the ticks on the primary Y-axis
# ax1.set_yticks(ratio_ticks)
# ax1.set_yticklabels(['{:.2g}'.format(tick) for tick in ratio_ticks])

# Create the secondary Y-axis
ax2 = ax1.twinx()
ax2.set_yscale('log')
ax2.set_ylim(ax1.get_ylim())

ax2.set_yticks(ratio_ticks)
ax2.set_yticklabels([str(tick) for tick in lifetime_ticks])
ax2.tick_params(axis='both', which='major', direction='in', length=tick_length, width=tick_width,
                labelsize=60, colors=color, pad=10)

ax2.set_ylabel('Lifetime (fs)', color=color, fontsize=60, labelpad=30)

# Invert the secondary Y-axis
# ax2.invert_yaxis()

# plt.title('Ratio $R_{Cu/Zn}$ and Lifetime vs Proton Energy $E_p$', fontsize=16)
# plt.show()
plt.savefig(r'D:\X\out\Bayesian_VS\Fig_PXCT_60Ga_Xratio.png')
plt.savefig(r'D:\X\out\Bayesian_VS\Fig_PXCT_60Ga_Xratio.eps')
