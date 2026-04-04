import math

# Alpha-particle energies (keV) extracted from Tables I and II of 2026ZhAA CP10868
# Filtered to include ONLY full-energy deposits. 
# Summed (footnote 'a') and Escaped (footnote 'b') events are excluded.
# Uncertainties: FWHM=30 keV. Authors report total uncertainty=SQRT((FWHM/2.355)^2/N+(ΔEcal)^2), where N is the number of events, and ΔEcal is the calibration uncertainty (not provided).

FWHM_KEV = 30.0
FWHM_TO_SIGMA = 2.355
MEAN_REPRODUCTION_TOLERANCE_KEV = 1.0

alpha_energy_lab_220np = {
    "220Np": [10083, 10091, 10073, 10093, 10065, 10075, 10074, 10071], # Final result: 10078(16)
    "216Pa": [7952, 7969, 7982, 7953, 7961, 7931], # Final result: 7958(9)
    "212Ac": [7367, 7314, 7377, 7367, 7380], # Final result: 7361(9)
    "208Fr": [6622, 6613, 6620, 6634, 6641, 6653] # Final result: 6632(6)
}
alpha_energy_lab_219np = {
    "219Np": [8948, 9014, 9032, 9066], # Final result: 9015(16)
    "215Pa": [8061, 8091, 8119, 8049, 8113], # Final result: 8087(12)
    "211Ac": [7467, 7473, 7401, 7454, 7490], # Final result: 7451(9)
    "207Fr": [6752, 6766, 6769] # Final result: 6762(8)
}

REPORTED_RESULTS = {
    "220Np": {"mean": 10078, "uncertainty": 16},
    "216Pa": {"mean": 7958, "uncertainty": 9},
    "212Ac": {"mean": 7361, "uncertainty": 9},
    "208Fr": {"mean": 6632, "uncertainty": 6},
    "219Np": {"mean": 9015, "uncertainty": 16},
    "215Pa": {"mean": 8087, "uncertainty": 12},
    "211Ac": {"mean": 7451, "uncertainty": 9},
    "207Fr": {"mean": 6762, "uncertainty": 8},
}

def calculate_mean_energy(energies):
    """Return the arithmetic mean alpha energy in keV."""
    count = len(energies)
    total_sum = sum(energies)
    return total_sum / count if count > 0 else 0.0


def calculate_counting_uncertainty(count, fwhm_kev=FWHM_KEV):
    """Return the event-count contribution to the uncertainty in keV."""
    if count <= 0:
        return 0.0
    return (fwhm_kev / FWHM_TO_SIGMA) / math.sqrt(count)


def calculate_total_uncertainty(count, delta_ecal_kev, fwhm_kev=FWHM_KEV):
    """Return the total uncertainty in keV for a chosen calibration uncertainty."""
    counting_uncertainty = calculate_counting_uncertainty(count, fwhm_kev=fwhm_kev)
    return math.sqrt(counting_uncertainty ** 2 + delta_ecal_kev ** 2)


def calculate_implied_delta_ecal(count, target_total_uncertainty, fwhm_kev=FWHM_KEV):
    """Return the calibration uncertainty needed to reproduce a reported total uncertainty."""
    counting_uncertainty = calculate_counting_uncertainty(count, fwhm_kev=fwhm_kev)
    remaining_variance = target_total_uncertainty ** 2 - counting_uncertainty ** 2
    if remaining_variance < 0:
        raise ValueError(
            "reported total uncertainty is smaller than the counting term from FWHM and N"
        )
    return math.sqrt(remaining_variance)


def format_mean_reproduction(mean_val, reported_mean, tolerance_kev=MEAN_REPRODUCTION_TOLERANCE_KEV):
    """Return a short statement describing whether the quoted mean is reproduced."""
    mean_offset = mean_val - reported_mean
    if abs(mean_offset) <= tolerance_kev:
        return "Eα reproduced within 1 keV"
    return f"Eα not reproduced within 1 keV (listed-event mean={mean_val:.2f} keV, offset={mean_offset:+.2f} keV)"


def analyze_nuclide_energies(nuclide, energies):
    """Return a concise report line for one nuclide."""
    count = len(energies)
    mean_val = calculate_mean_energy(energies)
    reported_result = REPORTED_RESULTS.get(nuclide)

    if reported_result is None:
        return f'"{nuclide}": N: {count}, no reported final result available.'

    reported_mean = reported_result["mean"]
    reported_uncertainty = reported_result["uncertainty"]
    implied_delta_ecal = calculate_implied_delta_ecal(count, reported_uncertainty)
    mean_reproduction = format_mean_reproduction(mean_val, reported_mean)

    return (
        f'"{nuclide}": Eα: {reported_mean}, ΔEα: {reported_uncertainty}, N: {count}, '
        f"{mean_reproduction}, ΔEα can be reproduced with ΔEcal={implied_delta_ecal:.2f} keV."
    )


def analyze_dataset(dataset_name, dataset):
    """Print concise report lines for every nuclide in one decay chain."""
    print(dataset_name)
    print("=" * len(dataset_name))
    for nuclide, energies in dataset.items():
        print(analyze_nuclide_energies(nuclide, energies))


def main():
    """Run the alpha-energy analysis for both decay chains."""
    print("Alpha-Energy Reproduction Summary")
    print("=" * 32)
    print(f"FWHM = {FWHM_KEV:.1f} keV")
    print("Model: total uncertainty = sqrt(((FWHM/2.355)^2)/N + (ΔEcal)^2)")
    print("Each line reports the quoted Eα, quoted ΔEα, event count N, mean-reproduction status, and implied ΔEcal.")
    print()

    analyze_dataset("220Np decay chain", alpha_energy_lab_220np)
    print()
    analyze_dataset("219Np decay chain", alpha_energy_lab_219np)


if __name__ == "__main__":
    main()