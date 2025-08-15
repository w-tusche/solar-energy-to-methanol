# -*- coding: utf-8 -*-
import os

import numpy as np
import pandas as pd
import scipy

import material_props as props
import photorespiration as Presp
import utils

# Constants
h = scipy.constants.h  # Planck's constant h = 6.626e-34 (m^2 kg / s) = (J*s)
c = scipy.constants.c  # Speed of light in vacuum c = 3e8 (m/s)
# Boltzmann constant k = 1.38e-23 (m^2 kg / s^2 K) = (J/K)
k = scipy.constants.k
q = scipy.constants.elementary_charge  # elementary charge
eV = scipy.constants.eV  # Electron Volt (J)

standard_enthalpy_of_formation_25C = props.standard_enthalpy_of_formation_25C


def import_spectrum():
    # Wavelength [nm], Solar spectrum ASTMG173 direct:  [W*m-2*nm-1], sq_l_object
    _, _, SQ_ASTMG = utils.create_shockley_queisser_limit_object_WL_and_irradiation()

    df_spectrum = pd.DataFrame(
        {
            "WL": SQ_ASTMG.WL_equally_spaced,
            "Power": SQ_ASTMG.AM15_WL,
            "Flux": SQ_ASTMG.AM15flux_WL,
        }
    )
    return df_spectrum


def photosynthetically_active_spectrum(
    df_spectrum, range_total_spectrum, range_photosynthesis_spectrum
):
    """Calcualte fraction of power and flux of a spectrum falling within a certain range

    Args:
        df_spectrum (pandas.DataFrame): df of spectrum with WL, Power and Flux
        range_total_spectrum (list): total range of the spectrum as reference
        range_photosynthesis_spectrum (list): specific range in spectrum to compare

    Returns:
        tuple of floats: power and flux of total spectrum, specific part and their fraction
    """

    # Creation of dataframes
    pow_per_WL = df_spectrum[["WL", "Power"]]
    flux_per_WL = df_spectrum[["WL", "Flux"]]

    # Calculations
    power_total_spectrum = utils.trapz_integration_of_area(
        pow_per_WL, *range_total_spectrum
    )  # [W m-2]
    flux_total_spectrum = utils.trapz_integration_of_area(
        flux_per_WL, *range_total_spectrum
    )  # [photons s-1 m-2]
    power_photosynthesis_spectrum = utils.trapz_integration_of_area(
        pow_per_WL, *range_photosynthesis_spectrum
    )  # [W m-2]
    flux_photosynthesis_spectrum = utils.trapz_integration_of_area(
        flux_per_WL, *range_photosynthesis_spectrum
    )  # [photons s-1 m-2]
    fraction_power_remaining = power_photosynthesis_spectrum / power_total_spectrum
    fraction_flow_remaining = flux_photosynthesis_spectrum / flux_total_spectrum

    # Compare http://astro1.panet.utoledo.edu/~relling2/teach/archives/4400.2013/20130124_PHYS4400_Lecture_5%20-%20Integrating%20the%20Solar%20Spectrum.pdf

    # About 1000 W/m2 in range 280 - 4000 nm
    # https://www.pveducation.org/pvcdrom/appendices/standard-solar-spectra
    # https://img.antpedia.com/standard/files/pdfs_ora/20210202/ASTM%20G173-03(2020).pdf
    print(
        f"Power in range from {range_total_spectrum[0]} to {range_total_spectrum[1]} nm: {power_total_spectrum:.2f} (W m-2)"
    )
    print(
        f"Power in range from {range_photosynthesis_spectrum[0]} to {range_photosynthesis_spectrum[1]} nm: {power_photosynthesis_spectrum:.2f} (W m-2)"
    )
    print(
        f"Power fraction photosynthesis of total range: {fraction_power_remaining:.4f}"
    )
    print(
        f"Flux in range from {range_total_spectrum[0]} to {range_total_spectrum[1]}: {flux_total_spectrum:.3e} (photons s-1 m-2)"
    )
    print(
        f"Flux in range from {range_photosynthesis_spectrum[0]} to {range_photosynthesis_spectrum[1]}: {flux_photosynthesis_spectrum:.3e} (photons s-1 m-2)"
    )
    print(f"Flux fraction photosynthesis of total range: {fraction_flow_remaining:.2f}")

    return fraction_power_remaining


def mcCree(df_spectrum, mcCree_min, mcCree_max):
    # Creation of dataframes
    pow_per_WL = df_spectrum[["WL", "Power"]]
    flux_per_WL = df_spectrum[["WL", "Flux"]]

    try:
        mcCree_absorptance = pd.read_csv(r"./data/raw/McCree_Absorptance.csv", sep=";")
    except:
        mcCree_absorptance = pd.read_csv(r"../data/raw/McCree_Absorptance.csv", sep=";")

    # Interpolate data to fit 0.5 nm steps from the rest of the calculation
    mcCree_spl = scipy.interpolate.CubicSpline(
        mcCree_absorptance["wavelength_nm"], mcCree_absorptance["Absorptance_percent"]
    )

    # New range to express the mcCree data on
    mcCree_WL = np.arange(mcCree_min, mcCree_max + 0.1, 0.5)  # [nm]
    mcCree_interp = mcCree_spl(mcCree_WL)  # [nm]

    ### Calculate losses
    df_mcCree = pd.DataFrame()
    df_mcCree = pow_per_WL.loc[
        (pow_per_WL["WL"] >= mcCree_min) & (pow_per_WL["WL"] <= mcCree_max)
    ].copy()
    df_mcCree["Flux"] = flux_per_WL.loc[
        (pow_per_WL["WL"] >= mcCree_min) & (pow_per_WL["WL"] <= mcCree_max)
    ]["Flux"]

    df_mcCree["Absorption"] = mcCree_interp
    df_mcCree["Power_absorbed"] = (
        df_mcCree["Power"] * df_mcCree["Absorption"]
    )  # [W m-2]
    df_mcCree["Flux_absorbed"] = (
        df_mcCree["Flux"] * df_mcCree["Absorption"]
    )  # [photons s-1 m-2]
    mcCree_efficiency = df_mcCree["Power_absorbed"].sum() / df_mcCree["Power"].sum()

    print(
        f'McCree "efficiency" in the range {mcCree_min} to {mcCree_max} is: {mcCree_efficiency:.3f}'
    )

    return mcCree_efficiency, df_mcCree, mcCree_absorptance


def photochemical_efficiency(df):
    #### Photochemical efficiency
    # Calculate energy absorbed in Photosynthetically active centrum I (BG = 680nm) and centrum II (BG = 700nm)
    # band gaps at 680 nm and 700 nm
    # 400nm < x <= 680 nm  is absorbed at band gap 680 nm
    # 680nm < x <= 700 nm is absorbed at band gap 700 nm

    # 1. Integrate photon flux from 400 to 680 nm; multiply by energy at 680nm
    # 2. Integrate photon flux from 680 to 700 nm; multiply by energy at 700nm
    # 3. Sum up

    wvl_PSI = 680.0  # [nm]
    wvl_PSII = 700.0  # [nm]
    wvl_lower_bound = 400.0  # [nm]
    photon_fluxes = []  # Store flux for PSI and PSII
    energies = []  # Store flux for PSI and PSII
    for wvl_bg in [wvl_PSI, wvl_PSII]:
        photon_flux = utils.trapz_integration_of_area(
            df[["WL", "Flux"]], wvl_lower_bound, wvl_bg
        )  # [photons s-1 m-2]
        photon_fluxes.append(photon_flux)
        bandgap_energy = utils.nm_to_eV(wvl_bg) * scipy.constants.eV
        energy = photon_flux * bandgap_energy
        energies.append(energy)
        wvl_lower_bound = wvl_bg

    # Total energy in photosynthetically active spectrum (400nm - 700nm)
    energy_tot_400_700 = utils.trapz_integration_of_area(
        df[["WL", "Power"]], 400.0, 700.0
    )

    photons_tot_400_700 = utils.trapz_integration_of_area(
        df[["WL", "Flux"]], 400.0, 700.0
    )

    # Caluclate efficency: (sum of energy per centrum)/total energy in spectrum slice
    efficiency_400_700 = sum(energies) / energy_tot_400_700

    # Avg. photon energy
    avg_energy = (energy_tot_400_700 / photons_tot_400_700) * scipy.constants.Avogadro

    # Calculate average photon energy
    avg_photon_energy = sum(energies) / sum(photon_fluxes)  # in (J)

    print(
        f"Photochemical efficiency (after McCree): Efficiency {efficiency_400_700:.4f}; avg. photon energy {avg_photon_energy/scipy.constants.eV:.2f} eV; in kJ/mol {avg_photon_energy * scipy.constants.Avogadro/1000:.2f}."
    )
    print(
        f"Mean energy of photon in PAR range (after McCree): {avg_energy/1000:.2f} in (kJ/mol)"
    )
    print(
        f"Difference mean PAR, mean absorbed: {avg_energy/1000 - (avg_photon_energy * scipy.constants.Avogadro/1000):.2f}"
    )

    return efficiency_400_700, avg_photon_energy


def carbohydrate_synthesis(avg_photon_energy_mcCree):
    # Carbohydrate Synthesis
    avg_photon_energy_mcCree_kJ_mol = (
        avg_photon_energy_mcCree / 1000
    ) * scipy.constants.Avogadro

    # https://en.wikipedia.org/wiki/Standard_enthalpy_of_formation
    # https://www.engineeringtoolbox.com/standard-state-enthalpy-formation-definition-value-Gibbs-free-energy-entropy-molar-heat-capacity-d_1978.html

    # delta H Photosynthesis
    # 6 CO2 + 12 H2O -> C6H12O6 + 6 O2 + 6H2O
    delta_H = (
        standard_enthalpy_of_formation_25C["C6H12O6(s)"]
        + 6 * standard_enthalpy_of_formation_25C["O2(g)"]
        + 6 * standard_enthalpy_of_formation_25C["H2O(l)"]
    ) - (
        6 * standard_enthalpy_of_formation_25C["CO2(g)"]
        + 12 * standard_enthalpy_of_formation_25C["H2O(l)"]
    )  # kJ/mol

    # Required photons:
    p_req = 48  # mol
    # energy_of_absorbed_centers = 174  # kJ/mol
    energy_of_absorbed_centers = avg_photon_energy_mcCree_kJ_mol  # kJ/mol

    energy_photons_absorbed_total = p_req * energy_of_absorbed_centers
    efficiency_carbohydrate_synthesis = delta_H / energy_photons_absorbed_total

    print(f"Avg. Photon Energy kJ/mol: {avg_photon_energy_mcCree_kJ_mol:.2f}")
    print(rf"Delta H = {delta_H:.2f} kJ/mol")
    print(
        f"Efficiency of carbohydrate synthesis: {efficiency_carbohydrate_synthesis:.4f}"
    )

    return efficiency_carbohydrate_synthesis


def photorespiration():
    Temperature_photorespiration_degC = 15  # °C
    o2_partial_pressure = 210  # O2 in air ~ 20.95 Vol.-%
    co2_ppm_in_air = props.y_CO2_air * 1e6  # CO2 in air ~ 420 ppm import
    factor_co2_concentration_at_leaf = 0.6
    co2_partial_pressure = (
        co2_ppm_in_air * factor_co2_concentration_at_leaf / 1000
    )  # ~ 0.252 at 420 ppm with 0.6 as factor

    PR = Presp.photorespiration(
        T=Temperature_photorespiration_degC,
        o2_partial_pressure=o2_partial_pressure,
        co2_partial_pressure=co2_partial_pressure,
    )
    d_pr = PR.d_pr

    return d_pr


def compute_photosynthesis_efficiency() -> pd.DataFrame:
    df_spectrum = import_spectrum()

    # Fraction of power and flux of the ASTM G173 1.5G spectrum falling within certain range
    range_total_spectrum = [280.0, 4000.0]
    range_photosynthesis_spectrum = [400.0, 700.0]
    fraction_power_remaining = photosynthetically_active_spectrum(
        df_spectrum, range_total_spectrum, range_photosynthesis_spectrum
    )

    # Calculation of mcCree efficiency
    # Definition of mcCree range
    mcCree_min = 400.0  # (nm)
    mcCree_max = 700.0  # (nm)
    mcCree_efficiency, df_mcCree, mcCree_absorptance = mcCree(
        df_spectrum, mcCree_min, mcCree_max
    )

    # Create Dataframe with McCree absorbed photons and energy for photosynthetically active spectrum
    df_photosynthetic_active_spectrum = pd.DataFrame(
        {
            "WL": df_mcCree["WL"],
            "Power": df_mcCree["Power_absorbed"],
            "Flux": df_mcCree["Flux_absorbed"],
        }
    )

    # Photochemical efficiency (incl. McCree -> using values after McCree)
    (
        efficiency_400_700_mcCree,
        avg_photon_energy_mcCree,
    ) = photochemical_efficiency(df_photosynthetic_active_spectrum)

    # Carbohydrate synthesis
    efficiency_carbohydrate_synthesis = carbohydrate_synthesis(avg_photon_energy_mcCree)

    # Photorespiration
    d_pr = photorespiration()

    # Respiration
    respiration_efficiency = 0.3

    # Initialize the efficiencies for the output table
    step_efficiencies = {
        "Photosynthetically_active_spectrum": fraction_power_remaining,
        "McCree_absorption": mcCree_efficiency,
        "Photochemical_efficiency": efficiency_400_700_mcCree,
        "Carbohydrate_synthesis": efficiency_carbohydrate_synthesis,
        "Photorespiration": 1 - d_pr,
        "Respiration": 1 - respiration_efficiency,
    }

    df_photosynthesis_efficiency = utils.table_creation(
        "Photosynthesis", step_efficiencies
    )

    return df_photosynthesis_efficiency


def main(output_dir="./data/results/"):
    print("\nPhotosynthesis efficiency calculation:\n")
    df_photosynthesis_efficiency = compute_photosynthesis_efficiency()
    # Write to new standard results path
    print(df_photosynthesis_efficiency)
    df_photosynthesis_efficiency.to_csv(
        os.path.join(output_dir, "photosynthesis_eff.csv"), sep=";"
    )


if __name__ == "__main__":
    main()
