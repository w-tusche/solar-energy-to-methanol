# -*- coding: utf-8 -*-
import os

import numpy as np
import pandas as pd

import utils


def create_blackbody_spectrum(T):
    """Create a table with a black body spectrum (Irradiation per Wavelength)
    Irradiation on earth

    Args:
        T (int,float): Temperature in K for the blackbody spectrum
    """

    # [m] Wavelengths from 4 nm to 4000 nm (took 4nm because of np.exp() function not working with 1)
    # np.linspace(4e-9, 4e-6, 3000)
    wavelengths_photovoltaic_nm = np.arange(start=4, stop=4000.1, step=0.5)
    wavelengths_photovoltaic = wavelengths_photovoltaic_nm * 1e-9

    # Calculate intensities
    spectral_radiance = np.array(
        list(
            map(
                lambda wvl: utils.planck_radiation_law(wvl, T), wavelengths_photovoltaic
            )
        )
    )
    spectral_irradiance = np.array(list(map(utils.irradiance_earth, spectral_radiance)))
    black_body_spectrum = pd.DataFrame(
        {
            "Wvlgth nm": wavelengths_photovoltaic_nm,
            "Etr W*m-2*nm-1": spectral_irradiance,
        }
    )

    return black_body_spectrum


def solar_cell_losses_calculation(SQ_ASTMG):
    """Calculate losses in solar cell and create a data frame with processes partial efficiency and the overall efficiency


    Returns:
        pd.DataFrame: _description_
    """

    # Calculate losses
    losses = SQ_ASTMG.E_loss(1.12)  # band gap energy 1.12 eV of Silicon
    transmission_loss = losses["Not Absorbed"].sum()
    thermalization_loss = losses["Thermalization Loss"].sum()
    extracttion_loss = losses[
        "Extraction Loss"
    ].sum()  # This is higher in the cell with updated detailed balance.
    available_energy = losses["Available"].sum()

    total_energy = (
        thermalization_loss + extracttion_loss + transmission_loss + available_energy
    )

    # Literature value: Schäfer 2018: Accurate Calculation of the Absorptance Enhances
    # Efficiency Limit of Crystalline Silicon Solar Cells With Lambertian Light Trapping
    # https://doi.org/10.1109/JPHOTOV.2018.2824024
    current_detailed_balance_limit = 0.2956

    # Calculate efficiencies
    transmission_step_efficiency = 1 - (transmission_loss / total_energy)
    thermalization_step_efficiency = 1 - (
        thermalization_loss / (total_energy - transmission_loss)
    )

    # Sum all further losses under the extraction losses using the current detailed balance limit
    extract_step_efficiency = current_detailed_balance_limit / (
        transmission_step_efficiency * thermalization_step_efficiency
    )

    step_efficiencies = {
        "Radiation_absorption": transmission_step_efficiency,  # Transmission step efficiency (absorbed part above band gap)
        "Thermalization": thermalization_step_efficiency,  # Thermalization step efficiency (above band gap thermalizes to band gap)
        "Power_extraction": extract_step_efficiency,  # Extraction step efficiency (further losses (Auger, Free Carrier Absorption etc.))
    }

    df_PV_efficiency = utils.table_creation("PV", step_efficiencies)

    return df_PV_efficiency


def main(output_dir="./data/results"):
    _, _, SQ_ASTMG = utils.create_shockley_queisser_limit_object_WL_and_irradiation()
    df_PV_efficiency = solar_cell_losses_calculation(SQ_ASTMG)
    print(f"Photovoltaic effciency {df_PV_efficiency}")
    df_PV_efficiency.to_csv(os.path.join(output_dir, "PV_eff.csv"), sep=";")

    # Create and store blackbody spectrum on earth with 5772 K sun temperature
    T = 5772  # Black body temperature (Kelvin)
    black_body_spectrum = create_blackbody_spectrum(T)
    # black_body_spectrum.to_csv(rf"./data/raw/solar_ref_{T}.csv", sep=";", index=None)


if __name__ == "__main__":
    main()
