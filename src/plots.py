# -*- coding: utf-8 -*-
import os

# import cmcrameri  # To register colormaps for matplotlib e.g. cmc.batlow
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sankeys
import utils


def plot_spectral_irridaiance(output_path):
    import photovoltaic

    plt.rcParams.update({"font.size": 14})
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = "Helvetica"
    plt.rcParams["text.usetex"] = True

    # blackbody spectrum on earth with 5772 K sun temperature
    T = 5772  # Black body temperature (Kelvin)
    # Create and store blackbody spectrum on earth with 5772 K sun temperature
    black_body_spectrum = photovoltaic.create_blackbody_spectrum(T)
    wavelengths_photovoltaic_nm = black_body_spectrum["Wvlgth nm"]
    spectral_irradiance = black_body_spectrum["Etr W*m-2*nm-1"]

    (
        WL_ASTMG,
        solar_per_nm_ASTMG,
        SQ_ASTMG,
    ) = utils.create_shockley_queisser_limit_object_WL_and_irradiation()

    fig, ax = plt.subplots()

    linewidth_plot = 0.7

    ax.plot(
        wavelengths_photovoltaic_nm,
        spectral_irradiance,
        color="k",
        linestyle="--",
        linewidth=linewidth_plot,
        label=f"Blackbody {T} K in " + r"($\mathrm{W~m^{-2}nm^{-1}}$)",
    )
    ax.plot(
        WL_ASTMG,
        solar_per_nm_ASTMG,
        color="b",
        linewidth=linewidth_plot,
        label=r"Spectral irradiance in ($\mathrm{W~m^{-2}nm^{-1}}$)",
    )
    ax2 = ax.twinx()
    p2 = ax2.plot(
        SQ_ASTMG.WL_equally_spaced,
        SQ_ASTMG.AM15flux_WL,
        color="red",
        linewidth=linewidth_plot,
        label=r"Photon flux in ($\mathrm{photons~s^{-1}m^{-2}nm^{-1}}$)",
    )

    ax.set_xlim(-0.1, 4000)
    ax.set_ylim(0, 1.8)
    ax2.set_ylim(0, 5e18)
    ax.set_xlabel(r"Wavelength in ($\mathrm{nm}$)")
    ax.set_ylabel(r"Spectral irradiance in ($\mathrm{W~m^{-2}nm^{-1}}$)")
    ax2.set_ylabel(r"Photon flux in ($\mathrm{photons~s^{-1}m^{-2}nm^{-1}}$)")

    # Create the rainbow background image
    y_min = 0
    y_max = 2
    # Visible range reported differently and changes over the lifetime. https://www.nature.com/articles/eye2015252
    img, extent = utils.prepare_rainbow_background(380, 750, y_min, y_max)

    ax.imshow(img, extent=extent, aspect="auto", alpha=0.5)

    # More info: https://matplotlib.org/stable/gallery/subplots_axes_and_figures/secondary_axis.html
    secax = ax.secondary_xaxis("top", functions=(utils.nm_to_eV, utils.eV_to_nm))
    secax.set_xlabel("Bandgap in (eV)")
    secax.set_xticks([3.5, 2.0, 1.4, 1, 0.7, 0.5, 0.4])
    # Fill area above the curve white!
    # ax.fill_between(wavelengths, spectrum, upper_limit, color="w")
    ax.grid()
    fig.legend(
        loc="upper right",
        bbox_to_anchor=(1, 1),
        bbox_transform=ax.transAxes,
        fontsize=12,
    )

    plt.savefig(os.path.join(output_path, "spectral_irradiance.png"), dpi=1000)


def plot_efficiency_over_band_gap_energy(output_path):
    _, _, SQ_ASTMG = utils.create_shockley_queisser_limit_object_WL_and_irradiation()

    fig, ax = plt.subplots()

    ax.plot(SQ_ASTMG.Es, SQ_ASTMG.PCE, color="black", label="Efficiency")
    ax.grid()
    ax.set_xlabel(r"Band gap in ($\mathrm{eV}$)")
    ax.set_ylabel(r"Power conversion efficiency in ($\%$)")
    ax.set_xlim(0.25, 4)
    ax.set_ylim(0, 40)
    ax.vlines(
        1.12,
        0,
        40,
        colors="red",
        linestyles=":",
        label=r"Band gap silicon $\rm{E}_g = 1.12 \mathrm{eV}$",
    )
    ax.legend()

    secax = ax.secondary_xaxis("top", functions=(utils.eV_to_nm, utils.nm_to_eV))
    secax.set_xticks([400, 500, 600, 800, 1000, 1300, 2000])
    secax.set_xlabel("Wavelength (nm)")

    fig.savefig(os.path.join(output_path, "efficiency_per_band_gap.png"), dpi=300)


def plot_theoretical_minimum_separation_energy(output_path):
    import dac

    dac.plot_figure_w_min_CO2(
        save=True,
        output_path=os.path.join(output_path, "minimum_separation_work_co2.png"),
    )


def plot_figure_w_min_CO2(output_path, save=False):
    import dac
    import material_props

    y_CO2_in_air = material_props.y_CO2_air
    M_CO2 = material_props.molar_mass_g_mol["CO2"]  # g/mol

    # Set global font size
    plt.rcParams.update({"font.size": 14})
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = "Helvetica"
    plt.rcParams["text.usetex"] = True

    fig, ax = plt.subplots()
    fig.set_size_inches(10, 5)

    x = np.linspace(0.000001, 0.99999, 1000)
    y_vals_per_mass = []
    y_vals_per_mol = []
    for y_CO2 in x:
        y_val_per_mass, y_val_per_mol = dac.minimum_work_of_gas_stream_separation(
            y_CO2, 0.99999, 0.99999, M_comp=M_CO2, T=298.15
        )
        y_vals_per_mass.append(y_val_per_mass)
        y_vals_per_mol.append(y_val_per_mol)

    y_vals = y_vals_per_mol

    ax.plot(x, y_vals, color="k", linewidth=2)

    # Calculate minimum work of air separation
    (
        co2_sep_air_per_mass,
        co2_sep_air_per_mol,
    ) = dac.minimum_work_of_gas_stream_separation(y_CO2_in_air, 0.99999, 0.99999)

    # Horizontal line to show energy requirement for air separation
    # Add filled area:
    fill_area = pd.DataFrame({"x": x, "y_bot": y_vals})
    fill_area = fill_area[(fill_area["x"] > 0.07) & (fill_area["x"] < 0.25)]
    ax.fill_between(
        x=fill_area["x"],
        y1=[co2_sep_air_per_mol] * len(fill_area["x"]),
        y2=fill_area["y_bot"],
        interpolate=True,
        color="grey",
        alpha=0.15,
    )
    fill_area = pd.DataFrame({"x": x, "y_bot": y_vals})
    fill_area = fill_area[(fill_area["x"] > 0.25) & (fill_area["x"] < 1)]
    ax.fill_between(
        x=fill_area["x"],
        y1=[co2_sep_air_per_mol] * len(fill_area["x"]),
        y2=fill_area["y_bot"],
        interpolate=True,
        color="black",
        alpha=0.15,
    )

    ax.hlines(
        co2_sep_air_per_mol, -0.1, 1.1, linestyles="dashed", color="r", linewidth=2
    )

    # Annotate air separation with text box
    ax.text(
        y_CO2_in_air + 0.06,
        co2_sep_air_per_mol + 3,
        "Air:\n"
        + str(
            round(y_CO2_in_air * 1e6),
        )
        + r"$~\mathrm{ppm_{CO_2}}$"
        + "\n"
        + str(round(co2_sep_air_per_mol, 2))
        + r" $\mathrm{kJ~mol_{CO_2}^{-1}}$"
        + "\n"
        + str(round(co2_sep_air_per_mass))
        + r" $\mathrm{kJ~kg_{CO_2}^{-1}}$",
        bbox=dict(facecolor="white", edgecolor="grey", boxstyle="round", pad=0.5),
        horizontalalignment="left",
        verticalalignment="bottom",
    )

    # Connect text box and air separation marker with arrow
    p1 = patches.FancyArrowPatch(
        (y_CO2_in_air, co2_sep_air_per_mol),
        (y_CO2_in_air + 0.052, co2_sep_air_per_mol + 2.35),
        arrowstyle="<-",
        mutation_scale=20,
        linewidth=1.5,
    )
    ax.add_patch(p1)

    x_val_for_arrow = 0.5
    y_mass_val_for_arrow, y_val_for_arrow = dac.minimum_work_of_gas_stream_separation(
        x_val_for_arrow, 0.99999, 0.99999
    )

    p2 = patches.FancyArrowPatch(
        (x_val_for_arrow, y_val_for_arrow),
        (x_val_for_arrow, co2_sep_air_per_mol),
        arrowstyle="<->",
        # linestyle="--",
        mutation_scale=20,
        linewidth=3,
        color="r",
    )
    ax.add_patch(p2)

    ax.text(
        x_val_for_arrow + 0.03,
        co2_sep_air_per_mol - 15.5,
        "Reduced separation\nenergy requirement\nfor biogas:\n"
        + str(
            round(x_val_for_arrow * 100),
        )
        + r"$~\mathrm{\%_{CO_2}}$"
        + "\n"
        + str(round(co2_sep_air_per_mol - y_val_for_arrow, 2))
        + r" $\mathrm{kJ~mol_{CO_2}^{-1}}$"
        + "\n"
        + str(round(co2_sep_air_per_mass - y_mass_val_for_arrow))
        + r" $\mathrm{kJ~kg_{CO_2}^{-1}}$",
        bbox=dict(facecolor="white", edgecolor="grey", boxstyle="round", pad=0.5),
        horizontalalignment="left",
        verticalalignment="bottom",
    )

    secax = ax.secondary_yaxis(
        "right", functions=(lambda x: x * (1000 / M_CO2), lambda y: y / (1000 / M_CO2))
    )
    secax.set_yticks([0, 100, 200, 300, 400, 500, 600])

    ax.set_xlim([-0.01, 1.01])
    ax.set_ylim([0, 35])
    ax.set_xlabel(r"CO$_2$ mole fraction ($x_{CO_2}$)")
    ax.set_ylabel(
        "Theoretical minimum work of total\nseparation in "
        + r"($\mathrm{kJ~mol_{CO_2}^{-1}}$)"
    )
    secax.set_ylabel(
        "Theoretical minimum work of total\n separation in "
        + r"($\mathrm{kJ~kg_{CO_2}^{-1}}$)"
    )
    ax.set_xticks(np.array(range(0, 11)) / 10)
    ax.grid()

    if save == True:
        fig.savefig(
            os.path.join(output_path, "biomass_energy_savings.png"),
            dpi=1000,
            bbox_inches="tight",
        )
    else:
        plt.show()


def plot_mcCree_absorption(output_path):
    import photosynthesis

    df_spectrum = photosynthesis.import_spectrum()
    # Calculation of mcCree efficiency
    # Definition of mcCree range
    mcCree_min = 400.0  # (nm)
    mcCree_max = 700.0  # (nm)
    mcCree_efficiency, df_mcCree, mcCree_absorptance = photosynthesis.mcCree(
        df_spectrum, mcCree_min, mcCree_max
    )

    # Plot mcCree efficiency
    # Compare: https://www.sciencedirect.com/science/article/pii/B9780323851527000185
    fig, ax = plt.subplots()
    ax.plot(df_mcCree["WL"], df_mcCree["Absorption"], label="Absorptance")
    ax.plot(
        mcCree_absorptance["wavelength_nm"],
        mcCree_absorptance["Absorptance_percent"],
        "bx",
    )
    ax.set_xlim(mcCree_min, mcCree_max)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"Wavelength in ($\mathrm{nm}$)")
    ax.set_ylabel("Absorption (single leaf)")
    ax.grid()

    fig.savefig(
        os.path.join(output_path, "mcCree_absorption.png"), dpi=300, bbox_inches="tight"
    )


def main(input_path="./data/results/", output_path="./img/"):
    ### Main paper
    plot_spectral_irridaiance(output_path)

    plot_theoretical_minimum_separation_energy(output_path)

    # Grey shading is done in inkscape
    plot_figure_w_min_CO2(output_path, save=True)

    ### Appendix
    plot_efficiency_over_band_gap_energy(output_path)

    # Label positioning finetuned afterwards
    sankeys.main(input_path, output_path)

    ### Additional
    plot_mcCree_absorption(output_path)


if __name__ == "__main__":
    main()
