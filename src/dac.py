# -*- coding: utf-8 -*-
import math

try:
    import material_props as material_props
except:
    import material_props as material_props

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants

y_CO2_in_air = material_props.y_CO2_air
M_CO2 = material_props.molar_mass_g_mol["CO2"]  # g/mol


def w_min_total_separation_one_component(mole_fraction_in_mixture, M=M_CO2, T=298.15):
    """
    Minimum work for 100% separation of a gas from a gas mixture
    assuming ideal gas (e.g., total removal of CO2 from air)
    Source: Socolow, R.; Desmond, M.; et al. (2011):
    Direct Air Capture of CO2 with Chemicals:
    A Technology Assessment for the APS Panel on Public Affairs
    https://api.semanticscholar.org/CorpusID:107482441

    Args:
        mole_fraction_in_mixture (float): mole fraction of gas to separate fully (e.g., CO2) in (0,1)
        M (float): molar mass (g/mol) of component to separate
        T (float): Temperature in Kelvin (defaults to 298.15)
    Returns:
        float: minimum work of separation per unit of mass in (kJ/kg_{CO_2})
    """

    assert (
        0 < mole_fraction_in_mixture < 1
    ), "Mole fraction has to be value between zero and one. (0,1)"

    y = mole_fraction_in_mixture  # fraction (btw.: 0 -> 1)
    R = scipy.constants.R  # 8.314 J/mole-K (gas constant)
    w_min_all = -((R * T) / (y * M)) * (y * math.log(y) + (1 - y) * math.log(1 - y))
    return w_min_all  # kJ/kg


def minimum_work_of_gas_stream_separation(y_A, y_B, CR, M_comp=44.0, T=298.15):
    """
    Minimum work required for separating a component ("comp") from a gas
    mixture for an isothermal and isobaric process (negative of difference in
    Gibbs free energy of separated final states (B,C) from the mixed initial
    state (A).

    Adjusted from:
    W.M. Budzianowski, Energy Efficient Solvents for CO2 Capture
    by Gas-Liquid Absorption, Green Energy and Technology,
    Assessment of Thermodynamic Efficiency of Carbon Dioxide
    Separation in Capture Plants by Using Gas-Liquid Absorption
    DOI 10.1007/978-3-319-47262-1_2
    (Formula 5)
                _______
               |       |----B----> enriched stream
     ----A---->|       |
               |_______|----C----> remainder

    Args:
        y_A (float): mole fraction of gas to separate in ingoing stream (A) (0,1)
        y_B (float): mole fraction of separated gas in enriched stream (B) (0,1)
        CR (float): capture fraction
        M_comp (float, optional): molar mass of component to separate. Defaults to 44 (g/mol) for CO2
        T (float, optional): Temperature in K. Defaults to 298.15.

    Returns:
        tuple: minimum work of separation from a gas mixture (w_min) per kg and mol (kJ/kg),(kJ/mol)
    """

    assert (
        (0 < y_A < 1) & (0 < y_B < 1) & (0 < CR < 1)
    ), "Mole fractions y_A, y_B and capture rate CR must be between 0 and 1."

    R = scipy.constants.R  # 8.314 J/mole-K (gas constant)

    n_A = 100_000.0  # Mole number in stream A

    # Derived variables
    x_A = n_A * y_A  # moles of component to separate in stream A
    x_B = CR * x_A  # moles of component to separate in stream B
    x_C = x_A - x_B  # moles of component to separate in stream C
    n_B = x_B / y_B  # total mole number in stream A
    n_C = n_A - n_B  # total mole number in stream A
    y_C = x_C / n_C  # molar fraction of component to separate in stream C

    # Vectors n and y
    n = np.array([n_A, n_B, n_C])
    y = np.array([y_A, y_B, y_C])

    # Compute the terms per stream
    terms = n * (y * np.log(y) + (1 - y) * np.log(1 - y))

    term = terms[1] + terms[2] - terms[0]
    w_min_spec_mol = (R * T * term) / (x_B)  # (J/mol)
    w_min_spec_mass = w_min_spec_mol / M_comp  # (kJ/kg)

    return (
        w_min_spec_mass,
        w_min_spec_mol / 1000,
    )  # per unit of mass(kJ/kg) and specific per mol (kJ/mol)


def plot_figure_w_min_CO2(save=False, output_path=None):
    # Set global font size
    plt.rcParams.update({"font.size": 14})
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = "Helvetica"
    plt.rcParams["text.usetex"] = True

    fig, ax = plt.subplots()
    fig.set_size_inches(10, 5)

    x = np.linspace(0.00001, 0.9999, 1000)
    y_vals_per_mass = []
    y_vals_per_mol = []
    for y_CO2 in x:
        y_val_per_mass, y_val_per_mol = minimum_work_of_gas_stream_separation(
            y_CO2, 0.99999, 0.99999, M_comp=44.0, T=298.15
        )
        y_vals_per_mass.append(y_val_per_mass)
        y_vals_per_mol.append(y_val_per_mol)

    y_vals = y_vals_per_mol

    ax.plot(x, y_vals, color="k", linewidth=2)

    # Calculate minimum work of air separation
    co2_sep_air_per_mass, co2_sep_air_per_mol = minimum_work_of_gas_stream_separation(
        y_CO2_in_air, 0.99999, 0.99999
    )

    # Annotate air separation with text box
    ax.text(
        y_CO2_in_air + 0.06,
        co2_sep_air_per_mol - 4,
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
        (y_CO2_in_air + 0.05, co2_sep_air_per_mol),
        arrowstyle="<-",
        mutation_scale=20,
        linewidth=1.5,
    )
    ax.add_patch(p1)

    secax = ax.secondary_yaxis(
        "right", functions=(lambda x: x * (1000 / M_CO2), lambda y: y / (1000 / M_CO2))
    )
    secax.set_yticks([0, 100, 200, 300, 400, 500, 600])

    ax.set_xlim([-0.01, 1.01])
    ax.set_ylim([0, 30])
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
        fig.savefig(output_path, dpi=1000, bbox_inches="tight")
    else:
        plt.show()


if __name__ == "__main__":
    plot_figure_w_min_CO2(
        save=False, output_path="./img/minimum_separation_work_co2.png"
    )
