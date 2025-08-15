# -*- coding: utf-8 -*-
import os
from typing import Tuple

import numpy as np
import pandas as pd
import scipy.constants
import scipy.optimize as optimize
import yaml

import material_props as props
import utils
from dac import minimum_work_of_gas_stream_separation as w_min_gas_sep

# Constants
R_GAS = scipy.constants.R  # J/(mol*K)
MOLAR_MASS = props.molar_mass_g_mol
M_CO2 = MOLAR_MASS["CO2"]  # g/mol (molar mass)
M_O2 = MOLAR_MASS["O2"]  # g/mol (molar mass)
E_OF = props.standard_enthalpy_of_formation_25C
HHV = props.HHV_kJ_g

O2_CONC_AIR = 0.20946  # O2 concentration in air (btw.: 0 -> 1)
O2_PURITY = 0.9999  # Purity of separated O2; pure
CO2_PURITY = 0.9999
CONVERSION_EFFICIENCY = 0.605  # Carnot efficiency at 520°C for electricity generation (60.5%; https://doi.org/10.1016/j.renene.2020.02.091)


def energy_demand_for_co2_separation(
    y_gas: float, y_co2: float, ccr: float
) -> Tuple[float, float]:
    """Calculate energy demand for CO2 separation

    Args:
        y_gas (float): Initial CO2 concentration in gas stream (0 to 1).
        y_co2 (float): Purity of separated CO2.
        ccr (float): Carbon capture rate or conversion rate.

    Returns:
        Tuple containing:
        - specific energy demand (kJ/kg CO2)
        - specific energy demand (kJ/mol CO2)
    """

    w_min_mass, w_min_mol = w_min_gas_sep(y_gas, y_co2, CR=ccr, M_comp=M_CO2, T=298.15)

    print(f"Specific energy demand of separation: {w_min_mol:.2f} J/mol CO2")

    return w_min_mass, w_min_mol


def energy_demand_for_o2_separation(o2_moles: float) -> float:
    """Calculate energy demand to separate given amount of O2 from air.

    Args:
        o2_moles (float): Moles of O2 to separate.

    Returns:
        Total energy demand (kJ) for separation process.
    """

    O2_CR = 0.9

    w_min_mass, w_min_mol = w_min_gas_sep(
        O2_CONC_AIR, O2_PURITY, CR=O2_CR, M_comp=M_O2, T=298.15
    )
    total_energy = o2_moles * w_min_mol
    return total_energy


def anaerobic_digestion() -> pd.DataFrame:
    """Calculate efficiencies for anaerobic digestion step."""

    # Reaction: C6H12O6 -> 3 CO2 + 3 CH4; calculate enthalpy change
    delta_h = 3 * E_OF["CO2(g)"] + 3 * E_OF["CH4(g)"] - E_OF["C6H12O6(s)"]

    # Reaction enthalpy loss (positive value)
    loss_reaction_enthalpy = -delta_h  # kJ/mol

    # Biomass regeneration loss (assumption: 7% of original biomass)
    loss_regeneration = HHV["C6H12O6"] * 0.07  # kJ/g

    total_loss = 1 - (
        (loss_regeneration / HHV["C6H12O6"])
        + (loss_reaction_enthalpy / MOLAR_MASS["C6H12O6"]) / HHV["C6H12O6"]
    )

    step_efficiencies = {
        "Biomass_formation_and_conversion": total_loss,
    }

    df_efficiencies = utils.table_creation("Anaerobic_digestion", step_efficiencies)

    return df_efficiencies


def biogas_reforming() -> Tuple[np.ndarray, float, float, float, float]:
    """Calculate stoichiometry and efficiency for biogas reforming.

    Returns:
        Tuple of (solution, reforming efficiency, energy_CH4, energy_CO, energy_H2)
    """
    # Biogas reforming
    # a CH4 + b CO2 + c O2 + d H2O -> x CO + y H2

    # Specifications:   - molar H2 to CO ratio = 2
    #                   - Sum enthalpy products - Sum enthalpy reactants = 0 (autothermal refomring)
    #                   - based on 1 mol CH4

    # C-balance: 1a + 1b = 1x
    # H-balance: 4a + 2d = 2y
    # O-balance: 2b + 2c + d = x
    # dH: -110.53x + 74.8a + 393.52b + 285.84d = 0
    # H2/CO: 2x = y
    # Set CH4: a = 1

    # System of equations for autothermal reforming based on 1 mol CH4
    KM = np.array(
        [
            [1, 1, 0, 0, -1, 0],
            [4, 0, 0, 2, 0, -2],
            [0, 2, 2, 1, -1, 0],
            [
                -E_OF["CH4(g)"],
                -E_OF["CO2(g)"],
                -E_OF["O2(g)"],
                -E_OF["H2O(l)"],
                E_OF["CO(g)"],
                E_OF["H2(g)"],
            ],
            [0, 0, 0, 0, -2, 1],
            [1, 0, 0, 0, 0, 0],
        ]
    )
    LV = np.array([0, 0, 0, 0, 0, 1])

    solution = np.linalg.solve(KM, LV)
    print(solution)

    # Report solution
    components = ["CH4(g)", "CO2(g)", "O2(g)", "H2O(l)", "CO(g)", "H2(g)"]
    for mol, comp in zip(solution, components):
        print(f"Per mole CH4 stoichiometric there are: {mol:.2f} moles of {comp}")

    energy_CH4 = HHV["CH4"] * MOLAR_MASS["CH4"] * solution[0]
    energy_CO = HHV["CO"] * MOLAR_MASS["CO"] * solution[4]
    energy_H2 = HHV["H2"] * MOLAR_MASS["H2(g)"] * solution[5]

    # Efficiency of pure reforming
    efficiency = (energy_CO + energy_H2) / energy_CH4

    return solution, efficiency, energy_CH4, energy_CO, energy_H2


def syngas_production_via_anaerobic_digestion_theory() -> (
    Tuple[pd.DataFrame, pd.DataFrame]
):
    print("\nSyngas Production via Anaerobic Digestion:\n")

    df_AD_efficiencies = anaerobic_digestion()

    (
        solution,
        eff_pure_reforming,
        energy_CH4,
        energy_CO,
        energy_H2,
    ) = biogas_reforming()

    moles_co2_anaerobic = solution[1]
    moles_o2_anaerobic = solution[2]

    # Energy demand of CO2 separation
    y_bg = 0.5  # CO2 concentration in biogas (btw.: 0 -> 1)
    y_co2 = CO2_PURITY  # Purity of separated CO2; approx. 1 to neglected CH4 losses

    # To have the right CO2 amount for further processing the carbon capture rate has to be found
    # Aim: remaining CO2 is the CO2 we can use downstream
    y_rest_ref = solution[1] / (
        1 + solution[1]
    )  # Control: CO2 share in Methane stream from Eq-System above
    n_BGges = 1  # arbitrary, relative value

    def calc_y_rest(CCR):
        n_CO2 = n_BGges * y_bg * CCR  # CO2 amount in CO2 stream
        n_CO2ges = n_CO2 / y_co2  # Moleflow CO2 stream
        n_Rest = n_BGges - n_CO2ges  # Moleflow Methane stream
        n_CO2Rest = n_BGges * y_bg - n_CO2  # CO2 flow in Methane stream
        return n_CO2Rest / n_Rest  # CO2 share in Methane stream

    def f(CCR):
        return calc_y_rest(CCR) - y_rest_ref

    # Bracket points
    a = 0.0
    b = 1.0

    # Verify bracket
    if f(a) * f(b) > 0:
        raise ValueError(f"No root bracketed in [{a}, {b}]")

    # Find root using brentq
    optimal_CCR = optimize.brentq(f, a, b, xtol=1e-11)

    print(f"Optimal CCR: {optimal_CCR:.4f}")
    print(f"Residual CO2 fraction at optimal CCR: {calc_y_rest(optimal_CCR):.4f}")

    CCR = optimal_CCR  # Based on amount of CO2 that has to be separated to reach y_Rest (given from above)

    w_min_co2_mass, w_min_co2_mol = energy_demand_for_co2_separation(y_bg, y_co2, CCR)

    w_min_bg_sep = (1 - moles_co2_anaerobic) * (w_min_co2_mol / CONVERSION_EFFICIENCY)
    print(f"Specific energy demand of separation: {w_min_co2_mass:.2f} kJ/kg CO2")
    print(f"Specific energy demand of separation: {w_min_bg_sep:.2f} kJ/mol CH4")

    # Efficiency of CO2 separation + reforming
    eff_ref_co2 = (energy_CO + energy_H2) / (energy_CH4 + w_min_bg_sep)
    eff_co2_sep = 1 - (eff_pure_reforming - eff_ref_co2)

    print(
        f"Efficiency pure reforming {eff_pure_reforming*100:.2f} % (Energy basis CH4| HHV)"
    )
    print(
        f"Efficiency including CO2 separation {eff_ref_co2*100:.2f} % (Energy basis CH4| HHV)"
    )

    w_min_o2_sep = energy_demand_for_o2_separation(moles_o2_anaerobic)
    print(f"Energy demand of O2 separation per mol CH4: {w_min_o2_sep:.2f} kJ/mol CH4")

    eff_ref_co2_o2 = (energy_CO + energy_H2) / (
        energy_CH4 + w_min_bg_sep + w_min_o2_sep / CONVERSION_EFFICIENCY
    )
    eff_o2_sep = 1 - (eff_ref_co2 - eff_ref_co2_o2)
    print(
        f"Efficiency including CO2 separation and O2 provision {eff_ref_co2_o2*100:.2f} % (Energy basis CH4| HHV)"
    )

    step_efficiencies = {
        "Reforming": eff_pure_reforming,
        "CO2_separation_and_conditioning": eff_co2_sep,
        "O2_separation": eff_o2_sep,
    }
    df_biogas_ref_efficiencies = utils.table_creation(
        "Biogas_reforming_and_conditioning", step_efficiencies
    )

    return df_AD_efficiencies, df_biogas_ref_efficiencies


def biomass_gasification() -> (
    Tuple[pd.DataFrame, np.ndarray, float, float, float, float]
):
    """Calculate stoichiometry and efficiency for biomass gasification."""
    # Biomass gasification

    # a C6H12O6 + b O2 + c H2O -> x CO + y H2 + z CO2

    # Specifications:   - molar H2 to CO ratio = 2
    #                   - Sum enthalpy products - Sum enthalpy reactants = 0
    #                   - based on 1 mol glucose

    # equation_1: 6a = x + z
    # equation_2: 12a + 2c = 2y
    # equation_3: 6a + 2b + c = x + 2z
    # equation_4: -110.53x - 393.52z + 1273.7a + 285.84c = 0
    # equation_5: 2x = y
    # equation_6: a = 1

    # Define the coefficient matrix A and the constant vector B
    A = np.array(
        [
            [6, 0, 0, -1, 0, -1],
            [12, 0, 2, 0, -2, 0],
            [6, 2, 1, -1, 0, -2],
            [
                -E_OF["C6H12O6(s)"],
                -E_OF["O2(g)"],
                -E_OF["H2O(l)"],
                E_OF["CO(g)"],
                E_OF["H2(g)"],
                E_OF["CO2(g)"],
            ],
            [0, 0, 0, -2, 1, 0],
            [1, 0, 0, 0, 0, 0],
        ]
    )
    B = np.array([0, 0, 0, 0, 0, 1])

    solution = np.linalg.solve(A, B)

    comps = ["C6H12O6(s)", "O2(g)", "H2O(l)", "CO(g)", "H2(g)", "CO2(g)"]
    for mol, comp in zip(solution, comps):
        print(f"Per mole glucose stoichiometric there are: {mol:.2f} moles of {comp}")

    energy_glucose = HHV["C6H12O6"] * MOLAR_MASS["C6H12O6"]
    energy_CO = HHV["CO"] * MOLAR_MASS["CO"] * solution[3]
    energy_H2 = HHV["H2"] * MOLAR_MASS["H2(g)"] * solution[4]

    # Efficiency of pure gasification
    efficiency = (energy_CO + energy_H2) / energy_glucose

    step_efficiencies = {
        "Autothermal_conversion_of_solids_into_gases": efficiency,
    }

    df_efficiencies = utils.table_creation("Gasification", step_efficiencies)

    return df_efficiencies, solution, efficiency, energy_glucose, energy_CO, energy_H2


def syngas_production_via_biomass_gasification_theory():
    print("\nSyngas Production via Gasification:\n")

    # Biomass gasification
    (
        df_efficiencies,
        solution,
        efficiency,
        energy_glucose,
        energy_CO,
        energy_H2,
    ) = biomass_gasification()

    moles_glucose = solution[0]
    moles_o2 = solution[1]
    moles_co = solution[3]
    moles_h2 = solution[4]
    moles_co2 = solution[5]

    total_moles_product = moles_co + moles_h2 + moles_co2

    y_gas = moles_co2 / total_moles_product

    print(f"CO2 mole fraction in syngas: {y_gas:.2f}")
    print(f"Total moles product gas: {total_moles_product:.2f} mol")

    # CO2 separation
    _, w_min_co2_mol = energy_demand_for_co2_separation(y_gas, CO2_PURITY, 0.9999)
    w_min_co2_sep = w_min_co2_mol * moles_co2
    print(
        f"Energy demand of CO2 separation per mol glucose: {w_min_co2_sep:.2f} kJ/mol"
    )

    w_min_o2_sep = energy_demand_for_o2_separation(moles_o2)

    eff_gas_co2 = (energy_CO + energy_H2) / (
        energy_glucose + w_min_co2_sep / CONVERSION_EFFICIENCY
    )
    eff_gas_co2_o2 = (energy_CO + energy_H2) / (
        energy_glucose
        + w_min_co2_sep / CONVERSION_EFFICIENCY
        + w_min_o2_sep / CONVERSION_EFFICIENCY
    )

    eff_co2_sep_gas = 1 - (efficiency - eff_gas_co2)
    eff_o2_prov_gas = 1 - (eff_gas_co2 - eff_gas_co2_o2)

    print(
        f"Efficiency pure gasification {efficiency*100:.2f} % (Energy basis glucose| HHV)"
    )
    print(
        f"Efficiency including CO2 separation {eff_gas_co2*100:.2f} % (Energy basis glucose| HHV)"
    )
    print(
        f"Efficiency including CO2 separation and O2 provision {eff_gas_co2_o2*100:.2f} % (Energy basis glucose| HHV)"
    )

    step_efficiencies = {
        "CO2_separation": eff_co2_sep_gas,
        "O2_provision": eff_o2_prov_gas,
    }
    df_gas_conditioning_efficiencies = utils.table_creation(
        "Gas_conditioning_efficiencies", step_efficiencies
    )

    # Carbon efficiency calculation
    carbon_separation_cost = (w_min_o2_sep + w_min_co2_sep) / (
        energy_glucose * CONVERSION_EFFICIENCY / 6
    )
    carbon_efficiency = moles_co / (6 * moles_glucose + carbon_separation_cost)
    print(f"Carbon efficiency (gasification): {carbon_efficiency*100:.2f} %")

    print(
        f"Energy demand for CO2 separation {w_min_co2_sep / CONVERSION_EFFICIENCY:.2f} kJ/mol Glucose"
    )
    print(
        f"Energy demand for O2 provision {w_min_o2_sep / CONVERSION_EFFICIENCY:.2f} kJ/mol Glucose"
    )

    return df_efficiencies, df_gas_conditioning_efficiencies


def efficiency_syngas_provision(
    w_min_DAC_spec_mol: float, efficiency_electrolysis: float
) -> float:
    """Calculate efficiency of syngas production via DAC and electrolysis.

    CO2 + 3 H2 -> CH3OH + H2O
    energy 1 mole CO2 + energy (HHV) 3 mole H2

    How much of the input (chemical) energy is still in the resulting gas stream?

    Args:
        w_min_DAC_spec_mol (float): Work for direct air capture (kJ/mol).
        efficiency_electrolysis (float): Efficiency of electrolysis (0 to 1).

    Returns:
        Overall efficiency (float).
    """
    if not (0 <= efficiency_electrolysis <= 1):
        raise ValueError("Electrolysis efficiency must be between 0 and 1.")

    energy_educts = 3 * props.HHV_kJ_mol["H2"]  # kJ/mol for 3 moles of H2
    energy_for_production = (
        w_min_DAC_spec_mol + energy_educts / efficiency_electrolysis
    )  # kJ/mol

    # DAC has no "Energy" efficiency (CO2 fully oxidized)
    # 3 HHV_H2 / (w_CO2 + 3 HHV_H2 * eta_H2)
    overall_efficiency = energy_educts / energy_for_production

    return overall_efficiency


def syngas_production_via_dac_and_electrolysis_theory_and_real() -> pd.DataFrame:
    print("\nSyngas Production via Direct Air Capture and Electrolysis:\n")

    # Minimum energy for CO2 separation (pure gas)
    w_min_mass, w_min_mol = w_min_gas_sep(
        props.y_CO2_air, 0.999999999, CR=0.9999999, M_comp=M_CO2, T=298.15
    )  # (kJ/kg) and specific per mol (kJ/mol)

    # Add fan energy (kJ/mol)
    fan_energy_kJ_per_mol = utils.kWh_per_t_to_kJ_per_mol(200, 44.0)
    w_min_DAC_spec_mol_opt = w_min_mol + fan_energy_kJ_per_mol

    eff_el_opt = 1.0  # Theoretical efficiency of electroylsis = 100 %

    # Real values from literature or yaml config (HHV)
    try:
        with open("./_variables.yml", "r") as vars_file:
            vars_dict = yaml.safe_load(vars_file)
    except FileNotFoundError:
        print(
            "Warning: _variables.yml not found. Using default electrolysis efficiencies."
        )
        vars_dict = {
            "eta_el_real_lower_value": 54,
            "eta_el_real_higher_value": 71,
        }

    eff_el_real_low = vars_dict["eta_el_real_lower_value"] / 100
    eff_el_real_high = vars_dict["eta_el_real_higher_value"] / 100

    # Present efficiency (values from the text!)
    w_min_DAC_spec_mol_real = 316.80  # kJ/mol (Fasihi et al.)

    syngas_prod_opt = efficiency_syngas_provision(w_min_DAC_spec_mol_opt, eff_el_opt)
    syngas_prod_real_low = efficiency_syngas_provision(
        w_min_DAC_spec_mol_real, eff_el_real_low
    )
    syngas_prod_real_high = efficiency_syngas_provision(
        w_min_DAC_spec_mol_real, eff_el_real_high
    )

    print(
        f"Optimal efficiency (100% EL, DAC theoretical min energy): {syngas_prod_opt*100:.2f} %"
    )
    print(f"Real efficiency lower bound (54% EL): {syngas_prod_real_low*100:.2f} %")
    print(f"Real efficiency upper bound (71% EL): {syngas_prod_real_high*100:.2f} %")

    df_efficiencies = utils.table_creation(
        "Provision_of_synthesis_gas", {"Electrolysis_and_DAC": syngas_prod_opt}
    )

    return df_efficiencies


def main(output_dir="./data/results/"):
    # Syngas provision via Anaerobic digestion
    df_AD_eff, df_Biogas_Ref_eff = syngas_production_via_anaerobic_digestion_theory()
    print(df_AD_eff)
    print(df_Biogas_Ref_eff)
    df_Biogas_Ref_eff.to_csv(
        os.path.join(output_dir, "Biogas_Ref_efficiencies.csv"), sep=";"
    )
    df_AD_eff.to_csv(os.path.join(output_dir, "AD_eff.csv"), sep=";")

    # Syngas provision via Gasification
    (
        df_Gasification_eff,
        df_Gas_conditioning_eff,
    ) = syngas_production_via_biomass_gasification_theory()
    print(df_Gasification_eff)
    print(df_Gas_conditioning_eff)
    df_Gasification_eff.to_csv(
        os.path.join(output_dir, "Gasification_eff.csv"), sep=";"
    )
    df_Gas_conditioning_eff.to_csv(
        os.path.join(output_dir, "Gas_conditioning_efficiencies.csv"), sep=";"
    )

    # Syngas provision via Direct Air Capture and Electrolysis
    df_Gas_Provision_eff_el = (
        syngas_production_via_dac_and_electrolysis_theory_and_real()
    )
    print(df_Gas_Provision_eff_el)
    df_Gas_Provision_eff_el.to_csv(
        os.path.join(output_dir, "Gas_Provision_efficiencies_el_eff.csv"), sep=";"
    )


if __name__ == "__main__":
    main()
