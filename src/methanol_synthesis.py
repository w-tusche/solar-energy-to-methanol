# -*- coding: utf-8 -*-
import os

import material_props as props
import utils

EoF = props.standard_enthalpy_of_formation_25C


def conventional_methanol_synthesis():
    # Biomass-based methanol synthesis

    # CO + 2 H2 -> CH3OH

    # Reaction enthalpy change \Delta H for conventional methanol synthesis (25°C) gaseous CH3OH
    dH_BM_MeOH = EoF["CH3OH(g)"] - EoF["CO(g)"]
    print(
        f"Conventional methanol synthesis: exothermic reaction with {dH_BM_MeOH:.2f} (kJ/mol_CH3OH) (CH3OH(g))"
    )
    # Reaction enthalpy change \Delta H for conventional methanol synthesis (25°C) liquid CH3OH
    dH_BM_MeOH = EoF["CH3OH(l)"] - EoF["CO(g)"]
    print(
        f"Conventional methanol synthesis: exothermic reaction with {dH_BM_MeOH:.2f} (kJ/mol_CH3OH) (CH3OH(l))"
    )

    HHV_MeOH = (2 * EoF["H2O(l)"] + EoF["CO2(g)"] - EoF["CH3OH(l)"]) * -1
    HHV_H2 = -EoF["H2O(l)"]
    HHV_CO = EoF["CO(g)"] - EoF["CO2(g)"]

    # Efficiency of biomass-based synthesis
    # Efficiency regarding liquid methanol as product
    Eff_BMSyn = HHV_MeOH / (HHV_CO + (2 * HHV_H2))
    print(f"Conventional methanol synthesis: efficiency {Eff_BMSyn*100:.2f} %")

    # Data-Frame Bio-based methanol synthesis

    step_efficiencies = {"Methanol_synthesis_from_CO_and_H2": Eff_BMSyn}

    df_Methanol_synthesis_efficiencies_BIO = utils.table_creation(
        "Methanol_synthesis", step_efficiencies
    )

    return df_Methanol_synthesis_efficiencies_BIO


def direct_methanol_synthesis():
    # Electricity-based methanol synthesis

    # CO2 + 3 H2 -> CH3OH + H2O

    # Reaction enthalpy change \Delta H for direct methanol synthesis (25°C) gaseous CH3OH and water
    dH_E_MeOH = EoF["CH3OH(g)"] + EoF["H2O(g)"] - EoF["CO2(g)"]
    print(
        f"Direct methanol synthesis: exothermic reaction with {dH_E_MeOH:.2f} (kJ/mol_CH3OH) (CH3OH(g) & H2O(g))"
    )
    # Reaction enthalpy change \Delta H for direct methanol synthesis (25°C) liquid CH3OH and water
    dH_E_MeOH = EoF["CH3OH(l)"] + EoF["H2O(l)"] - EoF["CO2(g)"]
    print(
        f"Direct methanol synthesis: exothermic reaction with {dH_E_MeOH:.2f} (kJ/mol_CH3OH) (CH3OH(l) & H2O(l))"
    )

    HHV_MeOH = (2 * EoF["H2O(l)"] + EoF["CO2(g)"] - EoF["CH3OH(l)"]) * -1
    HHV_H2 = -EoF["H2O(l)"]

    # Efficiency of electricity-based synthesis

    Eff_PSyn = HHV_MeOH / (3 * HHV_H2)

    minSep = 1.316  # kJ/mol from AspenPlus
    # Efficiency regarding liquid methanol as product
    Eff_PSyn_Sep = HHV_MeOH / (3 * HHV_H2 + minSep)
    print(f"Direct methanol synthesis: efficiency {Eff_PSyn_Sep*100:.2f} %")

    # Data-Frame Electricity-based methanol synthesis

    step_efficiencies = {"Methanol-synthesis": Eff_PSyn_Sep}

    df_Methanol_synthesis_efficiencies_Electricity = utils.table_creation(
        "Methanol_synthesis_efficiencies_Electricity", step_efficiencies
    )

    return df_Methanol_synthesis_efficiencies_Electricity


def main(output_dir="./data/results/"):
    print("\nMethanol synthesis efficiencies:\n")
    df_Methanol_synthesis_efficiencies_BIO = conventional_methanol_synthesis()
    print(df_Methanol_synthesis_efficiencies_BIO)
    df_Methanol_synthesis_efficiencies_BIO.to_csv(
        os.path.join(output_dir, "Methanol_synthesis_efficiencies_BIO.csv"), sep=";"
    )

    df_Methanol_synthesis_efficiencies_Electricity = direct_methanol_synthesis()
    print(df_Methanol_synthesis_efficiencies_Electricity)
    df_Methanol_synthesis_efficiencies_Electricity.to_csv(
        os.path.join(output_dir, "Methanol_synthesis_efficiencies_Electricity.csv"),
        sep=";",
    )


if __name__ == "__main__":
    main()
