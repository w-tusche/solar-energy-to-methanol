# -*- coding: utf-8 -*-
import os

import numpy as np
import pandas as pd
import yaml

try:
    with open("./_variables.yml", "r") as vars:
        vars_dict = yaml.safe_load(vars)
except:
    with open("../_variables.yml", "r") as vars:
        vars_dict = yaml.safe_load(vars)


def create_tables(input_dir, output_dir):
    def make_eff_table(df, lower, higher):
        df_effs = pd.DataFrame(
            columns=["Process", "Theoretical efficiency", "Present efficiency"]
        )
        for df_, low, high in zip(df, lower, higher):
            key = df_["Process"].iloc[0]
            val = round(df_["Overall_efficiency_in_%"].iloc[-1], 1)
            low_str = f"{low:.1f}"
            high_str = f"{high:.1f}"
            real_range = (
                f"{low_str}"
                if round(low) == round(high)
                else f"{low_str} to {high_str}"
            )
            df_effs.loc[len(df_effs.index)] = [key, val, real_range]
        df_effs["Process"] = df_effs["Process"].str.replace("_", " ")
        df_effs.set_index("Process", inplace=True)
        return df_effs

    df_photosynthesis = pd.read_csv(
        os.path.join(input_dir, "photosynthesis_eff.csv"), sep=";"
    )
    df_anaerobic_digestion = pd.read_csv(os.path.join(input_dir, "AD_eff.csv"), sep=";")
    df_gasification = pd.read_csv(
        os.path.join(input_dir, "Gasification_eff.csv"), sep=";"
    )
    df_Biogas_Ref = pd.read_csv(
        os.path.join(input_dir, "Biogas_Ref_efficiencies.csv"), sep=";"
    )
    df_gasification_conditioning = pd.read_csv(
        os.path.join(input_dir, "Gas_conditioning_efficiencies.csv"), sep=";"
    )
    df_methanol_synthesis_bio = pd.read_csv(
        os.path.join(input_dir, "Methanol_synthesis_efficiencies_BIO.csv"), sep=";"
    )
    df_PV = pd.read_csv(os.path.join(input_dir, "PV_eff.csv"), sep=";")
    df_synthesis_gas_provision = pd.read_csv(
        os.path.join(input_dir, "Gas_Provision_efficiencies_el_eff.csv"), sep=";"
    )
    df_methanol_synthesis_el = pd.read_csv(
        os.path.join(input_dir, "Methanol_synthesis_efficiencies_Electricity.csv"),
        sep=";",
    )

    df_bio_theo = [
        df_photosynthesis,
        df_anaerobic_digestion,
        df_Biogas_Ref,
        df_gasification,
        df_gasification_conditioning,
        df_methanol_synthesis_bio,
    ]
    bio_real_lower = [
        vars_dict["eta_photosynthesis_real"],
        vars_dict["eta_anaerobic_digestion_lower_value"],
        vars_dict["eta_ref_conditioning_lower_value"],
        vars_dict["eta_gasification_lower_value"],
        vars_dict["eta_gas_conditioning_lower_value"],
        vars_dict["eta_conv_meoh"],
    ]
    bio_real_higher = [
        vars_dict["eta_photosynthesis_real"],
        vars_dict["eta_anaerobic_digestion_higher_value"],
        vars_dict["eta_ref_conditioning_higher_value"],
        vars_dict["eta_gasification_higher_value"],
        vars_dict["eta_gas_conditioning_higher_value"],
        vars_dict["eta_conv_meoh"],
    ]
    df_el = [df_PV, df_synthesis_gas_provision, df_methanol_synthesis_el]
    el_real_lower = [
        vars_dict["eta_pv_real"],
        vars_dict["eta_syngas_prod_real_lower_value"],
        vars_dict["eta_direct_meoh"],
    ]
    el_real_higher = [
        vars_dict["eta_pv_real"],
        vars_dict["eta_syngas_prod_real_higher_value"],
        vars_dict["eta_direct_meoh"],
    ]

    df_bio = make_eff_table(df_bio_theo, bio_real_lower, bio_real_higher)
    df_bio.rename(
        index={"Methanol synthesis": "Conventional methanol synthesis"}, inplace=True
    )
    df_bio.rename(
        mapper={
            "Theoretical efficiency": "Theoretical efficiency in %",
            "Present efficiency": "Present efficiency in %",
        },
        axis=1,
        inplace=True,
    )

    print(df_bio)
    df_bio.to_csv(os.path.join(output_dir, "efficiencies-overall-bio.csv"))

    df_electro = make_eff_table(df_el, el_real_lower, el_real_higher)
    df_electro.rename(
        index={
            "Methanol synthesis efficiencies Electricity": "Direct methanol synthesis"
        },
        inplace=True,
    )
    df_electro.rename(
        mapper={
            "Theoretical efficiency": "Theoretical efficiency in %",
            "Present efficiency": "Present efficiency in %",
        },
        axis=1,
        inplace=True,
    )

    print(df_electro)
    df_electro.to_csv(os.path.join(output_dir, "efficiencies-overall-electro.csv"))

    # Combine table 1 (Biobased pathway)
    df_combined1 = pd.concat(
        [
            df_photosynthesis,
            df_anaerobic_digestion,
            df_Biogas_Ref,
            df_gasification,
            df_gasification_conditioning,
            df_methanol_synthesis_bio,
        ],
        ignore_index=True,
    )

    df_bio_step_efficiencies = df_combined1.round(2)
    df_bio_step_efficiencies.replace("_", " ", regex=True, inplace=True)
    df_bio_step_efficiencies.columns = df_bio_step_efficiencies.columns.str.replace(
        "_", " ", regex=True
    )
    df_bio_step_efficiencies = df_bio_step_efficiencies.loc[
        :, ["Process", "Partial efficiency in %"]
    ]  # kick out Overall efficiency
    df_bio_step_efficiencies.rename(
        index={
            "Methanol synthesis Methanol synthesis from CO and H2": "Conventional methanol synthesis"
        },
        inplace=True,
    )

    print(df_bio_step_efficiencies)
    df_bio_step_efficiencies.to_csv(
        os.path.join(output_dir, "step-efficiencies-overall-bio.csv"), index=None
    )

    # Combine table 2 (Electricity pathway)

    df_combined_2 = pd.concat(
        [
            df_PV,
            df_synthesis_gas_provision,
            df_methanol_synthesis_el,
        ],
        ignore_index=True,
    )

    # round table columns to 2 digits (columns need to be number format; not str or other)
    df_electro_step_efficiencies = df_combined_2.round(2)
    df_electro_step_efficiencies.replace("_", " ", regex=True, inplace=True)
    df_electro_step_efficiencies.columns = (
        df_electro_step_efficiencies.columns.str.replace("_", " ", regex=True)
    )
    df_electro_step_efficiencies = df_electro_step_efficiencies.loc[
        :, ["Process", "Partial efficiency in %"]
    ]  # kick out Overall efficiency
    df_electro_step_efficiencies.rename(
        index={
            "Methanol synthesis efficiencies Electricity Methanol-synthesis": "Direct methanol synthesis"
        },
        inplace=True,
    )

    print(df_electro_step_efficiencies)
    df_electro_step_efficiencies.to_csv(
        os.path.join(output_dir, "step-efficiencies-overall-electro.csv"), index=None
    )


def main(input_dir="./data/results/", output_dir="./data/results/"):
    create_tables(input_dir, output_dir)


if __name__ == "__main__":
    main()
