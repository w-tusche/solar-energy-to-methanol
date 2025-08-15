# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd


class photorespiration(object):
    def __init__(self, T=15, o2_partial_pressure=210, co2_partial_pressure=0.252):
        """
        T: temperature in °C
        o2_partial_pressure: O2 partial pressure [micromol mol-1]
        co2_partial_pressure: CO2 partial pressure [micromol mol-1]
        """
        self.T_c = T  # Temperature in Celsius [°C]
        self.T_k = T + 273.15  # Temperature in Kelvin [K]
        self.o2_partial_pressure = (
            o2_partial_pressure  # O2 partial pressure [micromol mol-1]
        )
        self.co2_partial_pressure = (
            co2_partial_pressure  # CO2 partial pressure [micromol mol-1]
        )
        self.__calc()  # Calculation upon initialization

    def __calc(self):
        self.Vc_max = self.Bernacchi(
            "Vc_max", self.T_k
        )  # Maximum carboxylation rate [micromol m-2 s-1]
        self.Vo_max = self.Bernacchi(
            "Vo_max", self.T_k
        )  # Maximum oxygenation rate [micromol m-2 s-1]
        self.Kc = self.Bernacchi(
            "Kc", self.T_k
        )  # Michaelis-Menten constant for CO2 [micromol mol-1]
        self.Ko = (
            self.Bernacchi("Ko", self.T_k) * 1000
        )  # Michaelis-Menten constant for O2 [micromol mol-1]

        # Rubisco specificity factor
        self.tau = self.specificity_factor(self.Vc_max, self.Vo_max, self.Kc, self.Ko)
        # Ratio of RuBP oxygenations to carboxylations
        self.phi = self.ratio_oxygenations_carboxylations(
            self.o2_partial_pressure, self.co2_partial_pressure, self.tau
        )
        # Energy losses due to photorespiration
        self.d_pr = self.photorespiration(self.phi)

    def Bernacchi(self, param, T):
        """Calcualte generic temperature responses of the six parameters Vc,max, Vo,max, Γ*, Ko, Kc and Rd
        Formula by Bernacchi et al.: https://doi.org/10.1111/j.1365-3040.2001.00668.x

        Args:
            param (str): Parameter to calculate. One in ['R_d','Vc_max','Vo_max','Tau_star','Kc','Ko']
            T (float): Temperature [K]

        Returns:
            float: value of parameter selected at specific temperature
        """

        assert (
            param in ["R_d", "Vc_max", "Vo_max", "Tau_star", "Kc", "Ko"]
        ) == True, (
            "Input needs to be one of 'R_d','Vc_max','Vo_max','Tau_star','Kc','Ko'"
        )

        # Constants for Bernacchi's formula
        bernacchi_table = {
            "unit": {
                "R_d": "(micromol m-2 s-1)",
                "Vc_max": "(micromol m-2 s-1)",
                "Vo_max": "(micromol m-2 s-1)",
                "Tau_star": "(micromol mol-1)",
                "Kc": "(micromol mol-11)",
                "Ko": "(mmol mol-1)",
            },
            "val_at_25degC": {
                "R_d": 1.0,
                "Vc_max": 1.0,
                "Vo_max": 1.0,  # As in table but seems to be wrong
                "Tau_star": 42.75,
                "Kc": 404.9,
                "Ko": 278.4,
            },
            "c_dimless": {
                "R_d": 18.72,
                "Vc_max": 26.35,
                "Vo_max": 22.98,
                "Tau_star": 19.02,
                "Kc": 38.05,
                "Ko": 20.3,
            },
            "d_H_a_kJmol_1": {
                "R_d": 46.39,
                "Vc_max": 65.33,
                "Vo_max": 60.11,
                "Tau_star": 37.83,
                "Kc": 79.43,
                "Ko": 36.38,
            },
        }
        R = 8.314  # Gas constant in [J/(molK)]
        vls = pd.DataFrame(bernacchi_table).loc[param]

        # Calculate parameter using Bernacchi's formula
        parameter = np.exp(vls["c_dimless"] - vls["d_H_a_kJmol_1"] / (((R * T) / 1000)))

        return parameter

    def specificity_factor(self, Vc_max, Vo_max, K_c, K_o):
        """Calculate Rubisco specificity factor
        specificity factor of Rubisco for CO2 where specificity is the ratio of the probabilities of carboxylation to oxygenation
        https://doi.org/10.1111/j.1365-3040.2001.00668.x

        Args:
            Vc_max (float): _description_ [µmol m-2 s-1]
            Vo_max (float): _description_ [µmol m-2 s-1]
            K_c (float): _description_ [µmol mol-1]
            K_o (float): _description_ [µmol mol-1]

        Returns:
            float: Rubisco specificity factor
        """

        tau = (Vc_max * K_o) / (K_c * Vo_max)
        return tau

    def ratio_oxygenations_carboxylations(self, inter_cell_O, inter_cell_co2, tau):
        """Calculation of ratio of RuBP oxygenation and carboxylations
        Zhu, Xin-Guang; Long, Stephen P.; Ort, Donald R. (2008):
        What is the maximum efficiency with which photosynthesis can convert solar energy into biomass?
        In: Current opinion in biotechnology 19 (2), S. 153-159. DOI: 10.1016/j.copbio.2008.02.004.

        Args:
            inter_cell_O (float): intercellular O2 concentration [µmol mol-1]
            inter_cell_co2 (float): intercellular CO2 concentration [µmol mol-1]
            tau (float): specificity factor of Rubisco for CO2 where specificity is the ratio of the probabilities of carboxylation to oxygenation


        Returns:
            float: phi (Ratio of RuBP oxygenation and carboxylations)
        """
        phi = inter_cell_O / (inter_cell_co2 * tau)
        return phi

    def photorespiration(self, phi):
        """Calculation of decrease in epsilon_c caused by photorespiration (d_pr)
        Zhu, Xin-Guang; Long, Stephen P.; Ort, Donald R. (2008):
        What is the maximum efficiency with which photosynthesis can convert solar energy into biomass?
        In: Current opinion in biotechnology 19 (2), S. 153-159. DOI: 10.1016/j.copbio.2008.02.004.

        Args:
            phi (float): ratio of oxygenations to carboxylations

        Returns:
            float: photorespiration
        """
        d_pr = 1 - ((3 * (1 - (0.5 * phi))) / (3 + (3.5 * phi)))
        return d_pr


if __name__ == "__main__":
    PR = photorespiration(T=15, o2_partial_pressure=210, co2_partial_pressure=0.252)
    print(PR.d_pr)
