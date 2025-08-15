# -*- coding: utf-8 -*-
standard_enthalpy_of_formation_25C = {
    "H2O(l)": -285.84,  # all data from DOI 10.1007/978-3-662-49924-5
    "H2O(g)": -241.83,
    "C6H12O6(s)": -1273.7,  # except glucose http://dx.doi.org/10.1016/j.jct.2012.11.031
    "CO2(g)": -393.52,
    "O2(g)": 0,
    "H2(g)": 0,
    "CO(g)": -110.53,
    "CH4(g)": -74.8,
    "CH3OH(l)": -238.57,
    "CH3OH(g)": -201.17,
}  # 298.15 K and 1 atm deltaHf0 (kJ/mol)

# HHV and LHV at 25°C
HHV_kJ_g = {"C6H12O6": 15.556, "CH4": 55.511, "CO": 10.103, "H2": 141.786}  # (kJ/g)
HHV_kJ_mol = {"H2": 285.8}  # (kJ/mol)
LHV_kJ_mol = {"H2": 241.8}  # (kJ/mol)
HHV_kWh_kg = {"H2": 39.39}  # (kWh/kg)
LHV_kWh_kg = {"H2": 33.33}  # (kWh/kg)

molar_mass_g_mol = {
    "CO2": 44.01,
    "O2": 31.999,
    "CH4": 16.04,
    "CO": 28.01,
    "C6H12O6": 180.156,
    "H2(g)": 2.016,
}  # g/mol (molar mass)

# Set values
y_CO2_air = 0.00042  # molar fraction of CO2 in air (420 ppm)
