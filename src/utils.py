# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import scipy


def create_shockley_queisser_limit_object_WL_and_irradiation():
    import SQLimit.SQlimit as SQL

    # Load ASTMG Data
    # data range: 280nm to 4000nm, 0.31eV to 4.42857 eV
    # https://www.nrel.gov/grid/solar-resource/spectra-am1.5.html
    try:
        ref_solar_ASTMG = pd.read_csv(
            r"./data/raw/astmg173.csv", sep=";"
        )  # (nm vs W m^-2 nm^-1)
    except:
        ref_solar_ASTMG = pd.read_csv(
            r"../data/raw/astmg173.csv", sep=";"
        )  # (nm vs W m^-2 nm^-1)

    # Extract Wavelength
    WL_ASTMG = ref_solar_ASTMG.iloc[:, 0]  # (nm)

    # Extract solar spectrum irradiation data from ASTMG173 (direct)
    solar_per_nm_ASTMG = ref_solar_ASTMG.iloc[:, 2]  # (W*m-2*nm-1)

    # Calculate the Shockley-Queisser Limit
    SQ_ASTMG = SQL.SQlim(WL=WL_ASTMG, solar_per_nm=solar_per_nm_ASTMG)

    return WL_ASTMG, solar_per_nm_ASTMG, SQ_ASTMG


def quantum_energy(wavelength):
    """
    Photon energy calculation.

    Args:
        wavelength (float): wavelength of photon in meters [m]

    Returns:
        float: energy of photon in joule [J]
    """
    h = scipy.constants.h  # Planck's constant [m^2 kg / s] = [J*s]
    c = scipy.constants.c  # Speed of light in vacuum [m/s]
    # Delta_E_q = h * v = (h*c)/lambda
    delta_Eq = (h * c) / wavelength
    return delta_Eq


def planck_radiation_law(wavelength, T):
    """Calculation of spectral radiance of a blackbody depending on wavelength and temperature.
    Compare to: https://en.wikipedia.org/wiki/Planck%27s_law B_lambda(lambda,T)

    Args:
        wavelength (float): wavelength in [m]
        T (int): temperature in [K]

    Returns:
        float: spectral radiance in [W/(m^-2 sr^-1 nm^-1)]
    """

    # https://en.wikipedia.org/wiki/Planck%27s_law
    h = scipy.constants.h  # Planck's constant h = 6.626e-34 (m^2 kg / s) = (J*s)
    c = scipy.constants.c  # Speed of light in vacuum c = 3e8 (m/s)
    k = scipy.constants.k  # Boltzmann constant k = 1.38e-23 [m^2 kg / s^2 K] = [J/K]

    # Spectral Radiance (Powerflow per unit of solid angle)
    a = 2 * h * c**2
    b = (h * c) / (wavelength * k * T)
    spectral_radiance = a / (
        wavelength**5 * (np.exp(b) - 1)
    )  # in [W/(m^-2 sr^-1 m^-1)] sr = Steradian
    spectral_radiance = spectral_radiance * 1e-9  # in [W/(m^-2 sr^-1 nm^-1)]

    return spectral_radiance


def irradiance_earth(spectral_radiance):
    """
    Calculates spectral irradiance (at earths atmosphere) from the spectral
    radiance of the sun.
    Compare: https://en.wikipedia.org/wiki/Planck%27s_law#/media/File:EffectiveTemperature_300dpi_e.png

    Args:
        spectral_radiance (float): spectral radiance of blackbody (sun) in [W/(m^-2 sr^-1 nm^-1)]

    Returns:
        float: spectral irradiance in [W/(m^-2 nm^-1)]
    """

    # Irradiance (integrated over solid angle)
    radius_sun = (
        695_700_000  # [m] https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html
    )
    mean_distance_sun_earth = 149.6e9  # [m] (= 1 Astronomical Unit (AU)) https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html

    # Account for 1/r^2 dependency of intensity on distance
    intensity_factor = (radius_sun / mean_distance_sun_earth) ** 2

    #  x pi because of Lambert's cosine law (integration over solid angle)
    spectral_irradiance = (
        spectral_radiance * intensity_factor * np.pi
    )  # Compare: https://en.wikipedia.org/wiki/Planck%27s_law#/media/File:EffectiveTemperature_300dpi_e.png

    return spectral_irradiance


def nm_to_eV(wavelength_nm):
    # hc/(\lambda [m] * eV)
    wavelength_nm = np.array(wavelength_nm, float)
    # Handling of wavelength_nm == 0
    near_zero = np.isclose(wavelength_nm, 0)
    wavelength_nm[near_zero] = np.inf
    wavelength_nm[~near_zero] = (scipy.constants.h * scipy.constants.c) / (
        wavelength_nm[~near_zero] * 1e-9 * scipy.constants.eV
    )
    return wavelength_nm


def eV_to_nm(eV):
    # hc/(E[eV] * eV* nm/m)
    eV = np.array(eV, float)
    # Handling of eV == 0
    near_zero = np.isclose(eV, 0)
    eV[near_zero] = np.inf
    eV[~near_zero] = (
        (scipy.constants.h * scipy.constants.c)
        / (eV[~near_zero] * scipy.constants.eV)
        * 1e9
    )
    return eV


def trapz_integration_of_area(df, x_min, x_max):
    """
    Calculate integral from numerica data in an dataframe. Dataframe needs to
    have the x values in the first column and the y values in the second
    column.

    Args:
        df (pd.Dataframe): Pandas Dataframe with x value in column 0 and y values in column 1
        x_min (float): minimum value to form lower bound for integration
        x_max (_type_): minimum value to form upper bound for integration

    Returns:
        float: Integral value
    """
    # Select area to integrade from the array
    assert x_min < x_max
    area = df[(df.iloc[:, 0] >= x_min) & (df.iloc[:, 0] <= x_max)]
    # Integration with trapez
    value = scipy.integrate.trapezoid(area.iloc[:, 1], x=area.iloc[:, 0])
    return value


def wavelength_to_rgb(wavelength, gamma=0.8, alpha_spectrum=1.0, alpha_outside=0.5):
    """
    From http://www.noah.org/wiki/Wavelength_to_RGB_in_Python
    This converts a given wavelength of light to an
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).

    Idea: https://stackoverflow.com/questions/44959955/matplotlib-color-under-curve-based-on-spectral-color"
    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    Additionally alpha value set to 0.5 outside range

    Visible range reported differently and changes over the lifetime. https://www.nature.com/articles/eye2015252
    """
    wavelength = float(wavelength)
    if wavelength >= 380 and wavelength < 750:
        A = alpha_spectrum
    else:
        A = alpha_outside
    if wavelength < 380:
        wavelength = 380.0
    if wavelength > 750:
        wavelength = 750.0
    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    return (R, G, B, A)


def prepare_rainbow_background(min_nm, max_nm, y_min, y_max):
    """
    Visible range reported differently and changes over the lifetime. https://www.nature.com/articles/eye2015252
    """
    wl = np.linspace(min_nm, max_nm, 500)
    colors = np.array([wavelength_to_rgb(w) for w in wl])

    # Create image with shape (height, width, RGBA)
    height = 20  # thickness of the background strip in pixels
    img = np.tile(colors[:, np.newaxis, :], (1, height, 1)).transpose(1, 0, 2)

    extent = (min_nm, max_nm, y_min, y_max)
    return img, extent


def table_creation(name, step_efficiencies_dict, overall_efficiency=100):
    """
    Create table with standard format for each process step with efficiency.

    Args:
        name (str): name of the overall process step (e.g. Photosynthesis)
        step_efficiencies_dict (dict): dict of step name, step efficiency value pairs
        overall_efficiency (float): Efficiency to start from. Defaults to 100

    Returns:
        pd.DataFrame: Table with standard format.
    """
    # Initialize the DataFrame
    df = pd.DataFrame(
        {
            "Partial_efficiency_in_%": {name: np.nan},
            "Overall_efficiency_in_%": {name: np.nan},
        }
    )
    df.rename_axis(index="Process", inplace=True)

    # Fill the DataFrame using a loop
    for process, step_efficiency in step_efficiencies_dict.items():
        partial_efficiency = step_efficiency * 100
        overall_efficiency *= step_efficiency
        df.loc[process] = [round(partial_efficiency, 6), round(overall_efficiency, 6)]

    return df


def kJ_per_mol_to_kJ_per_kg(kJ_per_mol, M):
    """
    Convert energy from kJ/mol to kJ/kg.

    Args:
        kJ_per_mol (float): Energy value in kJ/mol.
        M (float): Molar mass of the component in g/mol.

    Returns:
        float: Energy value in kJ/kg.
    """
    return kJ_per_mol * (1000 / M)


def kJ_per_kg_to_kWh_per_t(kJ_per_kg):
    """
    Convert energy from kJ/kg to kWh/ton.

    Args:
        kJ_per_kg (float): Energy value in kJ/kg.

    Returns:
        float: Energy value in kWh/ton.
    """
    return kJ_per_kg * (1000 / 3600)


def kJ_per_mol_to_kWh_per_t(kJ_per_mol, M):
    """
    Convert energy from kJ/mol to kWh/ton using the molar mass.

    Args:
        kJ_per_mol (float): Energy value in kJ/mol.
        M (float): Molar mass of the component in g/mol.

    Returns:
        float: Energy value in kWh/ton.
    """
    kJ_per_kg = kJ_per_mol_to_kJ_per_kg(kJ_per_mol, M)
    return kJ_per_kg_to_kWh_per_t(kJ_per_kg)


def kJ_per_kg_to_kJ_per_mol(kJ_per_kg, M):
    """
    Convert energy from kJ/kg to kJ/mol.

    Args:
        kJ_per_kg (float): Energy value in kJ/kg.
        M (float): Molar mass of the component in g/mol.

    Returns:
        float: Energy value in kJ/mol.
    """
    return kJ_per_kg * (M / 1000)


def kWh_per_t_to_kJ_per_kg(kWh_per_t):
    """
    Convert energy from kWh/ton to kJ/kg.

    Args:
        kWh_per_t (float): Energy value in kWh/ton.

    Returns:
        float: Energy value in kJ/kg.
    """
    return kWh_per_t * (3600 / 1000)


def kWh_per_t_to_kJ_per_mol(kWh_per_t, M):
    """
    Convert energy from kWh/ton to kJ/mol using the molar mass.

    Args:
        kWh_per_t (float): Energy value in kWh/ton.
        M (float): Molar mass of the component in g/mol.

    Returns:
        float: Energy value in kJ/mol.
    """
    kJ_per_kg = kWh_per_t_to_kJ_per_kg(kWh_per_t)
    return kJ_per_kg_to_kJ_per_mol(kJ_per_kg, M)
