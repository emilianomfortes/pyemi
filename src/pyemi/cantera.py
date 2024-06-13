from pyemi.math import grid_deriv
from functools import reduce
from pathlib import Path
import cantera as ct
import pandas as pd
import numpy as np
import json
import copy
import os


#
# Flame
#
def FreeFlame_good_convergence(
    gas,
    initial_grid,
    timestep=1e-5, 
    timestep_tries=[2, 5, 10, 20, 40, 50], 
    ratio=2, 
    slope=0.001, 
    curve=0.01, 
    max_grid_points=10000, 
    tol_ss = [1.0e-7, 1.0e-12], 
    tol_ts = [1.0e-7, 1.0e-12],
    max_jac_age=10,
    ):

    flame = ct.FreeFlame(gas, grid=initial_grid)
    flame.flame.set_steady_tolerances(default=tol_ss)
    flame.flame.set_transient_tolerances(default=tol_ts)
    flame.energy_enabled = True
    flame.set_max_jac_age(max_jac_age, max_jac_age)
    flame.set_max_grid_points(flame.flame, max_grid_points)
    flame.set_time_step(timestep, timestep_tries)
    flame.set_refine_criteria(ratio=ratio, slope=slope, curve=curve)

    return flame

#
# Flame properties
#
def laminar_flame_speed(v_u, v_b, rho_u, rho_b, df=None, stationary=False):

    if isinstance(df, pd.DataFrame):
        v_col = None
        rho_col = None

        potential_rho_cols = ["DENSI", "rho", "RHO", "rho(kg/m3)", "rho(kg/m^3)"]
        for col in potential_rho_cols:
            if col in list(df.columns):
                rho_col = col
        potential_vel_cols = ["u_x", "u(m/s)", "u"]
        for col in potential_vel_cols:
            if col in list(df.columns):
                vel_col = col

        v_u = df[vel_col].iloc[0]
        v_b = df[vel_col].iloc[-1]
        rho_u = df[rho_col].iloc[0]
        rho_b = df[rho_col].iloc[-1]

    """Calculates the laminar flame speed s_l for a 1d flame."""
    if stationary:
        # Flame front doesn't move
        flame_speed = v_b * rho_b / rho_u
    else:
        # There's an inlet of unburnt gas
        flame_speed = (v_b - v_u) / (rho_u / rho_b - 1)
    return flame_speed


def calc_progress_var_df(df_flamelet, Yc_species, Yc_weights):
    Yc_array = np.zeros_like(df_flamelet["u(m/s)"].to_numpy())
    for _species, _weights in zip(Yc_species, Yc_weights):
        Yc_array += _weights * df_flamelet[_species].to_numpy()
    return Yc_array


def thermal_flame_thickness(flame):
    """Calculates the thermal flame thickness for a 1d flame."""
    
    
    if isinstance(flame, pd.DataFrame):
        
        # Physical space
        if "x(m)" in flame.columns:
            print(f"Found x(m) in columns")
            x = flame["x(m)"].to_numpy()
            
        elif "x" in flame.columns:
            print(f"Found x in columns")
            x = flame["x"].to_numpy()
        
        # Temperature
        if "T(K)" in flame.columns:
            print(f"Found T(K) in columns")
            T = flame["T(K)"].to_numpy()
        elif "T" in flame.columns:
            print(f"Found T in columns")
            T = flame["T"].to_numpy()
        
        max_grad = np.max(
            grid_deriv(x, T)
        )
        T_u = T[0]
        T_b = T[-1]
        
    else:
        max_grad = np.max(grid_deriv(flame.grid, flame.T))
        T_u = flame.T[0]
        T_b = flame.T[-1]

    flame_thickness = (T_b - T_u) / max_grad
    return flame_thickness


def elemental_Z(mass_frac, element, gas):
    iel = gas.species_index(element)
    Wel = gas.atomic_weight(element)
    Zel = 0.0
    for k in gas.species_names:
        ik = gas.species_index(k)
        Zel = (
            Zel
            + gas.n_atoms(k, element)
            * mass_frac[ik]
            * Wel
            / gas.molecular_weights[ik]
        )
    return Zel


def bilger_beta(
    mass_frac, gas, dict_bilger={"C": 2.0, "S": 2.0, "H": 0.5, "O": -1.0}
):
    beta = 0.0
    for element in dict_bilger:
        if element in gas.species_names:
            beta = beta + dict_bilger[element] * elemental_Z(
                mass_frac, element, gas
            ) / gas.atomic_weight(element)
    return beta

def Z_from_flame_df(df_flame, mech, pressure, X_fuel, X_oxidizer, species_col_prefix="Y_", Z_weights={"C": 2.0, "S": 2.0, "H": 0.5, "O": -1.0}):

    T_in = df_flame["T"].iloc[0]
    gas_oxi = ct.Solution(mech)
    gas_oxi.TPX = T_in, pressure, X_oxidizer
    gas_fuel = ct.Solution(mech)
    gas_fuel.TPX = T_in, pressure, X_fuel
    Yk = df_flame[[f"Y_{s}" for s in gas_fuel.species_names]].to_numpy().T

    Z = bilger_Z(Yk, gas_fuel, gas_oxi, Z_weights)
    return Z

def bilger_Z(
    mass_frac,
    gas_fuel,
    gas_ox,
    dict_bilger={"C": 2.0, "S": 2.0, "H": 0.5, "O": -1.0},
):
    Z_mix = bilger_beta(mass_frac, gas_fuel, dict_bilger)
    # print(Z_mix)
    Z_fuel = bilger_beta(gas_fuel.Y, gas_fuel, dict_bilger)
    Z_ox = bilger_beta(gas_ox.Y, gas_ox, dict_bilger)
    Zbilger = (Z_mix - Z_ox) / (Z_fuel - Z_ox)
    return Zbilger


def get_species_fluxes(df, gas):
    """Calculates the species fluxes from the mass fractions and the
    velocity field. The species fluxes are returned in a pandas DataFrame
    with the same index as the input DataFrame."""

    df_aux = df.copy()

    for k in gas.species_names:
        ik = gas.species_index(k)
        df_aux[f"jstar_{k}"] = (
            -1.0
            * df_aux["rho"]
            * gas.molecular_weights[ik]
            * (df_aux["W"] ** (-1))
            * df_aux[f"Dk_{k}"]
            * grid_deriv(df_aux["x"], df_aux[f"Xk_{k}"])
        )

    for k in gas.species_names:
        ik = gas.species_index(k)
        df[f"j_{k}"] = df_aux[f"jstar_{k}"].copy()
        for j in gas.species_names:
            ij = gas.species_index(j)
            df[f"j_{k}"] -= df_aux[k] * df_aux[f"jstar_{j}"]

    return df


def getZfromflux(df, gas, dict_bilger, beta_fuel, beta_oxi):
    df["j_Z"] = 0.0
    for el, wel in dict_bilger.items():
        if el in gas.species_names:
            df[f"j_Z{el}"] = 0.0
            for k in gas.species_names:
                ik = gas.species_index(k)
                Wik = gas.molecular_weights[ik]
                akp = gas.n_atoms(k, el)
                df[f"j_Z{el}"] += (akp / Wik) * df[f"j_{k}"]
            df["j_Z"] += (wel / (beta_fuel - beta_oxi)) * df[f"j_Z{el}"]

    return df["j_Z"]


def heat_release_rate(data, mech):

    if isinstance(data, pd.DataFrame):
        _data = data.copy()
    elif isinstance(data, dict):
        _data = dict(data)

    gas = ct.Solution(mech)
    W = gas.molecular_weights
    gas.TP = 298.15, 1e5
    dHf0 = {}
    for k in gas.species_names:
        ik = gas.species_index(k)
        dHf0[k] = (
            298.15
            * ct.gas_constant
            / (W[ik])
            * gas.standard_enthalpies_RT[gas.species_index(k)]
        )
    _data["HeatRel"] = 0.0
    for k in gas.species_names:
        _data["HeatRel"] += _data[f"omega_{k}"] * dHf0[k]

    # Heat release should be positive, but sometimes
    # it is negative because of author's choice of
    # the sign convention for the reaction enthalpy
    # of formation. This is a hack to fix that.
    if isinstance(data, pd.DataFrame):
        _data["HeatRel"] = _data["HeatRel"].abs()
    elif isinstance(data, dict):
        _data["HeatRel"] = np.abs(_data["HeatRel"])
    
    return _data["HeatRel"]


#
# Mixture average tabulation
#
def get_prop_from_gas(prop, gas, gas_oxidizer, gas_fuel, YcWeights, mech):
    if prop == "Yc":
        Ylean = gas.Y
        prop_val = 0.0
        for k, alphak in YcWeights.items():
            ik = gas.species_index(k)
            prop_val += alphak * Ylean[ik]
    elif prop == "omegaYc":
        omegas = gas.net_production_rates
        prop_val = 0.0
        for k, alphak in YcWeights.items():
            ik = gas.species_index(k)
            prop_val += alphak * omegas[ik]
    elif prop == "omegaYc_Yc":
        omegas = gas.net_production_rates
        prop_val = 0.0
        for k, alphak in YcWeights.items():
            ik = gas.species_index(k)
            prop_val += alphak * omegas[ik] * omegas[ik]

    elif prop.startswith("Yk"):
        k = prop[3:]
        ik = gas.species_index(k)
        prop_val = gas.Y[ik]

    elif prop.startswith("Dk"):
        k = prop[3:]
        ik = gas.species_index(k)
        prop_val = gas.mix_diff_coeffs[ik]

    elif prop == "T(K)":
        prop_val = gas.T
    elif prop == "W":
        prop_val = gas.mean_molecular_weight
    elif prop == "cp":
        prop_val = gas.cp_mass
    elif prop == "Entha":
        prop_val = gas.enthalpy_mass
    elif prop == "u":
        prop_val = 0.0
    elif prop == "rho":
        prop_val = gas.density
    elif prop == "Z":
        prop_val = bilger_Z(gas.Y, gas_fuel, gas_oxidizer)
    elif prop == "mu":
        prop_val = gas.viscosity
    elif prop == "HeatRel":
        #
        # Need enthalpies of formation
        # 25C 1bar
        #
        gasheatrel = ct.Solution(mech)
        gasheatrel.TP = 298.15, 1e5
        dHf0 = {}
        for k in gas.species_names:
            ik = gas.species_index(k)
            dHf0[k] = (
                298.15
                * ct.gas_constant
                / (gas.molecular_weights[ik])
                * gas.standard_enthalpies_RT[gas.species_index(k)]
            )
        prop_val = 0.0
        for k in gas.species_names:
            ik = gas.species_index(k)
            prop_val -= gas.net_production_rates[ik] * dHf0[k]
    elif prop == "lambda":
        prop_val = gas.thermal_conductivity

    return prop_val


#
# PANDAS
#

def EqRat_from_gas_and_df(df, fuel, oxidizer, mech, pressure):
	
	# Iterate over df rows, create a gas object, set the mass fractions,
	# pressure and temperature, and compute the equivalence ratio.

	_gas_aux = ct.Solution(mech)
	species = _gas_aux.species_names
	df["eq_rat"] = 0.0

	# T column could be T(K) or just T
	if "T(K)" not in df.columns:
		df["T(K)"] = df["T"]

	for i, row in df.iterrows():
		_gas = ct.Solution(_gas_aux)
		Y = [row[f"{k}"] for k in species]
		_gas.TPY = row["T(K)"], pressure, Y
		df.loc[i, "eq_rat"] = _gas.equivalence_ratio(
			fuel, oxidizer
		)

	return df["eq_rat"]

def EqRat_from_Zbilger(ZBilger, beta_fuel, beta_oxi):
	AF_ratio_st_Bilger = -beta_fuel/beta_oxi
	ZBilger_st = 1.0 / (1.0 + AF_ratio_st_Bilger)
	S = 1.0 / ZBilger_st - 1.0

	phi = S*((ZBilger/(1-ZBilger)))
	return phi
