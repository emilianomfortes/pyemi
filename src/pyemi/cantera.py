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
# Flame properties
#
def laminar_flame_speed(v_u, v_b, rho_u, rho_b, stationary=False):
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
    for (_species,_weights) in zip(Yc_species, Yc_weights):
        Yc_array += _weights*df_flamelet[_species].to_numpy()
    return Yc_array


def thermal_flame_thickness(flame):

    if isinstance(flame, pd.DataFrame):
        max_grad = np.max(grid_deriv(flame["x(m)"].to_numpy(), flame["T(K)"].to_numpy()))
        T_u = flame["T(K)"].to_numpy()[0]
        T_b = flame["T(K)"].to_numpy()[-1]
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
		Zel = Zel + gas.n_atoms(k, element) * mass_frac[ik] * Wel / gas.molecular_weights[ik]
	return Zel

def bilger_beta(mass_frac, gas, dict_bilger = {"C": 2.0,"S": 2.0,"H": 0.5, "O": -1.0}):
	beta = 0.0
	for element in dict_bilger:
		if element in gas.species_names:
			beta = beta + dict_bilger[element] * elemental_Z(mass_frac, element, gas) / gas.atomic_weight(element)
	return beta

def bilger_Z(mass_frac, gas_fuel, gas_ox, dict_bilger = {"C": 2.0,"S": 2.0,"H": 0.5, "O": -1.0}):
    Z_mix = bilger_beta(mass_frac, gas_fuel, dict_bilger)
    #print(Z_mix)
    Z_fuel = bilger_beta(gas_fuel.Y, gas_fuel, dict_bilger)
    Z_ox = bilger_beta(gas_ox.Y, gas_ox, dict_bilger)
    Zbilger = (Z_mix - Z_ox) / (Z_fuel - Z_ox)
    return Zbilger

#
# Mixture average tabulation
#
def get_prop_from_gas(prop, gas, gas_oxidizer, gas_fuel, YcWeights, mech):
	if prop == "Yc":
		Ylean = gas.Y
		prop_val = 0.0
		for k, alphak in YcWeights.items():
			ik = gas.species_index(k)
			prop_val += alphak* Ylean[ik]
	elif prop == "omegaYc":
		omegas = gas.net_production_rates
		prop_val = 0.0
		for k, alphak in YcWeights.items():
			ik = gas.species_index(k)
			prop_val += alphak* omegas[ik]
	elif prop == "omegaYc_Yc":
		omegas = gas.net_production_rates
		prop_val = 0.0
		for k, alphak in YcWeights.items():
			ik = gas.species_index(k)
			prop_val += alphak* omegas[ik] * omegas[ik]

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
