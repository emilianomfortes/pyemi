from scipy import interpolate
from tqdm import tqdm
import cantera as ct
import pandas as pd
import numpy as np
import copy
import json
import os

#
# Alya + FlameGen
#
def read_alya_table(path_tab, table_type="laminar"):
    assert (table_type == "laminar") or (
        table_type == "prop"
    ), f"Valid inputs for type of table are 'laminar' or 'prop', got {table_type} instead."

    # Read file to get raw lines
    with open(path_tab, "r") as fprop:
        lines_prop = fprop.readlines()

    # Find last index before the important data begins
    final_iline = 0
    for iline, line in enumerate(lines_prop):
        if line.startswith("$"):
            final_iline = iline

    final_iline += 2
    print(final_iline)

    #df_tab = pd.read_csv(path_tab, skiprows=final_iline, delim_whitespace=True)
    df_tab = pd.read_csv(path_tab, skiprows=final_iline, sep='\s+')  # Updated to new pandas version
    return df_tab

def write_table_general(
    path, table_name, props_to_tab, controlvars_list, table_write, controlvars, chunk_data_per_line=1000000,
):

    # controlvarlist should be a list of [Z, C, I] ordered

    open(path / f"{table_name}.dat", "a")

    df = table_write[controlvars + props_to_tab]
    # Drop columns duplicated
    df = df.loc[:, ~df.columns.duplicated()]

    for col in df.columns:
        if df[col].dtype == "float64":
            df[col] = df[col].map(lambda x: "{:.8e} ".format(x))

    with open(path / f"{table_name}.dat", "w") as f:

        #
        # Write header
        #
        line = f"N_VALUES: "
        for icv, cv in enumerate(controlvars_list):
            if icv == 0:
                line += f"{len(cv)}"
            else:
                line += f", {len(cv)}"
        f.write(line)
        f.write(f"\n")
        f.write(f"N_VARIABLES: {len(props_to_tab)}\n")
        f.write(f"\n")
        f.write(f"N_CHUNK_DATA_PER_LINE: {chunk_data_per_line}\n")
        f.write(f"\n")
        f.write(f"\n")

        for icv, cv in enumerate(controlvars_list):
            f.write(f"{controlvars[icv]} SUBDIVISION:\n")
            cv_to_write = cv	
            for i in range(len(cv_to_write) // chunk_data_per_line + 1):
                cv_sub_write = cv_to_write[
                    i * chunk_data_per_line : (i + 1) * chunk_data_per_line
                ]
                string = ""
                for cc in cv_sub_write:
                    string += "{:.8e} ".format(cc) + " "
                string += "\n"
                f.write(string)

            f.write(f"\n")
            f.write(f"\n")

        for icv, cv in enumerate(controlvars):
            icv += 1
            f.write(f"$$     {icv}: {cv}     coordinate #{icv}\n")
        for iprop, prop in enumerate(props_to_tab):
            iprop += 1
            f.write(f"$$     {iprop+icv}: {prop}     variable #{iprop}\n")

        f.write(f"\n")
        string = ""
        for col in list(df.columns):
            string += col + "     "
        string += "\n"
        f.write(string)

        for i in tqdm(range(len(df))):
            string = ""
            for col in df.columns:
                string += df[col].iloc[i] + " "
            string += "\n"
            f.write(string)

def read_alya_table(path_tab, table_type="laminar"):
    assert (table_type == "laminar") or (
        table_type == "prop"
    ), f"Valid inputs for type of table are 'laminar' or 'prop', got {table_type} instead."

    # Read file to get raw lines
    with open(path_tab, "r") as fprop:
        lines_prop = fprop.readlines()

    # Find last index before the important data begins
    final_iline = 0
    for iline, line in enumerate(lines_prop):
        if line.startswith("$"):
            final_iline = iline

    final_iline += 2
    print(final_iline)

    #df_tab = pd.read_csv(path_tab, skiprows=final_iline, delim_whitespace=True)
    df_tab = pd.read_csv(path_tab, skiprows=final_iline, sep='\s+')  # Updated to new pandas version
    return df_tab

# 
# FlameGen functions robbed
#
def selectGridMethod(method):
    """Create numpy array based on different strategies.

    Parameters
    ----------
    method
        The input of this function. It can take various input types.

    Returns
    -------
    np.ndarray
        The resulting numpy array.

    Notes
    -----
    If the `method` is a list, it convets it to a numpy array.

    If the `method` is a dictionary, it is processed as follows:

    -  if `method["type"]` is `"linspace"`:

    .. math::
        min + (max - min) * [0..1]

    -  if `method["type"]` is `"power"`:

    .. math::
        min + (max - min) * [0..1]^{pow}

    -  if `method["type"]` is `"geom"`: Generates `method["n"]` elements with a geometrically growing step size with growth rate `method["r"]`.
    -  if `method["type"]` is `"attractor"`: Generates `method["n"]` elements concentrated around `method["loc"]` with a geometrically growing step size towards `method["max"]` and `method["min"]` with growth rate `method["r"]`. The number of elements on the two sides of `method["loc"]` is adjusted such, that teh step size is smooth at `method["loc"]`.
    -  if `method["type"]` is `"multipleAttractor"`: Based on a list of locations: `method["locList"]` and a list of manually defined section lengths: `method["sectionLengths"]` it uses the attractor function multiple times.
    """
    if isinstance(method, float):
        return np.array([method])

    elif isinstance(method, list):
        if isinstance(method[0], float) or isinstance(method[0], int):
            ret = [float(m) for m in method]
            return np.array(ret)

    elif isinstance(method, dict):
        if method["type"] == "linspace":
            if method["n"] == 1:
                return np.array([method["min"]])
            elif method["n"] == 2:
                return np.array([method["min"], method["max"]])
            else:
                return np.linspace(method["min"], method["max"], method["n"])

        elif method["type"] == "power":
            return method["min"] + (method["max"] - method["min"]) * np.power(
                np.linspace(0.0, 1.0, method["n"]), method["pow"]
            )

        elif method["type"] == "geom":
            r = method["r"]
            n = method["n"]
            if r != 1.0:
                dLeft = (r - 1) / (r ** (n - 1) - 1)
            else:
                dLeft = 1.0 / (n - 1)

            returnArray = [0]
            step = dLeft
            for i in range(n - 1):
                returnArray.append(returnArray[-1] + step)
                step *= r

            return method["min"] + (method["max"] - method["min"]) * np.array(
                returnArray
            )

        elif method["type"] == "attractor":
            r = method["r"]
            n = method["n"]
            nLow = 1
            done = False

            for i in range(n):
                partLow = selectGridMethod(
                    method={
                        "type": "geom",
                        "n": nLow + 1,
                        "min": method["min"],
                        "max": method["loc"],
                        "r": 1 / r,
                    }
                )
                partHigh = selectGridMethod(
                    method={
                        "type": "geom",
                        "n": n - nLow,
                        "min": method["loc"],
                        "max": method["max"],
                        "r": r,
                    }
                )
                stepLow = partLow[-1] - partLow[-2]
                stepHigh = partHigh[1] - partHigh[0]
                if stepLow > stepHigh:
                    nLow += 1
                else:
                    nLow -= 1
                    done = True

                nLow = max(nLow, 1)
                if done:
                    break

            return np.append(partLow[0:-1], partHigh)

        elif method["type"] == "multipleAttractor":
            locList = method["locList"]
            sectionLengths = method["sectionLengths"]
            nSections = len(locList)

            minList = []
            maxList = []

            for i in range(nSections):
                if i == 0:
                    minList.append(method["min"])
                else:
                    minList.append((locList[i - 1] + locList[i]) / 2)

                if i == nSections - 1:
                    maxList.append(method["max"])
                else:
                    maxList.append((locList[i] + locList[i + 1]) / 2)

            parts = np.array([])
            for i in range(nSections):
                p = selectGridMethod(
                    method={
                        "type": "attractor",
                        "n": sectionLengths[i],
                        "min": minList[i],
                        "max": maxList[i],
                        "r": method["r"],
                        "loc": locList[i],
                    }
                )
                parts = np.append(parts, p)

            parts = np.unique(parts)

            return parts
