from pathlib import Path
import numpy as np
import json
import yaml
import os

# 
# Directories
# 

def make_path_relative_to_dir(path, dir_name):
    # Try to find path above, if not return original path
    if not isinstance(path, Path):
        raise TypeError(f"path must be of type pathlib.Path, not {type(path)}")

    try:
        path_dir_name = find_dir_up_from_cwd(dir_name)
    except:
        print(f"Could not find relative path for dir_name={dir_name}")
        return path

    # Check if dir_name is in the parts of path
    path_parts = path.parts
    if dir_name not in path_parts:
        print(f"Could not find {dir_name} in path parts")
        return path
    else:
        # Find index of dir_name in path_parts
        i_dir_name = path_parts.index(dir_name)
        
        # Make relative path
        path_relative = Path(*path_parts[i_dir_name:])
        return path_relative


def find_dir_up_from_cwd(dir_name):
	"""Find the directory up from the current working directory"""
	cwd = Path.cwd()
	while True:
		if dir_name in [d.name for d in cwd.iterdir() if d.is_dir()]:
			print(f"Found directory {dir_name} from {cwd}, returning {cwd / dir_name}")
			return cwd / dir_name
		else:
			cwd = cwd.parent
			if cwd == Path.cwd().root:
				raise FileNotFoundError(f"Could not find directory {dir_name} from {Path.cwd()}")

def check_make_dir(path):
    # If folder doesn't exist, then create it.
    check_path = os.path.isdir(path)
    if not check_path:
        os.makedirs(path)
        print("created folder : ", path)

    else:
        print(path, "folder already exists.")

#
# Reading files
#

def read_json(path_json):
    with open(path_json) as f:
        dict_json = json.load(f)
    return dict_json


def read_yaml(path_file) -> dict:
    """Reads a yaml file as a dictionary"""
    with open(path_file, "r") as f:
        _file = yaml.safe_load(f)
    return _file

#
# Strings
#

def find_line(str_find, lines):
    for idx, line in enumerate(lines):
        if str_find in line:
            return idx


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
