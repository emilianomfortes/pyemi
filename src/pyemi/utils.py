from pathlib import Path
import numpy as np
import json
import yaml
import os


def check_make_dir(path):
    # If folder doesn't exist, then create it.
    check_path = os.path.isdir(path)
    if not check_path:
        os.makedirs(path)
        print("created folder : ", path)

    else:
        print(path, "folder already exists.")


def read_json(path_json):
    with open(path_json) as f:
        dict_json = json.load(f)
    return dict_json


def read_yaml(path_file) -> dict:
    """Reads a yaml file as a dictionary"""
    with open(path_file, "r") as f:
        _file = yaml.safe_load(f)
    return _file


def find_line(str_find, lines):
    for idx, line in enumerate(lines):
        if str_find in line:
            return idx


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
