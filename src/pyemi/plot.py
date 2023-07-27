import matplotlib.pyplot as plt
import pandas as pd

def plot_dimensions():

  dims = {}

  dims["in_to_pt"] = 72.27
  dims["in_to_cm"] = 2.54
  dims["cm_to_pt"] = 28.45

  dims["A4"] = {}
  dims["A4"]["in"] = {}
  dims["A4"]["pt"] = {}
  dims["A4"]["cm"] = {}
  dims["A4"]["in"]["width"] = 8.27
  dims["A4"]["in"]["height"] = 11.69
  dims["A4"]["pt"]["width"] = 595
  dims["A4"]["pt"]["height"] = 842
  dims["A4"]["cm"]["width"] = 21
  dims["A4"]["cm"]["height"] = 29.7

  return dims
