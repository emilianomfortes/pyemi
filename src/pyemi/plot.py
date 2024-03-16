from copy import deepcopy
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

def default_mpl_style():

  mpl.rcParams['text.usetex'] = True
  mpl.rcParams['font.family'] = 'serif'

  mpl.rcParams['font.size'] = 8.0

  #
  # Figure
  #
  mpl.rcParams['figure.titlesize'] = 8.0


  #
  # Axes
  #
  mpl.rcParams['axes.titlesize'] = 8.0

  #
  # Legend
  #
  mpl.rcParams['legend.fontsize'] = 8.0

  #
  # Lines
  #
  mpl.rcParams['lines.linewidth'] = 1.0
  mpl.rcParams['lines.markersize'] = 3.0

  #
  # Ticks
  #
  mpl.rcParams['xtick.labelsize'] = 8.0
  mpl.rcParams['ytick.labelsize'] = 8.0
  mpl.rcParams['xtick.direction'] = "inout"
  mpl.rcParams['ytick.direction'] = "inout" 

  rcParams = deepcopy(mpl.rcParams)
  return rcParams

def plot_dimensions():

  dims = {}

  dims["in_to_pt"] = 72.27
  dims["in_to_cm"] = 2.54
  dims["cm_to_pt"] = 28.45
  dims["pt_to_in"] = 1.0/72.27

  #
  # Page sizes
  #
  
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

  #
  # Journal LaTeX parameters
  #

  # Elsevier: https://www.elsevier.com/authors/policies-and-guidelines/artwork-and-media-instructions/artwork-sizing
  dims["Journal"] = {}
  dims["Journal"]["elsevier-twocolumn-A4"] = {}
  dims["Journal"]["elsevier-twocolumn-A4"]["pt"] = {}
  dims["Journal"]["elsevier-twocolumn-A4"]["pt"]["margin"] = 34.51605
  dims["Journal"]["elsevier-twocolumn-A4"]["pt"]["top"] = 52.74657

  dims["Journal"]["elsevier-twocolumn-A4"]["Figure"] = {}
  dims["Journal"]["elsevier-twocolumn-A4"]["Figure"]["Single column_width_pt"] = 255
  dims["Journal"]["elsevier-twocolumn-A4"]["Figure"]["1.5 column_width_pt"] = 297
  dims["Journal"]["elsevier-twocolumn-A4"]["Figure"]["Double column_width_pt"] = 539

  # International Combustion Symposium 40th edition
  # https://www.combustioninstitute.org/wp-content/uploads/2023/08/Instructions-to-Authors-for-Manuscript-Preparation.40thISOC.pdf
  journal = "40th International Symposium on Combustion"
  dims["Journal"][journal] = {}
  dims["Journal"][journal]["Figure"] = {}
  dims["Journal"][journal]["Figure"]["Single_width_in"] = 3.5416666667
  dims["Journal"][journal]["Figure"]["Double_width_in"] = 7.4861111111

  
  return dims
