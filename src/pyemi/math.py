from scipy.signal import savgol_filter
import pandas as pd
import numpy as np

def find_tangent_vector(data, col_n_x=None, col_n_y=None):

    if isinstance(data, pd.DataFrame):

        assert col_n_x is not None
        assert col_n_y is not None

        nx = data[col_n_x].to_numpy()
        ny = data[col_n_y].to_numpy()
        tx = np.zeros_like(nx)
        ty = np.zeros_like(ny)
        
    tx = ny
    ty = -nx

    return tx, ty

# -------
# Differentiation
# -------

def grid_deriv(grid, arr):
    if isinstance(grid, pd.Series):
        grid = grid.to_numpy()
    if isinstance(arr, pd.Series):
        arr = arr.to_numpy()
    
    gradi = np.zeros_like(grid)
    if len(grid) == 1:
        return gradi

    for k in range(len(grid)):
      tol = 1e-7
      if k == 0:
        h0 = np.inf
        h1 = grid[k + 1] - grid[k]
        gradi[k] = (arr[k + 1] - arr[k]) / h1

      elif k == len(grid) - 1:
        h0 = grid[k] - grid[k - 1]
        h1 = np.inf
        gradi[k] = (arr[k] - arr[k - 1]) / h0
      else:
        h0 = grid[k] - grid[k - 1]
        h1 = grid[k + 1] - grid[k]
        if (h0 < tol) or (h1 < tol):
            gradi[k] = gradi[k-1]
        else:
          gradi[k] = (
            (arr[k + 1] * h0 / (h1 * (h0 + h1)))
            + (arr[k] * (h1 - h0) / (h0 * h1))
            - (arr[k - 1] * h1 / (h0 * (h0 + h1)))
          )
    return gradi



def increasing_linspace(start, stop, num, base=10):
    """Returns a linspace with a base other than 10."""
    index = np.logspace(0, 1, num, base=base)
    return np.interp(index, [index.min(), index.max()], [start, stop])

def decreasing_linspace(start, stop, num, base=2):
    powers = np.linspace(0, 1, num)
    # The formula below ensures that the spacing decreases as the index increases
    spacings = (1 - base ** (-powers)) / (1 - base ** (-1))
    return start + spacings * (stop - start)
	
def interp_deriv(grid, arr, der=1):
	# https://en.wikipedia.org/wiki/Finite_difference_coefficient
	if isinstance(grid, pd.Series):
		grid = grid.to_numpy()
	if isinstance(arr, pd.Series):
		arr = arr.to_numpy()
	coeffs_1 = {
		-4: 1/280.0, 
		-3: -4/105.0,
		-2: 1/5.0, 
		-1: -4/5.0,
		0: 0, 
		1: 4/5.0, 
		2: -1/5.0, 
		3: 4/105.0,
		4:-1/280.0,
		}

	#arr = savgol_filter(arr, 10, 5)
	arr_interp = interp1d(grid, arr, kind="linear", fill_value="extrapolate")
	gradi = np.zeros_like(arr)

	if len(grid) == 1:
		return gradi

	for k in range(1, len(grid) - 1):
		dx_left = grid[k] - grid[k - 1]
		dx_right = grid[k + 1] - grid[k]
		h = min(dx_left, dx_right) / 4
		x = grid[k]
		for i in range(-4, 5):
			gradi[k] += coeffs_1[i] * arr_interp(x + i * h)
		gradi[k] /= h

	# Extreme left
	dx_right = grid[1] - grid[0]
	h = dx_right
	x = grid[0]
	for i in range(-4, 5):
		gradi[0] += coeffs_1[i] * arr_interp(x + i * h)
	gradi[0] /= h

	# Extreme right
	dx_left = grid[-1] - grid[-2]
	h = dx_left
	x = grid[-1]
	for i in range(-4, 5):
		gradi[-1] += coeffs_1[i] * arr_interp(x + i * h)
	gradi[-1] /= h

	return gradi

# ----------
# Integration
# ----------
def trapezoidal_integration(f, x):
    """
    Performs numerical integration of f(x) using the trapezoidal rule.
    Inputs:
        - f: function to integrate
        - x: array of uniformly spaced grid points
    Output:
        - F: array of integrated values, with the same shape as x
    """

    if isinstance(f, np.ndarray):
        f_eval = [f[i] for i in range(len(x))]
    else:
        f_eval = [f(x[i]) for i in range(len(x))]

    h = x[1] - x[0]  # grid spacing
    F = np.zeros_like(x)  # initialize integrated values array
    F[0] = 0  # initial condition
    
    for i in range(1, len(x)):
        F[i] = F[i-1] + h/2 * (f_eval[i-1] + f_eval[i])  # trapezoidal rule formula
    
    return F

# -------------
# Array scaling
# -------------
def arr_scale_linear(arr, minval, maxval):
	a = (maxval - minval) / (arr.max() - arr.min())
	b = maxval - a * arr.max()
	return a * arr + b
