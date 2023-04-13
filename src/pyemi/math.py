import numpy as np

def increasing_linspace(start, stop, num, base=10):
    """Returns a linspace with a base other than 10."""
    index = np.logspace(0, 1, num, base=base)
    return np.interp(index, [index.min(), index.max()], [start, stop])

def decreasing_linspace(start, stop, num, base=2):
    powers = np.linspace(0, 1, num)
    # The formula below ensures that the spacing decreases as the index increases
    spacings = (1 - base ** (-powers)) / (1 - base ** (-1))
    return start + spacings * (stop - start)

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