import numpy as np
import matplotlib.pyplot as plt
def lagrange_interpolation(x_data, y_data, x_eval):
    """
    Évalue le polynôme d'interpolation de Lagrange aux points x_eval.
    """
    x_data = np.asarray(x_data)
    y_data = np.asarray(y_data)
    n = len(x_data)
    result = np.zeros_like(x_eval, dtype=float)
    for i in range(n):
        # construction du i‑ème polynôme de base
        Li = np.ones_like(x_eval)
        for j in range(n):
            if j != i:
                Li *= (x_eval - x_data[j]) / (x_data[i] - x_data[j])
        result += y_data[i] * Li
    return result