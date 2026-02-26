import numpy as np
import matplotlib.pyplot as plt
from .difference_divisees import newton_divided_differences
def newton_interpolation(x_data, y_data, x_eval):
    """
    Évalue le polynôme d'interpolation de Newton (forme des différences divisées).
    """
    coef = newton_divided_differences(x_data, y_data)
    n = len(x_data)
    # Évaluation par l'algorithme de Horner
    result = coef[-1] * np.ones_like(x_eval)
    for i in range(n - 2, -1, -1):
        result = result * (x_eval - x_data[i]) + coef[i]
    return result