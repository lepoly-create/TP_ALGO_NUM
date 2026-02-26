import numpy as np
import matplotlib.pyplot as plt
def newton_divided_differences(x_data, y_data):
    """
    Calcule les différences divisées de Newton.
    Retourne un tableau des coefficients (les différences divisées diagonales).
    """
    n = len(x_data)
    coef = y_data.copy().astype(float)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / (x_data[i] - x_data[i - j])
    return coef