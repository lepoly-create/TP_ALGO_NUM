import numpy as np
import matplotlib.pyplot as plt
def least_squares_poly(x_data, y_data, deg):
    """
    Ajuste un polynôme de degré `deg` aux données par moindres carrés.
    Retourne les coefficients du polynôme (ordre croissant) et la matrice de Vandermonde.
    """
    x_data = np.asarray(x_data)
    y_data = np.asarray(y_data)
    # Construction de la matrice de Vandermonde
    V = np.vander(x_data, deg + 1, increasing=True)
    # Résolution des équations normales : V.T @ V @ c = V.T @ y
    # Utilisation de numpy.linalg.lstsq pour la stabilité
    coeffs, residuals, rank, s = np.linalg.lstsq(V, y_data, rcond=None)
    return coeffs, V