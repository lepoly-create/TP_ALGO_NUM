import numpy as np
import matplotlib.pyplot as plt
def cholesky(A):
    """
    Factorisation de Cholesky pour une matrice symétrique définie positive.
    Retourne la matrice triangulaire inférieure L telle que A = L @ L.T.
    """
    A = A.astype(float, copy=True)
    n = A.shape[0]
    L = np.zeros_like(A)
    for i in range(n):
        for j in range(i + 1):
            if i == j:
                s = A[i, i] - np.dot(L[i, :i], L[i, :i])
                if s <= 0:
                    raise ValueError("Matrice non définie positive")
                L[i, i] = np.sqrt(s)
            else:
                L[i, j] = (A[i, j] - np.dot(L[i, :j], L[j, :j])) / L[j, j]
    return L


def cholesky_solve(A, b):
    """
    Calcule la décomposition de Cholesky puis résout A x = b.
    Retourne (L, x) où A = L @ L.T.
    """
    try:
        L = cholesky(A)
    except ValueError:
        # Matrice non définie positive -> retourner None pour indiquer l'échec
        return None, None
    n = A.shape[0]
    # Substitution avant L y = b
    y = np.zeros(n)
    for i in range(n):
        y[i] = (b[i] - np.dot(L[i, :i], y[:i])) / L[i, i]
    # Substitution arrière L.T x = y
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - np.dot(L.T[i, i + 1:], x[i + 1:])) / L[i, i]
    return L, x