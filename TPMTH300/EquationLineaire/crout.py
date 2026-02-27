import numpy as np
import matplotlib.pyplot as plt
def crout_lu(A):
    """
    Factorisation LU de Crout : L triangulaire inférieure (diagonale quelconque),
    U triangulaire supérieure avec diagonale unité.
    Retourne (L, U).
    """
    A = A.astype(float, copy=True)
    n = A.shape[0]
    L = np.zeros_like(A)
    U = np.eye(n)
    for k in range(n):
        # Calcul de L[k:, k]
        for i in range(k, n):
            L[i, k] = A[i, k] - np.dot(L[i, :k], U[:k, k])
        # Calcul de U[k, k+1:]
        for j in range(k + 1, n):
            U[k, j] = (A[k, j] - np.dot(L[k, :k], U[:k, j])) / L[k, k]
    return L, U


def crout_solve(A, b):
    """
    Factorise A avec Crout puis résout A x = b.
    Retourne (L, U, x).
    """
    L, U = crout_lu(A)
    n = A.shape[0]
    # Substitution avant L y = b
    y = np.zeros(n)
    for i in range(n):
        y[i] = (b[i] - np.dot(L[i, :i], y[:i])) / L[i, i]
    # Substitution arrière U x = y
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - np.dot(U[i, i + 1:], x[i + 1:])) / U[i, i]
    return L, U, x
