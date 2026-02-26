import numpy as np
import matplotlib.pyplot as plt
def doolittle_lu(A):
    """
    Factorisation LU de Doolittle : L triangulaire inférieure (diagonale unité),
    U triangulaire supérieure.
    Retourne (L, U).
    """
    A = A.astype(float, copy=True)
    n = A.shape[0]
    L = np.eye(n)
    U = np.zeros_like(A)
    for k in range(n):
        # Calcul de U[k, k:]
        for j in range(k, n):
            U[k, j] = A[k, j] - np.dot(L[k, :k], U[:k, j])
        # Calcul de L[k+1:, k]
        for i in range(k + 1, n):
            L[i, k] = (A[i, k] - np.dot(L[i, :k], U[:k, k])) / U[k, k]
    return L, U