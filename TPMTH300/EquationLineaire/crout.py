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