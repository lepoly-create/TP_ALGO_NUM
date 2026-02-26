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