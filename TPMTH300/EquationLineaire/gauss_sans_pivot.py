import numpy as np
import matplotlib.pyplot as plt

def gauss_sans_pivot(A, b):
    """
    Résolution de Ax = b par élimination de Gauss sans pivot.
    Retourne la solution x.
    """
    A = A.astype(float, copy=True)
    b = b.astype(float, copy=True)
    n = len(b)
    # Étape de triangularisation
    for k in range(n - 1):
        if A[k, k] == 0:
            raise ValueError("Pivot nul rencontré. Essayez avec pivot partiel.")
        for i in range(k + 1, n):
            factor = A[i, k] / A[k, k]
            A[i, k:] -= factor * A[k, k:]
            b[i] -= factor * b[k]
    # Substitution arrière
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (b[i] - np.dot(A[i, i + 1:], x[i + 1:])) / A[i, i]
    return x

# if __name__ == "__main__":
#     # petit test rapide
#     A = [[2.0, 1.0, -1.0],
#          [ -3.0, -1.0, 2.0], 
#          [ -2.0, 1.0, 2.0]
#          ]
#     b = [8.0, -11.0, -3.0]
#     sol = gauss_sans_pivot(A, b)
#     print("Solution:", [f"{v:.2g}" for v in sol])
