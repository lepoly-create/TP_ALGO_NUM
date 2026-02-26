import numpy as np
import matplotlib.pyplot as plt
def gauss_pivot_total(A, b):
    """
    Élimination de Gauss avec pivot total (recherche dans toute la sous‑matrice).
    Les permutations de colonnes sont mémorisées pour réordonner la solution.
    """
    A = A.astype(float, copy=True)
    b = b.astype(float, copy=True)
    n = len(b)
    col_order = np.arange(n)   # ordre initial des colonnes
    for k in range(n - 1):
        # Recherche du pivot maximal dans la sous‑matrice (k..n-1, k..n-1)
        max_val = 0
        max_i, max_j = k, k
        for i in range(k, n):
            for j in range(k, n):
                if abs(A[i, j]) > max_val:
                    max_val = abs(A[i, j])
                    max_i, max_j = i, j
        if max_val == 0:
            raise ValueError("Matrice singulière")
        # Échange des lignes
        if max_i != k:
            A[[k, max_i]] = A[[max_i, k]]
            b[[k, max_i]] = b[[max_i, k]]
        # Échange des colonnes
        if max_j != k:
            A[:, [k, max_j]] = A[:, [max_j, k]]
            col_order[[k, max_j]] = col_order[[max_j, k]]
        # Élimination
        for i in range(k + 1, n):
            factor = A[i, k] / A[k, k]
            A[i, k:] -= factor * A[k, k:]
            b[i] -= factor * b[k]
    # Substitution arrière
    x_temp = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x_temp[i] = (b[i] - np.dot(A[i, i + 1:], x_temp[i + 1:])) / A[i, i]
    # Réordonner les composantes selon l'ordre initial des colonnes
    x = np.zeros(n)
    for i in range(n):
        x[col_order[i]] = x_temp[i]
    return x


if __name__ == "__main__":
    import numpy as _np

    # Exemple conçu pour nécessiter un pivot total (colonne 0 contient un zéro)
    A = _np.array([[0.0, 2.0, -1.0],
                   [1.0, -2.0, 3.0],
                   [2.0, 3.0, 1.0]])
    b = _np.array([3.0, -1.0, 4.0])

    x = gauss_pivot_total(A, b)
    print("Solution (pivot total):", [f"{v:.4g}" for v in x])