import numpy as np
import matplotlib.pyplot as plt
def gauss_jordan(A, b):
    """
    Méthode de Gauss‑Jordan : transforme la matrice augmentée en forme diagonale.
    Retourne la solution x.
    """
    A = A.astype(float, copy=True)
    b = b.astype(float, copy=True)
    n = len(b)
    M = np.column_stack([A, b])   # matrice augmentée
    for k in range(n):
        # Pivot partiel
        pivot_row = np.argmax(np.abs(M[k:, k])) + k
        if M[pivot_row, k] == 0:
            raise ValueError("Matrice singulière")
        if pivot_row != k:
            M[[k, pivot_row]] = M[[pivot_row, k]]
        # Normalisation de la ligne pivot
        M[k, k:] /= M[k, k]
        # Élimination dans toutes les autres lignes
        for i in range(n):
            if i != k:
                factor = M[i, k]
                M[i, k:] -= factor * M[k, k:]
    # La dernière colonne de M contient la solution
    return M[:, -1]


if __name__ == "__main__":
    import numpy as _np

    # Exemple pour tester Gauss-Jordan (nécessite permutations)
    A = _np.array([[0.0, 2.0, 1.0],
                   [1.0, -2.0, 3.0],
                   [2.0, 3.0, 1.0]])
    b = _np.array([3.0, -1.0, 4.0])

    x = gauss_jordan(A, b)
    print("Solution (Gauss-Jordan):", [f"{v:.4g}" for v in x])