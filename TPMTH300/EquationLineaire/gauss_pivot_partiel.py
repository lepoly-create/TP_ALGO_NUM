import numpy as np
import matplotlib.pyplot as plt
def gauss_pivot_partiel(A, b):
    """
    Élimination de Gauss avec pivot partiel (recherche du pivot dans la colonne).
    """
    A = A.astype(float, copy=True)
    b = b.astype(float, copy=True)
    n = len(b)
    for k in range(n - 1):
        # Recherche du pivot maximal en colonne k
        pivot_row = np.argmax(np.abs(A[k:, k])) + k
        if A[pivot_row, k] == 0:
            raise ValueError("Matrice singulière")
        # Échange des lignes
        if pivot_row != k:
            A[[k, pivot_row]] = A[[pivot_row, k]]
            b[[k, pivot_row]] = b[[pivot_row, k]]
        # Élimination
        for i in range(k + 1, n):
            factor = A[i, k] / A[k, k]
            A[i, k:] -= factor * A[k, k:]
            b[i] -= factor * b[k]
    # Substitution arrière
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (b[i] - np.dot(A[i, i + 1:], x[i + 1:])) / A[i, i]
    return x

if __name__ == "__main__":
    import numpy as _np

    # Exemple pour tester le pivot partiel : premier pivot très petit
    A = _np.array([[1e-12, 1.0, 2.0],
                   [1.0, 2.0, 3.0],
                   [4.0, 5.0, 6.0]])
    b = _np.array([3.0, 6.0, 15.0])

    sol = gauss_pivot_partiel(A, b)
    # Affichage compact
    print("Solution (pivot partiel):", [f"{float(v):.4g}" for v in sol])