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


if __name__ == "__main__":
    import numpy as np

    # Exemple pour vérifier la décomposition LU et la résolution A x = b
    A = np.array([[4.0, 3.0, 2.0],
                   [2.0, 1.0, -1.0],
                   [6.0, 5.0, 4.0]])
    b = np.array([9.0, 1.0, 15.0])

    L, U = doolittle_lu(A)
    print("L:")
    print( np.array2string(L, precision=4, suppress_small=True))
    print("U:")
    print(np.array2string(U, precision=4, suppress_small=True))

    # Résolution L y = b (substitution avant)
    n = A.shape[0]
    y = np.zeros(n)
    for i in range(n):
        s = 0.0
        for j in range(i):
            s += L[i, j] * y[j]
        y[i] = (b[i] - s) / L[i, i]

    # Résolution U x = y (substitution arrière)
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        s = 0.0
        for j in range(i + 1, n):
            s += U[i, j] * x[j]
        x[i] = (y[i] - s) / U[i, i]

    print("Solution (Doolittle):", [f"{v:.4g}" for v in x])