import numpy as np
import matplotlib.pyplot as plt
def jacobi(A, b, x0=None, tol=1e-10, max_iter=1000):
    """
    Méthode itérative de Jacobi.
    Retourne la solution approchée, le nombre d'itérations et la liste des erreurs.
    """
    A = A.astype(float)
    b = b.astype(float)
    n = len(b)
    if x0 is None:
        x = np.zeros(n)
    else:
        x = x0.copy()
    D = np.diag(A)
    if np.any(D == 0):
        raise ValueError("Élément diagonal nul")
    x_old = x.copy()
    errors = []
    for it in range(max_iter):
        x_new = np.zeros(n)
        for i in range(n):
            s = b[i] - np.dot(A[i, :i], x_old[:i]) - np.dot(A[i, i + 1:], x_old[i + 1:])
            x_new[i] = s / D[i]
        err = np.linalg.norm(x_new - x_old, np.inf)
        errors.append(err)
        if err < tol:
            return x_new, it + 1, errors
        x_old = x_new
    return x_old, max_iter, errors