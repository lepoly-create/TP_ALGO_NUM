import numpy as np
import matplotlib.pyplot as plt
def gauss_seidel(A, b, x0=None, tol=1e-10, max_iter=1000):
    """
    Méthode itérative de Gauss‑Seidel.
    Retourne la solution approchée, le nombre d'itérations et la liste des erreurs.
    """
    A = A.astype(float)
    b = b.astype(float)
    n = len(b)
    if x0 is None:
        x = np.zeros(n)
    else:
        x = x0.copy()
    if np.any(np.diag(A) == 0):
        raise ValueError("Élément diagonal nul")
    errors = []
    for it in range(max_iter):
        x_old = x.copy()
        for i in range(n):
            s = b[i] - np.dot(A[i, :i], x[:i]) - np.dot(A[i, i + 1:], x_old[i + 1:])
            x[i] = s / A[i, i]
        err = np.linalg.norm(x - x_old, np.inf)
        errors.append(err)
        if err < tol:
            return x, it + 1, errors
    return x, max_iter, errors