import numpy as np


def ResolutionSystTriSup(Taug):
    A = np.copy(Taug)
    n, m = np.shape(A)
    X = np.zeros(n)
    for k in range(n - 1, -1, -1):
        S = 0
        for j in range(k + 1, n):
            S = S + A[k, j] * X[j]
        X[k] = (A[k, -1] - S) / A[k, k]
    return X


def precision(A, X, B):
    n = len(B)
    B = np.reshape(B, (1, n))
    X = np.ravel(X)
    B = np.ravel(B)
    a = np.dot(A, X) - B
    return np.linalg.norm(a)