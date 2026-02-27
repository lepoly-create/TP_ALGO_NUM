import numpy as np

def thomas(a, b, c, d):
    n = len(b)
    # Copy to avoid modifying inputs
    ac = a.copy().astype(float)
    bc = b.copy().astype(float)
    cc = c.copy().astype(float)
    dc = d.copy().astype(float)
    # Forward sweep
    for i in range(1, n):
        m = ac[i - 1] / bc[i - 1]
        bc[i] = bc[i] - m * cc[i - 1]
        dc[i] = dc[i] - m * dc[i - 1]
    # Back substitution
    x = np.zeros(n)
    x[-1] = dc[-1] / bc[-1]
    for i in range(n - 2, -1, -1):
        x[i] = (dc[i] - cc[i] * x[i + 1]) / bc[i]
    return x


def thomas_solve(A, d):
    A = np.asarray(A, dtype=float)
    n = A.shape[0]
    # check tridiagonality (tolerance for floating point zeros)
    tol = 1e-12
    band = np.zeros_like(A)
    band += np.diag(np.diag(A))
    if n > 1:
        band += np.diag(np.diag(A, k=1), k=1)
        band += np.diag(np.diag(A, k=-1), k=-1)
    off_band = np.abs(A - band)
    if np.any(off_band > tol):
        # matrix is not tridiagonal -> do not attempt Thomas
        print("\nThomas: matrice non-tridiagonale — méthode non applicable.")
        return None
    b = np.diag(A).copy()
    a = np.diag(A, k=-1).copy()
    c = np.diag(A, k=1).copy()
    return thomas(a, b, c, np.asarray(d, dtype=float))


def thomas_lu(A):
    
    A = np.asarray(A, dtype=float)
    n = A.shape[0]
    # extract diagonals
    a = np.diag(A, k=-1).copy()  # sub-diagonal (length n-1)
    b = np.diag(A).copy()        # main diagonal (length n)
    c = np.diag(A, k=1).copy()   # super-diagonal (length n-1)
    # compute U and multipliers for L
    U = np.zeros((n, n), dtype=float)
    L = np.eye(n, dtype=float)
    U[0, 0] = b[0]
    if n > 1:
        U[0, 1] = c[0]
    for i in range(1, n):
        L[i, i - 1] = a[i - 1] / U[i - 1, i - 1]
        U[i, i] = b[i] - L[i, i - 1] * (c[i - 1] if i - 1 < len(c) else 0.0)
        if i < n - 1:
            U[i, i + 1] = c[i]
    return L, U
