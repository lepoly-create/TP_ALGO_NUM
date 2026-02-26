import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Méthodes directes de résolution de systèmes linéaires
# =============================================================================

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


def crout_lu(A):
    """
    Factorisation LU de Crout : L triangulaire inférieure (diagonale quelconque),
    U triangulaire supérieure avec diagonale unité.
    Retourne (L, U).
    """
    A = A.astype(float, copy=True)
    n = A.shape[0]
    L = np.zeros_like(A)
    U = np.eye(n)
    for k in range(n):
        # Calcul de L[k:, k]
        for i in range(k, n):
            L[i, k] = A[i, k] - np.dot(L[i, :k], U[:k, k])
        # Calcul de U[k, k+1:]
        for j in range(k + 1, n):
            U[k, j] = (A[k, j] - np.dot(L[k, :k], U[:k, j])) / L[k, k]
    return L, U


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


# =============================================================================
# Méthodes itératives
# =============================================================================

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


# =============================================================================
# Interpolation polynomiale
# =============================================================================

def lagrange_interpolation(x_data, y_data, x_eval):
    """
    Évalue le polynôme d'interpolation de Lagrange aux points x_eval.
    """
    x_data = np.asarray(x_data)
    y_data = np.asarray(y_data)
    n = len(x_data)
    result = np.zeros_like(x_eval, dtype=float)
    for i in range(n):
        # construction du i‑ème polynôme de base
        Li = np.ones_like(x_eval)
        for j in range(n):
            if j != i:
                Li *= (x_eval - x_data[j]) / (x_data[i] - x_data[j])
        result += y_data[i] * Li
    return result


def newton_divided_differences(x_data, y_data):
    """
    Calcule les différences divisées de Newton.
    Retourne un tableau des coefficients (les différences divisées diagonales).
    """
    n = len(x_data)
    coef = y_data.copy().astype(float)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / (x_data[i] - x_data[i - j])
    return coef


def newton_interpolation(x_data, y_data, x_eval):
    """
    Évalue le polynôme d'interpolation de Newton (forme des différences divisées).
    """
    coef = newton_divided_differences(x_data, y_data)
    n = len(x_data)
    # Évaluation par l'algorithme de Horner
    result = coef[-1] * np.ones_like(x_eval)
    for i in range(n - 2, -1, -1):
        result = result * (x_eval - x_data[i]) + coef[i]
    return result


# =============================================================================
# Approximation au sens des moindres carrés (polynomiale)
# =============================================================================

def least_squares_poly(x_data, y_data, deg):
    """
    Ajuste un polynôme de degré `deg` aux données par moindres carrés.
    Retourne les coefficients du polynôme (ordre croissant) et la matrice de Vandermonde.
    """
    x_data = np.asarray(x_data)
    y_data = np.asarray(y_data)
    # Construction de la matrice de Vandermonde
    V = np.vander(x_data, deg + 1, increasing=True)
    # Résolution des équations normales : V.T @ V @ c = V.T @ y
    # Utilisation de numpy.linalg.lstsq pour la stabilité
    coeffs, residuals, rank, s = np.linalg.lstsq(V, y_data, rcond=None)
    return coeffs, V


# =============================================================================
# Résolution numérique d'équations différentielles
# =============================================================================

def euler_explicite(f, y0, t_span, h):
    """
    Méthode d'Euler explicite.
    f : fonction f(t, y) du problème y' = f(t, y)
    y0 : condition initiale
    t_span : tuple (t0, tfinal)
    h : pas de temps
    Retourne les vecteurs t et y.
    """
    t0, tf = t_span
    n_steps = int(np.ceil((tf - t0) / h))
    t = np.linspace(t0, t0 + n_steps * h, n_steps + 1)
    y = np.zeros((n_steps + 1,) + np.shape(y0))
    y[0] = y0
    for i in range(n_steps):
        y[i + 1] = y[i] + h * f(t[i], y[i])
    return t, y


def runge_kutta_2(f, y0, t_span, h, method='midpoint'):
    """
    Méthode de Runge‑Kutta d'ordre 2.
    method : 'midpoint' (point milieu) ou 'heun' (Euler modifié)
    """
    t0, tf = t_span
    n_steps = int(np.ceil((tf - t0) / h))
    t = np.linspace(t0, t0 + n_steps * h, n_steps + 1)
    y = np.zeros((n_steps + 1,) + np.shape(y0))
    y[0] = y0
    for i in range(n_steps):
        k1 = f(t[i], y[i])
        if method == 'midpoint':
            k2 = f(t[i] + h / 2, y[i] + h / 2 * k1)
            y[i + 1] = y[i] + h * k2
        elif method == 'heun':
            k2 = f(t[i] + h, y[i] + h * k1)
            y[i + 1] = y[i] + h / 2 * (k1 + k2)
        else:
            raise ValueError("Méthode inconnue. Choisir 'midpoint' ou 'heun'.")
    return t, y


def runge_kutta_4(f, y0, t_span, h):
    """
    Méthode de Runge‑Kutta classique d'ordre 4.
    """
    t0, tf = t_span
    n_steps = int(np.ceil((tf - t0) / h))
    t = np.linspace(t0, t0 + n_steps * h, n_steps + 1)
    y = np.zeros((n_steps + 1,) + np.shape(y0))
    y[0] = y0
    for i in range(n_steps):
        k1 = f(t[i], y[i])
        k2 = f(t[i] + h / 2, y[i] + h / 2 * k1)
        k3 = f(t[i] + h / 2, y[i] + h / 2 * k2)
        k4 = f(t[i] + h, y[i] + h * k3)
        y[i + 1] = y[i] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return t, y


# =============================================================================
# Exemples d'utilisation
# =============================================================================

if __name__ == "__main__":
    print("=== Exemples d'utilisation des méthodes ===")
    print("\n1. Résolution de systèmes linéaires")
    A = np.array([[3, 2, -1], [2, -2, 4], [-1, 0.5, -1]])
    b = np.array([1, -2, 0])
    print("Système A x = b")
    print("A =\n", A)
    print("b =", b)

    x_gauss = gauss_sans_pivot(A, b)
    print("\nGauss sans pivot :", x_gauss)

    x_pp = gauss_pivot_partiel(A, b)
    print("Gauss pivot partiel :", x_pp)

    x_pt = gauss_pivot_total(A, b)
    print("Gauss pivot total :", x_pt)

    x_jordan = gauss_jordan(A, b)
    print("Gauss‑Jordan :", x_jordan)

    L, U = doolittle_lu(A)
    print("\nFactorisation Doolittle :")
    print("L =\n", L)
    print("U =\n", U)
    # Vérification : A ≈ L @ U
    print("L @ U =\n", L @ U)

    Lc, Uc = crout_lu(A)
    print("\nFactorisation Crout :")
    print("L =\n", Lc)
    print("U =\n", Uc)
    print("Lc @ Uc =\n", Lc @ Uc)

    # Pour Cholesky, il faut une matrice symétrique définie positive
    A_pos = np.array([[4, 2, 1], [2, 5, 3], [1, 3, 6]])
    L_chol = cholesky(A_pos)
    print("\nCholesky :")
    print("L =\n", L_chol)
    print("L @ L.T =\n", L_chol @ L_chol.T)

    print("\n2. Méthodes itératives")
    A_diag_dom = np.array([[4, 1, 0], [1, 4, 1], [0, 1, 4]])
    b_it = np.array([1, 2, 3])
    x0 = np.zeros(3)
    x_jac, it_jac, err_jac = jacobi(A_diag_dom, b_it, x0, tol=1e-6, max_iter=100)
    print("Jacobi :", x_jac, "en", it_jac, "itérations")
    x_gs, it_gs, err_gs = gauss_seidel(A_diag_dom, b_it, x0, tol=1e-6, max_iter=100)
    print("Gauss‑Seidel :", x_gs, "en", it_gs, "itérations")

    print("\n3. Interpolation")
    x_pts = np.array([0, 1, 2, 3, 4])
    y_pts = np.array([1, 2, 0, 5, 3])
    x_eval = np.linspace(0, 4, 100)
    y_lag = lagrange_interpolation(x_pts, y_pts, x_eval)
    y_new = newton_interpolation(x_pts, y_pts, x_eval)

    plt.figure()
    plt.plot(x_pts, y_pts, 'o', label='Données')
    plt.plot(x_eval, y_lag, '-', label='Lagrange')
    plt.plot(x_eval, y_new, '--', label='Newton')
    plt.legend()
    plt.title("Interpolation polynomiale")
    plt.show()

    print("\n4. Moindres carrés")
    x_noisy = np.linspace(0, 5, 20)
    y_true = 2 + 3 * x_noisy - 0.5 * x_noisy ** 2
    y_meas = y_true + 0.5 * np.random.randn(len(x_noisy))
    coeffs, V = least_squares_poly(x_noisy, y_meas, deg=2)
    print("Coefficients (ordre croissant) :", coeffs)
    y_fit = V @ coeffs
    plt.figure()
    plt.plot(x_noisy, y_meas, 'o', label='Mesures')
    plt.plot(x_noisy, y_fit, '-', label='Ajustement quadratique')
    plt.legend()
    plt.title("Moindres carrés")
    plt.show()

    print("\n5. Équations différentielles")
    def f_ode(t, y):
        return -2 * y + t + 1
    y0 = 1.0
    t_span = (0, 5)
    h = 0.1
    t, y_euler = euler_explicite(f_ode, y0, t_span, h)
    t, y_rk2 = runge_kutta_2(f_ode, y0, t_span, h, method='midpoint')
    t, y_rk4 = runge_kutta_4(f_ode, y0, t_span, h)

    plt.figure()
    plt.plot(t, y_euler, ':', label='Euler')
    plt.plot(t, y_rk2, '--', label='RK2')
    plt.plot(t, y_rk4, '-', label='RK4')
    plt.legend()
    plt.title("Résolution d'équation différentielle")
    plt.show()