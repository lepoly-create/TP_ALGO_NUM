import numpy as np
import matplotlib.pyplot as plt
from EquationLineaire.gauss_sans_pivot import gauss_sans_pivot
from EquationLineaire.gauss_pivot_partiel import gauss_pivot_partiel
from EquationLineaire.gauss_pivot_total import gauss_pivot_total
from EquationLineaire.gauss_jordan import gauss_jordan
from EquationLineaire.doolithe import doolittle_lu
from EquationLineaire.crout import crout_lu
from EquationLineaire.cholesky import cholesky
from ITERACTIVES.jacobi import jacobi
from ITERACTIVES.gauss_seidel import gauss_seidel
from INTERPOLATION.lagrange import lagrange_interpolation
from INTERPOLATION.newton_interpolation import newton_interpolation
from INTERPOLATION.moindres_carres import least_squares_poly
from EQUADIFF.euler_explicite import euler_explicite
from EQUADIFF.runge_kutta2 import runge_kutta_2
from EQUADIFF.runge_kutta4 import runge_kutta_4


# =============================================================================
# Exemples d'utilisation
# =============================================================================

if __name__ == "__main__":
    print("=== Exemples d'utilisation des méthodes ===")
    print("\n1. Résolution de systèmes linéaires")
    A = np.array([
        
        [3, 2, -1], 
        [2, -2, 4], 
        [-1, 0.5, -1]]
    )
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
    A_pos = np.array([
        [4, 2, 1], 
        [2, 5, 3], 
        [1, 3, 6]
    ])
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