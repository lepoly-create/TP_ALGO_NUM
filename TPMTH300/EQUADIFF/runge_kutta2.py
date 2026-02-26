import numpy as np
import matplotlib.pyplot as plt
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
