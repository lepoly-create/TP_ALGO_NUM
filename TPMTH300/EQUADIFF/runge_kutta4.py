import numpy as np
import matplotlib.pyplot as plt
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