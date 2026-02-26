import numpy as np
import matplotlib.pyplot as plt
def euler_explicite(f, y0, t_span, h):

    t0, tf = t_span
    n_steps = int(np.ceil((tf - t0) / h)) # arrondi vers le haut pour couvrir tout l'intervalle
    t = np.linspace(t0, t0 + n_steps * h, n_steps + 1)
    y = np.zeros((n_steps + 1,) + np.shape(y0))
    y[0] = y0
    for i in range(n_steps):
        y[i + 1] = y[i] + h * f(t[i], y[i])
    return t, y