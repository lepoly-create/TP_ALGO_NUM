import numpy as np
import matplotlib.pyplot as plt
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