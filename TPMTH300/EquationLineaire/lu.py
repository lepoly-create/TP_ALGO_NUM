#!/usr/bin/python
# -*- coding:utf-8 -*-

"""
    AMEGADJIN Komlan Josue - GENIE LOGICIEL S3
    MTH 1300: Méthode de Gauss pour la résolution de systèmes linéaires.
    Programme de résolution d'un système linéaire avec la méthode de la décomposition LU.
    École Polytechnique de Lomé - EPL
"""

# ---------------- Import des modules nécessaires
import numpy as np
import math
import matplotlib.pyplot as plt

try:
    from .common import precision
except ImportError:
    from common import precision


# ---------------- Partie décomposition LU
def DecompositionLU(A):
    # On récupère la taille de A
    n, m = np.shape(A)
    # On créé une copie de A
    U = np.copy(A)
    # On créé une matrice identité de taille n
    L = np.eye(n)
    if m != n:
        print("La matrice n'est pas carrée.")
        return(None)
    # On calcule les éléments de L et U
    for j in range(n):
        for i in range(j + 1, n):
            pivot = U[i, j] / U[j, j]
            U[i, :] = U[i, :] - pivot * U[j, :]
            L[i, j] = pivot
    # On renvoie les matrices L et U
    return(L, U)

def ResolutionLU(L, U, B):
    Y = []
    # On récupère la taille de B
    n, m = np.shape(B)
    # On résoud le sytème grâce à la décomposition L U
    for i in range(n):
        Y.append(B[i])
        for j in range(i):
            Y[i] = Y[i] - (L[i, j] * Y[j])
        Y[i] = Y[i] / L[i, i]
    X = np.zeros(n)
    for i in range(n, 0, - 1):
        X[i - 1] = (Y[i - 1] - np.dot(U[i - 1, i:], X[i:])) / U[i - 1, i - 1]
    # On renvoie la matrice solution X
    return(X)

def LU(A, B):
    # On décompose la matrice A en deux matrices L et U
    L, U = DecompositionLU(A)
    # On trouve la solution du système
    solution = ResolutionLU(L, U, B)
    # On renvoie la matrice solution du système
    return(solution)

