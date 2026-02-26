from math import *

import numpy as np
import matplotlib.pyplot as plt


def recherche_pivot(A , i):
    p = i #p est l'indice du pivot
    for k in range(i +1 , len(A)):
        if abs(A[k] [i]) > abs(A[p] [i]): #si on trouve un meilleur pivot
            p = k

    return p

def echange_ligne(A , i , j):
    for k in range (len(A[0])):
        A[i][k] , A[j][k] = A[j][k] ,A[i][k]

def transvection(A, i , j ,mu):
    for k in range(len(A[0])):
        A[i][k] += mu * A[j][k]


def dilatation(A , i , mu):
    for k in range(len(A[0])):
        A[i][k] *= mu

def élimination(A ,Y):
    n = len(A)
    for j in range(n):
        p = recherche_pivot(A, j)
        echange_ligne(A ,j, p)
        echange_ligne(Y, j, p)
        dilatation(Y ,j, 1 / A[j][j])
                # Y est modifié avant que A[j][i] ne soit modifié
        dilatation(A ,j , 1 /A[j][j])
        for i in range(n):
            if i!= j:
                transvection(Y ,j , i, -A[j][i])
                transvection(A ,j , i, -A[j][i])
