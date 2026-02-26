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

def fang_sheng(A ,Y):
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

## NUMPY

a = np.zeros((3 , 3))
b =np.array([
            [1 , 2, 3] ,
            [4, 5, 6],
            [7, 8, 9]
    ] , dtype= int)

a[0 , 0] = 20
a[2][1] = -1
print(a)
print(a.shape , a.dtype , a.size, b.dtype, sep=", ")
# la methode shape donne la taille(forme) du tableau alors que celle de size fournit le nbre total de cellules élémentaires

a = np.array([4, 7, 9])
b= np.array([[1, 2, 3] , [4, 5, 6] , [7, 8, 9]], dtype=int)

print(a.size , a.shape)
print(b.size , b.shape)


# on peut remplir un tableau de différentes façons
def f(i , j):
    return 5/ (i + j + 1)

H1 = np.array([[f(i , j) for j in range(5)] for i in range(5)])
H2 = np.fromfunction(f, (5, 5))

# échangeons le type de cellules

Hentier = H2.copy()
Hentier = Hentier.astype(int)

print(H2.dtype)
print(H2)
print(Hentier)


























