

"""
    AMEGADJIN KOMLAN JOSUE - GENIE LOGICIEL S3 
    MTH 1300: Méthode de Gauss pour la résolution de systèmes linéaires.
    Programme de résolution d'un système linéaire avec la méthode de Gauss.
    École Polytechnique de Lomé - EPL
"""

# ---------------- Import des modules nécessaires
import numpy as np
import math
import matplotlib.pyplot as plt
from time import process_time

try:
    from .common import ResolutionSystTriSup, precision
except ImportError:
    from common import ResolutionSystTriSup, precision


# ---------------- Partie pivot de Gauss

def ReductionGauss(Aaug):
    # On créé une copie de Aaug
    A = np.copy(Aaug)
    # On récupère la taille de A
    n, m = np.shape(A)
    # On calcule les nouveaux éléments de la matrice A
    for j in range(1, n):
        for i in range(j, n):
            pivot = A[i, j-1] / A[j-1, j-1]
            A[i, :] = A[i, :] - pivot * A[j-1, :]
    # On renvoie la matrice A, matrice réduite de Aaug
    return(A)

def Gauss(A, B):
    # On créé Aaug en ajoutant B à A
    stack = np.column_stack([A, B])
    # On réduit la matrice
    reduc = ReductionGauss(stack)
    # On trouve la solution du système
    solution = ResolutionSystTriSup(reduc)
    # On renvoie la matrice solution du système
    return(solution)


# ---------------- Graphiques
def graph():
    # Mise en page pour mettre les 3 graphiques
    plt.gcf().subplots_adjust(wspace = 0.5, hspace = 0.5)
    plt.subplot(3, 1, 1)

    # Demande de la taille de la matrice maximale à calculer
    taille = int(input("Taille max de la matrice souhaitée ? \n"))

    # ------ Temps de calcul
    print("------------Gauss------------")
    time_list_gauss = []
    for i in range(0, taille + 1):
        A = np.random.rand(i, i)
        B = np.random.rand(i, 1)
        Gauss(A, B)
        t = process_time()
        time_list_gauss.append(t)

    # Calculs de temps & création de la liste les contenant tous
    T_list_gauss = []
    for i in range(len(time_list_gauss)):
        if i == len(time_list_gauss)-1:
            T = time_list_gauss[-1] - time_list_gauss[0]
        elif i == len(time_list_gauss):
            None
        else:
            T = time_list_gauss[i + 1] - time_list_gauss[i]
        T_list_gauss.append(T)

    # Affichage des temps en console
    if taille > 100:
        for i in range(0, len(T_list_gauss) -1, 100):
            print("Le temps de calcul pour une matrice de taille ", i, "est de :", T_list_gauss[i], "secondes.")

    print("Le temps de calcul pour une matrice de taille ", taille, "est de :", T_list_gauss[-2], "secondes.")

    print("\nLe temps de calcul total est de", T_list_gauss[-1], "secondes")
    minutes = int(T_list_gauss[-1]//60)
    secondes = int(T_list_gauss[-1] % 60)
    if minutes == 1:
        print("Soit environ", minutes, "minute et", secondes, "secondes.")
    elif minutes > 1:
        print("Soit environ", minutes, "minutes et", secondes, "secondes.")

    # On supprime le temps total afin de pouvoir afficher les temps de calcul
    del(T_list_gauss[- 1])

    abscisse = []
    for i in range(0, taille):
        abscisse.append(i)

    #  -- Création de la courbe
    plt.plot(abscisse, T_list_gauss, color = 'b', label = 'Méthode de Gauss')

    # -- Affichage de la courbe
    # Graphique 1 : Temps / taille
    plt.title("Temps de calcul en fonction de la taille de la matrice")
    plt.ylabel('Temps de calcul en secondes')
    plt.xlabel('Taille de la matrice')
    plt.grid(True)
    plt.legend(loc = 'best')
    # Graphique 2 : Temps en échelle logarithmique / taille
    plt.subplot(3, 1, 2)
    plt.plot(abscisse, T_list_gauss, color = 'b', label = 'Méthode de Gauss')
    plt.title("Temps de calcul en fonction de la taille de la matrice \n Echelle logarithmique ")
    plt.ylabel('Temps de calcul en secondes (log)')
    plt.xlabel('Taille de la matrice')
    plt.yscale('log')
    plt.grid(True)
    plt.legend(loc = 'best')

    # ------ Erreur en fonction de la taille
    plt.subplot(3, 1, 3)
    print("\n------------Gauss------------")
    size = []
    erreur = []
    for i in range(0, taille + 1):
        A = np.random.rand(i, i)
        B = np.random.rand(i, 1)
        X = Gauss(A, B)
        result = precision(A, X, B)
        size.append(i)
        erreur.append(result)
    # Graphique 3 : Erreur / taille
    #  -- Création de la courbe
    plt.plot(size, erreur, color = 'b', label = 'Méthode de Gauss')
    # -- Affichage de la courbe
    plt.xlabel('Taille de la matrice')
    plt.ylabel('Erreur ||A X - B||')
    plt.title('Erreur en fonction de la taille de la matrice')
    plt.grid(True)
    plt.legend(loc = 'best')
    plt.show()

graph()
