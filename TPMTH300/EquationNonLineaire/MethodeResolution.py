import math
import sys

# Fonction pour la methode de dichotomie

def methodeDichotomie(fonction, borneInf, borneSup, tolerance, nombreMax_iterations):

    #evaluons la fonction à la borneInf

    valueBorneInf = fonction(borneInf)

    #si f(a) est proche de 0 , a est racine

    if (math.fabs(valueBorneInf) <= tolerance):
        return borneInf

    valueBorneSup = fonction(borneSup)

    if (math.fabs(valueBorneSup) <= tolerance):
        return borneSup

    #verifions si la fonction a une racine entre a et b

    if (valueBorneInf * valueBorneSup > 0.0):
        print("La racine n'est pas dans l'intervalle [", borneInf, ";" , borneSup , "]")
        sys.exit(0)


    #calcul du nombre max d'itérations

    nombreIterattions = int(math.ceil(math.log(math.fabs(borneSup - borneInf) / tolerance) / math.log(2.0)))

    #on éffectue la dicho
    for i in range(min(nombreIterattions + 1, nombreMax_iterations)):
        pointCenter = (borneInf + borneSup) * 0.5

        valuePointCenter = fonction(pointCenter)

        #affichons intervalle actuel

        print(f"Iteration {i + 1}: Intervalle actuel = [{borneInf} , {borneSup}]")

        #si f(c) est proche de 0
        if (valuePointCenter == 0.0 or (borneSup - borneInf) < tolerance ):
            print(f"Solution trouvée: x = {pointCenter} après {i + 1} itérations")
            return pointCenter
        #mettons à jour les bornes a ou b selon le signe de la fonction f(c)

        if (valuePointCenter * valueBorneSup < 0.0):
            borneInf = pointCenter
            valueBorneInf = valuePointCenter
        else:
            borneSup = pointCenter
            valueBorneSup = valuePointCenter
    #retournons la moyenne des bornes a et b comme approximations de la rzcine concernée
    return (borneInf + borneSup) * 0.5


#Fonction de la méthode secante

def methodeSecante(fonction, borneInf, borneSup, tolerance, nombreMax_iterations):

    #Evaluons la fonction à la borneInf
    valueBorneInf = fonction(borneInf)

    if (math.fabs(valueBorneInf) <= tolerance):
        return borneInf


    valueBorneSup = fonction(borneSup)

    if (math.fabs(valueBorneSup) <= tolerance):
        return borneSup

    #verifions que la re=acine est encardrée

    if (valueBorneInf * valueBorneSup > 0.0):
        print("La racine n'est pas encédré dans l'intervalle [", borneInf, ";", borneSup, "]")
        sys.exit(0)
    compteurIterations = 0

    # appliquons la methode de secante jusqu'à ce que la tolerance soit atteinte
    while (math.fabs(borneSup - borneInf) > tolerance and math.fabs(valueBorneSup) > tolerance and (compteurIterations < nombreMax_iterations)):
        compteurIterations += 1
        estimation = borneInf - valueBorneInf * (borneSup - borneInf) / (valueBorneSup - valueBorneInf)
        valueEstimation = fonction(estimation)

        #si f(estimation) est proche de 0

        if (math.fabs(valueEstimation) <= tolerance):
            return estimation

        #Mise à jour des bornes en fonction du signe de f(estimation)

        if (valueEstimation * valueBorneSup < 0.0):
            borneInf = estimation
            valueBorneInf = valueEstimation
        else:
            borneSup = estimation
            valueBorneSup = valueEstimation
    #retournons l'estimation à partir de la methode de la secante
    return (borneInf - valueBorneInf * (borneSup - borneInf) / (valueBorneSup - valueBorneInf))


#Fonction pour la methode de newton

def methodeNewtonRaphson(fonction, deriveeFonction, xInit, tolerance, nombreMax_iterations):
    compteurIterations = 0
    x = xInit
    valueX = fonction(x)

    # la méthode itère jusqu'à converger
    while ((math.fabs(valueX) > tolerance) and (compteurIterations < nombreMax_iterations)):
        valu_deriveX = deriveeFonction(x)

        # si dérivée est nulle
        if (valu_deriveX == 0):
            print("La dérivée est nulle. Donc la méthode de Newton échoue ")
            return None

        x = x - valueX / valu_deriveX
        valueX = fonction(x)
        compteurIterations += 1

    if (math.fabs(valueX) <= tolerance):
        return x
    else:
        print("Pas de convergence avec la méthode de Newton.")
        return None



#fonction of methode with point fix

def methodePointFix(phi, xInit, tolerance, nombreMax_iterations):
    compteurIterations = 0
    x = xInit

    # on itère jusqu'à ce que la différence entre x et phi(x) soit < tolérance
    while (compteurIterations < nombreMax_iterations):
        x_next = phi(x)
        if (math.fabs(x_next - x) <= tolerance):
            return x_next
        x = x_next
        compteurIterations += 1

    print("Pas de convergence avec la méthode du point fixe.\n")
    return None


# améliorons ensemble la methode de balayage pour newton

def balayageNewton(fonction, deriveeFonction, borneInf, borneSup, tolerance, nombreMax_iterations):

    valueInits = [borneInf , (borneInf + borneSup) / 2, borneSup]

    for valInit in valueInits:
        print(f"Essai de methodeNewton avec xInit = {valInit}")
        solution  = methodeNewtonRaphson(fonction, deriveeFonction, valInit, tolerance, nombreMax_iterations)

        if solution is not None:
            return solution
    return None































































