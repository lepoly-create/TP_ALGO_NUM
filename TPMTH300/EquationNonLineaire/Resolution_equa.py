# -*- coding:Utf-8 -*-
from sympy import symbols, sympify, lambdify

try:
    from TPALGONUM.TP.TPMTH300.EquationNonLineaire.MethodeResolution import methodeDichotomie, methodeSecante, methodeNewtonRaphson, methodePointFix
    from TPALGONUM.TP.TPMTH300.EquationNonLineaire.Accueil_menu import rectangle_bienvenu_complet, methodeDeResolution, saisirToleranceIterations, saisirFonction
except Exception:
    from MethodeResolution import methodeDichotomie, methodeSecante, methodeNewtonRaphson, methodePointFix
    from Accueil_menu import rectangle_bienvenu_complet, methodeDeResolution, saisirToleranceIterations, saisirFonction

x=  symbols('x')


def demanderBornes():
    while True:
        try:


            borneInf = float(input("Entrer la borne inf de l'intervalle: "))
            borneSup = float(input("Entrez la borne sup de l'intervalle: "))


            if borneInf < borneSup:
                return borneInf , borneSup
            else:
                print("Erreur: borne if doit être inf ")
        except ValueError:
            print("Erreur : veuillez entrez un nombre réel pour chaque borne")


def main():

    action = True #nous pourrons vérifier par celà si l'user veut commencer ou exit

    rectangle_bienvenu_complet()

    while action:
        methodeDeResolution()


        #assurons nous que l'user saisi un integer valide

        while True:
            try:
                choose = int(input("Entrez un number de la méthode (1 à 4); "))

                # validation correcte : 1 <= choose <= 4
                if 1 <= choose <= 4:
                    break
                else:
                    print("Le numéro doit être un entier compris entre 1 et 4")

            except ValueError:
                print("Erreur: veullez entrer un number integer. \n")


        #affichons la methode saisie

        if choose == 1:
            print("************METHODE DICHOTOMIE**************\n")
        elif choose ==2:
            print("************METHODE DE SECANTE**************\n")
        elif choose == 3:
            print("************METHODE NEWTON-RAPHSON**************\n")
        elif choose == 4:
            print("************METHODE POINT FIXE**************\n")


        #demandons à l'user de saisir la fonction à resoudre


        fonction_saisie = saisirFonction()

        try:

            fonction_expression = sympify(fonction_saisie)
            fonction = lambdify(x, fonction_expression, 'math') #permet de transformer sympy en fonction maths python

            print("Fonction interprétée : f(x) = ", fonction_expression)
        except Exception as e:
            print(f"Erreur lors de l'interpretation de la fonction : {e}")

            continue #demande à nouveau sinon

        # au cas où la méthode choisi ou opté exige les bornes

        if choose in [1 , 2]: #verifions si le choix est 1 ou 2
            borneInf , borneSup = demanderBornes()
            print(f"Votre intervalle est: [{borneInf} , {borneSup} ]")


        tolerance , nombreMax_iterations = saisirToleranceIterations()

        racine = None

        if choose == 1:
            # Exécution de la méthode de la dichotomie
            racine = methodeDichotomie(fonction, borneInf, borneSup, tolerance, nombreMax_iterations)
        elif choose == 2:
            # Exécution de la méthode de la sécante
            racine = methodeSecante(fonction, borneInf, borneSup, tolerance, nombreMax_iterations)
        elif choose  == 3:
            # Exécution de la méthode de Newton-Raphson
            derivee_fonction = lambdify(x, fonction_expression.diff(x), 'math')  # Calcul de la dérivée de la fonction
            x_initiale = float(input("Entrez une estimation initiale de la racine (exemple: 2) : Xo = "))
            racine = methodeNewtonRaphson(fonction, derivee_fonction, x_initiale, tolerance, nombreMax_iterations)
        elif choose == 4:
            # Exécution de la méthode du point fixe
            print("Saisissez la fonction phi(x) telle que x = phi(x) pour le point fixe.")
            phi_saisie = saisirFonction()
            try:
                phi_expression = sympify(phi_saisie)
                phi = lambdify(x, phi_expression, 'math')
            except Exception as e:
                print(f"Erreur lors de l'interprétation de phi(x): {e}")
                continue
            x_initiale = float(input("Entrez une estimation initiale de la racine (exemple: 2) : Xo = "))
            racine = methodePointFix(phi, x_initiale, tolerance, nombreMax_iterations)

        # Affichage de la racine trouvée
        if racine is not None:
            print(f"La racine trouvée est : {racine}")  # Afficher la racine si elle a été trouvée
        else:
            print("Aucune racine n'a été trouvée.")  # Si aucune racine n'a été trouvée

        # Demander à l'utilisateur s'il veut recommencer ou quitter
        choix_recommencer_Terminer = input("Souhaitez-vous recommencer ? (oui/non) : ").lower()
        if choix_recommencer_Terminer == "oui":
            action = True  # Continuer si l'utilisateur veut recommencer
        else:
            action = False  # Quitter si l'utilisateur ne veut pas recommencer

# Appel de la fonction principale pour démarrer le programme
if __name__ == "__main__":
    main()  # Lancer le programme



    





























