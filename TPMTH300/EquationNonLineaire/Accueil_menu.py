from  sympy import symbols, sympify, sin, cos, tan, log, sqrt
import re


def rectangle_bienvenu_complet():

    phrase = " BIENVENUE"
    largeur = len(phrase) + 6 #Largeur du rectangle
    hauteur = 5 #Hauteur du rectangle

    # Dessinons la première ligne du rectangle
    print('*' * largeur)
    #Ligne avec bordures internes et des espaces à l'intérieur
    print('*' + ' ' * (largeur -2) + '*')

    #Ligne du milieu avec texte "BIENVENUE" centré
    print('* *' + phrase.center(largeur - 6) + '* *')
    print('*' + ' ' * (largeur - 2) + '*')

    print('*' * largeur)

#fonction permettant de choisir la methode de résolution parmi tant
def methodeDeResolution():
    print("****************Choisissez la méthode de résolution ********************")
    print("                 1-Méthode de Dichotomie")
    print("                 2-Méthode de la sécante")
    print("                 3-Méthode de Newton-ramphson")
    print("                 4-Méthode du Point fixe\n")

#initialisons 'x' pour pouvoir l'utiliser comme symbole à la suite
x = symbols('x')

def saisirFonction():
    print("****************Menus de choix de fonction ********************\n")
    print("Choisissez le type de fonction que vous souhaitez entrer: \n")
    print("                 1-Polynôme (exemple: x**2 -3*x + 8) ")
    print("                 2-Trigonométrique (exemple: sin(x) ; cos(x) - x )")
    print("                 3-Logarithmique (exemple: ln(x), log(x)")
    print("                 4-Racine carrée (exemple: sqrt(x) + x/3")
    print("                 5-Rationnelle (exemple: (x+1) / (x-2) )")
    print("                 6-Fonction générale\n")

    # Dictionnaire contenant les types de fonction disponibles avec des exemples
    types_fonction = {
        "1": "Polynôme (exemple: x**2 -3*x + 8)",
        "2": "Trigonométrique (exemple: sin(x) ; cos(x) - x )",
        "3": "Logarithmique (exemple: ln(x), log(x)",
        "4": "Racine carrée (exemple: sqrt(x) + x/3",
        "5": "Rationnelle (exemple: (x+1) / (x-2) )",
        "6": "Fonction générale (ex: x**3 - sin(x) + log(x))"
    }

    #Boucle principale pour permettre à l'utilisateur de choisir le type de fonction

    while True:

        type_fonction = input("Entrer le numéro du type de fonction : ").strip()

        #Vérifions si le choix est valide

        if type_fonction not in types_fonction:
            print("Choix invalide. Veuillez entrer un nombre entre 1 et 6. \n")
            continue

        # Supprimons les espaces pour que tout soit clean

        fonction = input("Entrez votre fonction mathématique: ").replace(" ", "").lower()
        
        # Valider en s'appuyant directement sur SymPy (plus robuste que le regex)
        try:
            fonction_sympy = sympify(fonction)
            print(f"Votre fonction saisie '{fonction}' est correcte. \n")
            return fonction_sympy
        except Exception:
            print("Erreur: Impossible de parser la fonction, vérifiez votre saisie (utilisez x, ** pour les puissances, log pour le ln, sqrt).")


#Fonction pour demander la tolérance et le nombre d'itérations

def saisirToleranceIterations():
    while True:
        try:
            #demandons la tolérance en tenant compte des notations maths

            tolerance_input = input("Entrez la tolérance souhaité (exemple: 0.0001, 1.0e-5, 1e-6 : \n")
            #convertissons l'entrée en float

            tolerance = float(tolerance_input)

            nombreMax_iterations = int(input("Entrez le nombre max d'itérations (exemple: 100) : \n"))

            if tolerance > 0 and nombreMax_iterations > 0:
                print(f"Votre tolérance est : {tolerance}")
                print(f"Votre nombre max d'itérations est : {nombreMax_iterations} \n")

                return tolerance , nombreMax_iterations
            else:
                print("Erreur: tolérance et nombre d'itérations doivent être > 0")
        except ValueError:
            print("Erreur: veuillez entrez des valeurs numériques valides")




























