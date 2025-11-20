import math
def pointfixe(fonction, x0, precision= 1e-10, max_iteration = 100):
    x=x0
    for i in range(max_iteration):
        x_nouveau =fonction(x)
        
        if abs(x_nouveau - x) < precision :
            return x_nouveau
        
        if abs(x_nouveau)> 1e10:
            return None
        
        x = x_nouveau
    return "pas de convergence"


def fonction(x):
    return x**2 -1

resultat = pointfixe(fonction, 1)
print(resultat)