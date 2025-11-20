import math
        
    
def dichotomie(fonction, a,b, precision = 1e-10, max_iteration=100):
    if fonction(a)*fonction(b) >0:
        print("les images des deux points doivent être de signes contraires, pas de solution dans ")
        return None
    for iteration in range(max_iteration):
        m=(a+b)/2
        
        if abs(fonction(m) )< precision or (b-a)/2 < precision:
            return m
        
        if fonction(a) * fonction(m)<0:
            b=m
        else:
            a=m
            
    print("Nbre max d'iteration atteinte ")
    
    return m

def fonction(x):
    return (x-1)/(x+1)

racine = dichotomie(fonction, -2,2)
print(f" la racine pour cette fonction est : {racine}")
        
        
        
        
        
        
    

