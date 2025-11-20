import math

def secante(fonction, x0, x1, precision = 1e-10, max_iteration = 100):
    
    for i in range(max_iteration):
        f0 = fonction(x0)
        f1 = fonction(x1)
        
        if f1-f0 ==0:
            print("Division par 0 n'est pas valide : indéterminée")
            
            return None
        
        x2 = x1 -f1 * (x1 - x0) / (f1 - f0)
        
        if abs(x2 - x1) < precision:
            return x2
        
        x0,x1 = x1,  x2
        
    print("Itération atteitn")
        
    return x2
    
def fonction(x):
    return (x-1)/(x+1)
    
solution = secante(fonction, -2,2)

print(f"soluton:{solution}")