import math

def dicho(fonction, a,b, precision=1e-10,max_iteration=50 ):
    xm = (a + b)/2
    if(abs(b-a)/ 2*abs(xm)) < precision:
        print("La convergence est atteinte ")
        print(f"la racine atteinte est : {xm}")
        fonction(xm)
        
        print(f"x1 : {a} , f(x1) : {fonction(a)} ")
        print(f"x2 : {b} , f(x2) : {fonction(b)} ")
        print(f"xm : {xm} , f(xm) : {fonction(xm)} ")
        print("\n")
        
    if fonction(a) * fonction(xm) <0 :
        b = xm
    elif fonction(b) * fonction(xm) <0 :
        a = xm
    iteration = 0
    if(iteration >= max_iteration):
        print(f"Convergence non atteinte à  {max_iteration} itérations")
        
        iteration +=1
    else:   
        return dicho(fonction, a,b, precision, max_iteration -1)

def fonction(x):
    return math.exp(-x) -x

resultat =dicho(fonction, 0,1)

print(f"{resultat}")
        
        
    