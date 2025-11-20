import math
def newton(fonction, dfonction, x0, precision = 1e-10, max_iteration = 100):
    x=x0
    
    for i in range(max_iteration):
        fx= fonction(x)
        dfx = dfonction(x)
        
        if dfx ==0:
            print("la dérivée ne doit pas être nulle")
            return None
        
        x_nouveau = (x- fx)/ dfx
        
        if abs(x_nouveau -x)< precision:
            return x_nouveau
        
        x=x_nouveau
    
    print("maximum itérations atteint")
    
    return x



def fonction(x):
    return (x-1)/(x+1)

def dfonction(x):
    return (2/x**2 -1 )

solution = newton(fonction, dfonction, x0= 1)
solution = round(solution, 4)

print(f"racine trouvée est: {solution}")