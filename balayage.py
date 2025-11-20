import math
#import function


def balayage(fonction, a,b, h=0.001):
    x=a
    solutions = []
    while x<b:
        x_suivant = x + h
        
        if fonction(x) * fonction(x_suivant) <0:
            solutions.append((x, x_suivant))
            
        x= x_suivant
    return solutions

def fonction(x):

    return (x-1)/(x+1)
interval = balayage(fonction, -2,  2,h=0.1)

print(interval)