"""
Created on Sat Mar 20 15:25:34 2021
@author: grupo 7
"""

import numpy as np


def binf2dec(binf):
    numBin = np.array(binf)
    s = numBin[0]           #Obtenemos el signo del numero
    e = numBin[1:6]         #Obtenemos el exponente del numero
    m = numBin[6:16]        #Obtenemos la mantiza del numero
    
    exp = np.sum(e)         #Obtengo la cantidad de 1's del exponente
    man = np.sum(m)         #Obtengo la cantidad de 1's de la mantiza
    
    exp_all_ones = 5        #Si el exp son todos 1's tendriamos 5x1
    man_all_ones = 10       #Si el man son todos 1's tendriamos 10x1

    if((exp == 0) and (man == 0)):
        print("El numero ingresado es nulo")
        return 0
    elif((exp == exp_all_ones) and (man == 0)):
        if(s == 1):
            print("El numero ingresado es +inf")
            return np.Inf
        else:
            print("El numero ingresado es -inf")
            return np.NINF
    elif((exp == exp_all_ones) and (man != 0)):
        print("El numero ingresado no es un numero")
        return np.NAN
    elif((exp == 0) and (man != 0)):
        print("El numero ingresado es subnormal ")
        sesgo = 2**(exp_all_ones - 1) - 1
        mantiza = 0
        for i in np.arange(man_all_ones):
            mantiza = mantiza + m[i]*(2.0**-(i+1))
        val = ((-1)**s)*mantiza*(2.0**(1-sesgo))
        print("Su valor es:")
        print("{:.10f}".format(val))
        return val
    elif((exp > 0) and (exp < exp_all_ones)):
        print("El numero ingresado es normal ")
        sesgo = 2**(exp_all_ones - 1) - 1
        mantiza = 0
        exponte = 0
        for i in np.arange(exp_all_ones):
            exponte = exponte + e[i]*(2**(exp_all_ones - i - 1))
        for i in np.arange(man_all_ones):
            mantiza = mantiza + m[i]*(2.0**-(i+1))
        mantiza = mantiza + 1
        val = ((-1)**s)*mantiza*(2.0**(exponte-sesgo))
        print("Su valor es:")
        print("{:.10f}".format(val))
        return val
    
    
def test():
    a = np.array([0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0])
    b = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0])
    c = np.array([1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    d = np.array([0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1])
    x = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0])

    binf2dec(a)
    binf2dec(b)
    binf2dec(c)
    binf2dec(d)
    binf2dec(x)

