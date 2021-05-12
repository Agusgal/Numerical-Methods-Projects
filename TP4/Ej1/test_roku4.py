import numpy as np
import matplotlib.pyplot as plt
from Ej1.ej1 import roku4

# Defino ecuación conocida (circuito con capacitor y resistenci)
def dx(time, x):
    R = 1e3
    C = 1e-6
    W = 2.0 * np.pi * 1000
    A = 1.0
    T = 5*2*np.pi/W

    return ((A*np.cos(W*time)-x)/(R*C))

# Defino solucion conocida
def known(time):
    R = 1e3
    C = 1e-6
    W = 2.0 * np.pi * 1000
    A = 1.0

    x = -np.exp(-time/(R*C))+np.cos(W*time)+W*R*C*np.sin(W*time)
    x = (A/(1+(W*R*C)**2))*x
    return x

#Funcion de testeo
def testRoku4():

    x0 = np.zeros(1)
    W = 2.0 * np.pi * 1000
    T = 5 * 2 * np.pi / W
    h = T/10000

    t, xe = roku4(dx, 0, T, h, x0)
    xknown = known(t)

    fig, ax = plt.subplots()
    ax.plot(t, xknown, label='Conocida')
    ax.plot(t, xe[:, 0], label='Roku4')

    ax.legend()
    plt.show()


    fig, ax = plt.subplots()
    ax.plot(t, xe[:, 0] - xknown, label='Error')
    ax.legend()
    plt.show()

    # Falta agregar el tema del error con delta/2, delta

if __name__ == '__main__':
    testRoku4()
