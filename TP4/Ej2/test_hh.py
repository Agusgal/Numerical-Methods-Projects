from Ej1.ej1 import ruku4
from Ej2.ej2 import hhDer
import numpy as np
import matplotlib.pyplot as plt

def test_hh():
    x0 = np.zeros(2)
    x0[0] = 2
    x0[1] = 3

    T = 600
    h = T/10000

    t1, xout1 = ruku4(hhDer, 0, T, h, x0)
    t2, xout2 = ruku4(hhDer, 0, T, h/2, x0)


    estimacionS = abs(xout1[-1, 0] - xout2[-1, 0])/15
    estimacionP = abs(xout1[-1, 1] - xout2[-1, 1])/15

    print(estimacionS)
    print(estimacionP)


    fig, ax = plt.subplots()
    ax.plot(t1, xout1[:, 0], label='ds')
    ax.plot(t1, xout1[:, 1], label='dp')

    ax.legend()
    plt.show()




if __name__ == '__main__':
    test_hh()