import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def ruku4(f, t0, tf, deltaT, x0):

    # cantidad de puntos
    N = int((tf - t0)/deltaT)

    tk = np.linspace(t0, tf, N + 1)

    n = x0.shape[0]
    xk = np.zeros((N + 1, n))  # matriz del problema
    xk[0, :] = x0

    # Se itera sobre cada instante de tiempo
    for i in range(N):
        f1 = f(tk[i], xk[i, :])
        f2 = f(tk[i] + deltaT/2, xk[i, :] + ((deltaT/2) * f1))
        f3 = f(tk[i] + deltaT/2, xk[i, :] + ((deltaT/2) * f2))
        f4 = f(tk[i] + deltaT, xk[i, :] + deltaT * f3)

        xk[i + 1, :] = xk[i, :] + ((f1 + 2*f2 + 2*f3 + f4)/6) * deltaT

    return tk, xk


def other(time, x):
    u = x[0]
    v = x[1]

    du = v
    dv = -np.tanh(np.cos(time) * u * v) - u

    dx = np.zeros_like(x)
    dx[0] = du
    dx[1] = dv

    return dx


def eq(time, x):
    mu = 2

    u = x[0]
    v = x[1]

    du = v
    dv = 1 - u - mu*(u**4 - 1)*v

    dx = np.zeros_like(x)
    dx[0] = du
    dx[1] = dv

    return dx


def higginsselkov(time, x):
    v0 = 0.5

    s = x[0]
    p = x[1]

    ds = v0 - 0.23 * s * p ** 2
    dp = 0.23 * s * p ** 2 - 0.4 * p

    dx = np.zeros_like(x)
    dx[0] = ds
    dx[1] = dp

    return dx


def test():
    #testRoku4()
    #test_hs()
    #testEq()
    testOther()

def testOther():
    x0 = np.zeros(2)
    x0[0] = 1
    x0[1] = 1

    T = 50
    h = T / 10000

    t1, xout1 = ruku4(other, 0, T, h, x0)

    fig, ax = plt.subplots()
    ax.plot(t1, xout1[:, 0], label='u')
    ax.plot(t1, xout1[:, 1], label='v')
    ax.legend()
    plt.show()


def testEq():
    x0 = np.zeros(2)
    x0[0] = 1.5
    x0[1] = 1

    T = 75
    h = T / 10000

    t1, xout1 = ruku4(eq, 0, T, h, x0)

    fig, ax = plt.subplots()
    ax.plot(t1, xout1[:, 0], label='u')
    ax.plot(t1, xout1[:, 1], label='v')
    ax.legend()
    plt.show()

def test_hs():
    print('Comienza testeo algoritmo higginsselkov:\n')

    x0 = np.zeros(2)
    x0[0] = 2
    x0[1] = 3

    T = 600
    h = T / 10000

    t1, xout1 = ruku4(higginsselkov, 0, T, h, x0)
    t2, xout2 = ruku4(higginsselkov, 0, T, h / 2, x0)

    estimacionS = abs(xout1[-1, 0] - xout2[-1, 0]) / 15
    estimacionP = abs(xout1[-1, 1] - xout2[-1, 1]) / 15

    print('Estimación para S: ', estimacionS)
    print('Estimación para P: ', estimacionP)

    fig, ax = plt.subplots()
    ax.plot(t1, xout1[:, 0], label='s con v0 = 0.5')
    ax.plot(t1, xout1[:, 1], label='p con v0 = 0.5')

    ax.legend()
    plt.show()

    print('Aplicando algoritmo de scipy, por favor espere...\n')

    ##Comparacion con algoritmo de scipy
    s45 = solve_ivp(higginsselkov, [0, T], x0, method='RK45', t_eval=t1, rtol=1e-13, atol=1e-14)

    print('Listo!')

    fig, ax = plt.subplots()
    ax.plot(s45.t, s45.y[0], label='S en scipy')
    ax.plot(s45.t, s45.y[1], label='P en scipy')
    ax.legend()
    plt.show()

    fig, ax = plt.subplots()
    ax.plot(t1, xout1[:, 0], label='S calculada')
    ax.plot(s45.t, s45.y[0], label='S Scipy')
    ax.legend()
    plt.show()

    fig, ax = plt.subplots()
    ax.plot(t1, xout1[:, 1], label='P calculada')
    ax.plot(s45.t, s45.y[1], label='P Scipy')
    ax.legend()
    plt.show()

    ##Estos errores son respecto a lo que dice scipy
    lastknownS = s45.y[0, -1]
    lastknownP = s45.y[1, -1]
    errorS = abs(xout2[-1, 0] - lastknownS)
    errorP = abs(xout2[-1, 1] - lastknownP)

    porcentualS = (abs(errorS - estimacionS) / errorS * 100)
    porcentualP = (abs(errorP - estimacionP) / errorP * 100)

    print('Error porcentual para S (comparado con resultado de scipy): ', porcentualS)
    print('Error porcentual para P (comparado con resultado de scipy): ', porcentualP)


# Defino ecuación conocida (circuito con capacitor y resistencia)
def dx(time, x):
    R = 1e3
    C = 1e-6
    W = 2.0 * np.pi * 1000
    A = 1.0
    T = 5 * 2 * np.pi / W

    return (A*np.cos(W*time)-x)/(R*C)


# Defino solucion conocida
def known(time):
    R = 1e3
    C = 1e-6
    W = 2.0 * np.pi * 1000
    A = 1.0

    x = -np.exp(-time/(R*C))+np.cos(W*time)+W*R*C*np.sin(W*time)
    x = (A/(1+(W*R*C)**2))*x
    return x


# Funcion de testeo
def testRoku4():

    x0 = np.zeros(1)
    W = 2.0 * np.pi * 1000
    T = 5 * 2 * np.pi / W ##5 ciclos
    h = T/5750

    t, xe = ruku4(dx, 0, T, h, x0)
    xknown = known(t)

    fig, ax = plt.subplots()
    ax.plot(t, xknown, label='Conocida')
    ax.plot(t, xe[:, 0], label='Roku4')

    ax.legend()
    plt.show()

    fig, ax = plt.subplots(2)
    fig.suptitle('Comparación soluciones RK4/real')
    ax[0].plot(t, xknown, label='Conocida')
    ax[1].plot(t, xe[:, 0], label='RK4')

    ax[0].legend()
    ax[1].legend()
    plt.show()


    fig, ax = plt.subplots()
    ax.plot(t, xe[:, 0] - xknown, label='Error')
    ax.legend()
    plt.show()

    ## test contra scipy
    s45 = solve_ivp(dx, [0, T], x0, method='RK45', t_eval=t, rtol=1e-13, atol=1e-14)

    fig, ax = plt.subplots()
    ax.plot(s45.t, s45.y[0], label='scipy solve_ivp')
    ax.plot(s45.t, xe[:, 0], label='RK4 implementada')

    ax.legend()
    plt.show()


    lastknown = known(T)


    t2, xe2 = ruku4(dx, 0, T, h/2, x0)
    xknown2 = known(t2)

    error = abs(xe2[-1, 0] - lastknown)
    error1 = abs(xe[:, 0] - xknown)
    error2 = abs(xe2[:, 0] - xknown2)

    fig, ax = plt.subplots()
    ax.plot(t, error1, label='Error con paso h')
    ax.plot(t2, error2, label='Error con paso h/2')
    ax.legend()
    plt.show()


    estimacion = abs(xe[-1, 0] - xe2[-1, 0])/15

    print('Relación entre error con paso h y error con paso h/2 para solucion de circuito RC: ' , error1[-1]/error2[-1])

    print('estimacion de error en resolucion de circuito RC: ' , estimacion)
    print('error en resolucion de circuito RC: ' , error)

    porcentual = (abs(error - estimacion)/error*100)
    print('Error porcentual: ', porcentual , '\n')

    print('Realizando gráfico de error en función de h, por favor espere...\n')

    errorRK(T)


def errorRK(T):

    print('Aca esta el gráfico!')
    n = 5
    N = np.logspace(1, 5, n)
    h = T/N
    errorRK = np.zeros(n)
    x0 = np.zeros(1)
    x = known(T)

    for k in range(n):
        terror, xerror = ruku4(dx, 0, T, h[k], x0)
        errorRK[k] = abs(xerror[-1, 0] - x)

    fig, ax = plt.subplots()
    ax.loglog(h, errorRK, label='Ruku4')
    ax.loglog(h, (h/h[0])**4*errorRK[0], 'r--')
    ax.legend()
    plt.show()

if __name__ == '__main__':
    test()