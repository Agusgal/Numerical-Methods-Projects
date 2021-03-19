import numpy as np

a = np.array([0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0])
b = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0])
c = np.array([1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
d = np.array([0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1])


#Todavia tiene el error, TODO revisar error
def bin2dec(bin_num):
    sign_dec = (-1)**bin_num[0]
    significand_dec = 0

    exponent_dec = 0
    exponent_bin = np.copy(bin_num[1:6])

    if not exponent_bin.any():
        print("son todos 0s")
        exponent_dec = 2**-14
        #el significante es 0,...

    elif exponent_bin.all():
        print("son todos 1s")
        #Debe retornar infinito

    else:
        significand_dec += 1
        number = 0
        for i in range(0, exponent_bin.size):
            number += exponent_bin[i] * 2**(exponent_bin.size - 1 - i)

        exponent_dec = 2**(number - 15)

    mantissa_bin = bin_num[6:]

    var = 0
    for i in range(0, mantissa_bin.size):
        var += mantissa_bin[i] * 2 ** (mantissa_bin.size - 1 - i)

    significand_dec += (var/1024)

    print(sign_dec*exponent_dec*significand_dec)
    return sign_dec*exponent_dec*significand_dec













