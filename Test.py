import numpy as np
import math
from random import randrange

def RandomS(q, d):
    L = q - (2*d) - 1
    V0 = np.zeros(L)
    Seq0 = np.zeros(L)
    for i in range(0,(L/2)):
        V0[i] = d + i + 1
        V0[(L/2)+i]= q - V0[i]
    # Random permutation of the original sequence to check the second 
    # condition.
    T0 = L
    for i in range (0,(L/2)):
        cont = randrange(T0)
        Seq0[i] = V0[cont]
        Seq0[i+(L/2)] = (q - V0[cont])
        V0 = V0[V0 != Seq0[i]]
        V0 = V0[V0 != (q - Seq0[i])]
        T0 = T0 - 2
    return Seq0

print RandomS(19,2)