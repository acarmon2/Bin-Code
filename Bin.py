# Function to calculate the Bin generators based in article 1997
# The program produce a random search in the permutations to get
# the D vector with F properties

import numpy as np
import math
from random import randrange

def Bin(q, d, NT):
    #################################################################################################
    # Function to produce the random sequence with the first condition
    # of Bin construction.
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
   
    #################################################################################################
    # Function to check the second condition of Bin article, the elements of matrix D must be all 
    # differents. The First row is the permutation vector and the second third and so on are the 
    # D1(j), D2
    def Dtest(TCode, q):
        L = np.size(TCode)
        D = np.zeros((np.int_(L/2),L), dtype=np.int_)
        D[0][:] = TCode
        for i in range(1, (L/2)):
            for j in range(0,L):
                Acc = 0
                for k in range(j,(j+i+1)):
                    Val = np.mod(k, L)
                    Acc = Acc + D[0][Val]
                Temp = np.mod(Acc, q)
                D[i][j] = Temp
        # Check if all the elements in each row are differents, from i = 2 to i = L/2 and 
        # return boolean value
        cont = 0
        for i in range(1,(L/2)):
            if((np.size(np.unique(D[i][:])))==L):
                cont = cont + 1
        if(cont == ((L/2)-1)):
            flag = True
        else:
            flag = False
        return flag


    #################################################################################################
    # Function to produce the User codes using the Generator vector Gen0. The final matrix 
    def UCodes(Gen0, q):
        L = np.size(Gen0)
        UCodes = np.zeros((q,L), dtype=np.int_)
        # Get the generator adding the C0, C0+C1, C0+C1+C2, ....., this is the first row for the
        # total user code
        for i in range(0,L):
            Acc = 0
            for j in range(0,(i+1)):
                Acc = Acc + Gen0[j]
            UCodes[0][i] = np.mod(Acc,q)
        # Fill the total rows using the previous row
        for i in range(1,q):
            UCodes[i][:] = np.mod(((UCodes[i-1][:])+1),q)
        return UCodes
