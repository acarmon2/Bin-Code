# Function to calculate the Bin generators based in article 1997
# The program produce a random search in the permutations to get
# the D vector with F properties

import numpy as np
import math
from random import randrange
import time

#################################################################################################
# Function to produce the random sequence with the first condition
# of Bin construction.
def RandomS(q, d):
    L = q - (2*d) - 1
    V0 = np.zeros(L, dtype = np.int_)
    Seq0 = np.zeros(L, dtype = np.int_)
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
    D = np.zeros((np.int_(L/2),L), dtype = np.int_)
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
    UCodes = np.zeros((q,L), dtype = np.int_)
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

def CheckBin(U0, d):
    x0 = np.shape(U0)
    # rows and columns
    qB = np.int_(x0[0])
    LB = np.int_(x0[1])
    ############################################################################################
    # One coincidence Hamming condition
    cont1 = 0
    for i in range(0, (qB - 1)):
        for j in range((i + 1), qB):
            K0 = HammingCC(U0[i][:], U0[j][:])
            if(np.max(K0) <= 1):
                cont1 = cont1 + 1
            else:
                break
    # Count the total values where the one coincidence is fullfilment
    N1T = np.int_((qB*(qB - 1))/2)
    if(cont1 == N1T):
        flag1 = True
    else:
        flag1 = False
    ##############################################################################################
    # Minimum distance d between adajcent symbols
    contd = 0
    for i in range(0, qB):
        VTest = np.append(U0[i][:], U0[i][0])
        for j in range(0, LB):
            if((abs(VTest[j] - VTest[j + 1])) > d):
                contd = contd + 1
            else:
                break
    # Count the total values where the adjacent distance are major than a value d
    N1d = np.int_(LB*qB)
    if(cont1 == N1T):
        flagd = True
    else:
        flagd = False
    # return the or case beween both flags
    return (flag1 and flagd)

###################################################################################################
# Main function to calculate NT chances of random permutation to find differents generator vectors,
# To save the "possible" N options, the maximum case is 100 in this case
start = time.time()
q = 29
d = 8
NT = 10000
L0 = q - (2*d) - 1
Bcont = 0
BCodes = np.zeros((1000,q,L0), dtype = np.int_)
BGenerators = np.zeros((1000,L0), dtype = np.int_)
# Search the first family of codes
for i in range(0,NT):
    G0 = RandomS(q,d)
    G1 = Dtest(G0,q)
    if(G1):
        BCodes[0][:][:]=UCodes(G0,q)
        BGenerators[0][:] = G0
        Bcont = 1
        print "First Generator founded it"
        break
# Print if there are not generators
if(Bcont == 0):
    print "First trial not generator founded it"
# Second trial
for i in range(0,NT):
    G0 = RandomS(q,d)
    G1 = Dtest(G0,q)
    # Check of the generator is not repetitive with previous generators founded it yet
    cont = 0
    if(G1):
        for k in range(0, Bcont):
            if(np.array_equal(G0,BGenerators[k][:])):
                break
            else:
                cont = cont + 1
        if(cont == Bcont):
            BGenerators[Bcont][:] = G0
            BCodes[Bcont][:][:] = UCodes(G0,q)
            Bcont = Bcont + 1
    # Return the complete tensor with the codes or generator no found
    if(Bcont == 0):
        print "Second trial! Generator not founded it"
        Codes = False
    else:
        #Codes = BCodes[0:Bcont][:][:]
        Codes = BGenerators[0:Bcont][:]

print Bcont
print Codes
end = time.time()
print (end - start)

A = Codes[0][:]
B = Codes[1][:]
# Function to check the cross correlation Hamming function, it returns the maximum value of cross
# correlation function

def HammingCC(C0, C1):
    # Function to check Hamming function between two vectors
    def Hamming(A0, A1):
        Lh = np.size(A0)
        cont = 0
        for i in range(0, Lh):
            if(A0[i] == A1[i]):
                cont = cont + 1
        return cont
    # Main function to move the second vector 
    L0 = np.size(C0)
    CCVector0 = np.zeros(((2*L0)-1), dtype = np.int_)
    # Right movement of vector C1
    for i in range(0, L0):
        C2 = np.roll(C1, i)
        CCVector0[L0 + i - 1] = Hamming(C0, C2)
    # Left movement of vector C1
    for i in range(0, (L0 - 1)):
        C2 = np.roll(C1, -i)
        CCVector0[L0 - i - 1] = Hamming(C0, C2)
    # return the complete vector with delays
    return CCVector0

print A
print B
print HammingCC(A,B)

# Fathallah generator shown in the article
A0 = np.array([13, 9, 12, 11, 14, 10, 16, 20, 17, 18, 15, 19], dtype = np.int_)
A1 = UCodes(A0, 29)
print A1
print CheckBin(A1, d)