
from __future__ import division
from sympy import *
import numpy as np
import math
from functools import reduce
from operator import mul

#Modular inverse 
def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    print("-----")
    print(a)
    print(m)
    g, x, y = egcd(a, m)
    if g != 1:
        raise ZeroDivisionError('mod div error')
    else:
        return x % m

#Factorize over the FactorBase
def primeFactorsBase(n, factorBase):
    primes = []
    index = 0
    i = factorBase[index]
    
    while i <= n:

        if n % i == 0:
            primes.append(i)
            n = n // i 
        else:
            index += 1
            if index >= len(factorBase):
                return []
            else:
                i = factorBase[index]
    return primes
#Chinese Remainder Theorem

def crt(M, N):
    #Extended Euclidian Algorithm
    prod = reduce(mul, N, 1)
    result = 0

    for i in range(0, len(M)):
        pp = prod // N[i]
        result += M[i]*pp*modinv(pp,N[i])
        
    return(result % prod)

#LinearDependence(vector, matrix)
def linearDep(v, M):
    if(len(M) ==0):
        return False

    Y = list(M)
    Y.append(v)
    s = np.asarray(v)
    Y = np.asarray(Y)
    
   
    Z = np.asmatrix(Y)
   
    
    numRows, numCols = Z.shape
  
    
    r = np.linalg.matrix_rank(Z)

    if(r == min(numRows,len(v))):
        return False
    else:
        return True
    
def Gauss(N, mod):
    M = np.copy(N)
    n = len(M)
    print('------- GAUSS ---------')
    print(M)

    for q in range(0,n):
        for m in range(0,len(M[0])):
            M[q][m] = M[q][m] % mod

    for k in range(0,len(M[0])-1):
        maxindex = abs(M[k:,k]).argmax() + k

        #Swap
        if maxindex != k:
            M[[k,maxindex]]=M[[maxindex,k]]
        #Pivot to 1
        print(M)
        if(M[k][k] == 0):
                continue
        inv = modinv(M[k][k], mod)
        M[k] = M[k] * inv % mod
        
        #Row Below pivot
        for row in range(k+1,n):
            mult = M[row][k]
            print(M)
            M[row, k:] = (M[row, k:] - mult * M[k,k:]) % mod
            print(M)
    print("out:")
    print(M)
    
    return(M)
            
#Discrete Log    
def dilog(alpha, factorBase, p, beta):
    ReducedCoef = []
    TabCoef = []
    ProcessRow = []
    ProcessTable = []
    #k = np.random.randint(max(factorBase)+1,p-1)
    k = max(factorBase)+1
    C = []

    ######Part 1 Precompute
    while len(TabCoef) < len(factorBase) + 4 :

        n = alpha ** k % p
       
        F = primeFactorsBase(n,factorBase)
        
       
    ###### Part 2 Solve log_alpha of factorBase
        if( len(F) !=0 ):
        #empty list == false
            j = 0
            while len(F) !=0:
                if(F[0] == factorBase[j]):
                   ReducedCoef.append(F.count(F[0]))
                   F = F[F.count(F[0]):len(F)]
                else:
                   ReducedCoef.append(0)
                j = j+1

            while len(ReducedCoef) < len(factorBase):
                ReducedCoef.append(0)

                
            if  not(linearDep(ReducedCoef,TabCoef)) :
                
                TabCoef.append(ReducedCoef)
                #print(ReducedCoef)
                ProcessRow = list(ReducedCoef)
                ProcessRow.append(k)
                ProcessTable.append(ProcessRow)

            ReducedCoef = []
            ReducedF = []
 
        #k = np.random.randint(max(factorBase)+1,p-1)
        k = k + 1
        
    Syst = np.asarray(ProcessTable)
    n=len(Syst)

    #####Chinese Remainder Theorem
    o = 1
    while True:
        if((p-1) % (2**(o))) == 0:
            o = o + 1
        else:
            o = o - 1
            break
    print(o)

    print((p-1) // (2**o))

    Out = Gauss(Syst,(p-1) // (2**o))
    
    A = np.delete(Out,len(Syst[0])-1,axis=1)
    b = np.delete(Out,[range(0,len(Syst[0])-1)],axis=1)

    x = np.zeros(len(Syst[0]) - 1)


    for k in range(len(Syst[0])-2,-1,-1):
        dot = np.dot(A[k][k+1:],x[k+1:]) % ((p-1)//(2**o))
        
        # 1 = modinv(A[k,k],(p-1)//2))
        x[k] = ((b[k] - dot) * 1) % ((p-1)//(2**o)) 
    print(x)
    
    #####
    Out2 = Gauss(Syst, 2)
    print(Out2)

    A2 = np.delete(Out2,len(Syst[0])-1,axis=1)
    b2 = np.delete(Out2,[range(0,len(Syst[0])-1)],axis=1)
    print(A2)
    print(b2)

    y = np.zeros(len(Syst[0]) - 1)

    
    for k in range(len(Syst[0])-2,-1,-1):
        dot = np.dot(A2[k][k+1:],y[k+1:]) % 2
        if A2[k,k] == 0 :
            continue
        y[k] = ((b2[k] - dot) * 1) % 2
    print("*******")
    print(x)
    print(y)
    
    ######## PART 3: trial on s

    '''C = np.zeros(len(Syst[0]) - 1)
    for m in range(0,len(C)):
        C[m] = ((p-1) // (2**o) * (2**(o-1))*y[m] + (((p-1) // (2**o)) + 1) * x[m] ) % (p-1)
    '''
    C = np.zeros(len(Syst[0]) - 1)
    for m in range(0,len(C)):
        C[m] = crt([x[m],y[m]],[(p-1)//2,2])
    

  
        
    loop = True
    s = 0
    while loop:
        trial = (beta * alpha ** s) % p
        F = primeFactorsBase(trial,factorBase)

        #s found
        if(len(F) != 0):
            loop = False
            print(s)
            ReducedCoef = []
            j = 0
            while len(F) !=0:
                if(F[0] == factorBase[j]):
                   ReducedCoef.append(F.count(F[0]))
                   F = F[F.count(F[0]):len(F)]
                else:
                   ReducedCoef.append(0)
                j = j+1

            while len(ReducedCoef) < len(factorBase):
                ReducedCoef.append(0)

            result = -s
            for l in range(0,len(C)):
                result = (result + C[l]*ReducedCoef[l]) % (p-1)
            
        s = s+1

    return result



beta = 9451
#print(dilog(5,base, 14087, beta))
base = [2,3,5,7]
print(dilog(5,base, 10007, beta))

'''
base = [2,3,5,7,11,13,17,19,23]
p = 10930889
alpha = 2317547
beta = 5273437
print(dilog(alpha,base, 59407, beta))
'''



