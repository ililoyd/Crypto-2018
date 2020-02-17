#
# output:
# Answer is: 17102
#

import math
from functools import reduce
from operator import mul

def prime_factors(n):
    i = 2
    factors = []
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            factors.append(i)
    if n > 1:
        factors.append(n)
    return factors

def isPrime(x):
    a = 2
    while a <= math.ceil(x/2):
        if x % a == 0:
            return 0 #false
        a += 1
    return 1 #true

def pohligHellmanLoop(alpha, beta, n, q,c):
    p = int(n+1)
    q = int(q)

    A = []
    j = 0
    betaTab = [beta]
    while j <= c-1:
        i=0
        cont = 1
        xp1 = n /(q**(j+1))
        delta = pow(betaTab[j], int(xp1), p)

        #find i such that delta = ...
        while cont:
            xp2 = i*n // q

            if alpha ** xp2 % p == delta:
                cont = 0
                A.append(i)
                xp3 = i * (q**j)
      
                betaTab.append(betaTab[j] * modinv((alpha ** xp3), p) % p)
            i +=1
        j +=1


    return A

#Modular inverse 
def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m
    
#Chinese Remainder Theorem    
def crt(M, N):
    #Extended Euclidian Algorithm
    prod = reduce(mul, N, 1)
    result = 0

    for i in range(0, len(M)):
        pp = prod // N[i]
        result += M[i]*pp*modinv(pp,N[i])
        
    return(result % prod)
    
def pohligHellman(alpha, p, beta):
    n = p - 1
    #vectorized version (a0, ... ac-1)
    N = []  
    A = []
    F = []
    ReducedF = []
    
    F = prime_factors(n)
    #empty list == false
    while len(F) !=0:
        ReducedF.append((F[0],F.count(F[0])))
        F = F[F.count(F[0]):len(F)]
    
    for i in range(0,len(ReducedF)):
        q = ReducedF[i][0]
        c = ReducedF[i][1]
        a = pohligHellmanLoop(alpha, beta, n, q, c)
        m = 0
  
        for i in range(0,len(a)):
            m += a[i] * q**i
        A.append(m)
        N.append(q**c)

    return crt(A, N)

#print(pohligHellman(2,29,18))
#print(pohligHellman(6,8101,7531))
print(pohligHellman(10,31153,12611))
