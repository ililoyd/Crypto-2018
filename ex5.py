'''
Output:
31569


'''
def bbsCompute(i, s, p, q):
    n = p * q
    phi = (p - 1) * (q - 1)

    s_i = s ** ( 2 ** i % phi) % n

    return s_i


print(bbsCompute(10000, 17995, 227, 263))
