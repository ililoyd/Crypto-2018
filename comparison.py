'''
Output:
q: 15 - epsilon1: 0.25290131976368646 - epsilon2: 0.2499918702746372
q: 16 - epsilon1: 0.2836040052528501 - epsilon2: 0.28018937563793644
q: 17 - epsilon1: 0.3150076652965609 - epsilon2: 0.31106113346867104
q: 18 - epsilon1: 0.3469114178717896 - epsilon2: 0.34241291970846444
q: 19 - epsilon1: 0.37911852603153695 - epsilon2: 0.37405523755741676
q: 20 - epsilon1: 0.41143838358058027 - epsilon2: 0.40580512747932584
q: 21 - epsilon1: 0.443688335165206 - epsilon2: 0.4374878053458634
q: 22 - epsilon1: 0.4756953076625503 - epsilon2: 0.4689381107801478
q: 23 - epsilon1: 0.5072972343239857 - epsilon2: 0.5000017521827107
q: 24 - epsilon1: 0.538344257914529 - epsilon2: 0.5305363394090516
q: 25 - epsilon1: 0.568699703969464 - epsilon2: 0.5604121995072768
q: 26 - epsilon1: 0.598240820135939 - epsilon2: 0.5895129752132673
q: 27 - epsilon1: 0.6268592822632421 - epsilon2: 0.6177360099457612
q: 28 - epsilon1: 0.6544614723423995 - epsilon2: 0.6449925267628924
q: 29 - epsilon1: 0.6809685374777771 - epsilon2: 0.671207612066945
q: 30 - epsilon1: 0.7063162427192688 - epsilon2: 0.6963200177219244

'''


import math
def compare(M,minQ, maxQ):
    for q in range(minQ, maxQ + 1):
        prod = 1
        for i in range(1,q):
            prod *= (M - i)/M
        epsilon1 = 1 - prod
        epsilon2 = 1 - math.exp(-q * (q-1) / (2*M))
        print("q:", end=" ")
        print(q, end=" - ")
        print("epsilon1:", end=" ")
        print(epsilon1, end=" - ")
        print("epsilon2:", end=" ")
        print(epsilon2)
compare(365,15,30)
