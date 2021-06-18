"""
Fast non-dominated sorting procedure

K. Deb, A. Pratap, S. Agarwal and T. Meyarivan, "A fast and elitist 
multiobjective genetic algorithm: NSGA-II," in IEEE Transactions on 
Evolutionary Computation, vol. 6, no. 2, pp. 182-197, 2002.
"""

import numpy as np

def fastNonDominatedSort(evaluation):
    """Divides the population by fronts"""
    N = len(evaluation)
    S = []
    Fronts = []
    Front1 = []
    eta = np.zeros(N)
    ranks = np.zeros(N)
    for p in range(0,N):
        Sp = []
        for q in range(0,N):
            if (all(evaluation[p] <= evaluation[q]) and any(evaluation[p] < evaluation[q])):
                Sp.append(q)
            elif (all(evaluation[q] <= evaluation[p]) and any(evaluation[q] < evaluation[p])):
                eta[p] += 1
        S.append(Sp)
        if (eta[p] == 0):
            Front1.append(p)
            ranks[p] = 1
    Fronts.append(Front1)
    i = 0
    while (Fronts[i] != []):
        Q = []
        for p in Fronts[i]:
            for q in S[p]:
                eta[q] -= 1
                if (eta[q] == 0):
                    ranks[q] = i+2
                    Q.append(q)
        i += 1
        Fronts.append(Q)
    Fronts.pop(len(Fronts)-1)
    return Fronts
