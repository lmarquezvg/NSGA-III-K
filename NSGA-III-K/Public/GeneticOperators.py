"""
Generate offspring with simulated binary crossover and polynomial mutation.

K. Deb and R. Agrawal, "Simulated binary crossover for continuous search 
space," in Complex Systems, vol. 9, no. 2, pp. 115–148, 1995.

Y. Tian, R. Cheng, X. Zhang, and Y. Jin, "PlatEMO: A MATLAB platform for 
evolutionary multi-objective optimization," in IEEE Computational Intelligence 
Magazine, vol. 12, no. 4, pp. 73-87, 2017.
"""

import numpy as np

class population:
    """Creates a class for population"""
    def __init__(Pop, dec, obj):
        Pop.dec = dec
        Pop.obj = obj

def generateOffspring(P,N,lb,ub,evaluate):
    """Generates offspring population from mating pool"""
    n = np.shape(P.dec)[1]
    m = np.shape(P.obj)[1]
    Q = population(np.zeros([N,n]),np.zeros([N,m]))
    s = 0
    for i in range(0,N//2):
        parentA,parentB = P.dec[s],P.dec[s+1]
        offspringA,offspringB = simulatedBinaryCrossover(parentA, parentB, lb, ub)
        offspringA = polynomialMutation(offspringA, lb, ub)
        offspringB = polynomialMutation(offspringB, lb, ub)
        Q.dec[s] = offspringA
        Q.dec[s+1] = offspringB
        Q.obj[s] = evaluate(offspringA,m)
        Q.obj[s+1] = evaluate(offspringB,m)
        s += 2
    if (N%2 == 1):
        parentA,parentB = P.dec[s],P.dec[s+1]
        offspringA,offspringB = simulatedBinaryCrossover(parentA, parentB, lb, ub)
        offspring = np.copy(offspringA) if np.random.rand() > 0.5 else np.copy(offspringB)
        offspring = polynomialMutation(offspring, lb, ub)
        Q.dec[s] = offspring
        Q.obj[s] = evaluate(offspring,m)
    return Q

def simulatedBinaryCrossover(parentA, parentB, lb, ub):
    """Generates two offspring using simulated binary crossover"""
    n = len(parentA)
    pc, nc = 1, 30
    beta = np.zeros(n)
    mu = np.random.rand(n)
    beta[mu<=0.5] = (2*mu[mu<=0.5])**(1/(nc+1))
    beta[mu>0.5] = (2-2*mu[mu>0.5])**(-1/(nc+1))
    beta *= ((-1)**np.random.randint(0,2,n))
    beta[np.random.rand(n)<=0.5] = 1
    beta[np.matlib.repmat(np.random.rand()>pc,1,n)[0]] = 1
    offspringA = (parentA+parentB)/2+beta*(parentA-parentB)/2
    offspringB = (parentA+parentB)/2-beta*(parentA-parentB)/2
    offspringA = np.minimum(np.maximum(offspringA,lb),ub)
    offspringB = np.minimum(np.maximum(offspringB,lb),ub)
    return offspringA,offspringB

def polynomialMutation(individual, lb, ub):
    """Mutates an individual using polynomial mutation"""
    n = len(individual)
    pm, nm = 1/n, 20
    mutate = np.random.rand(n) <= pm
    mu = np.random.rand(n)
    temp = mutate & (mu <= 0.5)
    individual[temp] += (ub[temp]-lb[temp])*((2*mu[temp]+(1-2*mu[temp])*(1-(individual[temp]-lb[temp])/(ub[temp]-lb[temp]))**(nm+1))**(1/(nm+1))-1)
    temp = mutate & (mu > 0.5)
    individual[temp] += (ub[temp]-lb[temp])*(1-(2*(1-mu[temp])+2*(mu[temp]-0.5)*(1-(ub[temp]-individual[temp])/(ub[temp]-lb[temp]))**(nm+1))**(1/(nm+1)))
    individual = np.minimum(np.maximum(individual,lb),ub)
    return individual
