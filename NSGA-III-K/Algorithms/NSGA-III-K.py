"""
NSGA-III-K.
"""

import warnings
import numpy as np
from scipy.spatial import distance
from numpy.linalg import LinAlgError

from Public.UploadTestProblem import uploadTestProblem
from Public.RandomPopulation import randomPopulation
from Public.UniformPoints import uniformPoints
from Public.GeneticOperators import generateOffspring
from Public.FastNonDominatedSort import fastNonDominatedSort

def main(H1, H2, flag, problem, m, max_generations):
    """Runs main framework of NSGA-III-K"""
    parameters, evaluate = uploadTestProblem(problem)
    n, lb, ub = parameters(m)
    W = uniformPoints(H1, H2, m)
    N = len(W)
    P = randomPopulation(N, n, m, lb, ub, evaluate)
    zmin = np.min(P.obj, axis=0)
    extreme = None
    generations = 0
    while (generations < max_generations):
        M = matingSelection(P)
        Q = generateOffspring(M, N, lb, ub, evaluate)
        R = population(np.vstack((P.dec, Q.dec)), np.vstack((P.obj, Q.obj)))
        zmin = np.min(np.vstack(([zmin], R.obj)), axis=0)
        P, extreme = environmentalSelection(R, W, N, zmin, extreme, flag)
        generations += 1
    return P

class population:
    """Creates a class for population"""
    def __init__(Pop, dec, obj):
        Pop.dec = dec
        Pop.obj = obj

def matingSelection(P):
    """Selects random parent population"""
    N, n = np.shape(P.dec)
    N, m = np.shape(P.obj)
    N = N+1 if N%2 == 1 else N
    M = population(np.zeros([N,n]), np.zeros([N,m]))
    s = 0
    for i in range(0, N//2):
        indexA = np.random.randint(0, len(P.dec))
        while (True):
            indexB = np.random.randint(0, len(P.dec))
            if (indexA != indexB):
                break
        M.dec[s] = P.dec[indexA]
        M.dec[s+1] = P.dec[indexB]
        M.obj[s] = P.obj[indexA]
        M.obj[s+1] = P.obj[indexB]
        s += 2
    return M

def environmentalSelection(R, W, N, zmin, extreme, flag):
    """Returns population with best individuals"""
    unique = np.sort(np.unique(np.around(R.obj, 6), return_index=True, axis=0)[1])
    R.dec = R.dec[unique]
    R.obj = R.obj[unique]
    Fronts = fastNonDominatedSort(R.obj)
    zmax_front = np.max(R.obj[Fronts[0]], axis=0)
    zmax_pop = np.max(R.obj, axis=0)
    P, Fl = findCriticalFront(R, Fronts, N)
    N1 = len(P.dec)
    if (N1 < N):
        K = N-N1
        S = np.copy(Fl.obj) if N1 == 0 else np.vstack((P.obj, Fl.obj))        
        Sn, extreme = normalize(S, zmin, zmax_front, zmax_pop, extreme)
        pi, d = associate(Sn, W)
        if (N1 > 0):
            rho = np.array([np.count_nonzero(pi[0:N1] == i) for i in range(0, len(W))])
            P = niching(K, rho, pi, d, W, Fl, P)
        else:
            rho = np.array([np.count_nonzero(pi == i) for i in range(0, len(W))])
            P = nichingPairPotential(N, rho, pi, d, Fl, zmin, zmax_front, flag)
    return P, extreme

def findCriticalFront(R, Fronts, N):
    """Returns population in best fronts and population in critical front"""
    P = population([], [])
    for Front in Fronts:
        if (len(P.dec)+len(Front) > N):
            Fl = population(R.dec[Front], R.obj[Front])
            break
        else:
            P.dec = np.copy(R.dec[Front]) if len(P.dec) == 0 else np.vstack((P.dec, R.dec[Front]))
            P.obj = np.copy(R.obj[Front]) if len(P.obj) == 0 else np.vstack((P.obj, R.obj[Front]))
    return P,Fl

def normalize(S, zmin, zmax_front, zmax_pop, extreme):
    """Normalization procedure"""
    N, m = np.shape(S)
    weights = np.eye(m)
    weights[weights == 0] = 1e6
    Sadd = S
    if extreme is not None:
        Sadd = np.concatenate([extreme, Sadd], axis=0)
    Sprime = Sadd-zmin
    Sprime[Sprime < 1e-3] = 0
    Sasf = np.max(Sprime*weights[:,None,:], axis=2)
    I = np.argmin(Sasf, axis=1)
    extreme = Sadd[I, :]
    try:
        M = extreme-zmin
        b = np.ones(m)
        plane = np.linalg.solve(M, b)
        warnings.simplefilter("ignore")
        intercepts = 1/plane
        nadir = zmin+intercepts
        if not np.allclose(np.dot(M, plane), b) or np.any(intercepts <= 1e-6):
            raise LinAlgError()
    except LinAlgError:
        nadir = zmax_front
    b = nadir-zmin <= 1e-6
    nadir[b] = zmax_pop[b]
    denom = nadir-zmin
    denom[denom == 0] = 1e-12
    Sn = (S-zmin)/denom
    return Sn, extreme

def associate(S, W):
    """Returns closest reference point and its distance for each solution"""
    Snorm = np.sqrt(np.sum(S**2, axis=1))[:,np.newaxis]
    Cosine = 1-distance.cdist(S, W, 'cosine')
    Distance = np.matlib.repmat(Snorm, 1, len(W))*np.sqrt(1-Cosine**2)
    pi = np.argmin(Distance, axis=1)
    d = np.min(Distance, axis=1)
    return pi, d

def niching(K, rho, pi, d, W, Fl, P):
    """Selects K solutions from Fl to complete P using the niching procedure"""
    N1 = len(P.dec)
    Choose = np.zeros(len(Fl.dec), bool)
    Zchoose = np.ones(len(W), bool)
    while (np.sum(Choose) < K):
        temp = np.where(Zchoose)[0]
        Jmin = np.where(rho[temp]==min(rho[temp]))[0]
        j = temp[Jmin[np.random.randint(len(Jmin))]]
        I = np.where(np.bitwise_and(Choose==False, pi[N1:]==j))[0]
        if (len(I) != 0):
            s = np.argmin(d[N1+I]) if rho[j] == 0 else np.random.randint(len(I))
            Choose[I[s]] = True
            rho[j] += 1
        else:
            Zchoose[j] = False
    P.dec = np.copy(Fl.dec[Choose]) if N1 == 0 else np.vstack((P.dec, Fl.dec[Choose]))
    P.obj = np.copy(Fl.obj[Choose]) if N1 == 0 else np.vstack((P.obj, Fl.obj[Choose]))
    return P

def nichingPairPotential(N, rho, pi, d, Fl, zmin, zmax_front, flag):
    """Selects N solutions from Fl using niching procedure and a Pair-Potential Function"""
    Choose = np.zeros(len(Fl.dec),bool)
    J = np.where(rho>0)[0]
    for j in J:
        I = np.where(pi==j)[0]
        s = np.argmin(d[I])
        Choose[I[s]] = True
    selected = np.sum(Choose)
    if (selected == N):
        Fl.dec = Fl.dec[Choose]
        Fl.obj = Fl.obj[Choose]
    else:
        Fl.dec = np.vstack((Fl.dec[Choose],Fl.dec[~Choose]))
        Fl.obj = np.vstack((Fl.obj[Choose],Fl.obj[~Choose]))
        A = (Fl.obj-zmin)/(zmax_front-zmin)
        Diss = calculateDissimilarityMatrix(A, flag)
        Memo = np.sum(Diss,axis=1)
        while (len(Fl.dec) > N):
            C = np.abs(Memo)
            worst = np.argmax(C[selected:])+selected
            Diss = np.delete(Diss,worst,axis=0)
            Diss = np.delete(Diss,worst,axis=1)
            Memo = np.sum(Diss,axis=1)
            Fl.dec = np.delete(Fl.dec,worst,axis=0)
            Fl.obj = np.delete(Fl.obj,worst,axis=0)
    return population(Fl.dec, Fl.obj)

def calculateDissimilarityMatrix(A, flag):  
    if (flag == 1):
        m = np.shape(A)[1]
        Diss = calculateDissimilarityMatrixRSE(A, m-1)
    elif (flag == 2):
        Diss = calculateDissimilarityMatrixGAE(A, 512)
    elif (flag == 3):
        Diss = calculateDissimilarityMatrixCOU(A)
    elif (flag == 4):
        Diss = calculateDissimilarityMatrixPT(A, 5, 3, 0.02)
    elif (flag == 5):
        Diss = calculateDissimilarityMatrixMPT(A, 1, 25)
    elif (flag == 6):
        Diss = calculateDissimilarityMatrixGPT(A, 1, 1)
    elif (flag == 7):
        Diss = calculateDissimilarityMatrixKRA(A, 5, 3, 0.02)
    return Diss

def calculateDissimilarityMatrixRSE(A, s):
    """Returns dissimilarity matrix for Riesz s-energy"""
    N = len(A)
    Diss = np.zeros([N,N])
    for i in range(0,N-1):
        for j in range(i+1,N):
            d = np.sqrt(np.sum((A[i]-A[j])**2))
            Diss[i,j] = Diss[j,i] = 1/(d**s)
    return Diss

def calculateDissimilarityMatrixGAE(A, alpha):
    """Returns dissimilarity matrix for Gaussian alpha-energy"""
    N = len(A)
    Diss = np.zeros([N,N])
    for i in range(0,N-1):
        for j in range(i+1,N):
            d = np.sqrt(np.sum((A[i]-A[j])**2))
            Diss[i,j] = Diss[j,i] = np.e**(-alpha*(d**2))
    return Diss

def calculateDissimilarityMatrixCOU(A):
    """Returns dissimilarity matrix for Coulomb's law"""
    k = 1/(4*np.pi*8.854e-12)
    N = len(A)
    Diss = np.zeros([N,N])
    Anorm = np.sqrt(np.sum(A**2, axis=1))
    for i in range(0,N-1):
        for j in range(i+1,N):
            d = np.sqrt(np.sum((A[i]-A[j])**2))
            Diss[i,j] = Diss[j,i] = k*Anorm[i]*Anorm[j]/(d**2)
    return Diss

def calculateDissimilarityMatrixPT(A, V1, V2, alpha):
    """Returns dissimilarity matrix for Pösch-Teller Potential"""
    N = len(A)
    Diss = np.zeros([N,N])
    for i in range(0,N-1):
        for j in range(i+1,N):
            d = np.sqrt(np.sum((A[i]-A[j])**2))
            Diss[i,j] = Diss[j,i] = V1/(np.sin(alpha*d)**2)+V2/(np.cos(alpha*d)**2)
    return Diss

def calculateDissimilarityMatrixMPT(A, D, alpha):
    """Returns dissimilarity matrix for Modified Pösch-Teller Potential"""
    N = len(A)
    Diss = np.zeros([N,N])
    for i in range(0,N-1):
        for j in range(i+1,N):
            d = np.sqrt(np.sum((A[i]-A[j])**2))
            Diss[i,j] = Diss[j,i] = -D/(np.cosh(alpha*d)**2)
    return Diss

def calculateDissimilarityMatrixGPT(A, a, b):
    """Returns dissimilarity matrix for General form of Pösch-Teller Potential"""
    N = len(A)
    Diss = np.zeros([N,N])
    for i in range(0,N-1):
        for j in range(i+1,N):
            d = np.sqrt(np.sum((A[i]-A[j])**2))
            Diss[i,j] = Diss[j,i] = a/(1+np.cos(d))+b/(1-np.cos(d))
    return Diss

def calculateDissimilarityMatrixKRA(A, V1, V2, alpha):
    """Returns dissimilarity matrix for Kratzer Potential"""
    N = len(A)
    Diss = np.zeros([N,N])
    for i in range(0,N-1):
        for j in range(i+1,N):
            d = np.sqrt(np.sum((A[i]-A[j])**2))
            Diss[i,j] = Diss[j,i] = V1*(((d-(1/alpha))/d)**2)+V2
    return Diss
