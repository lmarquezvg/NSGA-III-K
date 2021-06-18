"""
Upload test problem.
"""

import importlib

def uploadTestProblem(problem):
    """Uploads parameters function and evaluate function for a given problem"""
    for i in range(1,8):
        if (problem == 'DTLZ'+str(i)):
            module = importlib.import_module('Problems.DTLZ.'+problem)
    for i in range(1,8):
        if (problem == 'DTLZ'+str(i)+'_MINUS'):
            module = importlib.import_module('Problems.DTLZ_MINUS.'+problem)
    return getattr(module,'parameters'),getattr(module,'evaluate')
