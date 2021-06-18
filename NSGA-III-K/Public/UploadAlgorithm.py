"""
Upload algorithm.
"""

import importlib

def uploadAlgorithm(algorithm):
    """Uploads main function for a given algorithm"""
    module = importlib.import_module('Algorithms.'+algorithm)
    return getattr(module,'main')
