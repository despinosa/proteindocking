from collections import namedtuple
from itertools import accumulate
from functools import reduce
import numpy as np
 
Cavity = namedtuple('Cavity', 'center min max')
 
def cavity_clustering(dataset):
    """Devuelve los puntos en dataset agrupados en clusters y su
    cavidad correspondiente"""

    results = []
    while dataset:
        center = dataset[0]
        min_ = np.array([np.inf, np.inf, np.inf], np.single)
        max_ = np.array([-np.inf, -np.inf, -np.inf], np.single)
        for n, sum_ in enumerate(accumulate(dataset)):
            if np.allclose(dataset[n], center,
                           rtol=0.1, atol=10): # item está cerca de centro
                center = sum_ / (n+1) # nuevo centro
            else: # item no está cerca de centro
                results.append((dataset[:n+1],
                                Cavity(center, min_, max_))) # nueva cavidad
                dataset = dataset[n:] # particionamos dataset
                break
        else:
            results.append((dataset[:n+1],
                            Cavity(center, min_, max_))) # ultima cavidad
            dataset = dataset[n+1:] # particionamos dataset
    return results
 
if __name__ == '__main__':
    from sys import stdin
    dataset = [np.array(line.split(maxsplit=3), np.single) for line in stdin]
    results = cavity_clustering(dataset)
    for item in results: print('\n{} : {}\n'.format(len(item[0]), item[1]))
