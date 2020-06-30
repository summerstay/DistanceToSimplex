import numpy as np
from scipy.spatial.distance import pdist


def simplex_content(new_simplex):
    '''
    note: to actually calculate the content, this result needs to be
    multiplied by a factor of (-1)/((-2)^simplex_size)((simplex_size)!)^2). But if you
    are just comparing simplices of the same size, this can be neglected.
    
    Matlab Implementation: Douglas Summers Stay
    Python Implementation: Snehesh Shrestha
    '''
    #P = sqdist(new_simplex, new_simplex)
    P = pdist(new_simplex, new_simplex)
    simplex_size = np.shape(new_simplex)[1]
    Phat = np.ones(simplex_size+1)
    Phat[1, 1] = 0
    Phat[2:, 2:] = P
    content = abs(np.linalg.det(Phat))

    return content
