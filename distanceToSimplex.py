import numpy as np
from sympy import Matrx
import random


def distanceToSimplex(point, S):
    '''
    An implementation of An Algorithm to Compute the Distance from a Point to a Simplex
    by Oleg Golubitsky, Vadim Mazalov and Stephen M. Watt

    Matlab Implementation: Douglas Summers Stay
    Python Implementation: Snehesh Shrestha
    '''

    n = np.shape(S)[0]

    if n == 1:
        distance = np.linalg.norm(point-S[0, :])
        projection = np.transpose(S[0, :])

    else:
        translatedS = []
        for t in range(0, n):
            translatedS[t] = S[t, :]-S[0, :]

        b = np.transpose(point-S[0, :])
        A = np.transpose(translatedS[1:])
        top = np.multiply(np.transpose(A), A)
        top = top+.00000001*random.random(np.shape(top))
        bottom = np.multiply(np.transpose(A), b)
        bottom = bottom+.00000001*random.random(np.shape(bottom))
        alpha = np.divide(bottom, top)
        _, pivcol = Matrx.rref(A)  # frref(A)
        A = A[:, pivcol]

        # P = A*inv(A'*A)*A'
        P = np.multiply(A, np.multiply(np.inv(np.multiply(np.transpose(A), A)), np.transpose(A)))
        pprime = np.multiply(P, b)

        posflag = 1
        for ii in range(0, np.shape(alpha)[0]):
            if alpha(ii) < 0:
                posflag = 0

        if sum(alpha) <= 1 and posflag == 1:
            # projection inside simplex
            distance = np.linalg.norm(pprime-b)
            projection = np.transpose(pprime+S[0, :])

        elif posflag == 0:
            Sprime = []
            Sprime[0] = S[0, :]
            count = 1
            for ii in range(0, np.shape(alpha)[0]):
                if alpha[ii] > 0:
                    Sprime[count] = S[ii+1, :]
                    count = count+1
            distance, projection = distanceToSimplex(point, Sprime)
        else:
            distance, projection = distanceToSimplex(point, S[1:, :])

    return distance, projection
