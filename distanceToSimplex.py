import numpy as np
from sympy import Matrix


def distanceToSimplex(point, S):
    '''
    An implementation of An Algorithm to Compute the Distance from a Point to a Simplex
    by Oleg Golubitsky, Vadim Mazalov and Stephen M. Watt

    Matlab Implementation: Douglas Summer Stay
    Pythong Implementation: Snehesh Shrestha
    '''

    n = np.shape(S)[0]

    if n == 1:
        distance = np.linalg.norm(point-S[0])
        projection = S[0].T

    else:
        translatedS = S
        for t in range(0, n):
            translatedS[t] = S[t]-S[0]

        b = (point-S[0]).T
        A = (translatedS[1:]).T

        top = np.dot(A.T, A)
        top = top+.00000001*np.random.rand(np.shape(top)[0], np.shape(top)[1])

        bottom = np.dot(A.T, b)
        bottom = bottom+.00000001*np.random.rand(np.shape(bottom)[0])

        #alpha = np.divide(bottom, top)
        alpha = np.linalg.solve(top, bottom)

        # _, pivcol = Matrix.rref(A)  # frref(A)
        __A = Matrix(A)
        _, pivcol = __A.rref()  # frref(A)
        A = A[:, pivcol]

        # P = A*inv(A'*A)*A'
        P = np.dot(A, np.dot(np.linalg.inv(np.dot(A.T, A)), A.T))
        pprime = np.dot(P, b)

        ############ TODO: CURRENT ISSUE ############
        posflag = 1
        for ii in range(0, np.shape(alpha)[0]):
            if alpha[ii] < 0:
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


# SANITY TEST
f = open('test_set1.csv', 'r')
lines = f.read().split('\n')
f.close()

svecs = []
for line in lines:
    if line.strip() == "":
        continue
    vec = line.split(',')[1].split(' ')[:-1]
    svecs.append(vec)

vecs = np.array(svecs).astype(float)

distance, projection = distanceToSimplex(vecs[0], vecs[1:])
print(distance)
