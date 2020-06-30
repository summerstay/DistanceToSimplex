import numpy as np
from sympy import Matrix


def distanceToSimplex(point, S):
    '''
    An implementation of An Algorithm to Compute the Distance from a Point to a Simplex
    by Oleg Golubitsky, Vadim Mazalov and Stephen M. Watt

    Matlab Implementation: Douglas Summers Stay
    Python Implementation: Snehesh Shrestha
    '''

    n = np.shape(S)[0]

    if n == 1:
        distance = np.linalg.norm(point-S[0])
        projection = np.matrix(S[0]).T

    else:
        translatedS = S-S[0]        

        b = np.matrix((point-S[0])).T
        A = np.matrix((translatedS[1:])).T

        top = np.dot(A.T, A)
        #top = top+.00000001*np.random.rand(np.shape(top)[0], np.shape(top)[1])

        bottom = np.dot(A.T, b)
        #bottom = bottom+.00000001*np.random.rand(np.shape(bottom)[0])

        #alpha = np.divide(bottom, top)
        alpha = np.linalg.solve(top, bottom)

        # _, pivcol = Matrix.rref(A)  # frref(A)
        __A = Matrix(A)
        _, pivcol = __A.rref()  # frref(A)
        A = A[:, pivcol]

        # P = A*inv(A'*A)*A'
        P = np.dot(A, np.dot(np.linalg.inv(np.dot(A.T, A)), A.T))
        pprime = np.dot(P, b)

        posflag = 1
        # for ii in range(0, np.shape(alpha)[0]):
        #     if alpha[ii] < 0:
        if any(alpha < 0):
            posflag = 0

        if sum(alpha) <= 1 and posflag == 1:
            # projection inside simplex
            distance = np.linalg.norm(pprime-b)
            projection = np.matrix(pprime+S[0]).T

        elif posflag == 0:
            Sprime = []

            #Sprime[0] = S[0]
            Sprime.append(S[0])

            for ii in range(0, np.shape(alpha)[0]):
                if alpha[ii] > 0:
                    Sprime.append(S[ii+1])
                    
            distance, projection = distanceToSimplex(point, np.array(Sprime))
        else:
            distance, projection = distanceToSimplex(point, np.array(S[1:]))

    return distance, projection


def test(svecs):
    vecs = np.array(svecs).astype(float)
    distance, projection = distanceToSimplex(vecs[0], vecs[1:])
    return distance, projection

def load_vecs(filename):
    f = open(filename, 'r')
    lines = f.read().split('\n')
    f.close()
    svecs = []
    for line in lines:
        if line.strip() == "":
            continue
        vec = line.split(',')[1].split(' ')[:-1]
        svecs.append(vec)

    return svecs

all_test_cases = []

# SANITY TEST
all_test_cases.append([[.5,.5,.5],[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
all_test_cases.append(load_vecs('test_set1.csv'))
all_test_cases.append(load_vecs('test_set2.csv'))

for sv in all_test_cases:
    distance, _ = test(sv)
    print(distance)
