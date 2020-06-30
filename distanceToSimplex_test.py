from distanceToSimplex import *

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