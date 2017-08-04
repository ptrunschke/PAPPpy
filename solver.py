from abc import ABCMeta, abstractproperty, abstractmethod

import numpy as np
norm = np.linalg.norm
argmin = np.argmin
solve = np.linalg.solve

import scipy.sparse as sp

class DiscreteFunction(object):
    def __init__(self, space):
        self.space = space
        self.coeffs = np.zeros(space.dim)

    def __setitem__(self, (simplex,vector), val):
        self.coeffs[self.space.idx(simplex, vector)] = val

    def __call__(self, x):
        t = self.space.triangulation.simplex4point(x) 
        return sum(e(x) for e in self.space.localBasis(t))

class LocalFunction(object):
    def __init__(self, simplex, space, idx): 
        """
        A restriction of the idx'th basis function onto simplex.
        """
        self.simplex = simplex
        self.space = space
        self.idx = idx # via self.space.idx(self.idx) bekommt man den globalen idx

class DiscreteSpace(object): 
    def localFunction(self, simplex, vector):
        # the local function on simplex for the given local basis vector
        return LocalFunction(simplex, self, vector)
    def assemble(self, product): # product takes two 
        ret = sp.dok_matrix((self.dim, self.dim))
        for i in range(self.dim):
            for j in self.indexPatch(i):
                for s in self.commonSimplices(i,j):
                    f1 = self.localFunction(s,i)
                    f2 = self.localFunction(s,j)
                    ret[i,j] += product(f1,f2)
        return ret.to_csr()

class P0Space(DiscreteSpace):
    def __init__(self, triangulation, imgDim=1):
        self.triangulation = triangulation
        self.localDim = imgDim # == 2 for P0(D, R2)
        self.dim = imgDim*len(triangulation.simplices)

    def idx(self, simplex, vector):
        "Returns the index of the basis function/the coefficient associated with the vector'th local basis function on the simplex simplex"
        v = simplex.vertices[1+vector] # v shall represent the vector pointing from vertices[0] to the vertex v
        return self.triangulation.vertices.index(v) #TODO: return a list of Vertex-objects (that contain their index) in .vertices 
        # ==> return simplex.vertices[vector].idx

    def project(self, fnc):
        ret = DiscreteFunction(self)
        for i,v in enumerate(self.triangulation.n4e):
            vs = [np.array(self.triangulation.c4n[n]) for n in v]
            center = sum(vs)/len(vs)
            #TODO: property of simplex...
            ret.coeffs[i] = fnc(center)
        return ret

    def indexPatch(self, idx):
        "Returns the indices for all basis functions that intersect with the domain of the idx'th"
        return {idx}
    def commonSimplices(self, i, j):
        """
        Gibt alle Simplizes zurueck, die in der Domain der Basisfunktion i UND in der der Basisfunktion j liegen.
        """
        if i!=j: return ()
        v = self.triangulation.n4e[i]
        s = Simplex(self.triangulation.c4n[m] for m in v)
        return (s,)
    #TODO: merge indexPatch and commonSimplices to return a list if pairs: (idx, [common simplices])

class P1Space(DiscreteSpace):
    def __init__(self, triangulation):
        self.triangulation = triangulation
        self.dim = len(triangulation.vertices)
        self.localDim = 3

    def idx(self, simplex, vector):
        v = simplex.vertices[:,vector].tolist() # vector is an int #TODO: das indizieren hier ist einfach sooooo falsch... verbessere Simplex!
        print "Vertex", v, "has Index", self.triangulation.vertices.index(v)
        return self.triangulation.vertices.index(v) # the index of v in the triangulation

    def project(self, fnc):
        ret = DiscreteFunction(self)
        for i,v in enumerate(self.triangulation.vertices):
            ret.coeffs[i] = fnc(v)
        return ret

    def indexPatch(self, idx):
        "Returns the indices for all basis functions that intersect with the domain of the idx'th"
        return set(sum((e for e in self.triangulation.n4e if idx in e), []))
    def commonSimplices(self, i, j):
        """
        Gibt alle Simplizes zurueck, die in der Domain der Basisfunktion i UND in der der Basisfunktion j liegen.
        """
        vs = (e for e in self.n4e if i in e and j in e)
        # Das funktioniert, da i und j auch die Indizes der Vertizes in der Triangulierung sind.
        ss = (Simplex(self.c4n[m] for m in v) for v in vs)
        return ss
    #TODO: merge indexPatch and commonSimplices to return a list if pairs: (idx, [common simplices])


class Simplex(object):
    def __init__(self, vertices):
        self.vertices = np.array(list(vertices)).T
        self.__cmat = self.vertices[:,1:] - self.vertices[:,:1]

    def __contains__(self, pnt, eps=1e-5):
        y = pnt-self.vertices[:,0]
        x = solve(self.__cmat, y) # solves __cmat.dot(x) == y
        # y has to be a convex combination of (0,0) and the vectors in __cmat
        return np.all(-eps<=x) and np.all(x<=1+eps) and np.sum(x)<=1+eps

    def __eq__(self, other, eps=1e-5):
        # ich gehe davon aus, dass die vertices in einem simplex nicht permutiert werden
        return norm(self.vertices-other.vertices)<eps

    def __repr__(self): return str(self.vertices)

class Triangulation(object):
    def __init__(self, c4n, n4e):
        self.c4n = c4n
        self.n4e = n4e

    @property
    def vertices(self): return self.c4n

    @property
    def simplices(self): return self.n4e #TODO: return Simplex(...)

    def getSimplex(self, i):  #TODO: deprectated
        v = self.n4e[i]
        s = Simplex(self.c4n[m] for m in v)
        return s

    def simplex4point(self, pnt):
        # compute the index of the closest vertex
        n = argmin(norm(v-pnt) for v in self.vertices)
        # pnt has to be in one of the simplices with vertex n
        vs = (e for e in self.n4e if n in e)
        ss = (Simplex(self.c4n[m] for m in v) for v in vs)
        return next(s for s in ss if pnt in s)

c4n = [
    [0, 0],
    [0, 1],
    [1, 0],
    [.5, .5]
]
n4e = [
    [0,1,3],
    [0,3,2]
]

T = Triangulation(c4n, n4e)
# print T.simplex4point((0,0))
# print T.simplex4point((0.5,0.1))

V = P1Space(T)
W = P0Space(T) 



f = V.project(lambda x: 1-norm(x)**2)
# print f.coeffs
g = W.project(lambda x: 1-x[0])
# print g.coeffs


def div(df): pass # DiscreteFunction(P1Space(T)) -> DiscreteFunction(P0Space(T))
def local_div(lf): # LocalFunction(P1Space(T)) -> (float,float) - the coefficients for the basis functions in P0Space(T,2)
    c = lf.space.idx(lf.simplex, lf.idx) # the coeff
    # c = lf.coeff
    s = lf.simplex.vertices.T #TODO: ueberarbeite Simplex (s.o.)
    c1 = c/norm(s[1]-s[0])
    c2 = c/norm(s[2]-s[0])
    if lf.idx==1:
        return (c1,0)
    elif lf.idx==2:
        return (0,c2)
    else:
        return (-c1,-c2)

p0 = V.localFunction(T.getSimplex(0), 0)
p1 = V.localFunction(T.getSimplex(0), 1)
p2 = V.localFunction(T.getSimplex(0), 2)

print local_div(p0)
print local_div(p1)
print local_div(p2)

def product(lf1, lf2): pass # takes two local functions and returns their inner product

# def P1product(df1, df2): # P1 scalar product
#     ret = 0
#     for i in range(df1..dim):
#         for j in self.indexPatch(i):
#             for s in self.commonSimplices(i,j):
#                 f1 = self.localFunction(s,i)
#                 f2 = self.localFunction(s,j)
#                 ret[i,j] += product(f1,f2)
