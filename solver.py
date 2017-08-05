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

    #TODO: you actually dont need that
    def __getitem__(self, (simplex,vector)): 
        return self.coeffs[self.space.idx(simplex, vector)]
    def __setitem__(self, (simplex,vector), val): 
        self.coeffs[self.space.idx(simplex, vector)] = val
    def __call__(self, x):
        t = self.space.triangulation.simplex4point(x) 
        return sum(e(x) for e in self.space.localBasis(t))

    def local(self, simplex):
        # the local function on simplex for the given local basis vector
        return LocalFunction(self, simplex, self.space.indices4simplex(simplex))

# class LocalFunction(object):
#     def __init__(self, simplex, space, idx): 
#         """
#         A restriction of the idx'th basis function onto simplex.
#         """ #TODO: That is not true anymore!!!
#         self.simplex = simplex
#         self.space = space
#         self.idx = idx # via self.space.idx(self.idx) bekommt man den globalen idx

class LocalFunction(object):
    def __init__(self, df, simplex, idcs): # idcs[i] ist der index des coeffs fuer die i-te Basisfunktion in der summe der baseselem ueber dem simplex
        self.df = df # the discrete function (parent)
        self.simplex = simplex
        self.idcs = idcs


class DiscreteSpace(object): 
    def unitLocal(self, simplex, idx):
        # print idx
        df = DiscreteFunction(self)
        df[simplex,idx] = 1
        return LocalFunction(df, simplex, self.indices4simplex(simplex))

    def assemble(self, product): # product takes two 
        ret = sp.dok_matrix((self.dim, self.dim))
        for i in range(self.dim):
            for j in self.indexPatch(i):
                for s in self.commonSimplices(i,j):
                    f1 = self.unitLocal(s,self.ridx(s,i)) # a restriction of the ith basis function on the simplex s
                    f2 = self.unitLocal(s,self.ridx(s,j))
                    ret[i,j] += product(f1,f2)
        return ret.tocsr()

class P0Space(DiscreteSpace):
    def __init__(self, triangulation, imgDim=1):
        self.triangulation = triangulation
        self.localDim = imgDim # == 2 for P0(D, R2)
        self.dim = imgDim*len(triangulation.simplices)
        super(P0Space, self).__init__()

    def indices4simplex(self, simplex):
        return [simplex.index] # [self.triangulation.simplices.index(simplex)]

    def idx(self, simplex, vector):
        "Returns the index of the basis function/the coefficient associated with the vector'th local basis function on the simplex simplex"
        return self.localDim*simplex.index + vector
        # v = simplex.vertices[1+vector].tolist() # v shall represent the vector pointing from vertices[0] to the vertex v
        # print v
        # i = self.triangulation.vertices.index(v) #TODO: return a list of Vertex-objects (that contain their index) in .vertices 
        # print i
        # assert vector in [0,1]
        # return i*self.localDim + vector
        # # ==> return simplex.vertices[vector].idx

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
        super(P1Space, self).__init__()

    def indices4simplex(self, simplex):
        return simplex.triangulation.n4e[simplex.index]
        # return [v.index for v in simplex.vertices]

    def idx(self, simplex, vector):
        # print simplex.index, vector
        v = simplex.vertices[vector].tolist()
        # print "Vertex", v, "has Index", self.triangulation.vertices.index(v)
        return self.triangulation.vertices.index(v) # the index of v in the triangulation

    def ridx(self, simplex, idx):
        v = self.triangulation.vertices[idx]
        ls = [vs.tolist() for vs in simplex.vertices]
        return ls.index(v)

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
        t = self.triangulation
        return (Simplex(t,k) for k,e in enumerate(t.n4e) if i in e and j in e)
        # vs = (e for e in self.triangulation.n4e if i in e and j in e)
        # # Das funktioniert, da i und j auch die Indizes der Vertizes in der Triangulierung sind.
        # ss = (Simplex(self.triangulation.c4n[m] for m in v) for v in vs)
        # return ss
    #TODO: merge indexPatch and commonSimplices to return a list if pairs: (idx, [common simplices])


class Simplex(object):
    def __init__(self, triangulation, index):
        self.triangulation = triangulation
        self.index = index

        ns = triangulation.n4e[index]
        self.vertices = [np.array(triangulation.c4n[n]) for n in ns]

        vs = np.array(self.vertices).T
        self.__cmat = vs[:,1:] - vs[:,:1]

    def __contains__(self, pnt, eps=1e-5):
        y = pnt-self.vertices[0]
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
        return Simplex(self, i)
        # v = self.n4e[i]
        # s = Simplex(self.c4n[m] for m in v)
        # return s

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


def local_div(lf): # LocalFunction(P1Space(T)) -> (float,float) - the coefficients for the basis functions in P0Space(T,2)
    r = 0
    # print lf.idcs
    for i,ci in enumerate(lf.idcs):
        c = lf.df.coeffs[ci]
        s = lf.simplex.vertices
        if i==0: c = -c
        if i in [0,1]:
            p = (s[2]-s[0])
            p = p/norm(p)
            v = s[1] - p*s[1].dot(p)
            c1 = c/norm(s[1]-s[0])
            r += np.array((c1,0))
        if i in [0,2]:
            p = (s[1]-s[0])
            p = p/norm(p)
            v = s[2] - p*s[2].dot(p)
            c2 = c/norm(v)
            r += np.array((0,c2))

    return r


# print T.getSimplex(0)
# p0 = V.unitLocal(T.getSimplex(0), 0)
# p1 = V.unitLocal(T.getSimplex(0), 1)
# p2 = V.unitLocal(T.getSimplex(0), 2)
# print local_div(p0)
# print local_div(p1)
# print local_div(p2)

def div(df): # DiscreteFunction(P1Space(T)) -> DiscreteFunction(P0Space(T))
    t = df.space.triangulation
    ret = DiscreteFunction(P0Space(t, 2))
    for si in range(len(t.simplices)):
        s = Simplex(t,si)
        ld = local_div(df.local(s))
        ret[s,0] = ld[0]
        ret[s,1] = ld[1]
    return ret

df = div(f)
# print df.coeffs

# note, that the sum,prod,diff,div of two discrete functions is given by the exact same operation on their coeff.-arrays

def P0product(f1, f2): # takes two P0 functions and returns their inner product
    xs = f1.coeffs.reshape(-1, f1.space.localDim)
    ys = f2.coeffs.reshape(-1, f2.space.localDim)
    r = 0
    for x,y in zip(xs,ys): r += x.dot(y)
    return r

local_a = lambda lu,lv: np.dot(local_div(lu), local_div(lv))
a = lambda u,v: P0product(div(u),div(v))
# a2 = lambda u,v: sum(local_a(u.local(s),v.local(s)) for s in u.space.triangulation.simplices)
# print a(g,g)
# print a2(g,g)

A = V.assemble(local_a)
print A.todense()

# def handleDBC(discretization, matrix, nodes, values): #ist es schlau nodes zu uebergeben?
#     """
#     Dirichlet boundary condition
#     nodes[i] hat values[i]
#     """
#     pass
