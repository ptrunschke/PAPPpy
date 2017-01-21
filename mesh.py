from collections import deque
import numpy as np
import matplotlib.pyplot as plt

def replace(ls, i, v):
    ret = ls[:]
    ret[i] = v
    return ret

class Refine(object):
    def __init__(self, simplex, mesh):
        self.simplex = simplex
        self.mesh = mesh
        simplex.untouched = False

    def run(self):
        s1,s2 = self.simplex.bisect()
        new_vertex = s1.vertices[s1.__refIdx__]
        self.mesh.refinementAgenda.append(SplitEdge(self.simplex.refEdge, new_vertex, self.mesh))
        self.mesh.removeSimplex(self.simplex)
        self.mesh.addSimplex(s1)
        self.mesh.addSimplex(s2)

class SplitEdge(object):
    def __init__(self, edge, vertex, mesh):
        self.__validateArgs__(edge, mesh)
        self.edge = edge
        self.mesh = mesh

    def run(self):
        v1s = {s for s in self.edge[0].simplices if s.untouched}
        v2s = {s for s in self.edge[1].simplices if s.untouched}
        unconformal = v1s & v2s # - {self.simplex} # only non-empty if a simplex still has vnew as a hanging node; self.simplex must have been in the agenda
        if not unconformal: return
        self.mesh.refinementAgenda.extend(map(lambda s: Refine(s, self.mesh), unconformal))
        self.mesh.refinementAgenda.append(self)

    def __validateArgs__(self, edge, mesh):
        if not all(v.mesh is mesh for v in edge): raise Exception()

class Mesh(object):
    def __init__(self, dim, adim):
        """
        dim : dimension of objects in Triangulation
        adim : dimension of surrounding space
        """
        # Mesh has two phases:
        #     0 - just initialized (only vertices can be added (vertices -> property))
        #     1 - all vertices added (after calling 'finishVertexAccumulation()' (does nothing but increment the private 'stage' integer - but good for the user - easy to find errors)
        #             only simplices can be added)

        self.__validateArgs__(dim, adim)

        self.dim = dim
        self.adim = adim
        
        self.__vertices__ = []
        self.__simplices__ = set()

        self.refinementAgenda = deque()

        self.__stage__ = 0

        # For this, the point needs to be uniquely represented! (two different simplices, sharing the same point, also point to the same pt)
        # ensure that a simplex can only be initialized with Vertex-objects - this itself shoud ensure uniqueness

        # commonSimplices(v1,v2) = v1.simplices & v2.simplices (the simplices sharing a point correspond to the edges going out from that point)
        # Simplex(vs=[v1,v2,...]):
        #     for vi in vs: vi.simplices.append(that simplex)

        # A Simplex-object may be annotated to contain information about its degrees of freedom (postpone)

    def read_c4n(self, c4n):
        cs = c4n.split('\n')[:-1]
        for c in cs:
            coords = map(float, c.split())
            Vertex(coords, self)

    def read_n4e(self, n4e):
        ns = n4e.split('\n')[:-1]
        for n in ns:
            vis = map(int, n.split())
            verts = [self.vertices[i] for i in vis]
            Simplex(verts, self)

    def write_c4n(self): raise NotImplementedError()
    def write_n4e(self): raise NotImplementedError()

    @property
    def vertices(self): return self.__vertices__

    def addVertex(self, vertex):
        # if self.__stage__>0: raise Exception()
        vertex.idx = len(self.__vertices__)
        self.__vertices__.append(vertex)

    @property
    def simplices(self): return self.__simplices__

    def addSimplex(self, simplex):
        # if self.__stage__>1: raise Exception()
        # simplex.idx = len(self.__simplices__)
        # self.__simplices__.append(simplex)
        self.__simplices__.add(simplex)
        # self.__stage__ = 1

    def removeSimplex(self, simplex):
        for v in simplex.vertices: #TODO: dont do this here
            v.simplices.remove(simplex)
        self.__simplices__.remove(simplex) 

    def __validateArgs__(self, dim, adim):
        if dim<=0: raise Exception()
        if adim<dim: raise Exception()

    def plot(self):
        if self.dim != 2: raise NotImplementedError()
        edges = set()
        for s in self.simplices:
            for i in range(len(s.vertices)):
                edges.add(tuple(sorted((s.vertices[i-1], s.vertices[i]))))
        for e in edges:
            plt.plot((e[0].coord[0], e[1].coord[0]), (e[0].coord[1], e[1].coord[1]), 'k-')
        for v in self.vertices:
            plt.plot((v.coord[0],), (v.coord[1],), 'ro')
        l = plt.xlim()
        sl = sum(l)*0.05
        plt.xlim(l[0]-sl, l[1]+sl)
        l = plt.ylim()
        sl = sum(l)*0.05
        plt.ylim(l[0]-sl, l[1]+sl)

    def refine(self):
        while self.refinementAgenda:
            # print self.refinementAgenda
            self.refinementAgenda.popleft().run()

class Vertex(object): 
    def __init__(self, coord, mesh):
        self.__validateArgs__(coord, mesh)

        self.coord = np.asarray(coord)
        self.simplices = set()

        self.mesh = mesh
        mesh.addVertex(self) #TODO: python make method only callable from some position

    def __validateArgs__(self, coord, mesh):
        if len(coord) != mesh.adim: raise Exception()

class Simplex(object): 
    def __init__(self, vertices, mesh):
        self.__validateArgs__(vertices, mesh)

        self.vertices = vertices
        for v in vertices: 
            v.simplices.add(self)

        self.mesh = mesh
        mesh.addSimplex(self)

        self.refinedEdges = dict()

        self.untouched = True # False when the simplex is sceduled for refinement
        self.__refIdx__ = 0 # the index of the vertex that was added most recently
    
    def __validateArgs__(self, vertices, mesh):
        if not all(isinstance(v, Vertex) for v in vertices): raise Exception()
        if len(vertices)!=mesh.dim+1: raise Exception()

    # @property
    # def refIdx(self): return self.__refIdx__ # the index of the most recently added vertex in self.vertices
    @property
    def refVertex(self): return self.vertices[self.__refIdx__]

    @property
    def refEdgeIdc(self): 
        d = len(self.vertices)
        return (self.__refIdx__+1)%d, (self.__refIdx__+2)%d
    @property
    def refEdge(self):
        i1,i2 = self.refEdgeIdc
        v1,v2 = self.vertices[i1], self.vertices[i2]
        return min(v1,v2), max(v1,v2)
    @refEdge.setter
    def refEdge(self, val):
        raise NotImplementedError()

    def mark(self):
        if self.untouched:
            self.mesh.refinementAgenda.append(Refine(self, self.mesh))
            # self.untouched = False # here or in Refine.init?
    
    @property
    def edges(self):
        return [tuple(sorted((self.vertices[i-1], self.vertices[i]))) for i in range(len(self.vertices))]

    def bisect(self):
        global bisections 
        bisections += 1

        # print("I'm not sure if the implementation is correct in nD.")
        if self.refEdge in self.refinedEdges:
            vn = self.refinedEdges[self.refEdge]
            del self.refinedEdges[self.refEdge]
        else:
            v1,v2 = self.refEdge
            vn = Vertex((v1.coord + v2.coord)/2, self.mesh)
            for s in v1.simplices & v2.simplices:
                s.refinedEdges[(v1,v2)] = vn

        i1,i2 = self.refEdgeIdc
        s1 = Simplex(replace(self.vertices, i1, vn), self.mesh)
        s1.__refIdx__ = i1
        s1.refinedEdges = {e:v for e,v in self.refinedEdges.items() if e in s1.edges}
        s2 = Simplex(replace(self.vertices, i2, vn), self.mesh)
        s2.refinedEdges = {e:v for e,v in self.refinedEdges.items() if e in s2.edges}
        s2.__refIdx__ = i2
        return s1,s2
bisections = 0

if __name__=='__main__':
    m = Mesh(2,2)
    m.read_c4n("""\
    0 0
    0 1
    1 0
    1 3
    """)
    m.read_n4e("""\
    0 1 2
    3 1 2
    """)

    # print("Uniform refinement ...")
    # for i in range(12):
    #     for s in m.simplices: s.mark() 
    #     m.refine()
    # print("Done.")
    # m.plot()
    # plt.show()

    from random import choice
    mark = lambda: choice(list(m.simplices)).mark()
    # mark = lambda: list(m.simplices)[0].mark()

    print("Random refinement ...")
    for i in range(1000):
        mark()
        m.refine()
    print("Done.")
    print("Marked 100 simplices, needed %d additional refinements."%(bisections-100))
    m.plot()
    plt.show()
