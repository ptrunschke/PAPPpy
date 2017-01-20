from collections import deque

def load(self, c4n_file, n4e_file): pass

class Refine(object):
    def __init__(self, simplex, mesh):
        simplex.untouched = False

    def exec(self):
        s1,s2 = self.simplex.bisect()
        new_vertex = s1.vertice[s1.refIdx]
        self.mesh.refinementAgenda.append(SplitEdge(self.simplex.refEdge, self.mesh))
        self.mesh.removeSimplex(self.simplex)
        self.mesh.simplices.extent((s1,s2))

class SplitEdge(object):
    def __init__(self, edge, mesh):
        self.__validateArgs__(edge, mesh)
        self.edge = edge
        self.mesh = mesh

    def exec(self):
        v1s = filter(Simplex.untouched, edge[0].simplices)
        v2s = filter(Simplex.untouched, edge[1].simplices)
        unconformal = v1s & v2s # - {self.simplex} # only non-empty if a simplex still has vnew as a hanging node; self.simplex must have been in the agenda
        # unconformal = filter(Simplex.untouched, unconformal)
        self.mesh.refinementAgenda.extent(map(Refine, unconformal))
        self.mesh.refinementAgenda.append(self)

    def __validateArgs__(self, edge, mesh):
        if not all(v.mesh is mesh for v in edge): raise Exception()

class Mesh(object):
    def __init__(self, dim, adim):
        """
        dim : dimension of objects in Triangulation
        adim : dimension of surrounding space
        """
        self.__validateArgs__(dim, adim)

        self.dim = dim
        self.adim = adim
        
        self.__vetices__ = []
        self.__simplices__ = set()

        self.refinementAgenda = decue()

        self.__stage__ = 0

        # # self.edges = ColSparseSymmetric() # edges[i,j] == edges[min(i,j), max(i,j)] - edges[i,:] is O(1), edges[i,j] takes O(ln(len(edges[i,:])))
        # you dont need the edges structure - it suffices to store all edges of a point in that point
        # For this, the point needs to be uniquely represented! (two different simplices, sharing the same point, also point to the same pt)
        # Two methods:
        #     readAFEMMesh - handles the condition internally
        #     addVertex - explicitly checks, if the point already exists.
        # Then:
        # commonSimplices(v1,v2) = v1.simplices & v2.simplices (the simplices sharing a point correspond to the edges going out from that point)
        # Simplex(vs=[v1,v2,...]):
        #     for vi in vs: vi.simplices.append(that simplex)
        # ensure that a simplex can only be initialized with Vertex-objects - this itself shoud ensure uniqueness

        # Mesh has three phases:
        #     0 - just initialized (only vertices can be added (vertices -> property))
        #     1 - all vertices added (after calling 'finishVertexAccumulation()' (does nothing but increment the private 'stage' integer - but good for the user - easy to find errors)
        #             only simplices can be added)
        #     2 - all simplices added (after calling 'finishSimplexAccumulation()' - only markings may be set by mark(vetex) (adds Refine(vertex) to agenda) and refinements executed by refine())

        # A Simplex-object may be annotated to contain information about its degrees of freedom (postpone)

    @property
    def vertices(self): return self.__vertices__

    def addVertex(self, vertex):
        if self.__stage__>0: raise Exception()
        vertex.idx = len(self.__vertices__)
        self.__vertices__.append(vertex)

    @property
    def simplices(self): return self.__simplices__

    def addSimplex(self, simplex):
        if self.__stage__>1: raise Exception()
        # simplex.idx = len(self.__simplices__)
        # self.__simplices__.append(simplex)
        self.__simplices__.add(simplex)
        self.__stage__ = 1

    def removeSimplex(self, simplex):
        # self.mesh.pop(self.simplex) # the simplices may need to carry their index in the mesh.simplices-list
        del self.__simplices__.remove(simplex) # this is stupid - then the idc of the other simplices are not correct anymore!

    def __validateArgs__(self, dim, adim):
        if dim<=0: raise Exception()
        if adim<dim: raise Exception()

    def plot(self):
        self.dim != 2: raise NotImplementedError()
        edges = {}
        for s in self.simplices:
            for i in range(len(s.vertices)):
                edges.add(sorted((s.vertices[i-1], s.vertices[i])))
        for e in edges:
            plt.plot((e[0].coords[0], e[1].coords[0]), (e[0].coords[1], e[1].coords[1]), 'k-')
        for v in self.vertices:
            plt.plot((v.coords[0],), (v.coords[1],), 'ro')
        plt.show()

class Vertex(object): 
    def __init__(self, coord, mesh):
        self.__validateArgs__(coord, mesh)

        self.coord = coord
        self.simplices = set()

        self.mesh = mesh
        mesh.addVertex(self) #TODO: python make method only callable from some position

    def __validateArgs__(self, coord, mesh):
        if len(coords) != mesh.adim: raise Exception()

class Simplex(object): 
    def __init__(self, vertices, mesh):
        self.__validateArgs__(vetices, mesh)

        self.vertices = vertices
        for v in vertices: 
            v.simplices.add(self)

        self.mesh = mesh
        mesh.addSimplex(self)

        self.untouched = True # False when the simplex is sceduled for refinement
        self.__refIdx__ = 0 # the index of the vertex that was added most recently
    
    def __validateArgs__(self, vertices, mesh):
        if not all(isinstance(v, Vertex) for v in vertices): raise Exception()
        if len(vertices)!=self.mesh.dim: raise Exception()

    # @property
    # def refIdx(self): return self.__refIdx__ # the index of the most recently added vertex in self.vertices
    @property
    def refVertex(self): return self.vertices[self.__refIdx__]

    @property
    def refEdge(self):
        raise NotImplementedError()
    @refEdge.setter
    def refEdge(self, val):
        raise NotImplementedError()

    def bisect(self):
        # s1,s2 = Simplex(...), Simplex(...)
        # s1.__refIdx = s2.__refIdx = ???
        raise NotImplementedError()
