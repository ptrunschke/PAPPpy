
def load(self, c4n_file, n4e_file): pass

class Refine(object):
    def __init__(self, simplex, mesh): pass
    def exec(self):
        v1,v2 = self.simplex.refEdge
        s1,s2 = self.simplex.bisect()
        vnew = s1.vertice[s1.__refIdx]
        self.mesh.__refinementAgenda.add(SplitEdge(v1,v2,vnew, self.mesh))

        self.mesh.pop(self.simplex) # the simplices may need to carry their index in the mesh.simplices-list
        self.mesh.simplices += [s1,s2]

class SplitEdge(object):
    def exec(self):
        unconformal = (v1.simplices & v2.simplices) - {self.simplex} # only non-empty if a simplex still has vnew as a hanging node
        self.mesh.__refinementAgenda.extent(map(Refine, unconformal))
        self.mesh.__refinementAgenda.append(self)

class Mesh(object):
    def __init__(self, dim, adim):
        """
        dim : dimension of objects in Triangulation
        adim : dimension of surrounding space
        """
        assert dim>0
        assert adim>=dim

        self.dim = dim
        self.adim = adim
        
        self.vetices = []
        self.simplices = []
        # self.marks = [] # indices of marked simplices

        self.__refinementAgenda = decue

        # self.edges = ColSparseSymmetric() # edges[i,j] == edges[min(i,j), max(i,j)] - edges[i,:] is O(1), edges[i,j] takes O(ln(len(edges[i,:])))
        you dont need the edges structure - it suffices to store all edges of a point in that point
        For this, the point needs to be uniquely represented! (two different simplices, sharing the same point, also point to the same pt)
        Two methods:
            readAFEMMesh - handles the condition internally
            addVertex - explicitly checks, if the point already exists.
        Then:
        commonSimplices(v1,v2) = v1.simplices & v2.simplices (the simplices sharing a point correspond to the edges going out from that point)
        Simplex(vs=[v1,v2,...]):
            for vi in vs: vi.simplices.append(that simplex)
        ensure that a simplex can only be initialized with Vertex-objects - this itself shoud ensure uniqueness

        Mesh has three phases:
            0 - just initialized (only vertices can be added (vertices -> property))
            1 - all vertices added (after calling 'finishVertexAccumulation()' (does nothing but increment the private 'stage' integer - but good for the user - easy to find errors)
                    only simplices can be added)
            2 - all simplices added (after calling 'finishSimplexAccumulation()' - only markings may be set by mark(vetex) (adds Refine(vertex) to agenda) and refinements executed by refine())

        A Simplex-object may be annotated to contain information about its degrees of freedom (postpone)

    # def addVertex(self, coord):
    #     check ...
    #     self.coords.append(coord)

    # def addSimplex(self, simplex):
    #     if isinstance(simplex, Simplex) and simplex.mesh is self:
    #         self.simplices.append(simplex)
    #     else: # simplex = (p1, p2, p3, ...)
    #         assert len(simplex)==self.dim
    #         assert len(simplex[0])==self.adim

    #         simplex = map(Coord, simplex)
    #         pidx = []
    #         for pt in simplex:
    #             pis.append(self.coords.index(pt))
    #             if pidx[-1]<0:
    #                 self.coords.append(pt)
    #                 pidx[-1] = len(self.coords)-1
            
    #         s = Simplex(self, pidx)
    #         for i in pidx
    #         self.simplices.append(s)



class Vertex(object): 
    def __init__(self, coord):
        self.coord = coord
        self.simplices = set()

class Simplex(object): 
    def __init__(self, vertices, mesh=None):
        self.vertices = vertices
        for v in vertices: 
            v.simplices.add(self)
        self.mesh = mesh
        self.__refIdx = 0 # the index of the vertex that was added most recently

    def bisect(self):
        s1,s2 = Simplex(...), Simplex(...)
        s1.__refIdx = s2.__refIdx = ???
