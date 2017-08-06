from collections import deque
import numpy as np

class Refine(object):
    def __init__(self, simplex, mesh):
        self.simplex = simplex
        self.mesh = mesh
    def run(self):
        self.simplex.bisect()
        # hier liegt das problem - du kannst hier eine seite bisecten, die nicht-konformitaet erzeugt
        # dann wird diese aber nicht aufgeloest, da kein SplitEdge zur agenda hinzugefuegt wird
        # self.mesh.refinementAgenda.append(SplitEdge(self.simplex.refEdge, self.mesh))

class SplitEdge(object):
    def __init__(self, edge, mesh):
        self.edge = edge
        self.mesh = mesh

    def run(self):
        if self.edge.simplices:
            for s in self.edge.simplices:
                # if s.refEdge is self.edge:
                # self.mesh.refinementAgenda.append(SplitEdge(s.refEdge, self.mesh))
                self.mesh.refinementAgenda.append(Refine(s,self.mesh))
            self.mesh.refinementAgenda.append(self)
        else: self.edge.destruct()

        # while self.edge.simplices:
        #     s = next(iter(self.edge.simplices))
        #     s.bisect() # a bisect is not a split of an edge!
        #     # so you may bisect a node but not 

class Mesh(object):
    def __init__(self, dim, adim):
        self.dim = dim # dimension of objects in Triangulation
        self.adim = adim # dimension of surrounding space
        
        self.vertices = []
        self.simplices = set()

        self.refinementAgenda = deque()
        self.refinements = 0

    def addVertex(self, vertex):
        vertex.idx = len(self.vertices)
        self.vertices.append(vertex)

    def addSimplex(self, simplex): self.simplices.add(simplex)
    def removeSimplex(self, simplex): self.simplices.remove(simplex) 

    def refine(self):
        while self.refinementAgenda:
            self.refinementAgenda.popleft().run()
            self.refinements += 1


class Vertex(object): 
    def __init__(self, coord, mesh):
        assert len(coord) == mesh.adim
        self.coord = np.asarray(coord)
        self.edges = set()
        self.mesh = mesh
        mesh.addVertex(self) #TODO: python make method only callable from some position

class Edge(object):
    def __new__(cls, vertices, mesh):
        edges = vertices[0].edges & vertices[1].edges
        if len(edges)==0: return object.__new__(cls)
        elif len(edges)==1: return edges.pop()
        else: raise Exception()

    def __init__(self, vertices, mesh):
        if hasattr(self, 'vertices'): return

        self.vertices = min(vertices), max(vertices) # sorted(vertices)
        self.vertices[0].edges.add(self)
        self.vertices[1].edges.add(self)
        self.mesh =  mesh

        self.marked = False
        self.simplices = set() # simplices that share this edge
        self.__center__ = None

    def mark(self):
        if not self.marked:
            self.mesh.refinementAgenda.append(SplitEdge(self, self.mesh))
            self.marked = True

    def removeSimplex(self, simplex):
        self.simplices.remove(simplex)
        # if not self.simplices: self.destruct() # dont do this - there may be simplices that have this edge but are not created (refined) from larger ones yet

    def destruct(self):
        self.vertices[0].edges.remove(self)
        self.vertices[1].edges.remove(self)

    @property
    def center(self): 
        if self.__center__ is None: self.__center__ = Vertex((self.vertices[0].coord + self.vertices[1].coord)/2, self.mesh)
        return self.__center__

class Simplex(object): 
    def __init__(self, vertices, mesh):
        assert all(isinstance(v, Vertex) for v in vertices)
        assert len(vertices) == mesh.dim+1

        self.vertices = vertices
        self.edges = set()
        for i in range(len(self.vertices)):
            e = Edge((self.vertices[i-1], self.vertices[i]), mesh)
            e.simplices.add(self)
            self.edges.add(e) # why would you need that? to remove the reference to self if self has been refined (bisect)

        self.mesh = mesh
        mesh.addSimplex(self)

        # self.marked = False
        self.__refIdx__ = 0 # the index of the vertex that was added most recently
    
    def contains(self, vertex):
        O = self.vertices[0].coord
        B = np.asarray([v.coord-O for v in self.vertices[1:]])
        V = B.dot(vertex.coord - O)
        return all(V >= 0) and np.linalg.norm(V,1) < 1+1e-10

    def destruct(self):
        for e in self.edges: 
            e.removeSimplex(self)
        self.mesh.removeSimplex(self)

    @property
    def refEdgeIdc(self): 
        d = len(self.vertices)
        return (self.__refIdx__+1)%d, (self.__refIdx__+2)%d

    @property
    def refEdge(self):
        i1,i2 = self.refEdgeIdc
        v1,v2 = self.vertices[i1], self.vertices[i2]
        return Edge((v1,v2), self.mesh)

    def mark(self):
        self.refEdge.mark()
        # if not self.marked:
        #     self.mesh.refinementAgenda.append(Refine(self, self.mesh))
        #     self.marked = True
    
    def bisect(self):
        def replaced(ls, i, v):
            ret = ls[:]
            ret[i] = v
            return ret

        new_vertex = self.refEdge.center
        s1 = Simplex(replaced(self.vertices, self.refEdgeIdc[0], new_vertex), self.mesh)
        s2 = Simplex(replaced(self.vertices, self.refEdgeIdc[1], new_vertex), self.mesh)
        s1.__refIdx__, s2.__refIdx__ = self.refEdgeIdc
        self.destruct()

        return s1,s2

