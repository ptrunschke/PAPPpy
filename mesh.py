from collections import deque
import numpy as np
import matplotlib.pyplot as plt


# class Refine(object):
#     def __init__(self, simplex, mesh):
#         self.simplex = simplex
#         self.mesh = mesh
#     def run(self):
#         s1,s2 = self.simplex.bisect()
#         new_vertex = self.simplex.refEdge.center # vertices[s1.__refIdx__]
#         # self.mesh.refinementAgenda.append(SplitEdge(self.simplex.refEdge, new_vertex, self.mesh)) # das ist unnoetig, falls dieses split-edge bereits in der agenda ist!
#         self.mesh.addSimplex(s1)
#         self.mesh.addSimplex(s2)

class SplitEdge(object):
    def __init__(self, edge, mesh):
        self.__validateArgs__(edge, mesh)
        self.edge = edge
        self.mesh = mesh

    def run(self):
        while self.edge.simplices:
            s = next(iter(self.edge.simplices))
            s.bisect()
        # self.edge.destruct() # has already been destructed by the destructor of the last simplex

        # if not self.edge.simplices: return
        # for s in self.edge.simplices: s.mark()
        # # self.mesh.refinementAgenda.extend(map(lambda s: Refine(s, self.mesh), self.edge.simplices))
        # self.mesh.refinementAgenda.append(self)

    def __validateArgs__(self, edge, mesh):
        if not all(v.mesh is mesh for v in edge): raise Exception()

class Mesh(object):
    def __init__(self, dim, adim):
        """
        dim : dimension of objects in Triangulation
        adim : dimension of surrounding space
        """
        # A Simplex-object may be annotated to contain information about its degrees of freedom (postpone)
        self.__validateArgs__(dim, adim)

        self.dim = dim
        self.adim = adim
        
        self.__vertices__ = []
        self.__simplices__ = set()

        self.refinementAgenda = deque()
        self.refinements = 0

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

    def write_c4n(self):
        ret = ""
        for v in self.vertices:
            ret += " ".join(map(str, v.coord)) + "\n"
        return ret

    def write_n4e(self): 
        ret = ""
        for s in self.simplices:
            ret += " ".join(str(v.idx) for v in s.vertices) + "\n"
        return ret

    @property
    def vertices(self): return self.__vertices__

    def addVertex(self, vertex):
        vertex.idx = len(self.__vertices__)
        self.__vertices__.append(vertex)

    @property
    def simplices(self): return self.__simplices__

    def addSimplex(self, simplex):
        self.__simplices__.add(simplex)

    def removeSimplex(self, simplex):
        self.__simplices__.remove(simplex) 

    def __validateArgs__(self, dim, adim):
        if dim<=0: raise Exception()
        if adim<dim: raise Exception()

    def plot(self, mark_nodes=True):
        if self.dim != 2: raise NotImplementedError()
        edges = set()
        for s in self.simplices:
            for i in range(len(s.vertices)):
                edges.add(tuple(sorted((s.vertices[i-1], s.vertices[i]))))
        for e in edges:
            plt.plot((e[0].coord[0], e[1].coord[0]), (e[0].coord[1], e[1].coord[1]), 'k-')
        if mark_nodes:
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
            self.refinements += 1


class Vertex(object): 
    def __init__(self, coord, mesh):
        self.__validateArgs__(coord, mesh)

        self.coord = np.asarray(coord)
        self.edges = set()

        self.mesh = mesh
        mesh.addVertex(self) #TODO: python make method only callable from some position

    def __validateArgs__(self, coord, mesh):
        if len(coord) != mesh.adim: raise Exception()


class Edge(object):
    def __new__(cls, vertices, mesh):
        edges = vertices[0].edges & vertices[1].edges
        if len(edges)==0: return object.__new__(cls)
        elif len(edges)==1: return edges.pop()
        else: raise Exception()

    def __init__(self, vertices, mesh):
        if hasattr(self, 'vertices'): return
        self.__validateArgs__(vertices, mesh)

        self.vertices = min(vertices), max(vertices)
        self[0].edges.add(self)
        self[1].edges.add(self)
        self.mesh =  mesh

        self.marked = False
        self.simplices = set() # simplices that share this edge
        self.__center__ = None

    def __validateArgs__(self, vertices, mesh):
        pass

    def mark(self):
        if not self.marked:
            self.mesh.refinementAgenda.append(SplitEdge(self, self.mesh))

    def __getitem__(self, idx): return self.vertices[idx]

    def removeSimplex(self, simplex):
        self.simplices.remove(simplex)
        if not self.simplices: self.destruct()

    def destruct(self):
        self[0].edges.remove(self)
        self[1].edges.remove(self)

    @property
    def center(self): 
        if self.__center__ is None: self.__center__ = Vertex((self.vertices[0].coord + self.vertices[1].coord)/2, self.mesh)
        return self.__center__

class Simplex(object): 
    def __init__(self, vertices, mesh):
        self.__validateArgs__(vertices, mesh)

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
    
    def __validateArgs__(self, vertices, mesh):
        if not all(isinstance(v, Vertex) for v in vertices): raise Exception()
        if len(vertices)!=mesh.dim+1: raise Exception()

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

if __name__=='__main__':

    m = Mesh(2,2)
    m.read_c4n("""\
    0 0
    0 1
    1 0
    """)
    m.read_n4e("""\
    0 1 2
    """)

    # print("Uniform refinement ...")
    # for i in range(10):
    #     for s in m.simplices: s.mark() 
    #     m.refine()
    # print("Done.")
    # m.plot()
    # plt.show()

    from random import choice
    mark = lambda: choice(list(m.simplices)).mark()
    # mark = lambda: list(m.simplices)[0].mark()

    print("This refinement produces hanging nodes - i dont know why")

    print("Random refinement ...")
    marks = 100
    for i in range(marks):
        mark()
        m.refine()
    print("Done.")
    print("Marked %d simplices, needed %d additional refinements."%(marks, m.refinements-marks))

    m.plot(False)
    plt.show()

    # [r]andomly [r]efined [r]eference [t]riangle
    # with open("rrrt.c4n","w") as f: f.write(m.write_c4n())
    # with open("rrrt.n4e","w") as f: f.write(m.write_n4e())

    # marksls = range(100,8000,200)
    # add_marks = []
    # for marks in marksls:
    #     for i in range(200):
    #         mark()
    #         m.refine()
    #     # m.refinements -= 200
    #     add_marks.append(m.refinements)

    # plt.plot(marksls, add_marks, 'o-')
    # plt.xlabel("markings")
    # plt.ylabel("refinements")
    # plt.show()
