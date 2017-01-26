from abc import ABCMeta, abstractmethod, abstractproperty
import numpy as np

# class Refine(object):
#     def __init__(self, simplex, mesh):
#         self.simplex = simplex
#         self.mesh = mesh
#         simplex.untouched = False
#     def run(self):
#         s1,s2 = self.simplex.bisect()
#         new_vertex = s1.vertices[s1.__refIdx__]
#         self.mesh.refinementAgenda.append(SplitEdge(self.simplex.refEdge, new_vertex, self.mesh))
#         self.mesh.removeSimplex(self.simplex)
#         self.mesh.addSimplex(s1)
#         self.mesh.addSimplex(s2)

# class SplitEdge(object):
#     def __init__(self, edge, vertex, mesh):
#         self.__validateArgs__(edge, mesh)
#         self.edge = edge
#         self.mesh = mesh
#     def run(self):
#         v1s = {s for s in self.edge[0].simplices if s.untouched}
#         v2s = {s for s in self.edge[1].simplices if s.untouched}
#         unconformal = v1s & v2s # - {self.simplex} # only non-empty if a simplex still has vnew as a hanging node; self.simplex must have been in the agenda
#         if not unconformal: return
#         self.mesh.refinementAgenda.extend(map(lambda s: Refine(s, self.mesh), unconformal))
#         self.mesh.refinementAgenda.append(self)
#     def __validateArgs__(self, edge, mesh):
#         if not all(v.mesh is mesh for v in edge): raise Exception()

class MeshABC(object):
    __metaclass__ = ABCMeta

    @abstractproperty
    def dim(self): pass
    @abstractproperty
    def adim(self): pass
    @abstractproperty
    def vertices(self): pass
    @abstractproperty
    def simplices(self): pass
    @abstractproperty
    def refinementAgenda(self): pass

    def read_c4n(self, c4n):
        cs = c4n.split('\n')[:-1]
        for c in cs:
            coords = map(float, c.split())
            self.addVertex(coords, self)

    def read_n4e(self, n4e):
        ns = n4e.split('\n')[:-1]
        for n in ns:
            vis = map(int, n.split())
            verts = [self.vertices[i] for i in vis]
            self.addSimplex(verts, self)

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

    @abstractmethod
    def addVertex(self, vertex): pass

    @abstractmethod
    def addSimplex(self, simplex): pass

    # @abstractmethod
    # def removeVertex(self, vertex): pass

    @abstractmethod
    def removeSimplex(self, simplex): pass

    def refine(self):
        while self.refinementAgenda:
            edge = self.refinementAgenda.popleft()
            done = True
            for s in edge.simplices: #TODO: rename
                if s.refEdge is edge:
                    s.bisect()
                else:
                    s.mark()
                    done = False
            if not done:
                # edge.mark()
                self.refinementAgenda.apprend(edge)

    print("also has attributes self.Vertex(the vertex class) and self.Simplex as well as a private class self.__Edge__")

def unitSimplex(dim):
    zero = np.zeros(dim)
    c4n = " ".join(map(str, zero)) + "\n"
    for i in range(dim):
        zero[i] = 1
        c4n += " ".join(map(str, zero)) + "\n"
        zero[i] = 0
    
    n4e = " ".join(map(str, range(dim+1))) + "\n"

    return c4n, n4e

class VertexABC(object): 
    __metaclass__ = ABCMeta

    @abstractproperty
    def coord(self): pass
    @abstractproperty
    def simplices(self): pass

class EdgeABC(object):
    @abstractmethod
    def bisect(self): pass # Simplex.bisect benoetigt nicht wirklich edge.center - nur edge.bisect()

class AugEdge(EdgeABC): #TODO: wirklich vererben?
    __metaclass__ = ABCMeta

    # @abstractmethod
    # def premarking_cb(self): pass #TODO: use this OR self.augmentation???

    # def mark(self):
    #     self.marked = True
    #     self.premarking_cb()
    #     self.mesh.refinementAgenda.append(Bisect(self))
    #     if self.postmarking_cb:
    #         self.postmarking_cb(self.augmentation)

    def bisect(self):
        if self.prebisection_cb(self.augmentation):
            ...
            self.postbisection_cb(self.augmentation)
        else:
            self.mesh.refinementAgenda.append(self)

Rename:
    mesh -> triangulation
    plot -> plot_mesh

class triangulation: contains mesh and simplices
class mesh: contains vertices and edges

oder: du defst deine triangulation nur durch simplizes und vertizes und 
dann muessen aber die augmentations die edges speichern (wenn optimiert wird)

def simplex_prebisection_cb(edge, simplices):
    for s in simplices:
        if s.refEdge is not edge:
            s.refine()
            #TODO: ignore the marking phase! there is no need to mark the edges before splitting them
            # (the only reason would be performence but since we dont know how often we need to split a single simplex for completion)
            # if one wants a marking phase one can artificially implement a "refinementAgenda" which then just splits the yet unsplitted objects one after another
            # actually, the refinementAgenda in mesh_slow is equivalent to the call stack that this recursive algorithm would create

AugEdge.prebisection_cb = simplex_prebisection_cb


und im prinzip braucht simplex auch seine knoten nicht zu kennen, sondern nur seine edges.
in der refinementAgenda werden dann die edges bisectet (BisectionCommand-Klasse) und die an ihnen haengenden simplices geupdated
(es wird getestet, ob eine edge bisected wurde und gegebenenfalls auch der simplex ersetzt) (UpdateCommand-Klasse)

das problem ist, dass man den simplex evtl nicht an dieser edge splitten moechte!
mache also folgendes:
    augmentiere die edge (-> AugEdge-Klasse) um eine Liste ihrer Simplizes
    diese AugEdge-Klasse erwartet fuer jede augmentierung zwei callbacks (default: None) (pre-mark und post-mark callbacks)
    die bei einem split ausgefuehrt werden.
    in diesem Fall fuer die pre-mark: lambda augmentation: map(lambda simplex: simplex.mark(), augmentation) # das muss nicht terminieren! -> lieber nur post mark!
    und simplex enthaelt eine liste seiner edges und eine refEdge property
    es markiert dann seine refEdge



# class Simplex(object): 
#     __metaclass__ = ABCMeta
#         self.vertices
#         self.marked
#         self.__refIdx__
#     @property
#     def refEdge(self): 
#     print("das bisherige verfahren i->(i+1, i+2) funzt nicht im 3D und hoeher: es gibt 6 Kanten aber nur 4 Vertices!")
#     print("es ist sowieso in vielen Faellen schlecht. -> speichere provate eine sortierte liste von edges")
#     print("Moeglichkeit 1: fuege die beim splitten entstehenden, neuen edges wie folgt hinzu - erst alle, die nicht die bisection edge sind, dann die bisection edge")
#     print("Moeglichkeit 2: sortiere die neuen edges nach ihrer laenge (lang->kurz) und extende die liste (dauert laenger ist aber bei unregelmaessigen meshes besser)")
#     def mark(self):
#     def bisect(self):

if __name__=='__main__':

    c4n_2d = """\
    0 0
    0 1
    1 0
    """
    n4e_2d = """\
    0 1 2
    """

    print(unitSimplex(2) == (c4n_2d,n4e_2d))

    print
    print("4D unit Simplex:")
    print("================")
    print
    us4 = unitSimplex(4)
    print("C4N:")
    print("----")
    print us4[0]
    print("N4E:")
    print("----")
    print us4[1]
