import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mesh import *

def read_c4n(mesh, c4n):
    cs = c4n.split('\n')[:-1]
    for c in cs:
        coords = map(float, c.split())
        Vertex(coords, mesh)

def read_n4e(mesh, n4e):
    ns = n4e.split('\n')[:-1]
    for n in ns:
        vis = map(int, n.split())
        verts = [mesh.vertices[i] for i in vis]
        Simplex(verts, mesh)

def write_c4n(mesh):
    ret = ""
    for v in mesh.vertices:
        ret += " ".join(map(str, v.coord)) + "\n"
    return ret

def write_n4e(mesh): 
    ret = ""
    for s in mesh.simplices:
        ret += " ".join(str(v.idx) for v in s.vertices) + "\n"
    return ret

def plotMesh(mesh, mark_nodes=True):
    if mesh.dim != 2: raise NotImplementedError()
    edges = set()
    for s in mesh.simplices:
        for i in range(len(s.vertices)):
            edges.add(tuple(sorted((s.vertices[i-1], s.vertices[i]))))
    for e in edges:
        plt.plot((e[0].coord[0], e[1].coord[0]), (e[0].coord[1], e[1].coord[1]), 'k-')
    if mark_nodes:
        for v in mesh.vertices:
            plt.plot((v.coord[0],), (v.coord[1],), 'ro')
    l = plt.xlim()
    sl = sum(l)*0.05
    plt.xlim(l[0]-sl, l[1]+sl)
    l = plt.ylim()
    sl = sum(l)*0.05
    plt.ylim(l[0]-sl, l[1]+sl)

def plotSimplex(simplex):
    if simplex.mesh.adim!=2: raise NotImplementedError()
    plt.gca().add_patch(Polygon(np.array([v.coord for v in s.vertices])))

def unconformingSimplices(self):
    nverts = lambda s: sum(1 for v in self.vertices if s.contains(v))
    return [s for s in self.simplices if nverts(s)>(self.dim+1)]

if __name__=='__main__':

    m = Mesh(2,2)
    read_c4n(m, """\
    0 0
    0 1
    1 0
    """)
    read_n4e(m, """\
    0 1 2
    """)

    # print("Uniform refinement ...")
    # for i in range(2):
    #     for s in m.simplices: s.mark() 
    #     m.refine()
    # print("Done.")

    for s in m.simplices: s.mark() 
    m.refine()
    plt.figure()
    plotMesh(m)

    for s in m.simplices: s.mark() 
    m.refine()
    u = unconformingSimplices(m)
    print("Unconforming simplices: %d"%len(u))
    plt.figure()
    plotMesh(m)
    if u: plotSimplex(u[0])

    plt.show()

    from random import choice
    # mark = lambda: choice(list(m.simplices)).mark()
    mark = lambda: list(m.simplices)[0].mark()

    # print("Random refinement ...")
    # marks = 100
    # for i in range(marks):
    #     mark()
    #     m.refine()
    # print("Done.")
    # print("Marked %d simplices, needed %d additional refinements."%(marks, m.refinements-marks))
    # print("Unconforming simplices: %d"%m.unconforming())

    # mark(); m.refine()
    # mark(); m.refine()
    # mark(); m.refine()
    # mark(); m.refine()
    # mark(); m.refine()
    # mark(); m.refine()
    # mark(); m.refine()

    # plt.figure()
    # m.plot(False)
    # mark(); m.refine()

    # plt.figure()
    # m.plot(False)
    # plt.show()

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
