import numpy
import tetmesh
import mesh_util

def ball(rad=1.0,cnp=numpy.zeros(3),refine=2,att=[]):
    nod=numpy.array([  0.,  0., 0.,  -1., 0., 0., \
        1., 0.,  0.,   0., -1., 0.,   0., 1., 0., \
        0., 0., -1.,   0.,  0., 1.] ).reshape((7,3)).T
    ele=numpy.array([1,2,4,7, 1,4,3,7, 1,3,5,7, 1,5,2,7, \
        6,2,4,1, 6,4,3,1, 6,3,5,1, 6,5,2,1],dtype=int).reshape(8,4).T
    tet=tetmesh.tetmesh(nod,ele,att=att)
    for i in range(refine):
        tet=tet.all_refine()
        tri=tet.bnd_trimesh()
        r=mesh_util.L2norm(tri.nod)
        tri.nod=tri.nod/numpy.tile(r,(3,1))
        tet.nod=tet.nod.T
        tet.nod[tri.nid2tetnod-1]=tri.nod.T
        tet.nod=tet.nod.T
    tet.nod=tet.nod*rad
    for i in range(3):
        tet.nod[i]=tet.nod[i]+cnp[i]
    return tet

def sphere(rad=1.0,cnp=numpy.zeros(3),refine=2,att=[]):
    tet=ball(rad,cnp=cnp,refine=refine)
    tri=tet.bnd_trimesh()
    tri.att=att
    return tri
