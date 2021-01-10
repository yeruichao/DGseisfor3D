import numpy
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import meshtools.tetmesh as tetmesh
import meshtools.trimesh as trimesh

xm=1.0;ym=1.0;zm=5.0
dh=1.0
ndh=numpy.ones(4)*1.0

x,y=numpy.meshgrid(numpy.arange(-xm,xm+dh,dh),
        numpy.arange(-ym,ym+dh,dh), sparse=False)
z=numpy.ones(x.shape)*0.0
xy=numpy.array([x.flatten(),y.flatten()]).T
Tri = Delaunay(xy)
print type(Tri.simplices.copy())
print Tri.simplices.copy().shape
tri = numpy.array(map(list, zip(*Tri.simplices.copy())))
xy=numpy.array([x.flatten(),y.flatten(),z.flatten()])
flag = trimesh.tri_facenormal(xy,tri+1)[2]<0
tri[1][flag],tri[2][flag]=tri[2][flag],tri[1][flag]
topsurf = trimesh.trimesh(nod=xy.copy(),ele=tri+1,elemap=False)
tet = tetmesh.tetmesh(nod=xy.copy(),ele=[])
topsurf.ele2tetnod=topsurf.ele.copy()
for i in range(1):
    tet,topsurf = tet.tile(1.0,direction='z',trisurf=topsurf)
#print topsurf.nod
#print topsurf.ele
#tet,topsurf = tet.tile_refineall(1.0,direction='z',trisurf=topsurf)
tet,topsurf = tet.tile_isolated(1.0,direction='z',trisurf=topsurf)
#print topsurf.nod
#print topsurf.ele
#print topsurf.ele2tetnod
#print topsurf.Nnod
print topsurf.Nele
#for i in range(2):
#    tet,topsurf = tet.tile(1.0,direction='z',trisurf=topsurf)
flag=(tet.nod[2][tet.ele[0]-1]
     +tet.nod[2][tet.ele[1]-1]
     +tet.nod[2][tet.ele[2]-1]
     +tet.nod[2][tet.ele[3]-1])/4>1.8
tet1=tetmesh.tetmesh(nod=tet.nod.copy(),ele=tet.ele[:,flag])

#plt.triplot(topsurf.nod[0], topsurf.nod[1], topsurf.ele-1)
#plt.plot(topsurf.nod[0], topsurf.nod[1], 'o')
#plt.show()

#plt.triplot(tet.nod[0], tet.nod[1], topsurf.ele2tetnod-1)
#plt.plot(topsurf.nod[0], topsurf.nod[1], 'o')
#plt.show()

tet.quality()
tet1.output_vtu('rupt_test')

