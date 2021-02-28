import numpy
import sys
sys.path.insert(0,'../../')
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import meshtools.tetmesh as tetmesh
import meshtools.trimesh as trimesh
import meshtools.mesh_util as mesh_util
import meshtools.DG_util as DG_util
import meshtools.vis_util as vis_util

def boxcar(x,W,w):
    X=numpy.absolute(x)
    F=numpy.zeros(x.shape)
    F[X<=W]=1.
    flag=numpy.logical_and(X>W, X<W+w)
    F[flag]=0.5*(1.+numpy.tanh(w/(X[flag]-W-w)+w/(X[flag]-W)))
    return F

## generate fault block volume mesh

xm=6.0;ym=6.0
dh=2.0
N=int(xm/dh)
fname='rupt_N%d'%N

x,y=numpy.meshgrid(numpy.arange(-xm,xm+dh,dh),
        numpy.arange(-ym,ym+dh,dh), sparse=False)
z=numpy.ones(x.shape)*0.0
xy=numpy.array([x.flatten(),y.flatten()]).T
Tri = Delaunay(xy)
#tri = numpy.array(map(list, zip(*Tri.simplices.copy())))
tri = numpy.array(Tri.simplices.copy()).T
print (tri.shape)
xy=numpy.array([x.flatten(),y.flatten(),z.flatten()])
flag = trimesh.tri_facenormal(xy,tri+1)[2]<0
tri[1][flag],tri[2][flag]=tri[2][flag],tri[1][flag]
topsurf = trimesh.trimesh(nod=xy.copy(),ele=tri+1,elemap=False)
tet = tetmesh.tetmesh(nod=xy.copy(),ele=[])
topsurf.ele2tetnod=topsurf.ele.copy()

for i in range(max(1,N-2)):
    tet,topsurf = tet.tile(dh,direction='z',trisurf=topsurf,att=1)
tet,topsurf = tet.tile_refineall(dh,direction='z',trisurf=topsurf,att=1)
#tet,topsurf = tet.tile_refineall(dh/2,direction='z',trisurf=topsurf,att=1)
#tet,topsurf = tet.tile_isolated(dh/2,direction='z',trisurf=topsurf,att=1)
tet,topsurf = tet.tile_isolated(dh,direction='z',trisurf=topsurf,att=1)

## mirror reflect block volume mesh and glue together

tet.nod[2]=tet.nod[2].max()-tet.nod[2]
for i in range(3):
    print (tet.nod[i].min(),tet.nod[i].max())
nod=tet.nod.copy()
ele=tet.ele.copy()
tet.ele[[0,1],:]=tet.ele[[1,0],:]
att=-tet.att.copy()
top_tri=topsurf.ele2tetnod.copy()
nod[2]=-nod[2]
ele=mesh_util.vtx_merge(ele, top_tri, top_tri+tet.Nnod)
ele=numpy.hstack((ele,tet.ele+tet.Nnod))
nod=numpy.hstack((nod,tet.nod))
att=numpy.hstack((att,tet.att))
Tet = tetmesh.tetmesh(nod=nod,ele=ele,att=att).pickeles().elemap()
Rupt= Tet.int_trimesh()

## output volume mesh

Tet.quality()
Tet.output_tetgen(fname)
Rupt_tet=Tet.pickeles(Rupt.F2T[1]-1)

TetBnd=Tet.bnd_trimesh()
cnp=TetBnd.centers()
ID=numpy.logical_or(cnp[0]>5.9,cnp[2]>5.9)
ID=numpy.logical_or(ID,cnp[1]<-5.9)
TetBnd=TetBnd.pickeles(ID)
## output rupture surface mesh

f = open(fname+'.face', 'w')
f.write('%d\n'%Rupt.Nele)
for i in range(Rupt.Nele):
    f.write('%d\t%d\t%d\t%d\t%d\t%d\n'
        %(tuple(Rupt.ele2tetnod.T[i])+tuple(Rupt.F2T.T[i])+(i+1,)))
f.close()

## output rupture coefficient

a=0.008;
da=0.008;
b=0.012;
v0=1e-9;
f0=0.6;
L=0.02e-3;
Wx = 4.0; # half-length of fault
wx = 2.0; # width of transition region in x-direction
Wy = 6.0; # half-width of fault
wy = 2.0; # width of transition region in y-direction
x0 = 0; # hypocenter, km
y0 = 2.0; # hypocenter, km
dVt0=1e-15;
porder=3

#for i in range(3):
#    print (Rupt.nod[i].min(), Rupt.nod[i].max())

p2n=DG_util.porder2nodal()
x,y,z=p2n.tri_nodal_coord(Rupt,porder)
x=x.T.flatten()
y=y.T.flatten()
z=z.T.flatten()
B = boxcar(x-x0,Wx,wx)*boxcar(y-y0,Wy,wy)
a = a*numpy.ones(x.shape)+da*(1.-B)

f = open(fname+'.rupt.p%d'%porder, 'wb')
f.write(numpy.int32(Rupt.Nele))
data=numpy.arange(Rupt.Nele)+1
f.write(data.astype(numpy.int32))
f.write(x.astype(numpy.float64))
f.write(y.astype(numpy.float64))
f.write(z.astype(numpy.float64))
f.write(a.astype(numpy.float64))
data=b*numpy.ones(x.shape)
f.write(data.astype(numpy.float64))
data=L*numpy.ones(x.shape)
f.write(data.astype(numpy.float64))
data=v0*numpy.ones(x.shape)
f.write(data.astype(numpy.float64))
data=f0*numpy.ones(x.shape)
f.write(data.astype(numpy.float64))
data=dVt0*numpy.ones(x.shape)
f.write(data.astype(numpy.float64))
f.close()

subelemesh=DG_util.subele_mesh(porder,Rupt.Nele,dim=2)

# Visualize

Tet.nod[1]=Tet.nod[1]-6.0
Rupt_tet.nod[1]=Rupt_tet.nod[1]-6.0
TetBnd.nod[1]=TetBnd.nod[1]-6.0
Rupt.nod[1]=Rupt.nod[1]-6.0
y=y-6.0

Tet.output_vtu(fname)
Rupt_tet.output_vtu(fname+'Rupt_tet')
TetBnd.output_vtu(fname+'TetBnd')
Rupt.output_vtu(fname+'.face')

fh=vis_util.vtuxml_head(fname+'_p%d.rupt'%porder,x.shape[0],subelemesh.shape[1],3,3,PointScals=['a'])
vis_util.vtuxml_point(fh,numpy.array([x,y,z]))
vis_util.vtuxml_cell(fh,subelemesh)
vis_util.vtuxml_var(fh,a)
vis_util.vtuxml_end(fh)


#print subelemesh

#plt.triplot(x, y, subelemesh.T-1)
#plt.plot(x, y, 'o')
#plt.show()

