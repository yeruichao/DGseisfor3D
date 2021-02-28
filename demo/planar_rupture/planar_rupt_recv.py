import numpy
import sys
sys.path.insert(0,'../../')
import meshtools.draw_util as draw_util
import meshtools.trimesh as trimesh
import meshtools.mesh_util as mesh_util

recvx=numpy.arange(-4.0,5.0,2.0)
recvy=recvx-6.0
recvx,recvy=numpy.meshgrid(recvx,recvy,sparse=False)
recvx=recvx.flatten()
recvy=recvy.flatten()
rad=0.2
tri=draw_util.sphere(rad=rad,cnp=numpy.array([recvx[0],recvy[0],0.0]))
for i in range(1,recvx.shape[0]):
    tri1=draw_util.sphere(rad=rad,
            cnp=numpy.array([recvx[i],recvy[i],0.0]))
    tri=tri.join(tri1)
tri.output_vtu('rupt_recvs')
