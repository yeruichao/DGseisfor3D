import numpy
from .mesh_util import pickeles,ele_center,det3,L2norm,mesh_join
from .vis_util import vtuxml_head,vtuxml_point,vtuxml_cell,\
        vtuxml_var,vtuxml_end

def tri_facenormal(nod,ele):
    # computes outer normal direction of triangular mesh
    # normal direction follows right-hand-rule
    # nod: vertices coordinates, dim[3,Nnod]
    # ele: element to vertices, dim[3,Nele]
    # return: unit normal vectors, dim[3,Nele]
    a1=nod.T[ele[1]-1]-nod.T[ele[0]-1]
    a2=nod.T[ele[2]-1]-nod.T[ele[0]-1]
    n=numpy.cross(a1,a2).T
    a=L2norm(n)
    return n/numpy.tile(a,(3,1))

def tri2edge(ele):
    # obtain edge information of a triangular mesh
    # ele: element to vertices, dim[3,Nele]
    # edge: edge to vertices, dim[2,Nedge], small vertex id ahead
    # T2E: triangle element to edge, dim[3,Nele], 
    #      edges in the order [[2,3],[1,3],[1,2]]
    # E2T: edge to triangle element, dim[2,Nedge], small ele id ahead
    #      if and edge is on boundary, the first entry is -1
    # neigh: triangle to triangle, dim[3,Nele]
    # E2V: edge to opposite vetices that are belong to a same triangle,
    #      corresponding to E2T, dim[2,Nedge]
    Nele=ele.shape[1]
    l=numpy.array([2,3,1,3,1,2],dtype=int)-1
    edge=ele[l].T.reshape((Nele*3,2)).T
    tid=numpy.arange(Nele)+1
    tid=numpy.tile(tid,(3,1)).T.flatten()
    ii=edge[0]<edge[1]
    tmp=edge[0][ii]
    edge[0][ii]=edge[1][ii]
    edge[1][ii]=tmp
    edge,T2E = numpy.unique(edge,return_inverse=True,axis=1)
    Nedge=edge.shape[1]
    Tid=T2E.argsort()
    eid=T2E[Tid]
    Vid=ele.T.flatten()[Tid]
    Tid=tid[Tid]
    Eid,tmp=numpy.unique(eid,return_index=True)
    E2T=-numpy.ones((2,Nedge),dtype=int)
    E2V=-numpy.ones((2,Nedge),dtype=int)
    E2T[0][Eid]=Tid[tmp]
    E2V[0][Eid]=Vid[tmp]
    ID=numpy.ones(Nele*3, dtype=bool)
    ID[tmp]=False
    E2T[1][eid[ID]]=Tid[ID]
    E2V[1][eid[ID]]=Vid[ID]
    T2E.shape=(Nele,3)
    T2E=T2E.T
    tid=numpy.arange(Nele)+1
    neigh=E2T.T[T2E.flatten()].sum(axis=1) \
            -numpy.hstack((tid,tid,tid))
    neigh.shape=(3,Nele)
    ID=E2T[0]>E2T[1]
    tmp=E2T[0][ID]
    E2T[0][ID]=E2T[1][ID]
    E2T[1][ID]=tmp
    tmp=E2V[0][ID]
    E2V[0][ID]=E2V[1][ID]
    E2V[1][ID]=tmp
    T2E=T2E+1
    return edge,T2E,E2T,neigh,E2V

class trimesh:

    def __init__(self,nod,ele,nei=[],att=[],elemap=False):
        self.nod=nod
        self.ele=ele
        self.nei=nei
        self.att=att
        self.Nnod=nod.shape[1]
        self.Nele=ele.shape[1]
        if elemap:
            edge,T2E,E2T,neigh,E2V=tri2edge(ele)
            self.edge=edge
            self.T2E=T2E
            self.E2T=E2T
            self.E2V=E2V
            if nei==[]:
                self.nei=neigh

    def pickeles(self,eid=[]):
        ele,nid=pickeles(self.ele.copy(),eid=eid)
        nod=self.nod[:,nid]
        if eid!=[] and self.att!=[]:
            att=self.att[eid]
        else:
            att=self.att
        return trimesh(nod=nod,ele=ele,att=att)

    def areas(self):
        if self.nod.shape[0] == 2:
            tmp=numpy.ones(self.Nele)
            return det3(\
                self.nod[0][self.ele[0]-1],\
                self.nod[0][self.ele[1]-1],\
                self.nod[0][self.ele[2]-1],\
                self.nod[1][self.ele[0]-1],\
                self.nod[1][self.ele[1]-1],\
                self.nod[1][self.ele[2]-1],\
                tmp,tmp,tmp)
        else:
            p=self.nod.T
            a1=p[self.ele[0]-1]-p[self.ele[1]-1]
            a2=p[self.ele[0]-1]-p[self.ele[2]-1]
            a=numpy.cross(a1,a2).T
            return L2norm(a)

    def centers(self):
        return ele_center(self.nod,self.ele)

    def output_vtu(self,filename):
        if self.att != []:
            fh=vtuxml_head(filename,self.Nnod,self.Nele,\
                self.nod.shape[0],self.ele.shape[0],CellScals=['att'])
        else:
            fh=vtuxml_head(filename,self.Nnod,self.Nele,\
                self.nod.shape[0],self.ele.shape[0])
        vtuxml_point(fh,self.nod)
        vtuxml_cell(fh,self.ele)
        if self.att != []:
            vtuxml_var(fh,self.att)
        vtuxml_end(fh)

    def facenormal(self):
        nod=self.nod
        ele=self.ele
        if nod.shape[0]==2:
            nod=numpy.vstack((nod,numpy.zeros(self.Nnod)))
        return tri_facenormal(nod,ele)

    def join(self,tri):
        nod,ele,att=mesh_join(self,tri)
        return trimesh(nod=nod,ele=ele,att=att)

def read_triangle_file(basename,elemap=False,dim=2):

    with open(basename+'.ele') as f:
        Nele,Ndim,Natt = [int(x) for x in next(f).split()]
        a = []
        for iele in range(Nele):
            a.append([int(x) for x in next(f).split()])
    ele=numpy.array(a,dtype=int).reshape((Nele,4+Natt)).T
    att=[]
    if Natt > 0:
        att=ele[4:]
    ele=ele[1:4]

    with open(basename+'.node') as f:
        Nnod,Ndim,Natt,Kbnd = [int(x) for x in next(f).split()]
        a = []
        for iele in range(Nele):
            a.append([x for x in next(f).split()])
    nod=numpy.array(a,dtype=numpy.float32).\
            reshape((Nele,3+Natt+Kbnd)).T[1:3]

    nei=[]
    if basename+'.neigh' in locals():
        with open(basename+'.neigh') as f:
            Nele,Nnei = [int(x) for x in next(f).split()]
            a = []
            for iele in range(Nele):
                a.append([int(x) for x in next(f).split()])
        nei=numpy.array(a,dtype=int).reshape((Nele,4)).T[1:4]

    if dim==3:
        nod=numpy.vstack((nod,numpy.zeros(Nnod)))

    tri=trimesh(nod,ele,nei,att,elemap)
    return tri

