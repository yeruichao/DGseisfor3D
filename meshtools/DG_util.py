import numpy
from scipy import interpolate
#import scipy.interpolate.interp2d as interp2d
#import scipy.interpolate.RegularGridInterpolator as interp3d
from p2n_init import p2n_init

def tetsubele(N):
    if N==1:
        tt=numpy.array([numpy.arange(4)])+1
        return tt.T
    if N==2: 
        Nsele = N*(N+1)*(N+2)/6 + (N-1)*N*(N+1)/6*4
    elif N>=3:
        Nsele = N*(N+1)*(N+2)/6 + (N-1)*N*(N+1)/6*4 + (N-2)*(N-1)*N/6
    tt=numpy.zeros((4,Nsele),dtype=int)
    M=numpy.zeros((N+1,N+1,N+1),dtype=int)
    sk=0
    for k in range(N+1):
        for j in range(N+1-k):
            for i in range(N+1-k-j):
                M[i,j,k]=sk;sk=sk+1
    sk=0
    for k in range(N):
        for j in range(N-k):
            for i in range(N-k-j):
                tt[0,sk]=M[i  ,j  ,k  ]
                tt[1,sk]=M[i+1,j  ,k  ]
                tt[2,sk]=M[i  ,j+1,k  ]
                tt[3,sk]=M[i  ,j  ,k+1]
                sk=sk+1
    for k in range(N-1):
        for j in range(N-1-k):
            for i in range(N-1-k-j):
                t1=M[i+1,j,k];t2=M[i+1,j+1,k];t3=M[i,j+1,k  ]
                t4=M[i,j,k+1];t5=M[i+1,j,k+1];t6=M[i,j+1,k+1]
                tt[0,sk]=t4;tt[1,sk]=t2;tt[2,sk]=t5;tt[3,sk]=t1
                sk=sk+1
                tt[0,sk]=t4;tt[1,sk]=t2;tt[2,sk]=t1;tt[3,sk]=t3
                sk=sk+1
                tt[0,sk]=t4;tt[1,sk]=t2;tt[2,sk]=t3;tt[3,sk]=t6
                sk=sk+1
                tt[0,sk]=t4;tt[1,sk]=t2;tt[2,sk]=t6;tt[3,sk]=t5
                sk=sk+1
    if N==2:
        return tt
    for k in range(N-2):
        for j in range(N-2-k):
            for i in range(N-2-k-j):
                tt[0,sk]=M[i,j+1,k+1  ]
                tt[1,sk]=M[i+1,j+1,k+1]
                tt[2,sk]=M[i+1,j,k+1  ]
                tt[3,sk]=M[i+1,j+1,k  ]
                sk=sk+1
    return tt+1

def trisubele(N):
    if N==1:
       tt=numpy.array([numpy.arange(3)])+1
       return tt.T
    Nsele = N*(N+1)/2 + (N-1)*N/2
    tt=numpy.zeros((3,Nsele),dtype=int)
    M=numpy.zeros((N+1,N+1),dtype=int)
    sk=0
    for k in range(N+1):
        for j in range(N+1-k):
            M[j,k]=sk;sk=sk+1
    sk=0
    for k in range(N):
        for j in range(N-k):
            tt[0,sk]=M[j,k  ]
            tt[1,sk]=M[j+1,k]
            tt[2,sk]=M[j,k+1]
            sk=sk+1
    for k in range(N-1):
        for j in range(N-1-k):
            tt[0,sk]=M[j+1,k  ]
            tt[1,sk]=M[j+1,k+1]
            tt[2,sk]=M[j,k+1  ]
            sk=sk+1
    return tt+1

def subele_mesh(porder,Nele,dim=3):
    if dim==2:
        Np=(porder+1)*(porder+2)/2
        tt=trisubele(porder)
    elif dim==3:
        Np=(porder+1)*(porder+2)*(porder+3)/6
        tt=tetsubele(porder)
    Nvtx,Nsele=tt.shape
    tt=numpy.array([tt.T.flatten()])
    ele=numpy.dot(numpy.ones((Nele,1),dtype=int),tt)\
       +numpy.dot(numpy.ones((Nvtx*Nsele,1),dtype=int),\
       numpy.array([numpy.arange(Nele)])).T*Np
    ele = ele.reshape((Nsele*Nele,Nvtx)).T
    return ele

class porder2nodal:

    def __init__(self):
        self.p2n = p2n_init()

    def reftet_nodal_coord(self,porder):
        r=self.p2n[porder-1]['r']
        s=self.p2n[porder-1]['s']
        t=self.p2n[porder-1]['t']
        return r,s,t

    def tetmesh_nodal_coord(self,nod,ele,porder):
        r,s,t=self.reftet_nodal_coord(porder)
        x = 0.5*(
            -numpy.dot((1.+r+s+t),numpy.array([nod[0][ele[0]-1]]))\
            +numpy.dot((1.+r)    ,numpy.array([nod[0][ele[1]-1]]))\
            +numpy.dot((1.+s)    ,numpy.array([nod[0][ele[2]-1]]))\
            +numpy.dot((1.+t)    ,numpy.array([nod[0][ele[3]-1]])))
        y = 0.5*(
            -numpy.dot((1.+r+s+t),numpy.array([nod[1][ele[0]-1]]))\
            +numpy.dot((1.+r)    ,numpy.array([nod[1][ele[1]-1]]))\
            +numpy.dot((1.+s)    ,numpy.array([nod[1][ele[2]-1]]))\
            +numpy.dot((1.+t)    ,numpy.array([nod[1][ele[3]-1]])))
        z = 0.5*(
            -numpy.dot((1.+r+s+t),numpy.array([nod[2][ele[0]-1]]))\
            +numpy.dot((1.+r)    ,numpy.array([nod[2][ele[1]-1]]))\
            +numpy.dot((1.+s)    ,numpy.array([nod[2][ele[2]-1]]))\
            +numpy.dot((1.+t)    ,numpy.array([nod[2][ele[3]-1]])))
        return x,y,z

    def tri_nodal_coord(self,tri,porder):
        nod=tri.nod
        Nnod=tri.Nnod
        ele=tri.ele
        Nele=tri.Nele
        if nod.shape[0]==2:
            nod=numpy.vstack((nod,numpy.zeros(Nnod)))
        vtx=nod[:,ele[0]-1]+tri.facenormal()
        ele=numpy.vstack((ele,numpy.arange(Nele)+Nnod+1))
        nod=numpy.hstack((nod,vtx))
        Nfp=(porder+1)*(porder+2)/2
        x,y,z=self.tetmesh_nodal_coord(nod,ele,porder)
        return x[:Nfp,:],y[:Nfp,:],z[:Nfp,:]

    def tet_nodal_coord(self,tet,porder):
        return self.tetmesh_nodal_coord(tet.nod,tet.ele,porder)

    def print_all(self):
        numpy.set_printoptions(precision=20)
        print 'def p2n_init():'
        print '    from numpy import *'
        print '    p2n=[]'
        for porder in range(1,11):
            r=self.p2n[porder-1]['r']
            s=self.p2n[porder-1]['s']
            t=self.p2n[porder-1]['t']
            x=self.p2n[porder-1]['x']
            y=self.p2n[porder-1]['y']
            z=self.p2n[porder-1]['z']
            print '    r=',numpy.array_repr(r.flatten()).replace('\n', '')
            print '    s=',numpy.array_repr(s.flatten()).replace('\n', '')
            print '    t=',numpy.array_repr(t.flatten()).replace('\n', '')
            print '    x=',numpy.array_repr(x.flatten()).replace('\n', '')
            print '    y=',numpy.array_repr(y.flatten()).replace('\n', '')
            print '    z=',numpy.array_repr(z.flatten()).replace('\n', '')
            print '    p2n.append({"r":r,"s":s,"t":t,"x":x,"y":y,"z":z})'

    def verify(self):
        x0=self.p2n['x']
        y0=self.p2n['y']
        z0=self.p2n['z']
        nod0=numpy.hstack((x0,y0,z0)).T
        ele0=numpy.array([[1],[2],[3],[4]])
        err=0
        for porder in range(1,11):
            x_ref=self.p2n[porder-1]['x']
            y_ref=self.p2n[porder-1]['y']
            z_ref=self.p2n[porder-1]['z']
            x,y,z=self.tetmesh_nodal_coord(nod0,ele0,porder)
            err=err\
                +numpy.linalg.norm(x_ref-x)**2\
                +numpy.linalg.norm(y_ref-y)**2\
                +numpy.linalg.norm(z_ref-z)**2
        print numpy.sqrt(err)

class cartesian_map:

    def __init__(self,data,d0=1.0,d1=1.0,d2=1.0,
            o0=0.0,o1=0.0,o2=0.0):
        # data coordinates are arranged as [x1,x0] or [x1,x0,x2]
        self.data = data
        self.dim=len(data.shape)
        if self.dim==2:
            n0,n1=data.shape
            x0=numpy.arange(n0)*d0+o0
            x1=numpy.arange(n1)*d1+o1
            self.intpf = interpolate.RegularGridInterpolator(
                    (x0,x1),data)
            self.xm=[x0.min(),x0.max()]
            self.ym=[x1.min(),x1.max()]
        elif self.dim==3:
            n0,n1,n2=data.shape
            x0=numpy.arange(n0)*d0+o0
            x1=numpy.arange(n1)*d1+o1
            x2=numpy.arange(n2)*d2+o2
            self.intpf = interpolate.RegularGridInterpolator(
                    (x0,x1,x2),data)
            self.xm=[x0.min(),x0.max()]
            self.ym=[x1.min(),x1.max()]
            self.zm=[x2.min(),x2.max()]

    def intp(self,p):
        if self.dim == 2:
            if p.shape[0] == 2:
                print p[0,:],p[1,:]
                return self.intpf(p.T)
            else:
                return self.intpf(p)
        elif self.dim == 3:
            if p.shape[0] == 3:
                return self.intpf(p.T)
            else:
                return self.intpf(p)
