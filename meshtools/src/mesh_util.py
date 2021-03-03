import numpy

def pickeles(ele,eid=[]):
    if eid!=[]:
        ele=ele[:,eid]
    Nnod=ele.max()
    flag=numpy.zeros(Nnod,dtype=int)
    flag[ele.flatten()-1]=1
    nid=flag.nonzero()[0]
    Nnod=len(nid)
    flag[nid]=numpy.arange(Nnod)+1 
    ele=flag[ele-1]
    return ele,nid

def vtx_merge(ele,map1,map2):
    Nnod=ele.max()
    map0=numpy.arange(Nnod)+1
    map0[map1-1]=map2
    ele=map0[ele-1]
    return ele

def ele_center(nod,ele):
    Nvtx,Nele=ele.shape
    Ndim=nod.shape[0]
    ctp=numpy.zeros((Nele,Ndim))
    nod=nod.T
    for i in range(Nvtx):
        ctp = ctp + nod[ele[i]-1]
    ctp = ctp / float(Nvtx)
    return ctp.T

def det2(a11,a12,a21,a22):
    return a11*a22-a12*a21

def det3(a11,a12,a13,a21,a22,a23,a31,a32,a33):
    return a11*a22*a33+a12*a23*a31+a21*a32*a13 \
          -a13*a22*a31-a12*a21*a33-a11*a23*a32

def det4(a11,a12,a13,a14, \
         a21,a22,a23,a24, \
         a31,a32,a33,a34, \
         a41,a42,a43,a44):
    return a11*det3(a22,a23,a24,a32,a33,a34,a42,a43,a44) \
          -a12*det3(a21,a23,a24,a31,a33,a34,a41,a43,a44) \
          +a13*det3(a21,a22,a24,a31,a32,a34,a41,a42,a44) \
          -a14*det3(a21,a22,a23,a31,a32,a33,a41,a42,a43)

def L2norm(a):
    Ndim,Nele=a.shape
    b=numpy.zeros(Nele)
    for i in range(Ndim):
        b=b+a[i]*a[i]
    return numpy.sqrt(b)

def dotprod(a,b):
    Ndim,Nele=a.shape
    c=numpy.zeros(Nele)
    for i in range(Ndim):
        c=c+a[i]*b[i]
    return c

def cart2sph(nod,center=[]):
    if center!=[]:
        nod=nod-numpy.tile(center,(nod.shape[1],1)).T
    azi = numpy.arctan2(nod[1],nod[0])
    elv = numpy.arctan2(nod[2],numpy.sqrt(nod[1]**2 + nod[0]**2))
    r = L2norm(nod)
    return azi, elv, r

def sph2cart(azi,elv,r,center=[]):
    x = r * numpy.cos(elv) * numpy.cos(azi)
    y = r * numpy.cos(elv) * numpy.sin(azi)
    z = r * numpy.sin(elv)
    if center==[]:
        return numpy.vstack((x,y,z))
    else:
        return numpy.vstack((x,y,z))\
                + numpy.tile(center,(x.size,1)).T

def mesh_join(mesh1,mesh2):
    nod=numpy.hstack((mesh1.nod,mesh2.nod))
    ele=numpy.hstack((mesh1.ele,mesh2.ele+mesh1.Nnod))
    if mesh1.att!=[] and mesh2.att!=[]:
        att=numpy.hstack((mesh1.att,mesh2.att))
    else:
        att=[]
    return nod,ele,att
