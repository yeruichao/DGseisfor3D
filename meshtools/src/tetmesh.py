import numpy
from .trimesh import trimesh,tri_facenormal,tri2edge
from .mesh_util import pickeles,ele_center,det4,dotprod,L2norm,\
        sph2cart,cart2sph,mesh_join
from .vis_util import vtuxml_head,vtuxml_point,vtuxml_cell,\
        vtuxml_var,vtuxml_end

def tet2face(tet):
    Nele=tet.shape[1]
    l=numpy.array([2,3,4,1,4,3,1,2,4,1,3,2],dtype=int)-1
    face=tet[l].T.reshape((Nele*4,3)).T
    tid=numpy.arange(Nele)+1
    tid=numpy.tile(tid,(4,1)).T.flatten()
    ii=face[0]>face[1]
    tmp=face[0][ii]
    face[0][ii]=face[1][ii]
    face[1][ii]=tmp
    ii=face[1]>face[2]
    tmp=face[1][ii]
    face[1][ii]=face[2][ii]
    face[2][ii]=tmp
    ii=face[0]>face[1]
    tmp=face[0][ii]
    face[0][ii]=face[1][ii]
    face[1][ii]=tmp
    face,T2F = numpy.unique(face,return_inverse=True,axis=1)
    Nface=face.shape[1]
    Tid=T2F.argsort()
    fid=T2F[Tid]
    Vid=tet.T.flatten()[Tid]
    Tid=tid[Tid]
    Fid,tmp=numpy.unique(fid,return_index=True)
    F2T=-numpy.ones((2,Nface),dtype=int)
    F2V=-numpy.ones((2,Nface),dtype=int)
    F2T[0][Fid]=Tid[tmp]
    F2V[0][Fid]=Vid[tmp]
    ID=numpy.ones(Nele*4, dtype=bool)
    ID[tmp]=False
    F2T[1][fid[ID]]=Tid[ID]
    F2V[1][fid[ID]]=Vid[ID]
    T2F.shape=(Nele,4)
    T2F=T2F.T
    tid=numpy.arange(Nele)+1
    neigh=F2T.T[T2F.flatten()].sum(axis=1) \
            -numpy.hstack((tid,tid,tid,tid))
    neigh.shape=(4,Nele)
    ID=F2T[0]>F2T[1]
    tmp=F2T[0][ID]
    F2T[0][ID]=F2T[1][ID]
    F2T[1][ID]=tmp
    tmp=F2V[0][ID]
    F2V[0][ID]=F2V[1][ID]
    F2V[1][ID]=tmp
    T2F=T2F+1
    return face,T2F,F2T,neigh,F2V

def tet2edge(tet):
    Nele=tet.shape[1]
    l=numpy.array([1,2,1,3,1,4,2,3,2,4,3,4],dtype=int)-1
    edge=tet[l].T.reshape((Nele*6,2)).T
    ii=edge[0]>edge[1]
    tmp=edge[0][ii]
    edge[0][ii]=edge[1][ii]
    edge[1][ii]=tmp
    edge,T2E = numpy.unique(edge,return_inverse=True,axis=1)
    T2E.shape=(Nele,6)
    T2E=T2E.T+1
    return edge,T2E

class tetmesh:

    def __init__(self,nod,ele,nei=[],att=[],elemap=False):
        self.nod=nod
        self.ele=ele
        self.Nnod=nod.shape[1]
        if self.ele == []: 
            self.Nele=0
        else:
            self.Nele=ele.shape[1]
        if self.Nele==0 or att==[]:
            self.att=[]
        elif type(att) is numpy.ndarray:
            self.att=att
        else:
            self.att=numpy.ones(self.Nele)*att
        self.nei=nei
        if elemap:
            face,T2F,F2T,neigh,F2V=tet2face(ele)
            self.face=face
            self.T2F=T2F
            self.F2T=F2T
            self.F2V=F2V
            if nei==[]:
                self.nei=neigh

    def elemap(self):
        face,T2F,F2T,neigh,F2V=tet2face(self.ele)
        self.face=face
        self.T2F=T2F
        self.F2T=F2T
        self.F2V=F2V
        self.nei=neigh
        return self

    def sample(self,ID):
        ele,nid=pickeles(self.ele,ID)
        nod=self.nod.T[nid].T
        at = []
        if self.att != []:
            at = self.att[ID]
        tet=tetmesh(nod,ele,att=at)
        return tet

    def bnd_trimesh(self):
        if not hasattr(self, 'F2T'):
            face,T2F,F2T,neigh,F2V=tet2face(self.ele)
            self.F2T=F2T
            self.face=face
            self.F2V=F2V
        ID=self.F2T[0]<0
        ele2tetnod=self.face[:,ID]
        F2V=self.F2V[:,ID]
        F2T=self.F2T[:,ID]
        tet1=tetmesh(self.nod,numpy.vstack((ele2tetnod,F2V[1])))
        ID1=tet1.volumes()>0.0
        del tet1
        ele,nid=pickeles(self.face,ID)
        tmp=ele[[1,0],:];ele[:2,ID1]=tmp[:,ID1]
        tmp=ele2tetnod[[1,0],:];ele2tetnod[:2,ID1]=tmp[:,ID1]
        nod=self.nod[:,nid]
        tri=trimesh(nod,ele)
        tri.ele2tetnod=ele2tetnod
        tri.nid2tetnod=nid+1
        tri.F2V=F2V
        tri.F2V=F2T
        return tri

    def int_trimesh(self,bnd=False):
        if self.att==[]:
            print ('No interface found')
            return []
        if not hasattr(self, 'F2T'):
            face,T2F,F2T,neigh,F2V=tet2face(self.ele)
            self.F2T=F2T
            self.F2V=F2V
            self.face=face
        Nface=self.face.shape[1]
        fAt=-numpy.ones((2,Nface),dtype=int)
        ID=self.F2T>0
        fAt[ID]=self.att[self.F2T[ID]-1]-self.att.min()+1
        ID = fAt[0]!=fAt[1] 
        if not bnd:
            ID = numpy.logical_and(ID,fAt[0]>=0)
        ele,nid=pickeles(self.face,eid=ID)
        ele2tetnod=self.face[:,ID]
        nod=self.nod[:,nid]
        F2V=self.F2V[:,ID]
        F2T=self.F2T[:,ID]
        tet1=tetmesh(self.nod,numpy.vstack((ele2tetnod,F2V[1])))
        ID1=tet1.volumes()>0.0
        del tet1
        tmp=ele[[1,0],:];ele[:2,ID1]=tmp[:,ID1]
        tmp=ele2tetnod[[1,0],:];ele2tetnod[:2,ID1]=tmp[:,ID1]
        tri=trimesh(nod,ele)
        tri.ele2tetnod=ele2tetnod
        tri.nid2tetnod=nid+1
        tri.F2V=F2V
        tri.F2T=F2T
        return tri

    def pickeles(self,eid=[]):
        ele,nid=pickeles(self.ele.copy(),eid=eid)
        nod=self.nod[:,nid]
        if eid!=[] and self.att!=[]:
            att=self.att[eid]
        else:
            att=self.att
        return tetmesh(nod=nod,ele=ele,att=att)

    def centers(self):
        return ele_center(self.nod,self.ele)

    def volumes(self):
        tmp=numpy.ones(self.Nele)
        return -det4(\
                self.nod[0][self.ele[0]-1],\
                self.nod[0][self.ele[1]-1],\
                self.nod[0][self.ele[2]-1],\
                self.nod[0][self.ele[3]-1],\
                self.nod[1][self.ele[0]-1],\
                self.nod[1][self.ele[1]-1],\
                self.nod[1][self.ele[2]-1],\
                self.nod[1][self.ele[3]-1],\
                self.nod[2][self.ele[0]-1],\
                self.nod[2][self.ele[1]-1],\
                self.nod[2][self.ele[2]-1],\
                self.nod[2][self.ele[3]-1],\
                tmp,tmp,tmp,tmp)

    def inscrib_radius(self,vol=[]):
        if vol == []:
            vol=self.volumes()
        a1=trimesh(self.nod,self.ele[[0,2,1]]).areas()
        a2=trimesh(self.nod,self.ele[[0,1,3]]).areas()
        a3=trimesh(self.nod,self.ele[[0,3,2]]).areas()
        a4=trimesh(self.nod,self.ele[[1,2,3]]).areas()
        return vol/(a1+a2+a3+a4)/3.0

    def subscrib_radius(self,vol=[]):
        if vol == []:
            vol=self.volumes()
        a=self.nod.T[self.ele[1]-1]-self.nod.T[self.ele[0]-1]
        b=self.nod.T[self.ele[2]-1]-self.nod.T[self.ele[0]-1]
        c=self.nod.T[self.ele[3]-1]-self.nod.T[self.ele[0]-1]
        r = numpy.tile(dotprod(a.T,a.T),(3,1))\
                *numpy.cross(b,c).T + \
            numpy.tile(dotprod(b.T,b.T),(3,1))\
                *numpy.cross(c,a).T + \
            numpy.tile(dotprod(c.T,c.T),(3,1))\
                *numpy.cross(a,b).T 
        return L2norm(r)/vol/12.0

    def quality(self,vols=[]):
        if vols==[]:
            vols=self.volumes()
        r=self.inscrib_radius(vols)
        R=self.subscrib_radius(vols)
        l=numpy.array([1,3,2],dtype=int)-1
        n1=tri_facenormal(self.nod,self.ele[l])
        l=numpy.array([1,2,4],dtype=int)-1
        n2=tri_facenormal(self.nod,self.ele[l])
        l=numpy.array([1,4,3],dtype=int)-1
        n3=tri_facenormal(self.nod,self.ele[l])
        l=numpy.array([2,3,4],dtype=int)-1
        n4=tri_facenormal(self.nod,self.ele[l])
        th=numpy.pi-numpy.vstack(( \
                numpy.arccos(dotprod(n1,n2)),\
                numpy.arccos(dotprod(n1,n2)),\
                numpy.arccos(dotprod(n1,n2)),\
                numpy.arccos(dotprod(n1,n2)),\
                numpy.arccos(dotprod(n1,n2)),\
                numpy.arccos(dotprod(n1,n2)) ))
        ang=th.min(axis=1)
        Ang=th.max(axis=1)
        k=R/r
        print ("ele number = ", self.Nele               )
        print ("vtx number = ", self.Nnod               )
        print ("max volume = ", vols.max()              )
        print ("min volume = ", vols.min()              )
        print ("max insc r = ", r.max()                 )
        print ("min insc r = ", r.min()                 )
        print ("max R/r    = ", k.max()                 )
        print ("min R/r    = ", k.min()                 )
        print ("max Angle  = ", Ang.max()/numpy.pi*180.0)
        print ("min Angle  = ", ang.min()/numpy.pi*180.0)

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

    def output_tetgen(self,filename):
        nod=self.nod.T
        f = open(filename+'.nod', 'w')
        f.write('%d\t%d\t%d\t%d\n'%(self.Nnod,3,0,0))
        for i in range(self.Nnod):
            f.write('%d\t%f\t%f\t%f\n'%((i+1,)+tuple(nod[i])))
        f.write('# Generated by tetmesh.py')
        f.close()
        del nod
        ele=self.ele.T
        f = open(filename+'.ele', 'w')
        if self.att!=[]:
            f.write('%d\t%d\t%d\n'%(self.Nele,4,1))
            att=self.att-self.att.min()
            for i in range(self.Nele):
                f.write('%d\t%d\t%d\t%d\t%d\t%d\n'
                        %((i+1,)+tuple(ele[i])+(att[i],)))
        else:
            f.write('%d\t%d\t%d\n'%(self.Nele,4,0))
            for i in range(self.Nele):
                f.write('%d\t%d\t%d\t%d\t%d\n'
                        %((i+1,)+tuple(ele[i])))
        f.write('# Generated by tetmesh.py')
        f.close()
        del ele
        if self.nei!=[]:
            nei=self.nei.T
            f = open(filename+'.neigh', 'w')
            f.write('%d\t%d\n'%(self.Nele,4))
            for i in range(self.Nele):
                f.write('%d\t%d\t%d\t%d\t%d\n'
                        %((i+1,)+tuple(nei[i])))
            f.write('# Generated by tetmesh.py')
            f.close()

    def all_refine(self):
        T1=numpy.array([ \
             1,5,6,7,   5,2,8,9,   3,6,8,10,   7,4,9,10,\
             6,9,7,5,   6,9,5,8,   6,9,8,10,   6,9,10,7,\
             ],dtype=int)
        T2=numpy.array([ \
             1,5,6,7,   5,2,8,9,   3,6,8,10,   7,4,9,10,\
             10,5,7,6,  10,5,6,8,  10,5,8,9,   10,5,9,7,\
             ],dtype=int)
        T3=numpy.array([ \
             1,5,6,7,   5,2,8,9,   3,6,8,10,   7,4,9,10,\
             7,8,5,6,   7,8,6,10,  7,8,10,9,   7,8,9,5,\
             ],dtype=int)
        edge,T2E=tet2edge(self.ele)
        nod=(self.nod.T[edge[0]-1]+self.nod.T[edge[1]-1])/2.0
        nod=numpy.hstack((self.nod,nod.T))
        Nele=self.Nele*8
        ele=numpy.vstack((self.ele,T2E+self.Nnod))
        d1=L2norm(nod.T[ele[5]-1].T-nod.T[ele[8]-1].T)
        d2=L2norm(nod.T[ele[9]-1].T-nod.T[ele[4]-1].T)
        d3=L2norm(nod.T[ele[6]-1].T-nod.T[ele[7]-1].T)
        id1 = numpy.logical_and(d1<=d2, d1<=d3)
        id2 = numpy.logical_and(d2< d1, d2<=d3)
        id3 = numpy.logical_not(numpy.logical_or(id1, id2))
        ele = numpy.hstack((ele.T[id1].T[T1-1], \
                ele.T[id2].T[T2-1], ele.T[id3].T[T3-1])\
                ).T.reshape((Nele,4)).T
        a=[]
        if self.att != []:
            a=numpy.hstack((\
                    numpy.tile(self.att[id1],(8,1)).T.flatten(),\
                    numpy.tile(self.att[id2],(8,1)).T.flatten(),\
                    numpy.tile(self.att[id3],(8,1)).T.flatten(),\
                    )).T.flatten()
        return tetmesh(nod,ele,att=a)

    def tile(self,ndh,nid=[],trisurf=[],direction='r',center=[],att=[]):
        '''
        # By default spherically tile along r direction, 
        #   with center at (0,0,0)
        # If tiling along x,y,z direction, trisurf must be provided
        #   containing the components: 
        #   ele, nod, Nele, Nnod, ele2tetnod, nid2tetnod
        '''
        if trisurf==[]:
            trisurf=self.bnd_trimesh()
        Nfnod=trisurf.Nnod
        Nfele=trisurf.Nele
        levels=triface_tile_sample(trisurf.ele,nid)
        Nlevel=levels.max()
        if nid==[]:
            nod_add=trisurf.nod.copy()
        else:
            if len(nid) != Nfnod:
                nod_add=numpy.zeros(Nfnod,dtype=bool)
                nod_add[nid]=True
                nid=nod_add
            nod_add=trisurf.nod.T[nid].T
            ndh=ndh[nid]
        if direction=='r':
            azi,elv,r=cart2sph(nod_add,center)
            nod_add=sph2cart(azi,elv,r+ndh,center)
        elif direction=='x':
            nod_add[0]=nod_add[0]+ndh
        elif direction=='y':
            nod_add[1]=nod_add[1]+ndh
        elif direction=='z':
            nod_add[2]=nod_add[2]+ndh
        if nid==[]:
            nod=nod_add
            nid2tetnod=numpy.arange(Nfnod)+1+self.Nnod
            t_top1=trisurf.ele+self.Nnod
        else:
            nod=trisurf.nod.T
            nod[nid]=nod_add.T
            nod=nod.T
            Nnod_add=nod_add.shape[1]
            nid2tetnod=trisurf.nid2tetnod
            nid2tetnod[nid]=numpy.arange(Nnod_add)+1+self.Nnod
            t_top1=nid2tetnod[trisurf.ele-1]
        t_new=[]
        t_top=trisurf.ele2tetnod
        for ilevel in range(1,Nlevel+1):
            it_add = levels==ilevel
            if nid != []:
                it_add = numpy.logical_and(it_add, nid)
            it_add = it_add[trisurf.ele-1]
            for i in range(3):
                itadd=it_add[i]
                if any(itadd):
                    t_new.append(numpy.vstack(\
                            (t_top.T[itadd].T,t_top1[i][itadd])).T)
            t_top[it_add]=t_top1[it_add]
        t_new=numpy.vstack(tuple(t_new)).T
        Nt_new=t_new.shape[1]
        if self.Nele > 0:
            t_new=numpy.hstack((self.ele,t_new))
        nod_add=numpy.hstack((self.nod,nod_add))
        trisurf1=trimesh(nod,trisurf.ele)
        trisurf1.ele2tetnod=t_top
        trisurf1.nid2tetnod=nid2tetnod
        if att!=[]:
            if self.att == []:
                att = numpy.ones(Nt_new)*att
            else:
                att=numpy.hstack((self.att,numpy.ones(Nt_new)*att))
        return tetmesh(nod_add,t_new,att=att),trisurf1

    def tile_refineall(self,ndh,trisurf=[],direction='r',center=[],\
            att=[]):
        '''
        # By default spherically tile along r direction, 
        #   with center at (0,0,0)
        # If tiling along x,y,z direction, trisurf must be provided
        #   containing the components: 
        #   ele, nod, Nele, Nnod, ele2tetnod, nid2tetnod
        # 1 ----- 2
        #  \     /
        #   \   /
        #    \ /
        #     3
        #
        #    10
        #
        # 4 --9-- 5
        #  \     /
        #   8   7
        #    \ /
        #     6
        '''
        if trisurf==[]:
            trisurf=self.bnd_trimesh()
        Nfnod=trisurf.Nnod
        Nfele=trisurf.Nele
        tt=numpy.array([ 
             1,4,9,8,     2,5,7,9,     3,6,8,7,    10,1,3,2,
            10,7,8,9,    10,1,9,8,    10,2,7,9,    10,3,8,7,
            10,1,2,9,    10,2,3,7,    10,3,1,8,
            ],dtype=int)
        ttt=numpy.array([ 
            4,9,8,   9,5,7,   8,7,6,   8,9,7,
            ],dtype=int)
        nod_add=trisurf.nod.copy()
        if direction=='r':
            azi,elv,r=cart2sph(nod_add,center)
            nod_add=sph2cart(azi,elv,r+ndh,center)
        elif direction=='x':
            nod_add[0]=nod_add[0]+ndh
        elif direction=='y':
            nod_add[1]=nod_add[1]+ndh
        elif direction=='z':
            nod_add[2]=nod_add[2]+ndh
        nod_adt=nod_add-trisurf.nod
        edge,T2E,E2T,neigh,E2V=tri2edge(trisurf.ele)
        Nedg=edge.shape[1]
        nod_edg=(nod_add.T[edge[0]-1]+nod_add.T[edge[1]-1])/2.0
        nod_cnt=(nod_add.T[trisurf.ele[0]-1] \
               + nod_add.T[trisurf.ele[1]-1] \
               + nod_add.T[trisurf.ele[2]-1] \
               - nod_adt.T[trisurf.ele[0]-1]/2.0 \
               - nod_adt.T[trisurf.ele[1]-1]/2.0 \
               - nod_adt.T[trisurf.ele[2]-1]/2.0 \
               )/3.0
        nod_top=numpy.vstack((nod_add.T,nod_edg)).T
        nod_new=numpy.vstack((self.nod.T,nod_top.T,nod_cnt))
        t_add=numpy.vstack((trisurf.ele2tetnod,\
                trisurf.ele+self.Nnod,\
                T2E+self.Nnod+Nfnod,\
                numpy.arange(Nfele)+1+self.Nnod+Nfnod+Nedg))
        T_add=t_add[tt-1].T.reshape((11*Nfele,4))
        T_top=t_add[ttt-1].T.reshape((4*Nfele,3))
        T_ref=T_top-self.Nnod
        Nt_new=T_add.shape[0]
        if self.Nele > 0:
            t_new=numpy.hstack((self.ele,T_add.T))
        if att!=[]:
            if self.att == []:
                att = numpy.ones(Nt_new)*att
            else:
                att=numpy.hstack((self.att,numpy.ones(Nt_new)*att))
        tet=tetmesh(nod_new.T,t_new,att=att)
        tri=trimesh(nod_top,T_ref.T)
        tri.ele2tetnod=T_top.T
        tri.nid2tetnod=numpy.arange(Nfnod+Nedg)+1*self.Nnod
        return tet,tri

    def tile_isolated(self,ndh,trisurf=[],direction='r',center=[],\
            att=[]):
        '''
        # By default spherically tile along r direction, 
        #   with center at (0,0,0)
        # If tiling along x,y,z direction, trisurf must be provided
        #   containing the components: 
        #   ele, nod, Nele, Nnod, ele2tetnod, nid2tetnod
        # 1 ----- 2
        #  \     /
        #   \   /
        #    \ /
        #     3
        #
        # ____ 9 ___
        # \        /
        #  \  10  /
        #   8    7
        #    \  /
        #     \/
        #
        # 4 ----- 5
        #  \     /
        #   \   /
        #    \ /
        #     6
        '''
        if trisurf==[]:
            trisurf=self.bnd_trimesh()
        Nfnod=trisurf.Nnod
        Nfele=trisurf.Nele
        tt=numpy.array([ 
            1,2,3,10,    5,4,6,10,    2,1,9,10,    5,2,9,10,\
            4,5,9,10,    1,4,9,10,    1,3,8,10,    4,1,8,10,\
            3,6,8,10,    6,4,8,10,    3,2,7,10,    2,5,7,10,\
            6,3,7,10,    5,6,7,10,\
            ],dtype=int)
        nod_add=trisurf.nod.copy()
        if direction=='r':
            azi,elv,r=cart2sph(nod_add,center)
            nod_add=sph2cart(azi,elv,r+ndh,center)
        elif direction=='x':
            nod_add[0]=nod_add[0]+ndh
        elif direction=='y':
            nod_add[1]=nod_add[1]+ndh
        elif direction=='z':
            nod_add[2]=nod_add[2]+ndh
        nod_adt=nod_add-trisurf.nod
        edge,T2E,E2T,neigh,E2V=tri2edge(trisurf.ele)
        Nedg=edge.shape[1]
        nod_edg=(nod_add.T[edge[0]-1] \
                +nod_add.T[edge[1]-1] \
                +trisurf.nod.T[edge[0]-1] \
                +trisurf.nod.T[edge[1]-1] \
                )/4.0
        nod_cnt=(nod_add.T[trisurf.ele[0]-1] \
               + nod_add.T[trisurf.ele[1]-1] \
               + nod_add.T[trisurf.ele[2]-1] \
               - nod_adt.T[trisurf.ele[0]-1]/2.0 \
               - nod_adt.T[trisurf.ele[1]-1]/2.0 \
               - nod_adt.T[trisurf.ele[2]-1]/2.0 \
               )/3.0
        nod_new=numpy.vstack((self.nod.T,nod_add.T,nod_edg,nod_cnt))
        t_add=numpy.vstack((trisurf.ele2tetnod,\
                trisurf.ele+self.Nnod,\
                T2E+self.Nnod+Nfnod,\
                numpy.arange(Nfele)+1+self.Nnod+Nfnod+Nedg))
        T_add=t_add[tt-1].T.reshape((14*Nfele,4))
        T_top=trisurf.ele+self.Nnod
        if self.Nele > 0:
            t_new=numpy.hstack((self.ele,T_add.T))
        else:
            t_new=T_add.T
        if att!=[]:
            Nt_new=T_add.shape[0]
            if self.att == []:
                att = numpy.ones(Nt_new)*att
            else:
                att=numpy.hstack((self.att,numpy.ones(Nt_new)*att))
        else:
            if (self.att == [] and self.Nele == 0) or \
                    self.att != []:
                if trisurf.att != []:
                    att=numpy.dot(numpy.ones((14,1)),
                            numpy.array([trisurf.att])).T.flatten()
                if self.att != []:
                    att=numpy.hstack((self.att,att))
        tet=tetmesh(nod_new.T,t_new,att=att)
        trisurf.nod=nod_add
        trisurf.ele2tetnod=T_top
        trisurf.nid2tetnod=numpy.arange(Nfnod)+1*self.Nnod
        return tet,trisurf

    def join(self,tet):
        nod,ele,att=mesh_join(self,tet)
        return tetmesh(nod=nod,ele=ele,att=att)

def read_tetgen_file(basename,elemap=False):

    with open(basename+'.ele') as f:
        Nele,Ndim,Natt = [int(x) for x in next(f).split()]
        a = []
        for iele in range(Nele):
            a.append([int(x) for x in next(f).split()])
    ele=numpy.array(a,dtype=int).reshape((Nele,5+Natt)).T
    att=[]
    if Natt > 0:
        att=ele[5:]
    ele=ele[1:5]

    with open(basename+'.node') as f:
        Nnod,Ndim,Natt,Kbnd = [int(x) for x in next(f).split()]
        a = []
        for iele in range(Nele):
            a.append([x for x in next(f).split()])
    nod=numpy.array(a,dtype=numpy.float32).\
            reshape((Nele,4+Natt+Kbnd)).T[1:4]

    nei=[]
    if basename+'.neigh' in locals():
        with open(basename+'.neigh') as f:
            Nele,Nnei = [int(x) for x in next(f).split()]
            a = []
            for iele in range(Nele):
                a.append([int(x) for x in next(f).split()])
        nei=numpy.array(a,dtype=int).reshape((Nele,5)).T[1:5]

    tet=tetmesh(nod,ele,nei,att,elemap)
    return tet

def triface_tile_sample(tri,nid=[]):
    Ntri=tri.shape[1]
    Nnod=tri.max()
    flag=-numpy.ones(Nnod,dtype=int)
    if nid!=[]:
        flag=flag*2
        flag[nid]=flag[nid]+1
    level=0
    ic=numpy.random.rand(Ntri)
    ic=ic.argsort()
    while True:
        if all(flag != -1):
            break
        level=level+1
        for itri in range(Ntri):
            i=ic[itri]
            l=0
            for j in range(3):
                if flag[tri[j][i]-1]==level:
                    l=l+1
            if l == 0:
                for j in range(3):
                    inod=tri[j][i]-1
                    if flag[inod]==-1:
                        flag[inod]=level
                        for k in range(j+1,3):
                            inod=tri[k][i]-1
                            if flag[inod]==-1:
                                flag[inod]=0
                        break
            elif l==1:
                for j in range(3):
                    inod=tri[j][i]-1
                    if flag[inod]==-1:
                        flag[inod]=0
            else:
                k=True
                for j in range(3):
                    inod=tri[j][i]-1
                    if flag[inod]==-1:
                        flag[inod]=0
                    if flag[inod]==level:
                        if k:
                            k=False
                        else:
                            flag[inod]=0
        flag[flag==0]=-1 
    print ('level=',level)
    return flag


