!*******************************************************************!
!*  This module build up the entire mesh geometry and connections. *!
!*                                                                 *!
!*  Author: Ruichao Ye    Email: Ruichao.ye@gmail.com              *!
!*******************************************************************!

!--------------------------------------------------------------------
module geometry_mod
!--------------------------------------------------------------------

    use para_mod,     only : def_pOrder,pNp,Nfp,&
                             verbose,logfid,withrupt
    use datatype_mod, only : rkind,PI,&
                             tetmesh_geometry,matrices,&
                             pointer_vec_int,ini_vector
    use jacobi_mod,   only : JacobiGL,Vandermonde3D,&
                             GradVandermonde3D,Vandermonde2D,&
                             GradVandermonde2D,Vandermonde1D
    use meshfile_mod, only : xmax, xmin, ymax, ymin, zmax, zmin,&
                             glob_Nface, triface, F2T, Fatt

    implicit none
    public :: Build_geometry,rupt_perm,surf_perm,rupt_perm_info

    real(kind=rkind),parameter :: TOL=1d-10
    integer,allocatable :: rupt_perm(:,:,:),surf_perm(:,:)
    integer,allocatable :: rupt_perm_info(:,:)

!--------------------------------------------------------------------
    contains
!--------------------------------------------------------------------

subroutine Build_geometry(Nnode,Nele,nodes,tet,neigh, &
        pOrder,pNp,Nfp,Nsele,mesh,matrix)

!    input:
    integer     :: Nnode, Nele
    ! number of nodes and number of elements
    integer     :: tet(4,Nele)
    ! tetrahedral-vertices map
    integer     :: neigh(4,Nele)
    ! tetrahedral-tetrahedral map
    real(kind=rkind)    :: nodes(3,Nnode)
    ! vertices coordinates
    integer         :: pOrder,pNp,Nfp,Nsele
    ! predefined element order
!       output:
    type(tetmesh_geometry)  :: mesh
    ! tetrahedral elements geometry
    type(matrices)      :: matrix
    ! discreted matrices
!       pointers:
    integer,pointer :: vmapM(:),vmapP(:)
    real(kind=rkind),pointer :: x(:),y(:),z(:)
    real(kind=rkind),pointer :: nx(:,:),ny(:,:),nz(:,:)
    real(kind=rkind),pointer :: Jac(:,:),invJac(:,:)
    real(kind=rkind),pointer :: detJ(:)
    real(kind=rkind),pointer :: Fscale(:,:)
!       auxilary:
    integer :: i
    integer :: j,k,l,order,face,pMp,Mfp,&
           head,tail,head1,tail1,tetface(3)
    real(kind=rkind) :: a(4),b(4),c(4)
    real(kind=rkind) :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
    real(kind=rkind) :: r(pNp),s(pNp),t(pNp)
    integer,pointer :: Fmask(:,:)
    integer,pointer :: Fperm(:)
    integer,allocatable :: subtet(:,:)
    integer,allocatable :: VtoE(:),EtoF(:,:)
    integer,parameter :: FtoV(4,3)=reshape((/1,1,2,1,3,2,3,4,2,4,4,3/),shape(FtoV))
    ! Face 1 -> Nodes 1 3 2
    ! Face 2 -> Nodes 1 2 4
    ! Face 3 -> Nodes 2 3 4
    ! Face 4 -> Nodes 1 4 3


    allocate(matrix%Fmask(Nfp,4));Fmask=>matrix%Fmask
    call ini_vector(pNp*Nele,mesh%coord)
    x=>mesh%coord%x
    y=>mesh%coord%y
    z=>mesh%coord%z
    allocate(mesh%Jac   (    Nele,9));     Jac=>mesh%Jac
    allocate(mesh%invJac(    Nele,9));  invJac=>mesh%invJac
    allocate(mesh%detJ  (    Nele  ));    detJ=>mesh%detJ
    allocate(mesh%Fscale(  4,Nele  ));  Fscale=>mesh%Fscale
    allocate(mesh%nx    (  4,Nele  ));      nx=>mesh%nx
    allocate(mesh%ny    (  4,Nele  ));      ny=>mesh%ny
    allocate(mesh%nz    (  4,Nele  ));      nz=>mesh%nz
    allocate(mesh%vmapM (Nfp*4)     );   vmapM=>mesh%vmapM
    allocate(mesh%vmapP (Nfp*Nele*4));   vmapP=>mesh%vmapP
    allocate(VtoE(Nele*pNp))
    allocate(EtoF(4,Nele))

    mesh%xmax=xmax;mesh%xmin=xmin
    mesh%ymax=ymax;mesh%ymin=ymin
    mesh%zmax=zmax;mesh%zmin=zmin

    call form_ref_matrices(pOrder,pNp,Nfp,r,s,t,Fmask,matrix)
    EtoF=0
    do i=1,Nele,1
        do j=1,pNp
            x((i-1)*pNp+j)=0.5d0*(&
            -(1.0d0+r(j)+s(j)+t(j))*nodes(1,tet(1,i))&
                  +(1.0d0+r(j))*nodes(1,tet(2,i))&
                  +(1.0d0+s(j))*nodes(1,tet(3,i))&
                  +(1.0d0+t(j))*nodes(1,tet(4,i)))
                    y((i-1)*pNp+j)=0.5d0*(&
            -(1.0d0+r(j)+s(j)+t(j))*nodes(2,tet(1,i))&
                  +(1.0d0+r(j))*nodes(2,tet(2,i))&
                  +(1.0d0+s(j))*nodes(2,tet(3,i))&
                  +(1.0d0+t(j))*nodes(2,tet(4,i)))
                    z((i-1)*pNp+j)=0.5d0*(&
            -(1.0d0+r(j)+s(j)+t(j))*nodes(3,tet(1,i))&
                  +(1.0d0+r(j))*nodes(3,tet(2,i))&
                  +(1.0d0+s(j))*nodes(3,tet(3,i))&
                  +(1.0d0+t(j))*nodes(3,tet(4,i)))
            VtoE((i-1)*pNp+j)=i
        enddo
! dx/dr=[xr,yr,zr;xs,ys,zs;xt,yt,zt]
! dr/dx=[rx;sx;tx,ry;sy;ty,rz;sz;tz]
! Jac   =[xr,xs,xt,yr,ys,yt,zr,zs,zt]
! invJac=[rx,sx,tx,ry,sy,ty,rz,sz,tz]
        Jac(i,1:3)=(nodes(1,tet(2:4,i))-nodes(1,tet(1,i)))/2d0
        Jac(i,4:6)=(nodes(2,tet(2:4,i))-nodes(2,tet(1,i)))/2d0
        Jac(i,7:9)=(nodes(3,tet(2:4,i))-nodes(3,tet(1,i)))/2d0
        detJ(i)=Jac(i,1)*Jac(i,5)*Jac(i,9)&
               +Jac(i,3)*Jac(i,4)*Jac(i,8)&
               +Jac(i,2)*Jac(i,6)*Jac(i,7)&
               -Jac(i,3)*Jac(i,5)*Jac(i,7)&
               -Jac(i,2)*Jac(i,4)*Jac(i,9)&
               -Jac(i,1)*Jac(i,6)*Jac(i,8)

        x1=nodes(1,tet(1,i));y1=nodes(2,tet(1,i));z1=nodes(3,tet(1,i))
        x2=nodes(1,tet(2,i));y2=nodes(2,tet(2,i));z2=nodes(3,tet(2,i))
        x3=nodes(1,tet(3,i));y3=nodes(2,tet(3,i));z3=nodes(3,tet(3,i))
        x4=nodes(1,tet(4,i));y4=nodes(2,tet(4,i));z4=nodes(3,tet(4,i))

        invJac(i,1)=(y1*(z3-z4)+y3*(z4-z1)+y4*(z1-z3))/4.0d0/detJ(i)
        invJac(i,2)=(y1*(z4-z2)+y2*(z1-z4)+y4*(z2-z1))/4.0d0/detJ(i)
        invJac(i,3)=(y1*(z2-z3)+y2*(z3-z1)+y3*(z1-z2))/4.0d0/detJ(i)
        invJac(i,4)=(x1*(z4-z3)+x3*(z1-z4)+x4*(z3-z1))/4.0d0/detJ(i)
        invJac(i,5)=(x1*(z2-z4)+x2*(z4-z1)+x4*(z1-z2))/4.0d0/detJ(i)
        invJac(i,6)=(x1*(z3-z2)+x2*(z1-z3)+x3*(z2-z1))/4.0d0/detJ(i)
        invJac(i,7)=(x1*(y3-y4)+x3*(y4-y1)+x4*(y1-y3))/4.0d0/detJ(i)
        invJac(i,8)=(x1*(y4-y2)+x2*(y1-y4)+x4*(y2-y1))/4.0d0/detJ(i)
        invJac(i,9)=(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/4.0d0/detJ(i)

        a=(nodes(2,tet(FtoV(:,2),i))-nodes(2,tet(FtoV(:,1),i)))&
         *(nodes(3,tet(FtoV(:,3),i))-nodes(3,tet(FtoV(:,1),i)))&
         -(nodes(2,tet(FtoV(:,3),i))-nodes(2,tet(FtoV(:,1),i)))&
         *(nodes(3,tet(FtoV(:,2),i))-nodes(3,tet(FtoV(:,1),i)))
        b=(nodes(1,tet(FtoV(:,3),i))-nodes(1,tet(FtoV(:,1),i)))&
         *(nodes(3,tet(FtoV(:,2),i))-nodes(3,tet(FtoV(:,1),i)))&
         -(nodes(1,tet(FtoV(:,2),i))-nodes(1,tet(FtoV(:,1),i)))&
         *(nodes(3,tet(FtoV(:,3),i))-nodes(3,tet(FtoV(:,1),i)))
        c=(nodes(1,tet(FtoV(:,2),i))-nodes(1,tet(FtoV(:,1),i)))&
         *(nodes(2,tet(FtoV(:,3),i))-nodes(2,tet(FtoV(:,1),i)))&
         -(nodes(2,tet(FtoV(:,2),i))-nodes(2,tet(FtoV(:,1),i)))&
         *(nodes(1,tet(FtoV(:,3),i))-nodes(1,tet(FtoV(:,1),i)))

        Fscale(:,i)=sqrt(a**2+b**2+c**2)
        nx(:,i)=a/Fscale(:,i)
        ny(:,i)=b/Fscale(:,i)
        nz(:,i)=c/Fscale(:,i)
        do j=1,4
            if(neigh(j,i).le.0)then
                EtoF(j,i)=j
            else
                do k=1,4
                    if(neigh(k,neigh(j,i)).eq.i)then
                        EtoF(j,i)=k;exit
                    endif
                enddo
            endif
        enddo
        Fscale(:,i)=Fscale(:,i)/detJ(i)/4d0
    enddo
    call Build_maps3D(tet,neigh,EtoF,Fmask,vmapM,vmapP,Nele,&
        pOrder,pNp,Nfp)
    allocate(subtet(Nsele,4))

    deallocate(VtoE)
    deallocate(EtoF)

    allocate(mesh%fbctype(2,4,mesh%Nele));mesh%fbctype=0
    if(withrupt)then
        allocate(Fperm(Nfp))
        allocate(rupt_perm(Nfp,2,glob_Nface))
        allocate(rupt_perm_info(4,glob_Nface));rupt_perm_info=0
        do i=1,glob_Nface
            do j=1,2
                k=F2T(j,i)
                if(k .gt. 0)then
                    face=face_permutation(triface(:,i),tet(:,k),FtoV)
                    if(face .gt. 0)then
                        if(Fatt(i).gt.0)then ! rupture
                            mesh%fbctype(1,face,k)=2
                            tetface=tet(FtoV(face,:),k)
                            call facenodal_perm(pOrder+1,Nfp,face,&
                                triface(:,i),tetface,&
                                Fmask(:,1),Fperm,l)
                            if(l.gt.0)then
                                mesh%fbctype(2,face,k)=Fatt(i)
                                rupt_perm(:,1,Fatt(i))=Fperm
                                rupt_perm_info(1,Fatt(i))=l
                                rupt_perm_info(3,Fatt(i))=face
                            elseif(l.lt.0)then
                                mesh%fbctype(2,face,k)=-Fatt(i)
                                rupt_perm(:,2,Fatt(i))=Fperm
                                rupt_perm_info(2,Fatt(i))=l
                                rupt_perm_info(4,Fatt(i))=face
                            endif
!                        else ! external boundary source
!                            mesh%fbctype(1,face,k)=3
!                            mesh%fbctype(2,face,k)=i
                        endif
!                        print*,'face=',i,', tet=',j,', face=',face
                    else
                        print*,triface(:,i),tet(:,k)
!                        stop
                    endif
                endif
            enddo
        enddo
        deallocate(Fperm)
    endif

end subroutine Build_geometry

subroutine form_ref_matrices(pOrder,pNp,Nfp,r,s,t,Fmask,matrix)
    integer, intent(in) :: pOrder,pNp,Nfp
    real(kind=rkind),intent(out) :: r(pNp),s(pNp),t(pNp)
    integer, intent(out) :: Fmask(Nfp,4)
    type(matrices) :: matrix
!       pointers:
    real(kind=rkind),pointer :: Dr(:,:),Ds(:,:),Dt(:,:)
    real(kind=rkind),pointer :: Drw(:,:),Dsw(:,:),Dtw(:,:)
    real(kind=rkind),pointer :: V3D(:,:),iV3D(:,:),M3D(:,:),iM3D(:,:)
    real(kind=rkind),pointer :: V2D(:,:),iV2D(:,:),M2D(:,:),iM2D(:,:)
    real(kind=rkind),pointer :: LIFT(:,:),LIFT_2D(:,:)
    real(kind=rkind),pointer :: Dr2D(:,:),Ds2D(:,:)
!       auxilary memory:
    real(kind=rkind) :: V3Dr(pNp,pNp),V3Ds(pNp,pNp),V3Dt(pNp,pNp)
    real(kind=rkind) :: V2Dr(Nfp,Nfp),V2Ds(Nfp,Nfp)
    real(kind=rkind) :: V3DF(Nfp,Nfp,4),iV3DF(Nfp,Nfp,4)
    real(kind=rkind) :: tmp(pNp,pNp),tmp2d(Nfp,Nfp)
    integer :: face,info,ipiv(pNp)
    integer :: n
!       Allocations:
    allocate(matrix% V3D(pNp,pNp)  );   V3D=>matrix% V3D
    allocate(matrix%iV3D(pNp,pNp)  );  iV3D=>matrix%iV3D
    allocate(matrix% M3D(pNp,pNp)  );   M3D=>matrix% M3D
    allocate(matrix%iM3D(pNp,pNp)  );  iM3D=>matrix%iM3D
    allocate(matrix%  Dr(pNp,pNp)  );    Dr=>matrix%Dr
    allocate(matrix%  Ds(pNp,pNp)  );    Ds=>matrix%Ds
    allocate(matrix%  Dt(pNp,pNp)  );    Dt=>matrix%Dt
    allocate(matrix% Drw(pNp,pNp)  );   Drw=>matrix%Drw
    allocate(matrix% Dsw(pNp,pNp)  );   Dsw=>matrix%Dsw
    allocate(matrix% Dtw(pNp,pNp)  );   Dtw=>matrix%Dtw
    allocate(matrix% V2D(Nfp,Nfp)  );   V2D=>matrix% V2D
    allocate(matrix%iV2D(Nfp,Nfp)  );  iV2D=>matrix%iV2D
    allocate(matrix% M2D(Nfp,Nfp)  );   M2D=>matrix% M2D
    allocate(matrix%iM2D(Nfp,Nfp)  );  iM2D=>matrix%iM2D
    allocate(matrix%LIFT(pNp,4*Nfp));  LIFT=>matrix%LIFT
    allocate(matrix%LIFT2D(Nfp,3*(def_pOrder+1)));
                                    LIFT_2D=>matrix%LIFT2D
    allocate(matrix%Dr2D(Nfp,Nfp)  );  Dr2D=>matrix%Dr2D
    allocate(matrix%Ds2D(Nfp,Nfp)  );  Ds2D=>matrix%Ds2D

    call blend_nodes(r,s,t,pOrder,pNp)

    call Vandermonde3D(V3D,pOrder,r,s,t,pNp)

    iV3D=0d0;tmp=V3D
    do n=1,pNp,1
        iV3D(n,n)=1d0;
    enddo
    call dgesv(pNp, pNp, tmp, pNp, ipiv, iV3D, pNp, info)

    M3D=matmul(transpose(iV3D),iV3D)
    iM3D=matmul(V3D,transpose(V3D))

    call GradVandermonde3D(V3Dr,V3Ds,V3Dt,pOrder,r,s,t,pNp)

    Dr=transpose(V3Dr);tmp=transpose(V3D)
    call dgesv(pNp, pNp, tmp, pNp, ipiv, Dr, pNp, info)
    Drw=-matmul(matmul(iM3D,Dr),M3D)
    Dr=transpose(Dr)

    Ds=transpose(V3Ds);tmp=transpose(V3D)
    call dgesv(pNp, pNp, tmp, pNp, ipiv, Ds, pNp, info)
    Dsw=-matmul(matmul(iM3D,Ds),M3D)
    Ds=transpose(Ds)

    Dt=transpose(V3Dt);tmp=transpose(V3D)
    call dgesv(pNp, pNp, tmp, pNp, ipiv, Dt, pNp, info)
    Dtw=-matmul(matmul(iM3D,Dt),M3D)
    Dt=transpose(Dt)

    call Lift3D(LIFT,iV3DF,r,s,t,Fmask,V3D,pOrder,pNp,Nfp)

    call Vandermonde2D(V2D,pOrder,r,s,Nfp)
    iV2D=0d0;tmp2d=V2D
    do n=1,Nfp,1
        iV2D(n,n)=1d0;
    enddo
    call dgesv(Nfp,Nfp,tmp2d,Nfp,ipiv,iV2D,Nfp,info)

    M2D=matmul(transpose(iV2D),iV2D)
    iM2D=matmul(V2D,transpose(V2D))

    call GradVandermonde2D(V2Dr,V2Ds,pOrder,r,s,Nfp)
    Dr2D=transpose(V2Dr);tmp2d=transpose(V2D)
    call dgesv(Nfp, Nfp, tmp2d, Nfp, ipiv, Dr2D, Nfp, info)
    Dr2D=transpose(Dr2D)
    Ds2D=transpose(V2Ds);tmp2d=transpose(V2D)
    call dgesv(Nfp, Nfp, tmp2d, Nfp, ipiv, Ds2D, Nfp, info)
    Ds2D=transpose(Ds2D)

    call Lift2D(LIFT_2D,r,s,V2D,Nfp,def_pOrder+1)

end subroutine form_ref_matrices

subroutine blend_nodes(x,y,z,pOrder,pNp)

    integer :: pOrder,pNp
    real(kind=rkind) :: x(pNp),y(pNp),z(pNp)
    integer :: i,j,k,sk
    real(kind=rkind), allocatable :: r(:), s(:), t(:)
    real(kind=rkind) :: alphastore(15),alpha
    real(kind=rkind),allocatable :: L1(:),L2(:),L3(:),L4(:)
    real(kind=rkind),allocatable :: La(:),Lb(:),Lc(:),Ld(:),&
                    XYZ(:,:),shift(:,:)
    real(kind=rkind),allocatable :: warp1(:),warp2(:),blend(:),denom(:)

    real(kind=rkind) :: v1(3),v2(3),v3(3),v4(3),t1(3,4),t2(3,4)
    real(kind=rkind) :: rtmp
    data alphastore /0d0,0d0,0d0,0.1002d0,1.1332d0,1.5608d0, &
            1.3413d0,1.2577d0,1.1603d0,1.10153d0,0.6080d0,&
            0.4523d0,0.8856d0,0.8717d0,0.9655d0/
    v1(1)=-1d0;v1(2)=-1d0/sqrt(3d0);v1(3)=-1d0/sqrt(6d0)
    v2(1)= 1d0;v2(2)=-1d0/sqrt(3d0);v2(3)=-1d0/sqrt(6d0)
    v3(1)= 0d0;v3(2)= 2d0/sqrt(3d0);v3(3)=-1d0/sqrt(6d0)
    v4(1)= 0d0;v4(2)= 0d0          ;v4(3)= 3d0/sqrt(6d0)

    allocate(r(pNp));allocate(s(pNp));allocate(t(pNp));
    allocate(L1(pNp));allocate(L2(pNp))
    allocate(L3(pNp));allocate(L4(pNp))
    allocate(La(pNp));allocate(Lb(pNp))
    allocate(Lc(pNp));allocate(Ld(pNp))
    allocate(warp1(pNp));allocate(warp2(pNp))
    allocate(blend(pNp));allocate(denom(pNp))
    allocate(XYZ(3,pNp));allocate(shift(3,pNp))

    if(pOrder .le. 15)then
        alpha=alphastore(pOrder)
    else
        alpha=1.0d0
    endif

    sk=1
    do i=1,pOrder+1,1
        do j=1,pOrder+2-i,1
            do k=1,pOrder+3-i-j,1
                r(sk) = -1d0 + dble(k-1)*2d0/dble(pOrder)
                s(sk) = -1d0 + dble(j-1)*2d0/dble(pOrder)
                t(sk) = -1d0 + dble(i-1)*2d0/dble(pOrder)
                sk = sk+1
            enddo
        enddo
    enddo

    L1 = (1d0+t)/2d0; L2 = (1d0+s)/2d0
    L3 = -(1d0+r+s+t)/2d0; L4 = (1d0+r)/2d0

    t1(:,1) = v2-v1;        
    t1(:,2) = v2-v1;
    t1(:,3) = v3-v2;        
    t1(:,4) = v3-v1;   
    t2(:,1) = v3-0.5d0*(v1+v2); 
    t2(:,2) = v4-0.5d0*(v1+v2);
    t2(:,3) = v4-0.5d0*(v2+v3); 
    t2(:,4) = v4-0.5d0*(v1+v3);  

! normalize tangents
    do i=1,4,1
        rtmp=sqrt(sum((t1(:,i)*t1(:,i))))
        if(rtmp.gt.TOL)t1(:,i)=t1(:,i)/rtmp
        rtmp=sqrt(sum((t2(:,i)*t2(:,i))))
        if(rtmp.gt.TOL)t2(:,i)=t2(:,i)/rtmp
    enddo

    do i=1,3,1
        do j=1,pNp,1
            XYZ(i,j)=L3(j)*v1(i)+L4(j)*v2(i) &
                +L2(j)*v3(i)+L1(j)*v4(i)
        enddo
    enddo
    shift=0.0d0

    do i=1,4,1
        if(i .eq. 1)then
            La = L1; Lb = L2; Lc = L3; Ld = L4
        elseif(i .eq. 2)then
            La = L2; Lb = L1; Lc = L3; Ld = L4
        elseif(i .eq. 3)then
            La = L3; Lb = L1; Lc = L4; Ld = L2
        else
            La = L4; Lb = L1; Lc = L3; Ld = L2
        endif

        call evalshift(warp1,warp2,alpha,Lb,Lc,Ld,pOrder,pNp)

        blend=Lb*Lc*Ld
        denom=(Lb+0.5d0*La)*(Lc+0.5d0*La)*(Ld+0.5d0*La)
        do j=1,pNp,1
            if(denom(j)>TOL)then
            blend(j)=(1d0+(alpha*La(j))**2)*blend(j)/denom(j)
            endif
        enddo

        do j=1,pNp,1
            shift(:,j)=shift(:,j)+&
            blend(j)*warp1(j)*t1(:,i)+blend(j)*warp2(j)*t2(:,i)
        enddo

        do j=1,pNp,1
            if((La(j).lt.TOL) .and. (.not. ((Lb(j).gt.TOL).and.&
                (Lc(j).gt.TOL).and.(Ld(j).gt.TOL))))then
                shift(:,j)=warp1(j)*t1(:,i)+warp2(j)*t2(:,i)
            endif
        enddo
    enddo
    XYZ=XYZ+shift;
    r=XYZ(1,:);s=XYZ(2,:);t=XYZ(3,:)
    x=(v1(1)*(v4(2)*v3(3)-v3(2)*v4(3))&
      +v3(1)*(v1(2)*v4(3)-v4(2)*v1(3))&
      +v4(1)*(v3(2)*v1(3)-v1(2)*v3(3))&
   +(v1(2)*(v3(3)-v4(3))+v3(2)*(v4(3)-v1(3))+v4(2)*(v1(3)-v3(3)))*r &
   +(v1(1)*(v4(3)-v3(3))+v3(1)*(v1(3)-v4(3))+v4(1)*(v3(3)-v1(3)))*s &
   +(v1(1)*(v3(2)-v4(2))+v3(1)*(v4(2)-v1(2))+v4(1)*(v1(2)-v3(2)))*t)&
     /4.0d0/sqrt(2.0d0);
    y=(v1(2)*(v4(1)*v2(3)-v2(1)*v4(3))&
      +v2(2)*(v1(1)*v4(3)-v4(1)*v1(3))&
      +v4(2)*(v2(1)*v1(3)-v1(1)*v2(3))&
   +(v1(2)*(v4(3)-v2(3))+v2(2)*(v1(3)-v4(3))+v4(2)*(v2(3)-v1(3)))*r &
   +(v1(1)*(v2(3)-v4(3))+v2(1)*(v4(3)-v1(3))+v4(1)*(v1(3)-v2(3)))*s &
   +(v1(1)*(v4(2)-v2(2))+v2(1)*(v1(2)-v4(2))+v4(1)*(v2(2)-v1(2)))*t)&
     /4.0d0/sqrt(2.0d0);
    z=(v1(3)*(v3(1)*v2(2)-v2(1)*v3(2))&
      +v2(3)*(v1(1)*v3(2)-v3(1)*v1(2))&
      +v3(3)*(v2(1)*v1(2)-v1(1)*v2(2))&
   +(v1(2)*(v2(3)-v3(3))+v2(2)*(v3(3)-v1(3))+v3(2)*(v1(3)-v2(3)))*r &
   +(v1(1)*(v3(3)-v2(3))+v2(1)*(v1(3)-v3(3))+v3(1)*(v2(3)-v1(3)))*s &
   +(v1(1)*(v2(2)-v3(2))+v2(1)*(v3(2)-v1(2))+v3(1)*(v1(2)-v2(2)))*t)&
     /4.0d0/sqrt(2.0d0);

    x=x*2d0-1.0d0;y=y*2d0-1.0d0;z=z*2d0-1.0d0

    deallocate(XYZ)
    deallocate(La);deallocate(Lb);deallocate(Lc);deallocate(Ld)
    deallocate(L1);deallocate(L2);deallocate(L3);deallocate(L4)
    deallocate(r);deallocate(s);deallocate(t)
    deallocate(warp1);deallocate(warp2)
    deallocate(shift);deallocate(blend);deallocate(denom)

end subroutine blend_nodes

subroutine evalshift(w1,w2,alpha,L1,L2,L3,pOrder,pNp)
    integer :: pOrder,pNp
    real(kind=rkind) :: w1(pNp),w2(pNp),w3(pNp),alpha
    real(kind=rkind) :: L1(pNp),L2(pNp),L3(pNp)
    real(kind=rkind) :: gaussX(pOrder+1),&
                        warp1(pNp),warp2(pNp),warp3(pNp)

    call JacobiGL(gaussX,0.0D0,0.0D0,pOrder)
    gaussX=-gaussX
    w1=L3-L2;w2=L1-L3;w3=L2-L1
    call evalwarp(warp1,gaussX,w1,pOrder,pNp)
    call evalwarp(warp2,gaussX,w2,pOrder,pNp)
    call evalwarp(warp3,gaussX,w3,pOrder,pNp)

    warp1=L2*L3*4d0*warp1*(1d0+(alpha*L1)**2)
    warp2=L1*L3*4d0*warp2*(1d0+(alpha*L2)**2)
    warp3=L1*L2*4d0*warp3*(1d0+(alpha*L3)**2)
    w1=warp1-0.5d0*(warp2+warp3)
    w2=sqrt(3.0d0)/2d0*(warp2-warp3)

end subroutine evalshift

subroutine evalwarp(warp,xnodes,xout,pOrder,pNp)
    integer :: pOrder,pNp
    real(kind=rkind) :: warp(pNp),xnodes(pOrder+1),xout(pNp)
    real(kind=rkind),allocatable :: xeq(:),d(:)
    integer i,j
    warp=0d0;
    allocate(xeq(pOrder+1));allocate(d(pNp))
    do i=1,pOrder+1,1
        xeq(i)=-1d0+2d0*dble(pOrder+1-i)/dble(pOrder)
    enddo
    do i=1,pOrder+1,1
        d=xnodes(i)-xeq(i)
        do j=2,pOrder,1
            if(i .ne. j)then
                d=d*(xout-xeq(j))/(xeq(i)-xeq(j))
            endif
        enddo

        if(i.ne.1)then
            d=-d/(xeq(i)-xeq(1))
        endif
        if(i.ne.pOrder+1)then
            d=d/(xeq(i)-xeq(pOrder+1))
        endif
        warp=warp+d
    enddo

end subroutine evalwarp

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

subroutine Lift3D(lift,V3DF,r,s,t,Fmask,V3D,pOrder,pNp,Nfp)
    integer      :: pOrder,pNp,Nfp
    integer      :: N,xdim,Fmask(Nfp,4),face,i,j,k,l
    real(kind=rkind) :: lift(pNp,4*Nfp),r(pNp),s(pNp),t(pNp)
    real(kind=rkind) :: V3D(pNp,pNp),V3DF(Nfp,Nfp,4)
    real(kind=rkind) :: Emat(pNp,4*Nfp),Vface(Nfp,Nfp),&
                tmp(Nfp,Nfp),faceR(Nfp),faceS(Nfp),&
                eye(Nfp,Nfp),ap((Nfp*(Nfp+1))/2)
    integer :: IPIV(Nfp)

    Fmask=0
    i=1;j=1;k=1;l=1
    do n=1,pNp,1
        if(abs(1.0+t(n)).le.TOL)then
            Fmask(i,1)=n;i=i+1
        endif
        if(abs(1.0+s(n)).le.TOL)then
            Fmask(j,2)=n;j=j+1
        endif
        if(abs(1.0+r(n)+s(n)+t(n)).le.TOL)then
            Fmask(k,3)=n;k=k+1
        endif
        if(abs(1.0+r(n)).le.TOL)then
            Fmask(l,4)=n;l=l+1
        endif
    enddo

    Emat=0.0d0
    do face=1,4
        eye=0.0d0
        do i=1,Nfp
            eye(i,i)=1.0d0
        enddo
        if(face.eq.1)then
            faceR=r(Fmask(:,1));faceS=s(Fmask(:,1))
        elseif(face.eq.2)then
            faceR=r(Fmask(:,2));faceS=t(Fmask(:,2))
        elseif(face.eq.3)then
            faceR=s(Fmask(:,3));faceS=t(Fmask(:,3))
        else
            faceR=s(Fmask(:,4));faceS=t(Fmask(:,4))
        endif
        call Vandermonde2D(Vface,pOrder,faceR,faceS,Nfp)
        V3DF(:,:,face)=Vface
        tmp=transpose(Vface)
        tmp=matmul(Vface,tmp)
        do i=1,Nfp
            do j=i,Nfp
                ap(i+((j-1)*j)/2)=tmp(i,j)
            enddo
        enddo
        call dspsv('U',Nfp,Nfp,ap,IPIV,eye,Nfp,i)
        Emat(Fmask(:,face),(face-1)*Nfp+1:face*Nfp)=&
        Emat(Fmask(:,face),(face-1)*Nfp+1:face*Nfp)+eye
    enddo
    LIFT=matmul(V3D,matmul(transpose(V3D),Emat))
end subroutine Lift3D

subroutine Build_maps3D(tet,neigh,EtoF,Fmask,vmapM,vmapP,NNele,&
        pOrder,pNp,Nfp)
    integer :: NNele
    integer :: pOrder,pNp,Nfp
    integer :: tet(4,NNele),neigh(4,NNele),&
                   EtoF(4,NNele)
    integer ::Fmask(Nfp,4),FtoV(4,3)
    integer :: vmapM(Nfp*4),vmapP(Nfp*NNele*4)
    integer :: i,j,k
    integer :: N
    integer, allocatable :: tri(:,:),tri1(:),tri2(:),tri3(:)
    integer, allocatable :: ttri1(:),ttri2(:),ttri3(:)
    data((FtoV(i,j),i=1,4),j=1,3)/1,1,2,1,2,2,3,3,3,4,4,4/
    character (len=50) :: filename
    logical alive

    N=pOrder+1
    allocate(tri(N,N))
    k=1
    do i=1,N
        do j=1,N+1-i
            tri(j,i)=k
            k=k+1
        enddo
    enddo
    k=(N*(N+1))/2
    allocate(tri1(k));allocate(tri2(k));allocate(tri3(k))
    allocate(ttri1(k));allocate(ttri2(k));allocate(ttri3(k))
    tri1=0;tri2=0;tri3=0;ttri1=0;ttri2=0;ttri3=0;k=1
    do i=1,N
        do j=1,N+1-i
            tri1(tri(i,j))=k
            tri3(tri(j,N+2-i-j))=k
            tri2(tri(N+2-i-j,i))=k
            ttri1(tri(j,i))=k
            ttri3(tri(N+2-i-j,j))=k
            ttri2(tri(i,N+2-i-j))=k
            k=k+1
        enddo
    enddo
    do j=1,4
        k=(j-1)*Nfp
        vmapM(k+1:k+Nfp)=Fmask(:,j)
    enddo

    do i=1,NNele
        do j=1,4
            k=(i-1)*4*Nfp+(j-1)*Nfp
            ! print*,k
            if(neigh(j,i).gt.0)then
            ! not self-mapping
                if(j+EtoF(j,i).eq.5.or.j.eq.EtoF(j,i))then
                    if    ( tet(FtoV(j,1),i).eq.&
                            tet((FtoV(EtoF(j,i),1)),neigh(j,i)))then
                    ! ith ele's jth fac's 1st node == 
                    ! ith ele's jth neigh's mapping fac's 1st node
                        vmapP(k+1:k+NfP)=&
                            (neigh(j,i)-1)*pNp+Fmask(tri1,EtoF(j,i))
                    elseif( tet(FtoV(j,1),i).eq.&
                            tet((FtoV(EtoF(j,i),2)),neigh(j,i)))then
                    ! ith ele's jth fac's 1st node == 
                    ! ith ele's jth neigh's mapping fac's 2nd node
                        vmapP(k+1:k+Nfp)=&
                            (neigh(j,i)-1)*pNp+Fmask(tri2,EtoF(j,i))
                    elseif( tet(FtoV(j,1),i).eq.&
                            tet((FtoV(EtoF(j,i),3)),neigh(j,i)))then
                    ! ith ele's jth fac's 1st node == 
                    ! ith ele's jth neigh's mapping fac's 3rd node
                        vmapP(k+1:k+Nfp)=&
                            (neigh(j,i)-1)*pNp+Fmask(tri3,EtoF(j,i))
                    endif
                else
                    if(     tet(FtoV(j,1),i).eq.&
                            tet((FtoV(EtoF(j,i),1)),neigh(j,i)))then
                    ! ith ele's jth fac's 1st node == 
                    ! ith ele's jth neigh's mapping fac's 1st node
                        vmapP(k+1:k+NfP)=&
                            (neigh(j,i)-1)*pNp+Fmask(ttri1,EtoF(j,i))
                    elseif( tet(FtoV(j,1),i).eq.&
                            tet((FtoV(EtoF(j,i),2)),neigh(j,i)))then
                    ! ith ele's jth fac's 1st node == 
                    ! ith ele's jth neigh's mapping fac's 2nd node
                        vmapP(k+1:k+Nfp)=&
                            (neigh(j,i)-1)*pNp+Fmask(ttri2,EtoF(j,i))
                    elseif( tet(FtoV(j,1),i).eq.&
                            tet((FtoV(EtoF(j,i),3)),neigh(j,i)))then
                    ! ith ele's jth fac's 1st node == 
                    ! ith ele's jth neigh's mapping fac's 3rd node
                        vmapP(k+1:k+Nfp)=&
                            (neigh(j,i)-1)*pNp+Fmask(ttri3,EtoF(j,i))
                    endif
                endif
            else
            ! self-mapping
                vmapP(k+1:k+Nfp)=(i-1)*pNp+Fmask(:,j)
                ! vmapP(k+1:k+Nfp)=-1
            endif
        enddo
    enddo
    deallocate( tri1, tri2, tri3)
    deallocate(ttri1,ttri2,ttri3)

end subroutine Build_maps3D

subroutine Lift2D(lift,r,s,V2D,pNp,Nfp)
    integer          :: pNp,Nfp
    integer          :: xdim,Fmask(Nfp,3),face,i,j,k
    real(kind=rkind) :: V2D(pNp,pNp)
    real(kind=rkind) :: lift(pNp,3*Nfp),r(pNp),s(pNp)
    real(kind=rkind) :: Emat(pNp,3*Nfp),tmp(Nfp,Nfp)
    real(kind=rkind) :: faceR(Nfp),Vface(Nfp,Nfp),eye(Nfp,Nfp)
    real(kind=rkind) :: ap((Nfp*(Nfp+1))/2)
    integer :: IPIV(Nfp)

    do i=1,Nfp
        Fmask(i,1)=i
        do j=1,Nfp+1-i
            if(j.eq. 1)Fmask(Nfp+1-i,3)=k
            if(j.eq. Nfp+1-i)Fmask(i,2)=k
            k=k+1
        enddo
    enddo
    Emat=0.0d0
    do face=1,3
        eye=0.0d0
        do i=1,Nfp
            eye(i,i)=1.0d0
        enddo
        if(face.eq.1)faceR=r(Fmask(:,1))
        if(face.eq.2)faceR=r(Fmask(:,2))
        if(face.eq.3)faceR=s(Fmask(:,3))
        call Vandermonde1D(Vface,def_pOrder,faceR,Nfp)
        tmp=transpose(Vface)
        tmp=matmul(Vface,tmp)
        do i=1,Nfp
            do j=i,Nfp
                ap(i+((j-1)*j)/2)=tmp(i,j)
            enddo
        enddo
        call dspsv('U',Nfp,Nfp,ap,IPIV,eye,Nfp,i)
        Emat(Fmask(:,face),(face-1)*Nfp+1:face*Nfp)=&
        Emat(Fmask(:,face),(face-1)*Nfp+1:face*Nfp)+eye
    enddo
    LIFT=matmul(V2D,matmul(transpose(V2D),Emat))
end subroutine Lift2D

subroutine Build_maps2D(tri,neigh,EtoF,Fmask,vmapM,vmapP,Nele,&
        pOrder,pNp,Nfp)
    integer :: Nele
    integer :: pOrder,pNp,Nfp
    integer :: tri(3,Nele),neigh(3,Nele),EtoF(3,Nele)
    integer ::Fmask(Nfp,3),side(Nfp)
    integer :: vmapM(Nfp*Nele*3),vmapP(Nfp*Nele*3)
    integer i,j,k,N

    N=def_pOrder+1
    do i=1,Nfp
        side(i)=Nfp+1-i
    enddo
    do i=1,Nele
        do j=1,3
            k=(i-1)*3*Nfp+(j-1)*Nfp
            vmapM(k+1:k+Nfp)=(i-1)*pNp+Fmask(:,j)
        enddo
    enddo
    do i=1,Nele
        do j=1,3
            k=(i-1)*3*Nfp+(j-1)*Nfp
            if(neigh(j,i).gt.0)then
            ! not self-mapping
                if(EToF(j,i).eq.1)then
    vmapP(k+1:k+NfP)=(neigh(j,i)-1)*pNp+Fmask(side,EtoF(j,i))
                elseif(EtoF(j,i).eq.2)then
    vmapP(k+1:k+Nfp)=(neigh(j,i)-1)*pNp+Fmask(side,EtoF(j,i))
                elseif(EtoF(j,i).eq.3)then
    vmapP(k+1:k+Nfp)=(neigh(j,i)-1)*pNp+Fmask(side,EtoF(j,i))
                endif
            else
            ! self-mapping
                vmapP(k+1:k+Nfp)=(i-1)*pNp+Fmask(:,j)
            endif
        enddo
    enddo

end subroutine Build_maps2D

subroutine tet_subelement(N,tt,Nsele)

    ! Input:
    !   N: pOrder
    !   Nsele: number of subelements per element
    ! if(N == 1) Nsele = N*(N+1)*(N+2)/6
    ! if(N == 2) Nsele = N*(N+1)*(N+2)/6 + (N-1)*N*(N+1)/6*4
    ! if(N >= 3) Nsele = N*(N+1)*(N+2)/6 + (N-1)*N*(N+1)/6*4 
    !                  + (N-2)*(N-1)*N/6
    ! Output:
    !   tt(Nsele)
    integer,intent(in)  :: N, Nsele
    integer,intent(out) :: tt(Nsele,4)
    integer :: t1,t2,t3,t4,t5,t6
    integer,allocatable :: M(:,:,:)
    integer :: sk,i,j,k
    
    if(N.eq.1)then
        tt(1,1)=1;tt(1,2)=2;tt(1,3)=3;tt(1,4)=4
        return
    endif
    
    allocate(M(N+1,N+1,N+1))
    sk=0
    do k=1,N+1
        do j=1,N+2-k
            do i=1,N+3-k-j
                sk=sk+1
                M(i,j,k)=sk
            enddo
        enddo
    enddo
    sk=0
    do k=1,N
        do j=1,N+1-k
            do i=1,N+2-k-j
                sk=sk+1
                tt(sk,1)=M(i,j,k)
                tt(sk,2)=M(i+1,j,k)
                tt(sk,3)=M(i,j+1,k)
                tt(sk,4)=M(i,j,k+1)
            enddo
        enddo
    enddo
    do k=1,N-1
        do j=1,N-k
            do i=1,N+1-k-j
                t1=M(i+1,j,k);t2=M(i+1,j+1,k);t3=M(i,j+1,k)
                t4=M(i,j,k+1);t5=M(i+1,j,k+1);t6=M(i,j+1,k+1)
                sk=sk+1;
                tt(sk,1)=t4;tt(sk,2)=t2;tt(sk,3)=t5;tt(sk,4)=t1
                sk=sk+1;
                tt(sk,1)=t4;tt(sk,2)=t2;tt(sk,3)=t1;tt(sk,4)=t3
                sk=sk+1;
                tt(sk,1)=t4;tt(sk,2)=t2;tt(sk,3)=t3;tt(sk,4)=t6
                sk=sk+1;
                tt(sk,1)=t4;tt(sk,2)=t2;tt(sk,3)=t6;tt(sk,4)=t5
            enddo
        enddo
    enddo
    if(N.eq.2)then
        return
    endif
    do k=1,N-2
        do j=1,N-k-1
            do i=1,N-k-j
                sk=sk+1;
                tt(sk,1)=M(i,j+1,k+1)
                tt(sk,2)=M(i+1,j+1,k+1)
                tt(sk,3)=M(i+1,j,k+1)
                tt(sk,4)=M(i+1,j+1,k)
            enddo
        enddo
    enddo

end subroutine tet_subelement

function face_permutation(face,tet,FtoV)
        integer :: face(3),tet(4),FtoV(4,3)
        integer :: face_permutation,i,f(3)
        face_permutation=-1
        do i=1,4
            f=tet(FtoV(i,:))
            if(    sum(face).eq.   sum(f) .and. &
                maxval(face).eq.maxval(f) .and. &
                minval(face).eq.minval(f) )then
                face_permutation=i
                exit
            endif
        enddo
end function face_permutation

subroutine facenodal_perm(N,Nfp,f,surftri,facetri,Fmask,map,kk)
        integer,intent(in) :: N,Nfp,surftri(3),facetri(3),Fmask(Nfp),f
        integer,intent(out) :: map(Nfp),kk
        integer :: tri(N,N),tri1(Nfp),tri2(Nfp),tri3(Nfp)
        integer :: ttri1(Nfp),ttri2(Nfp),ttri3(Nfp)
        integer :: i,j,k

        k=0
        do i=1,N
            do j=1,N+1-i
                k=k+1
                tri(j,i)=k
            enddo
        enddo
        k=0
        do i=1,N
            do j=1,N+1-i
                k=k+1
                tri1(tri(i,j))=k
                tri3(tri(j,N+2-i-j))=k
                tri2(tri(N+2-i-j,i))=k
            enddo
        enddo

        kk=0
        if(f.eq.2 .or. f.eq.3)then
            if(facetri(1).eq.surftri(1))then
! element face's 1st node == surface facet's 1st node
                if(facetri(2).eq.surftri(3))then
                        ! positive map
                        map=tri1
                        kk=1
                else
                        ! negative map
                        map=tri1(tri1)
                        kk=-1
                endif
            elseif(facetri(1).eq.surftri(2))then
! element face's 1st node == surface facet's 2nd node
                if(facetri(2).eq.surftri(1))then
                        ! positive map
                        map=tri2
                        kk=2
                else
                        ! negative map
                        map=tri2(tri1)
                        kk=-2
                endif
            elseif(facetri(1).eq.surftri(3))then
! element face's 1st node == surface facet's 3rd node
                if(facetri(2).eq.surftri(2))then
                        ! positive map
                        map=tri3
                        kk=3
                else
                        ! negative map
                        map=tri3(tri1)
                        kk=-3
                endif
            endif
        else
            if(facetri(1).eq.surftri(1))then
! element face's 1st node == surface facet's 1st node
                if(facetri(2).eq.surftri(3))then
                        ! positive map
                        map=tri1(tri1)
                        kk=1
                else
                        ! negative map
                        map=tri1
                        kk=-1
                endif
            elseif(facetri(1).eq.surftri(2))then
! element face's 1st node == surface facet's 2nd node
                if(facetri(2).eq.surftri(1))then
                        ! positive map
                        map=tri2(tri1)
                        kk=2
                else
                        ! negative map
                        map=tri2
                        kk=-2
                endif
            elseif(facetri(1).eq.surftri(3))then
! element face's 1st node == surface facet's 3rd node
                if(facetri(2).eq.surftri(2))then
                        ! positive map
                        map=tri3(tri1)
                        kk=3
                else
                        ! negative map
                        map=tri3
                        kk=-3
                endif
            endif
        endif

end subroutine facenodal_perm

end module geometry_mod

