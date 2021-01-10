!*******************************************************************!
!*  This module obtains semidiscretized wave equation              *!
!*                                                                 *!
!*  Author: Ruichao Ye    Email: Ruichao.ye@gmail.com              *!
!*******************************************************************!

!--------------------------------------------------------------------
module grav_mod
!--------------------------------------------------------------------

    use datatype_mod, only : rkind,matrices,sources,rupture,&
                             tetmesh_geometry,tetmesh_domain,&
                             tetmesh_material,&
                             vector_array,tensor_array,wavefield,&
                             auxilary_array,reset_vector,&
                             reset_tensor,&
                             reset_wavefield,reset_auxilary_array,&
                             BLK_axpy,DiagMM,Array_permute
    use para_mod,     only : def_pOrder

!--------------------------------------------------------------------

    implicit none

!--------------------------------------------------------------------
    contains
!--------------------------------------------------------------------

subroutine form_grav_blocMatrix(pNp,Dr,Ds,Dt,Drw,Dsw,Dtw,iJac,&
        T0xx,T0yy,T0zz,T0yz,T0xz,T0xy,gx,gy,gz,rho,Mass,invM,A)
! input:
    integer,intent(in):: pNp
    real(kind=rkind),intent(in) :: Dr(pNp,pNp),Drw(pNp,pNp)
    real(kind=rkind),intent(in) :: Ds(pNp,pNp),Dsw(pNp,pNp)
    real(kind=rkind),intent(in) :: Dt(pNp,pNp),Dtw(pNp,pNp)
    real(kind=rkind),intent(in) :: Mass(pNp,pNp),invM(pNp,pNp)
    real(kind=rkind),intent(in) :: iJac(9)
    real(kind=rkind),intent(in) :: T0xx(pNp),T0yy(pNp),T0zz(pNp)
    real(kind=rkind),intent(in) :: T0yz(pNp),T0xz(pNp),T0xy(pNp)
    real(kind=rkind),intent(in) :: gx(pNp),gy(pNp),gz(pNp),rho(pNp)
! output
    real(kind=rkind),intent(out) :: A(pNp*3,pNp*3)
! auxilary:
    real(kind=rkind),target :: DD(pNp,pNp,3),DDw(pNp,pNp,3)
    real(kind=rkind) :: g0(pNp,3),Dg0(pNp,3,3)
    real(kind=rkind) :: irT0(pNp,3,3),irDivT0(pNp,3),irDT0(pNp,3,3,3)
    real(kind=rkind),target :: tmp(pNp,pNp)
    real(kind=rkind) :: tmp1(pNp,pNp),tmp2(pNp,pNp),tmp3(pNp,pNp)
    integer :: i,j,l,m
    real(kind=rkind),parameter :: delta(3,3) = &
        reshape((/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/),&
		shape(delta))
    real(kind=rkind),pointer:: Dx(:,:),Dy(:,:),Dz(:,:),Dtmp(:)
    Dtmp=>tmp(1:pNp**2:pNp+1,1)
    Dx=>DD(:,:,1);Dy=>DD(:,:,2);Dz=>DD(:,:,3)

    DD(:,:,1)=iJac(1)*Dr+iJac(2)*Ds+iJac(3)*Dt
    DD(:,:,2)=iJac(4)*Dr+iJac(5)*Ds+iJac(6)*Dt
    DD(:,:,3)=iJac(7)*Dr+iJac(8)*Ds+iJac(9)*Dt

    DDw(:,:,1)=iJac(1)*Drw+iJac(2)*Dsw+iJac(3)*Dtw
    DDw(:,:,2)=iJac(4)*Drw+iJac(5)*Dsw+iJac(6)*Dtw
    DDw(:,:,3)=iJac(7)*Drw+iJac(8)*Dsw+iJac(9)*Dtw

!print*,maxval(abs(DD)),maxval(abs(DDW))
!stop
    g0(:,1)=gx;g0(:,2)=gy;g0(:,3)=gz

    Dg0(:,1,1)=matmul(Dx,gx)
    Dg0(:,2,2)=matmul(Dy,gy)
    Dg0(:,3,3)=matmul(Dz,gz)
    Dg0(:,2,3)=(matmul(Dy,gz)+matmul(Dz,gy))/2d0
    Dg0(:,1,3)=(matmul(Dx,gz)+matmul(Dz,gx))/2d0
    Dg0(:,1,2)=(matmul(Dx,gy)+matmul(Dy,gx))/2d0
    Dg0(:,3,2)=Dg0(:,2,3)
    Dg0(:,3,1)=Dg0(:,1,3)
    Dg0(:,2,1)=Dg0(:,1,2)

    irT0(:,1,1)=T0xx/rho;irT0(:,2,3)=T0yz/rho;irT0(:,3,2)=irT0(:,2,3)
    irT0(:,2,2)=T0yy/rho;irT0(:,1,3)=T0xz/rho;irT0(:,3,1)=irT0(:,1,3)
    irT0(:,3,3)=T0zz/rho;irT0(:,1,2)=T0xy/rho;irT0(:,2,1)=irT0(:,1,2)

    irDT0(:,1,1,1)=matmul(Dx,T0xx)/rho
    irDT0(:,2,1,1)=matmul(Dy,T0xx)/rho
    irDT0(:,3,1,1)=matmul(Dz,T0xx)/rho
    irDT0(:,1,2,2)=matmul(Dx,T0yy)/rho
    irDT0(:,2,2,2)=matmul(Dy,T0yy)/rho
    irDT0(:,3,2,2)=matmul(Dz,T0yy)/rho
    irDT0(:,1,3,3)=matmul(Dx,T0zz)/rho
    irDT0(:,2,3,3)=matmul(Dy,T0zz)/rho
    irDT0(:,3,3,3)=matmul(Dz,T0zz)/rho
    irDT0(:,1,2,3)=matmul(Dx,T0yz)/rho; irDT0(:,1,3,2)=irDT0(:,1,2,3)
    irDT0(:,1,1,3)=matmul(Dx,T0xz)/rho; irDT0(:,1,3,1)=irDT0(:,1,1,3)
    irDT0(:,1,1,2)=matmul(Dx,T0xy)/rho; irDT0(:,1,2,1)=irDT0(:,1,1,2)
    irDT0(:,2,2,3)=matmul(Dy,T0yz)/rho; irDT0(:,2,3,2)=irDT0(:,2,2,3)
    irDT0(:,2,1,3)=matmul(Dy,T0xz)/rho; irDT0(:,2,3,1)=irDT0(:,2,1,3)
    irDT0(:,2,1,2)=matmul(Dy,T0xy)/rho; irDT0(:,2,2,1)=irDT0(:,2,1,2)
    irDT0(:,3,2,3)=matmul(Dz,T0yz)/rho; irDT0(:,3,3,2)=irDT0(:,3,2,3)
    irDT0(:,3,1,3)=matmul(Dz,T0xz)/rho; irDT0(:,3,3,1)=irDT0(:,3,1,3)
    irDT0(:,3,1,2)=matmul(Dz,T0xy)/rho; irDT0(:,3,2,1)=irDT0(:,3,1,2)

    irDivT0(:,1)=irDT0(:,1,1,1)+irDT0(:,2,2,1)+irDT0(:,3,3,1)
    irDivT0(:,2)=irDT0(:,1,1,2)+irDT0(:,2,2,2)+irDT0(:,3,3,2)
    irDivT0(:,3)=irDT0(:,1,1,3)+irDT0(:,2,2,3)+irDT0(:,3,3,3)

!print*,irT0

!!!! math formulation
!!!! term  1*: - D^T_l M irT_lm D_m U_i
!!!! term  2*: + D^T_l M irT_ij D_l U_j
!!!! term  3 : - D^T_i M irDdT_j U_j
!!!! term  4*: - irDdT_i M D_j U_j
!!!! term  5 : + D^T_l M irDT_jil U_j
!!!! term  6*: + irDT_ijl M D_l U_j
!!!! term  7 : - 2 M Dg_ij U_j
!!!! term  8 : - D^T_j M g_i U_j
!!!! term  9*: - g_j M D_i U_j
!!!! term 10 : + D^T_i M g_j U_j
!!!! term 11*: + g_i M D_j U_j

!!!! term 1
    tmp=0d0;tmp1=0d0
    do l=1,3
    do m=1,3
        Dtmp=irT0(:,l,m)
        tmp3=matmul(tmp,Mass)
        tmp3=(tmp3+transpose(tmp3))/2d0
        tmp3=matmul(transpose(DDw(:,:,l)),tmp3)
        tmp3=matmul(tmp3,DDw(:,:,m))
        tmp1=tmp1+tmp3
    enddo
    enddo
    do i=1,3
    do j=1,3
        tmp2=0d0
!!!!! term 1
!        if(i.eq.j)then
!            tmp2=-tmp1
!        endif
!!!!! term 2
!        Dtmp=irT0(:,i,j)
!        tmp3=matmul(tmp,Mass)
!        tmp3=(tmp3+transpose(tmp3))/2d0
!        do l=1,3
!            tmp2=tmp2+matmul(&
!                matmul(transpose(DDw(:,:,l)),tmp3),DDw(:,:,l))
!        enddo
!!!! term 3
        tmp3=matmul(transpose(DDw(:,:,i)),Mass)
        Dtmp=irDivT0(:,j)
        tmp3=matmul(tmp3,tmp)
        tmp2=tmp2-tmp3
!!!!! term 4
!        Dtmp=irDivT0(:,i)
!        tmp3=matmul(tmp,Mass)
!        tmp3=matmul(tmp3,DDw(:,:,j))
!        tmp2=tmp2-tmp3
!!!! term 5
        do l=1,3
            tmp3=matmul(transpose(DDw(:,:,l)),Mass)
            Dtmp=irDT0(:,j,i,l)
            tmp3=matmul(tmp3,tmp)
            tmp2=tmp2+tmp3
        enddo
!!!!! term 6
!        do l=1,3
!            Dtmp=irDT0(:,i,j,l)
!            tmp3=matmul(tmp,Mass)
!            tmp3=matmul(tmp3,DDw(:,:,l))
!            tmp2=tmp2+tmp3
!        enddo
!!!! term 7
        Dtmp=Dg0(:,i,j)
        tmp3=matmul(tmp,Mass)
        tmp3=tmp3+transpose(tmp3)
        tmp2=tmp2-tmp3
!!!! term 8
        tmp3=matmul(transpose(DDw(:,:,j)),Mass)
        Dtmp=g0(:,i)
        tmp3=matmul(tmp3,tmp)
        tmp2=tmp2-tmp3
!!!!! term 9
!        Dtmp=g0(:,j)
!        tmp3=matmul(tmp,Mass)
!        tmp3=matmul(tmp3,DDw(:,:,i))
!        tmp2=tmp2-tmp3
!!! term 10
        tmp3=matmul(transpose(DDw(:,:,i)),Mass)
        Dtmp=g0(:,j)
        tmp3=matmul(tmp3,tmp)
        tmp2=tmp2+tmp3
!!!!! term 11
!        Dtmp=g0(:,i)
!        tmp3=matmul(tmp,Mass)
!        tmp3=matmul(tmp3,DDw(:,:,j))
!        tmp2=tmp2+tmp3
!!!! assign to matrix
        A((i-1)*pNp+1:i*pNp,(j-1)*pNp+1:j*pNp)=0.5d0*matmul(invM,tmp2)
!        A((i-1)*pNp+1:i*pNp,(j-1)*pNp+1:j*pNp)=tmp2
    enddo
    enddo

!    ! debugging
!    print*,sum((A-transpose(A))**2)
!    print*,A
!    stop
!!!    ! end debugging

end subroutine form_grav_blocMatrix

subroutine grav_blocMatrix(pNp,mesh,matrix,T0,g0,rho)

    integer                :: pNp,pNpm1
    type(tetmesh_geometry) :: mesh
    type(matrices)         :: matrix
    type(tensor_array)     :: T0
    type(vector_array)     :: g0
    real(kind=rkind)       :: rho(*)
    real(kind=rkind) :: sig0(pNp),T0xx(pNp),T0yy(pNp),T0zz(pNp)
    real(kind=rkind) :: M3D(pNp,pNp)
    integer :: i,head,tail
    M3D=matrix%M3D
!    pNpm1=def_pOrder*(def_pOrder+1)*(def_pOrder+2)/6
!    M3D=matmul(&
!        transpose(matrix%iV3D(1:pNpm1,:)),matrix%iV3D(1:pNpm1,:))
    matrix%G0=0d0
    do i=1,mesh%Nhele
        head=(i-1)*pNp+1;tail=i*pNp
        T0xx=T0%xx(head:tail)
        T0yy=T0%yy(head:tail)
        T0zz=T0%zz(head:tail)
!        sig0=(T0xx+T0yy+T0zz)/3d0
!        T0xx=T0xx-sig0;T0yy=T0yy-sig0;T0zz=T0zz-sig0
        call form_grav_blocMatrix(pNp,&
            matrix%Dr,matrix%Ds,matrix%Dt,&
            matrix%Drw,matrix%Dsw,matrix%Dtw,&
            mesh%invJac(i,1:9),T0xx,T0yy,T0zz,&
            T0%yz(head:tail),&
            T0%xz(head:tail),&
            T0%xy(head:tail),&
            g0%x(head:tail),&
            g0%y(head:tail),&
            g0%z(head:tail),&
            rho(head:tail),M3D,matrix%iM3D,&
            matrix%G0(:,:,i))
    enddo

end subroutine grav_blocMatrix

subroutine apply_grav_blocMV(pNp,Nele,hVx,hVy,hVz,Ux,Uy,Uz,G0)
    integer :: pNp,Nele,i,head,tail
    real(kind=rkind) :: G0(pNp,pNp,Nele),&
        hVx(pNp*Nele),hVy(pNp*Nele),hVz(pNp*Nele),&
         Ux(pNp*Nele), Uy(pNp*Nele), Uz(pNp*Nele)
    real(kind=rkind) :: U(pNp*3),V(pNp*3)
    do i=1,Nele
        head=(i-1)*pNp+1;tail=i*pNp
        U(      1:  pNp)=Ux(head:tail)
        U(  pNp+1:2*pNp)=Uy(head:tail)
        U(2*pNp+1:3*pNp)=Uz(head:tail)
        call DGEMV('n',3*pNp,3*pNp,1d0,G0(1:,1,i),3*pNp,U,1,&
            0d0,V,1)
        hVx(head:tail)=hVx(head:tail)+V(      1:  pNp)
        hVy(head:tail)=hVy(head:tail)+V(  pNp+1:2*pNp)
        hVz(head:tail)=hVz(head:tail)+V(2*pNp+1:3*pNp)
    enddo
end subroutine apply_grav_blocMV

subroutine form_rupt_grav_blocMatrix(Np,pNp,Dr,Ds,Dt,iJac,&
        irJL,nx,ny,nz,T0xx,T0yy,T0zz,T0yz,T0xz,T0xy,A)
! input:
    integer,intent(in):: Np,pNp
    real(kind=rkind),intent(in) :: Dr(pNp,pNp)
    real(kind=rkind),intent(in) :: Ds(pNp,pNp)
    real(kind=rkind),intent(in) :: Dt(pNp,pNp)
    real(kind=rkind),intent(in) :: iJac(9),nx,ny,nz
    real(kind=rkind),intent(in) :: irJL(Np,pNp)
    real(kind=rkind),intent(in) :: T0xx(pNp),T0yy(pNp),T0zz(pNp)
    real(kind=rkind),intent(in) :: T0yz(pNp),T0xz(pNp),T0xy(pNp)
! output
    real(kind=rkind),intent(inout) :: A(Np*3,pNp*3)
! auxilary:
    real(kind=rkind),target :: DD(pNp,pNp,3)
    real(kind=rkind) :: T0(pNp,3,3),n(3)
    real(kind=rkind),target :: tmp(pNp,pNp)
    real(kind=rkind) :: tmp1(pNp,pNp),tmp2(pNp,pNp),tmp3(pNp,pNp)
    integer :: i,j,k,l,m
    real(kind=rkind),parameter :: delta(3,3) = &
        reshape((/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/),&
		shape(delta))
    real(kind=rkind),pointer:: Dx(:,:),Dy(:,:),Dz(:,:),Dtmp(:)
    Dtmp=>tmp(1:pNp**2:pNp+1,1)
    Dx=>DD(:,:,1);Dy=>DD(:,:,2);Dz=>DD(:,:,3)

    DD(:,:,1)=iJac(1)*Dr+iJac(2)*Ds+iJac(3)*Dt
    DD(:,:,2)=iJac(4)*Dr+iJac(5)*Ds+iJac(6)*Dt
    DD(:,:,3)=iJac(7)*Dr+iJac(8)*Ds+iJac(9)*Dt

    T0(:,1,1)=T0xx;T0(:,2,3)=T0yz;T0(:,3,2)=T0yz
    T0(:,2,2)=T0yy;T0(:,1,3)=T0xz;T0(:,3,1)=T0xz
    T0(:,3,3)=T0zz;T0(:,1,2)=T0xy;T0(:,2,1)=T0xy

    n(1)=nx;n(2)=ny;n(3)=nz

!!!! math formulation
!!!! term 1:   0.5 N_k T_ik D_j U_j
!!!! term 2: + 0.5 N_i T_jk D_k U_j
!!!! term 3: + 0.5 N_k T_kl D_l U_i
!!!! term 4: - 0.5 N_k T_ij D_k U_j
!!!! term 5: - 0.5 N_j T_ik D_k U_j
!!!! term 6: - 0.5 N_k T_jk D_i U_j
!!!! term 7: - D^S_j ( n_k T_ik U_j )
!!!! 

!!!! term 3
    tmp=0d0;tmp1=0d0
    do l=1,3
    do k=1,3
        Dtmp=n(k)*T0(:,k,l)
        tmp1=tmp1+0.5d0*matmul(tmp,DD(:,:,l))
    enddo
    enddo
    do i=1,3
    do j=1,3
!!!! term 1
        Dtmp=0d0
        do k=1,3
            Dtmp=Dtmp+n(k)*T0(:,k,i)
        enddo
        tmp2=0.5d0*matmul(tmp,DD(:,:,j))
!!!! term 2,4,5
        do k=1,3
            Dtmp=n(i)*T0(:,k,j)-n(k)*T0(:,i,j)-n(j)*T0(:,i,k)
            tmp2=tmp2+0.5d0*matmul(tmp,DD(:,:,k))
        enddo
!!!! term 3
        if(i.eq.j)then
            tmp2=tmp2+tmp1
        endif
!!!! term 6
        Dtmp=0d0
        do k=1,3
            Dtmp=Dtmp+n(k)*T0(:,k,j)
            tmp2=tmp2-0.5d0*matmul(tmp,DD(:,:,i))
        enddo
!!!! term 7
        tmp3=DD(:,:,j)
        do l=1,3
            tmp3=tmp3-n(j)*n(l)*DD(:,:,l)
        enddo
        do k=1,3
            Dtmp=n(k)*T0(:,k,i)
            tmp2=tmp2-matmul(tmp3,tmp)
        enddo
!!!! asign to matrix
!        print*,i,j,&
!        maxval(abs(A((i-1)*Np+1:i*Np,(j-1)*pNp+1:j*pNp))),&
!        maxval(abs(matmul(irJL,tmp2)))
        A((i-1)*Np+1:i*Np,(j-1)*pNp+1:j*pNp)=&
        A((i-1)*Np+1:i*Np,(j-1)*pNp+1:j*pNp)-matmul(irJL,tmp2)
    enddo
    enddo
!    stop

end subroutine form_rupt_grav_blocMatrix

subroutine rupt_grav_blocMatrix(pNp,Nfp,mesh,rupt,matrix,T0)
    integer                :: pNp,Nfp,iface,itet,itri,head,tail
    type(tetmesh_geometry) :: mesh
    type(rupture)          :: rupt
    type(matrices)         :: matrix
    type(tensor_array)     :: T0
    real(kind=rkind) :: irJL(pNp,pNp)
    real(kind=rkind),target  :: tmp(Nfp,Nfp)
    real(kind=rkind),pointer :: Dtmp(:)
    integer :: map(Nfp)
    real(kind=rkind) :: tmp1,tmp2
    Dtmp=>tmp(1:Nfp**2:Nfp+1,1)
    tmp=0d0
    do iface=1,rupt%Nhface
        rupt%M(iface)%g0m=0d0
        rupt%M(iface)%g0p=0d0
        itet=rupt%T2E(1,iface);itri=rupt%T2E(3,iface)
        head=(itet-1)*pNp+1;tail=itet*pNp
        call form_rupt_grav_blocMatrix(Nfp,pNp,&
            matrix%Dr,matrix%Ds,matrix%Dt,&
            mesh%invJac(itet,1:9),rupt%M(iface)%Lm,&
            mesh%nx(itri,itet),mesh%ny(itri,itet),mesh%nz(itri,itet),&
            T0%xx(head:tail),&
            T0%yy(head:tail),&
            T0%zz(head:tail),&
            T0%yz(head:tail),&
            T0%xz(head:tail),&
            T0%xy(head:tail),&
            rupt%M(iface)%g0m)
!        tmp1=mesh%nz(itri,itet)
        itet=rupt%T2E(2,iface);itri=rupt%T2E(4,iface)
        head=(itet-1)*pNp+1
        call form_rupt_grav_blocMatrix(Nfp,pNp,&
            matrix%Dr,matrix%Ds,matrix%Dt,&
            mesh%invJac(itet,1:9),rupt%M(iface)%Lp,&
            mesh%nx(itri,itet),mesh%ny(itri,itet),mesh%nz(itri,itet),&
            T0%xx(head:tail),&
            T0%yy(head:tail),&
            T0%zz(head:tail),&
            T0%yz(head:tail),&
            T0%xz(head:tail),&
            T0%xy(head:tail),&
            rupt%M(iface)%g0p)
!        tmp2=mesh%nz(itri,itet)
!        print*,tmp1,tmp2
!        stop
    enddo
end subroutine rupt_grav_blocMatrix

!--------------------------------------------------------------------

end module grav_mod

