!*******************************************************************!
!*  This module obtains semidiscretized wave equation              *!
!*                                                                 *!
!*  Author: Ruichao Ye    Email: Ruichao.ye@gmail.com              *!
!*******************************************************************!

!--------------------------------------------------------------------
module solve_mod
!--------------------------------------------------------------------

    use para_mod,     only : def_pOrder,pNp,Nfp,Nfp4,&
                             rpml,convtest,withrupt,rupt_gamma
    use datatype_mod, only : rkind,matrices,sources,rupture,&
                             tetmesh_geometry,tetmesh_domain,&
                             tetmesh_material,&
                             vector_array,tensor_array,wavefield,&
                             auxilary_array,reset_vector,&
                             reset_tensor,vector_axpy,tensor_axpy,&
                             reset_wavefield,reset_auxilary_array,&
                             BLK_axpy,DiagMM,Array_permute
    use ruptsolve_mod,only : rupt_on,rupt_alpha_V,rupt_alpha_E

!--------------------------------------------------------------------

    implicit none

    private :: strain_stress0,strain_stress1,&
               strain_stress2,strain_stress3

    public :: bilinear_Elastic3D,strain_stress,penalty_flux

!--------------------------------------------------------------------
    contains
!--------------------------------------------------------------------

subroutine bilinear_Elastic3D(matrix,mesh,subdomain,material,srcs,&
        rhsW,divW,VS_W,localtime)
! inout:
    type(matrices)         :: matrix
    type(tetmesh_geometry) :: mesh
    type(tetmesh_Domain)   :: subdomain
    type(tetmesh_material) :: material
    type(sources)          :: srcs
    type(wavefield)        :: rhsW,VS_W
    type(auxilary_array)   :: divW
    real(kind=rkind),intent(in) :: localtime
! auxilary:
    integer :: head,tail,fhead,fh,offset,i,j,k,iele,iface,ii
    integer :: ierr,Ndim,Ndof,NdofE,NdofV
    logical :: symm
    real(kind=rkind),pointer :: invJ(:,:)

!!!!!!!!!!!!!!!! compute volume terms !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    NdofV=mesh%Nele*3

    symm=rhsW%E%symmetric
    if(symm)then
        NdofE=mesh%Nele*6
        Ndof=mesh%Nele*9
    else
        NdofE=mesh%Nele*9
        Ndof=mesh%Nele*12
    endif

    Ndim=mesh%Nhele
    invJ=>mesh%invJac

!    call reset_wavefield(rhsW)
!    if(rupt_gamma.gt.0)then
!        call vector_AXPY(Ndim*pNp,-rupt_gamma,VS_W%V,rhsW%V)
!    endif

    call DGEMM('n','n',pNp,Ndof,pNp,1d0,&
        matrix%Dr,pNp,VS_W%array,pNp,0d0,divW%array,pNp)

    call DiagMM('r',Ndim,pNp,1d0,invJ(:,1),divW%V%x,rhsW%E%xx,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,1),divW%V%y,rhsW%E%xy,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,1),divW%V%z,rhsW%E%xz,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,4),divW%V%x,rhsW%E%yx,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,4),divW%V%y,rhsW%E%yy,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,4),divW%V%z,rhsW%E%yz,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,7),divW%V%x,rhsW%E%zx,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,7),divW%V%y,rhsW%E%zy,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,7),divW%V%z,rhsW%E%zz,pNp,pNp)
                                                                     
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,1),divW%S%xx,rhsW%V%x,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,1),divW%S%xy,rhsW%V%y,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,1),divW%S%xz,rhsW%V%z,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,4),divW%S%yx,rhsW%V%x,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,4),divW%S%yy,rhsW%V%y,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,4),divW%S%yz,rhsW%V%z,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,7),divW%S%zx,rhsW%V%x,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,7),divW%S%zy,rhsW%V%y,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,7),divW%S%zz,rhsW%V%z,pNp,pNp)
                                                                     
    call DGEMM('n','n',pNp,Ndof,pNp,1d0,&                            
        matrix%Ds,pNp,VS_W%array,pNp,0d0,divW%array,pNp)             
                                                                     
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,2),divW%V%x,rhsW%E%xx,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,2),divW%V%y,rhsW%E%xy,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,2),divW%V%z,rhsW%E%xz,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,5),divW%V%x,rhsW%E%yx,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,5),divW%V%y,rhsW%E%yy,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,5),divW%V%z,rhsW%E%yz,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,8),divW%V%x,rhsW%E%zx,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,8),divW%V%y,rhsW%E%zy,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,8),divW%V%z,rhsW%E%zz,pNp,pNp)
                                                                     
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,2),divW%S%xx,rhsW%V%x,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,2),divW%S%xy,rhsW%V%y,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,2),divW%S%xz,rhsW%V%z,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,5),divW%S%yx,rhsW%V%x,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,5),divW%S%yy,rhsW%V%y,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,5),divW%S%yz,rhsW%V%z,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,8),divW%S%zx,rhsW%V%x,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,8),divW%S%zy,rhsW%V%y,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,8),divW%S%zz,rhsW%V%z,pNp,pNp)
                                                                     
    call DGEMM('n','n',pNp,Ndof,pNp,1d0,&                            
        matrix%Dt,pNp,VS_W%array,pNp,0d0,divW%array,pNp)             
                                                                     
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,3),divW%V%x,rhsW%E%xx,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,3),divW%V%y,rhsW%E%xy,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,3),divW%V%z,rhsW%E%xz,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,6),divW%V%x,rhsW%E%yx,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,6),divW%V%y,rhsW%E%yy,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,6),divW%V%z,rhsW%E%yz,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,9),divW%V%x,rhsW%E%zx,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,9),divW%V%y,rhsW%E%zy,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,9),divW%V%z,rhsW%E%zz,pNp,pNp)
                                                                     
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,3),divW%S%xx,rhsW%V%x,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,3),divW%S%xy,rhsW%V%y,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,3),divW%S%xz,rhsW%V%z,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,6),divW%S%yx,rhsW%V%x,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,6),divW%S%yy,rhsW%V%y,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,6),divW%S%yz,rhsW%V%z,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,9),divW%S%zx,rhsW%V%x,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,9),divW%S%zy,rhsW%V%y,pNp,pNp)
    call DiagMM('r',Ndim,pNp,1d0,invJ(:,9),divW%S%zz,rhsW%V%z,pNp,pNp)

    if(symm)call DSCAL(mesh%Nele*pNp*3,0.5d0,rhsW%E%yz,1)

    if(.true. .and. rupt_on .and. mesh%rupt%Nhface.gt.0)then
        do j=1,mesh%rupt%Nhface
            do k=1,2
                i=mesh%rupt%T2E(k,j)
                head=(i-1)*pNp+1;tail=i*pNp
                rhsW%V%x(head:tail)=0d0
                rhsW%V%y(head:tail)=0d0
                rhsW%V%z(head:tail)=0d0
                rhsW%E%xx(head:tail)=0d0
                rhsW%E%yx(head:tail)=0d0
                rhsW%E%zx(head:tail)=0d0
                rhsW%E%xy(head:tail)=0d0
                rhsW%E%yy(head:tail)=0d0
                rhsW%E%zy(head:tail)=0d0
                rhsW%E%xz(head:tail)=0d0
                rhsW%E%yz(head:tail)=0d0
                rhsW%E%zz(head:tail)=0d0
            enddo
        enddo
!open(1111,file='mat.dat',status='replace')
!write(1111,*)rhsW%V%array
!close(1111)
!stop
    endif

!!!!!!!!!!!!!!!! compute flux !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call Array_permute(pNp,4*Nfp,Ndof,VS_W%array,divW%array,mesh%vmapM)

    if(symm)then
        do i=1,mesh%Nhele
            ! Further than first arrival
            if(localtime.lt.srcs%t_direct(i))cycle

            head=(i-1)*pNp+1;offset=(i-1)*Nfp*4
            do j=1,4
                fhead=(j-1)*Nfp+1
                fh=fhead+offset
                call penalty_flux(Nfp,&
                    mesh%fbctype(:,j,i),material%k_media(i),&
                    mesh%nx(j,i),mesh%ny(j,i),mesh%nz(j,i),&
                    mesh%Fscale(j,i),mesh%vmapP(fh:),&
                    VS_W%V%x        ,VS_W%V%y        ,VS_W%V%z        ,&
                    VS_W%S%xx       ,VS_W%S%yy       ,VS_W%S%zz       ,&
                    VS_W%S%yz       ,VS_W%S%xz       ,VS_W%S%xy       ,&
                    divW%fV%x (fh:) ,divW%fV%y (fh:) ,divW%fV%z (fh:) ,&
                    divW%fE%xx(fh:) ,divW%fE%yy(fh:) ,divW%fE%zz(fh:) ,&
                    divW%fE%yz(fh:) ,divW%fE%xz(fh:) ,divW%fE%xy(fh:) ,&
                    mesh)
            enddo
        enddo
    else
        do i=1,mesh%Nhele
            ! Further than first arrival
            if(localtime.lt.srcs%t_direct(i))cycle

            head=(i-1)*pNp+1;offset=(i-1)*Nfp*4
            do j=1,4
                fhead=(j-1)*Nfp+1
                fh=fhead+offset
                call penalty_flux_asym(Nfp,&
                    mesh%fbctype(:,j,i),material%k_media(i),&
                    mesh%nx(j,i),mesh%ny(j,i),mesh%nz(j,i),&
                    mesh%Fscale(j,i),mesh%vmapP(fh:),&
                    VS_W%V%x        ,VS_W%V%y        ,VS_W%V%z        ,&
                    VS_W%S%xx       ,VS_W%S%yy       ,VS_W%S%zz       ,&
                    VS_W%S%yz       ,VS_W%S%xz       ,VS_W%S%xy       ,&
                    VS_W%S%zy       ,VS_W%S%zx       ,VS_W%S%yx       ,&
                    divW%fV%x (fh:) ,divW%fV%y (fh:) ,divW%fV%z (fh:) ,&
                    divW%fE%xx(fh:) ,divW%fE%yy(fh:) ,divW%fE%zz(fh:) ,&
                    divW%fE%yz(fh:) ,divW%fE%xz(fh:) ,divW%fE%xy(fh:) ,&
                    divW%fE%zy(fh:) ,divW%fE%zx(fh:) ,divW%fE%yx(fh:) ,&
                    mesh)
            enddo
        enddo
    endif

    call DGEMM('n','n',pNp,Ndof,4*Nfp,1d0,&
        matrix%LIFT,pNp,divW%array,4*Nfp,1d0,rhsW%array,pNp)

    rhsW%V%x=rhsW%V%x/material%rho
    rhsW%V%y=rhsW%V%y/material%rho
    rhsW%V%z=rhsW%V%z/material%rho

end subroutine bilinear_Elastic3D

!--------------------------------------------------------------------

subroutine strain_stress(pNp,Nele,mat,E,S)
    type(tensor_array) :: E,S
    integer,intent(in) :: pNp,Nele
    type(tetmesh_material) :: mat
    logical :: symm
    symm=E%symmetric

    if(mat%job.eq.0)then
        call strain_stress0(Nele*pNp,&
            mat%C(:,1),&
            E%xx,E%yy,E%zz,S%xx,S%yy,S%zz)
    elseif(mat%job.eq.1)then
        if(symm)then
            call strain_stress1(Nele*pNp,&
                mat%C(:,1),mat%C(:,2),&
                E%xx,E%yy,E%zz,E%yz,E%xz,E%xy,&
                S%xx,S%yy,S%zz,S%yz,S%xz,S%xy)
        else
            call strain_stress_asym1(Nele*pNp,&
                mat%C(:,1),mat%C(:,2),&
                E%xx,E%yy,E%zz,E%yz,E%xz,E%xy,E%zy,E%zx,E%yx,&
                S%xx,S%yy,S%zz,S%yz,S%xz,S%xy,S%zy,S%zx,S%yx)
        endif
    elseif(mat%job.eq.2)then
        if(symm)then
            call strain_stress2(Nele*pNp,&
                mat%C(1:,1),mat%C(1:,2),mat%C(1:,3),&
                mat%C(1:,4),mat%C(1:,5),mat%C(1:,6),&
                mat%C(1:,7),mat%C(1:,8),mat%C(1:,9),&
                E%xx,E%yy,E%zz,E%yz,E%xz,E%xy,&
                S%xx,S%yy,S%zz,S%yz,S%xz,S%xy)
        else
            call strain_stress_asym2(Nele*pNp,&
                mat%C(1:,1),mat%C(1:,2),mat%C(1:,3),&
                mat%C(1:,4),mat%C(1:,5),mat%C(1:,6),&
                mat%C(1:,7),mat%C(1:,8),mat%C(1:,9),&
                E%xx,E%yy,E%zz,E%yz,E%xz,E%xy,E%zy,E%zx,E%yx,&
                S%xx,S%yy,S%zz,S%yz,S%xz,S%xy,S%zy,S%zx,S%yx)
        endif
    elseif(mat%job.eq.3)then
        if(symm)then
            call strain_stress3(Nele*pNp,&
                mat%C(1:, 1),mat%C(1:, 2),mat%C(1:, 3),&
                mat%C(1:, 4),mat%C(1:, 5),mat%C(1:, 6),&
                mat%C(1:, 7),mat%C(1:, 8),mat%C(1:, 9),&
                mat%C(1:,10),mat%C(1:,11),mat%C(1:,12),&
                mat%C(1:,13),mat%C(1:,14),mat%C(1:,15),&
                mat%C(1:,16),mat%C(1:,17),mat%C(1:,18),&
                mat%C(1:,19),mat%C(1:,20),mat%C(1:,21),&
                E%xx,E%yy,E%zz,E%yz,E%xz,E%xy,&
                S%xx,S%yy,S%zz,S%yz,S%xz,S%xy)
        else
            call strain_stress_asym3(Nele*pNp,&
                mat%C(1:, 1),mat%C(1:, 2),mat%C(1:, 3),&
                mat%C(1:, 4),mat%C(1:, 5),mat%C(1:, 6),&
                mat%C(1:, 7),mat%C(1:, 8),mat%C(1:, 9),&
                mat%C(1:,10),mat%C(1:,11),mat%C(1:,12),&
                mat%C(1:,13),mat%C(1:,14),mat%C(1:,15),&
                mat%C(1:,16),mat%C(1:,17),mat%C(1:,18),&
                mat%C(1:,19),mat%C(1:,20),mat%C(1:,21),&
                E%xx,E%yy,E%zz,E%yz,E%xz,E%xy,E%zy,E%zx,E%yx,&
                S%xx,S%yy,S%zz,S%yz,S%xz,S%xy,S%zy,S%zx,S%yx)
        endif
    endif

end subroutine strain_stress

!--------------------------------------------------------------------

subroutine penalty_flux(&
        Nfp,fbctype,k_media,nx,ny,nz,Fscale,mapP,&
         v1 ,v2 ,v3 ,S11 ,S22 ,S33 ,S23 ,S13 ,S12,&
        fv1,fv2,fv3,fE11,fE22,fE33,fE23,fE13,fE12,&
        mesh)
! input
    integer,intent(in) :: Nfp,k_media,fbctype(2)
    real(kind=rkind),intent(in) :: nx,ny,nz,Fscale
    type(tetmesh_geometry) :: mesh
    real(kind=rkind),intent(in) :: v1(1),v2(1),v3(1),&
            S11(1),S22(1),S33(1),S23(1),S13(1),S12(1)
    integer,intent(in) :: mapP(Nfp)
! output
    real(kind=rkind) ::  fv1(Nfp), fv2(Nfp), fv3(Nfp),&
                        fE11(Nfp),fE22(Nfp),fE33(Nfp),&
                        fE23(Nfp),fE13(Nfp),fE12(Nfp)
! auxilary
    real(kind=rkind) :: &
               S11p(Nfp), S22p(Nfp), S33p(Nfp),&
               S23p(Nfp), S13p(Nfp), S12p(Nfp),&
                dv1(Nfp),  dv2(Nfp),  dv3(Nfp),&
               dSn1(Nfp), dSn2(Nfp), dSn3(Nfp),&
               Sn1m(Nfp), Sn2m(Nfp), Sn3m(Nfp),&
               ndSn(Nfp),  ndv(Nfp)
    real(kind=rkind),parameter :: alpha=0.5d0
    integer :: j,fhead,ftail

    Sn1m = fE11*nx + fE12*ny + fE13*nz 
    Sn2m = fE12*nx + fE22*ny + fE23*nz 
    Sn3m = fE13*nx + fE23*ny + fE33*nz 
    if(fbctype(1).eq.0)then 
        dv1 = v1(mapP) - fv1 
        dv2 = v2(mapP) - fv2 
        dv3 = v3(mapP) - fv3
        S11p = S11(mapP); S22p = S22(mapP); S33p = S33(mapP)
        S23p = S23(mapP); S13p = S13(mapP); S12p = S12(mapP)
        dSn1 = S11p*nx + S12p*ny + S13p*nz - Sn1m 
        dSn2 = S12p*nx + S22p*ny + S23p*nz - Sn2m 
        dSn3 = S13p*nx + S23p*ny + S33p*nz - Sn3m
        if(k_media.ge.1)then ! solid media
            if(fbctype(2).eq.1)then ! solid-fluid interface
                ndv = dv1*nx+dv2*ny+dv3*nz 
                dv1 = ndv*nx
                dv2 = ndv*ny
                dv3 = ndv*nz
                ndSn=dSn1*nx+dSn2*ny+dSn3*nz
                dSn1 = 2d0*dSn1-ndSn*nx
                dSn2 = 2d0*dSn2-ndSn*ny
                dSn3 = 2d0*dSn3-ndSn*nz
            endif
        else ! fluid media
            if(fbctype(2).eq.1)then ! solid-fluid interface
                ndv = dv1*nx+dv2*ny+dv3*nz 
                dv1 = ndv*nx
                dv2 = ndv*ny
                dv3 = ndv*nz
                ndSn=dSn1*nx+dSn2*ny+dSn3*nz
                dSn1 = ndSn*nx
                dSn2 = ndSn*ny
                dSn3 = ndSn*nz
            endif
        endif
        fE11 = 0.5d0* dv1 + alpha*dSn1 + fv1 
        fE22 = 0.5d0* dv2 + alpha*dSn2 + fv2 
        fE33 = 0.5d0* dv3 + alpha*dSn3 + fv3 
         fv1 = 0.5d0*dSn1 + alpha* dv1 + Sn1m 
         fv2 = 0.5d0*dSn2 + alpha* dv2 + Sn2m 
         fv3 = 0.5d0*dSn3 + alpha* dv3 + Sn3m 
    elseif(fbctype(1).eq.1)then
        if(fbctype(2).eq.1)then ! Free surface boundary
            dSn1=-Sn1m; dv1=dSn1
            dSn2=-Sn2m; dv2=dSn2
            dSn3=-Sn3m; dv3=dSn3
        elseif(fbctype(2).eq.2)then ! Absrobing boundary
            dv1=-Sn1m-fv1; dSn1=0d0
            dv2=-Sn2m-fv2; dSn2=0d0
            dv3=-Sn3m-fv3; dSn3=0d0
        endif
        if(k_media.ge.1)then ! solid media
            fE11 = 0.5d0* dv1 + alpha*dSn1 + fv1 
            fE22 = 0.5d0* dv2 + alpha*dSn2 + fv2 
            fE33 = 0.5d0* dv3 + alpha*dSn3 + fv3 
             fv1 = 0.5d0*dSn1 + alpha* dv1 + Sn1m 
             fv2 = 0.5d0*dSn2 + alpha* dv2 + Sn2m 
             fv3 = 0.5d0*dSn3 + alpha* dv3 + Sn3m 
        else ! fluid media
             dSn1=dSn1*nx+dSn2*ny+dSn3*nz
              ndv= dv1*nx+ dv2*ny+ dv3*nz
             ndSn=0.5d0*dSn1+alpha*ndv
              ndv=0.5d0* ndv+alpha*dSn1
            fE11 =  ndv*nx + fv1 
            fE22 =  ndv*ny + fv2 
            fE33 =  ndv*nz + fv3 
             fv1 = ndSn*nx + Sn1m 
             fv2 = ndSn*ny + Sn2m 
             fv3 = ndSn*nz + Sn3m 
        endif
    elseif(fbctype(1).eq.2)then ! rupture
        fE11 = 0d0
        fE22 = 0d0
        fE33 = 0d0
         fv1 = 0d0
         fv2 = 0d0
         fv3 = 0d0
    elseif(fbctype(1).eq.3)then ! external Neumann surface source
        j=fbctype(2)
        fhead=(j-1)*Nfp+1;ftail=j*Nfp
        dSn1=mesh%surf%Sn%x(fhead:ftail)-Sn1m;dv1=0d0
        dSn2=mesh%surf%Sn%y(fhead:ftail)-Sn2m;dv2=0d0
        dSn3=mesh%surf%Sn%z(fhead:ftail)-Sn3m;dv3=0d0
        fE11 = 0.5d0* dv1 + alpha*dSn1 + fv1 
        fE22 = 0.5d0* dv2 + alpha*dSn2 + fv2 
        fE33 = 0.5d0* dv3 + alpha*dSn3 + fv3 
         fv1 = 0.5d0*dSn1 + alpha* dv1 + Sn1m 
         fv2 = 0.5d0*dSn2 + alpha* dv2 + Sn2m 
         fv3 = 0.5d0*dSn3 + alpha* dv3 + Sn3m 
    endif

    fv1  = Fscale * fv1
    fv2  = Fscale * fv2
    fv3  = Fscale * fv3
    fE11 = Fscale * fE11
    fE22 = Fscale * fE22
    fE33 = Fscale * fE33
    fE23 = 0.5d0 *( ny*fE33 + nz*fE22 )
    fE13 = 0.5d0 *( nx*fE33 + nz*fE11 )
    fE12 = 0.5d0 *( nx*fE22 + ny*fE11 )
    fE11 = fE11*nx
    fE22 = fE22*ny
    fE33 = fE33*nz

end subroutine penalty_flux

subroutine penalty_flux_asym(&
        Nfp,fbctype,k_media,nx,ny,nz,Fscale,mapP,&
         v1 ,v2 ,v3 ,S11 ,S22 ,S33 ,S23 ,S13 ,S12 ,S32 ,S31 ,S21,&
        fv1,fv2,fv3,fE11,fE22,fE33,fE23,fE13,fE12,fE32,fE31,fE21,&
        mesh)
! input
    integer,intent(in) :: Nfp,k_media,fbctype(2)
    real(kind=rkind),intent(in) :: nx,ny,nz,Fscale
    type(tetmesh_geometry) :: mesh
    real(kind=rkind),intent(in) :: v1(1),v2(1),v3(1),&
            S11(1),S22(1),S33(1),&
            S23(1),S13(1),S12(1),&
            S32(1),S31(1),S21(1)
    integer,intent(in) :: mapP(Nfp)
! output
    real(kind=rkind) ::  fv1(Nfp), fv2(Nfp), fv3(Nfp),&
                        fE11(Nfp),fE22(Nfp),fE33(Nfp),&
                        fE23(Nfp),fE13(Nfp),fE12(Nfp),&
                        fE32(Nfp),fE31(Nfp),fE21(Nfp)
! auxilary
    real(kind=rkind) :: &
                dv1(Nfp),  dv2(Nfp),  dv3(Nfp),&
               dSn1(Nfp), dSn2(Nfp), dSn3(Nfp),&
               Sn1m(Nfp), Sn2m(Nfp), Sn3m(Nfp),&
               ndSn(Nfp),  ndv(Nfp)
    real(kind=rkind),parameter :: alpha=0.1d0
    integer :: i,j,k,fhead,ftail,flag,map(Nfp)

    Sn1m = fE11*nx + fE21*ny + fE31*nz 
    Sn2m = fE12*nx + fE22*ny + fE32*nz 
    Sn3m = fE13*nx + fE23*ny + fE33*nz 
    if(fbctype(1).eq.0)then 
        dv1 = v1(mapP) - fv1 
        dv2 = v2(mapP) - fv2 
        dv3 = v3(mapP) - fv3
        dSn1 = S11(mapP)*nx + S21(mapP)*ny + S31(mapP)*nz - Sn1m 
        dSn2 = S12(mapP)*nx + S22(mapP)*ny + S32(mapP)*nz - Sn2m 
        dSn3 = S13(mapP)*nx + S23(mapP)*ny + S33(mapP)*nz - Sn3m
        if(k_media.ge.1)then ! solid media
            if(fbctype(2).eq.1)then ! solid-fluid interface
                ndv = dv1*nx+dv2*ny+dv3*nz 
                dv1 = ndv*nx
                dv2 = ndv*ny
                dv3 = ndv*nz
            endif
            fE11 = 0.5d0* dv1 + alpha*dSn1 + fv1 
            fE12 = 0.5d0* dv2 + alpha*dSn2 + fv2 
            fE13 = 0.5d0* dv3 + alpha*dSn3 + fv3 
             fv1 = 0.5d0*dSn1 + alpha* dv1 + Sn1m 
             fv2 = 0.5d0*dSn2 + alpha* dv2 + Sn2m 
             fv3 = 0.5d0*dSn3 + alpha* dv3 + Sn3m 
        else ! fluid media
             dSn1=dSn1*nx+dSn2*ny+dSn3*nz
              ndv= dv1*nx+ dv2*ny+ dv3*nz
             ndSn=0.5d0*dSn1+alpha*ndv
              ndv=0.5d0* ndv+alpha*dSn1
            fE11 =  ndv*nx + fv1 
            fE12 =  ndv*ny + fv2 
            fE13 =  ndv*nz + fv3 
             fv1 = ndSn*nx + Sn1m 
             fv2 = ndSn*ny + Sn2m 
             fv3 = ndSn*nz + Sn3m 
        endif
    elseif(fbctype(1).eq.1)then
        if(fbctype(2).eq.1)then ! Free surface boundary
            dSn1=-Sn1m; dv1=dSn1
            dSn2=-Sn2m; dv2=dSn2
            dSn3=-Sn3m; dv3=dSn3
        elseif(fbctype(2).eq.2)then ! Absrobing boundary
            dv1=-Sn1m-fv1; dSn1=0d0
            dv2=-Sn2m-fv2; dSn2=0d0
            dv3=-Sn3m-fv3; dSn3=0d0
        endif
        if(k_media.ge.1)then ! solid media
            fE11 = 0.5d0* dv1 + alpha*dSn1 + fv1 
            fE12 = 0.5d0* dv2 + alpha*dSn2 + fv2 
            fE13 = 0.5d0* dv3 + alpha*dSn3 + fv3 
             fv1 = 0.5d0*dSn1 + alpha* dv1 + Sn1m 
             fv2 = 0.5d0*dSn2 + alpha* dv2 + Sn2m 
             fv3 = 0.5d0*dSn3 + alpha* dv3 + Sn3m 
        else ! fluid media
             dSn1=dSn1*nx+dSn2*ny+dSn3*nz
              ndv= dv1*nx+ dv2*ny+ dv3*nz
             ndSn=0.5d0*dSn1+alpha*ndv
              ndv=0.5d0* ndv+alpha*dSn1
            fE11 =  ndv*nx + fv1 
            fE12 =  ndv*ny + fv2 
            fE13 =  ndv*nz + fv3 
             fv1 = ndSn*nx + Sn1m 
             fv2 = ndSn*ny + Sn2m 
             fv3 = ndSn*nz + Sn3m 
        endif
    elseif(fbctype(1).eq.2)then ! rupture
        if(rupt_on)then
            if(convtest)then
                k=fbctype(2)
                if(k.gt.0)then
                    map=mesh%rupt%perm(:,1,k)
                    flag=1
                else
                    k=-k
                    map=mesh%rupt%perm(:,2,k)
                    flag=-1
                endif
                do i=1,Nfp
                    j=(k-1)*Nfp+map(i)
                    dv1(i) = mesh%rupt%Vt%x(j)
                    dv2(i) = mesh%rupt%Vt%y(j)
                    dv3(i) = mesh%rupt%Vt%z(j)
                enddo
                fE11 = flag * 0.5d0*dv1 
                fE12 = flag * 0.5d0*dv2 
                fE13 = flag * 0.5d0*dv3 
                 fv1 = rupt_alpha_V * dv1 
                 fv2 = rupt_alpha_V * dv2 
                 fv3 = rupt_alpha_V * dv3 
            else
                fE11 = 0d0
                fE12 = 0d0
                fE13 = 0d0
                 fv1 = 0d0
                 fv2 = 0d0
                 fv3 = 0d0
            endif
         else
            dv1 = v1(mapP) - fv1 
            dv2 = v2(mapP) - fv2 
            dv3 = v3(mapP) - fv3
            dSn1 = S11(mapP)*nx + S21(mapP)*ny + S31(mapP)*nz - Sn1m 
            dSn2 = S12(mapP)*nx + S22(mapP)*ny + S32(mapP)*nz - Sn2m 
            dSn3 = S13(mapP)*nx + S23(mapP)*ny + S33(mapP)*nz - Sn3m
            fE11 = 0.5d0* dv1 + rupt_alpha_E*dSn1 + fv1 
            fE12 = 0.5d0* dv2 + rupt_alpha_E*dSn2 + fv2 
            fE13 = 0.5d0* dv3 + rupt_alpha_E*dSn3 + fv3 
             fv1 = 0.5d0*dSn1 + rupt_alpha_V*dv1 + Sn1m 
             fv2 = 0.5d0*dSn2 + rupt_alpha_V*dv2 + Sn2m 
             fv3 = 0.5d0*dSn3 + rupt_alpha_V*dv3 + Sn3m 
         endif
    elseif(fbctype(1).eq.3)then ! external Neumann surface source
        j=fbctype(2)
        fhead=(j-1)*Nfp+1;ftail=j*Nfp
        dSn1=mesh%surf%Sn%x(fhead:ftail)-Sn1m;dv1=0d0
        dSn2=mesh%surf%Sn%y(fhead:ftail)-Sn2m;dv2=0d0
        dSn3=mesh%surf%Sn%z(fhead:ftail)-Sn3m;dv3=0d0
        fE11 = 0.5d0* dv1 + alpha*dSn1 + fv1 
        fE12 = 0.5d0* dv2 + alpha*dSn2 + fv2 
        fE13 = 0.5d0* dv3 + alpha*dSn3 + fv3 
         fv1 = 0.5d0*dSn1 + alpha* dv1 + Sn1m 
         fv2 = 0.5d0*dSn2 + alpha* dv2 + Sn2m 
         fv3 = 0.5d0*dSn3 + alpha* dv3 + Sn3m 
!print*,dSn1
    endif

    fv1  = Fscale * fv1
    fv2  = Fscale * fv2
    fv3  = Fscale * fv3
    fE11 = Fscale * fE11
    fE12 = Fscale * fE12
    fE13 = Fscale * fE13
    fE21 = ny*fE11
    fE22 = ny*fE12
    fE23 = ny*fE13
    fE31 = nz*fE11
    fE32 = nz*fE12
    fE33 = nz*fE13
    fE11 = nx*fE11
    fE12 = nx*fE12
    fE13 = nx*fE13

end subroutine penalty_flux_asym

!--------------------------------------------------------------------

subroutine strain_stress0(Ndim,lambda,&
        E11,E22,E33,S11,S22,S33)
    integer,intent(in) :: Ndim
    real(kind=rkind),intent(in)  :: lambda(Ndim)
    real(kind=rkind),intent(in)  :: E11(Ndim),E22(Ndim),E33(Ndim)
    real(kind=rkind),intent(out) :: S11(Ndim),S22(Ndim),S33(Ndim)
    S33=lambda*(E11+E22+E33)
    S11=S33
    S22=S33
end subroutine strain_stress0

subroutine strain_stress1(Ndim,lambda,muX2,&
        E11,E22,E33,E23,E13,E12,S11,S22,S33,S23,S13,S12)
    integer,intent(in) :: Ndim
    real(kind=rkind),intent(in)  :: lambda(Ndim),muX2(Ndim)
    real(kind=rkind),intent(in)  :: E11(Ndim),E22(Ndim),E33(Ndim)
    real(kind=rkind),intent(in)  :: E23(Ndim),E13(Ndim),E12(Ndim)
    real(kind=rkind),intent(out) :: S11(Ndim),S22(Ndim),S33(Ndim)
    real(kind=rkind),intent(out) :: S23(Ndim),S13(Ndim),S12(Ndim)
    S33=lambda*(E11+E22+E33)
    S11=S33+muX2*E11
    S22=S33+muX2*E22
    S33=S33+muX2*E33
    S23=muX2*E23
    S13=muX2*E13
    S12=muX2*E12
end subroutine strain_stress1

subroutine strain_stress2(Ndim,C11,C22,C33,&
        C44,C55,C66,C23,C13,C12,&
        E11,E22,E33,E23,E13,E12,S11,S22,S33,S23,S13,S12)
    integer,intent(in) :: Ndim
    real(kind=rkind),intent(in)  :: &
        C11(Ndim),C22(Ndim),C33(Ndim),&
        C23(Ndim),C13(Ndim),C12(Ndim),&
        C44(Ndim),C55(Ndim),C66(Ndim)
    real(kind=rkind),intent(in)  :: E11(Ndim),E22(Ndim),E33(Ndim)
    real(kind=rkind),intent(in)  :: E23(Ndim),E13(Ndim),E12(Ndim)
    real(kind=rkind),intent(out) :: S11(Ndim),S22(Ndim),S33(Ndim)
    real(kind=rkind),intent(out) :: S23(Ndim),S13(Ndim),S12(Ndim)
    S11=C11*E11+C12*E22+C13*E33
    S22=C12*E11+C22*E22+C23*E33
    S33=C13*E11+C23*E22+C33*E33
    S23=2d0*C44*E23
    S13=2d0*C55*E13
    S12=2d0*C66*E12
end subroutine strain_stress2

subroutine strain_stress3(Ndim,C11,C22,C33,&
        C44,C55,C66,C23,C13,C12,&
        C14,C15,C16,C24,C25,C26,&
        C34,C35,C36,C45,C46,C56,&
        E11,E22,E33,E23,E13,E12,S11,S22,S33,S23,S13,S12)
    integer,intent(in) :: Ndim
    real(kind=rkind),intent(in)  :: &
        C11(Ndim),C22(Ndim),C33(Ndim),&
        C23(Ndim),C13(Ndim),C12(Ndim),&
        C44(Ndim),C55(Ndim),C66(Ndim),&
        C14(Ndim),C15(Ndim),C16(Ndim),&
        C24(Ndim),C25(Ndim),C26(Ndim),&
        C34(Ndim),C35(Ndim),C36(Ndim),&
        C45(Ndim),C46(Ndim),C56(Ndim)
    real(kind=rkind),intent(in)  :: E11(Ndim),E22(Ndim),E33(Ndim)
    real(kind=rkind),intent(in)  :: E23(Ndim),E13(Ndim),E12(Ndim)
    real(kind=rkind),intent(out) :: S11(Ndim),S22(Ndim),S33(Ndim)
    real(kind=rkind),intent(out) :: S23(Ndim),S13(Ndim),S12(Ndim)
    S11=C11*E11+C12*E22+C13*E33+2d0*C14*E23+2d0*C15*E13+2d0*C16*E12
    S22=C12*E11+C22*E22+C23*E33+2d0*C24*E23+2d0*C25*E13+2d0*C26*E12
    S33=C13*E11+C23*E22+C33*E33+2d0*C34*E23+2d0*C35*E13+2d0*C36*E12
    S23=C14*E11+C24*E22+C34*E33+2d0*C44*E23+2d0*C45*E13+2d0*C46*E12
    S13=C15*E11+C25*E22+C35*E33+2d0*C45*E23+2d0*C55*E13+2d0*C56*E12
    S12=C16*E11+C26*E22+C36*E33+2d0*C46*E23+2d0*C55*E13+2d0*C66*E12
end subroutine strain_stress3

subroutine strain_stress_asym1(Ndim,lambda,muX2,&
        E11,E22,E33,E23,E13,E12,E32,E31,E21,&
        S11,S22,S33,S23,S13,S12,S32,S31,S21)
    integer,intent(in) :: Ndim
    real(kind=rkind),intent(in)  :: lambda(Ndim),muX2(Ndim)
    real(kind=rkind),intent(in)  :: E11(Ndim),E22(Ndim),E33(Ndim)
    real(kind=rkind),intent(in)  :: E23(Ndim),E13(Ndim),E12(Ndim)
    real(kind=rkind),intent(in)  :: E32(Ndim),E31(Ndim),E21(Ndim)
    real(kind=rkind),intent(out) :: S11(Ndim),S22(Ndim),S33(Ndim)
    real(kind=rkind),intent(out) :: S23(Ndim),S13(Ndim),S12(Ndim)
    real(kind=rkind),intent(out) :: S32(Ndim),S31(Ndim),S21(Ndim)
    real(kind=rkind) :: trS,ttS
    integer :: i
    do i=1,Ndim
        trS=lambda(i)*(E11(i)+E22(i)+E33(i))
        S11(i)=S11(i)+trS+muX2(i)*E11(i)
        S22(i)=S22(i)+trS+muX2(i)*E22(i)
        S33(i)=S33(i)+trS+muX2(i)*E33(i)
        ttS=muX2(i)*(E23(i)+E32(i))*0.5d0
        S23(i)=S23(i)+ttS
        S32(i)=S32(i)+ttS
        ttS=muX2(i)*(E13(i)+E31(i))*0.5d0
        S13(i)=S13(i)+ttS
        S31(i)=S31(i)+ttS
        ttS=muX2(i)*(E12(i)+E21(i))*0.5d0
        S12(i)=S12(i)+ttS
        S21(i)=S21(i)+ttS
    enddo
!    S33=lambda*(E11+E22+E33)
!    S11=S33+muX2*E11
!    S22=S33+muX2*E22
!    S33=S33+muX2*E33
!    S23=muX2*(E23+E32)*0.5d0
!    S13=muX2*(E13+E31)*0.5d0
!    S12=muX2*(E12+E21)*0.5d0
!    S32=S23
!    S31=S13
!    S21=S12
end subroutine strain_stress_asym1

subroutine strain_stress_asym2(Ndim,C11,C22,C33,&
        C44,C55,C66,C23,C13,C12,&
        E11,E22,E33,E23,E13,E12,E32,E31,E21,&
        S11,S22,S33,S23,S13,S12,S32,S31,S21)
    integer,intent(in) :: Ndim
    real(kind=rkind),intent(in)  :: &
        C11(Ndim),C22(Ndim),C33(Ndim),&
        C23(Ndim),C13(Ndim),C12(Ndim),&
        C44(Ndim),C55(Ndim),C66(Ndim)
    real(kind=rkind),intent(in)  :: E11(Ndim),E22(Ndim),E33(Ndim)
    real(kind=rkind),intent(in)  :: E23(Ndim),E13(Ndim),E12(Ndim)
    real(kind=rkind),intent(in)  :: E32(Ndim),E31(Ndim),E21(Ndim)
    real(kind=rkind),intent(out) :: S11(Ndim),S22(Ndim),S33(Ndim)
    real(kind=rkind),intent(out) :: S23(Ndim),S13(Ndim),S12(Ndim)
    real(kind=rkind),intent(out) :: S32(Ndim),S31(Ndim),S21(Ndim)
    S11=C11*E11+C12*E22+C13*E33
    S22=C12*E11+C22*E22+C23*E33
    S33=C13*E11+C23*E22+C33*E33
    S23=C44*(E23+E32)
    S13=C55*(E13+E31)
    S12=C66*(E12+E21)
    S32=S23
    S31=S13
    S21=S12
end subroutine strain_stress_asym2

subroutine strain_stress_asym3(Ndim,C11,C22,C33,&
        C44,C55,C66,C23,C13,C12,&
        C14,C15,C16,C24,C25,C26,&
        C34,C35,C36,C45,C46,C56,&
        E11,E22,E33,E23,E13,E12,E32,E31,E21,&
        S11,S22,S33,S23,S13,S12,S32,S31,S21)
    integer,intent(in) :: Ndim
    real(kind=rkind),intent(in)  :: &
        C11(Ndim),C22(Ndim),C33(Ndim),&
        C23(Ndim),C13(Ndim),C12(Ndim),&
        C44(Ndim),C55(Ndim),C66(Ndim),&
        C14(Ndim),C15(Ndim),C16(Ndim),&
        C24(Ndim),C25(Ndim),C26(Ndim),&
        C34(Ndim),C35(Ndim),C36(Ndim),&
        C45(Ndim),C46(Ndim),C56(Ndim)
    real(kind=rkind),intent(in)  :: E11(Ndim),E22(Ndim),E33(Ndim)
    real(kind=rkind),intent(in)  :: E23(Ndim),E13(Ndim),E12(Ndim)
    real(kind=rkind),intent(in)  :: E32(Ndim),E31(Ndim),E21(Ndim)
    real(kind=rkind),intent(out) :: S11(Ndim),S22(Ndim),S33(Ndim)
    real(kind=rkind),intent(out) :: S23(Ndim),S13(Ndim),S12(Ndim)
    real(kind=rkind),intent(out) :: S32(Ndim),S31(Ndim),S21(Ndim)
    S11=C11*E11+C12*E22+C13*E33+C14*(E32+E23)+C15*(E31+E13)+C16*(E21+E12)
    S22=C12*E11+C22*E22+C23*E33+C24*(E32+E23)+C25*(E31+E13)+C26*(E21+E12)
    S33=C13*E11+C23*E22+C33*E33+C34*(E32+E23)+C35*(E31+E13)+C36*(E21+E12)
    S23=C14*E11+C24*E22+C34*E33+C44*(E32+E23)+C45*(E31+E13)+C46*(E21+E12)
    S13=C15*E11+C25*E22+C35*E33+C45*(E32+E23)+C55*(E31+E13)+C56*(E21+E12)
    S12=C16*E11+C26*E22+C36*E33+C46*(E32+E23)+C55*(E31+E13)+C66*(E21+E12)
    S32=S23
    S31=S13
    S21=S12
end subroutine strain_stress_asym3

subroutine stress_TPK1(Ndim,T11,T22,T33,T23,T13,T12,&
        E11,E22,E33,E23,E13,E12,E32,E31,E21,&
        S11,S22,S33,S23,S13,S12,S32,S31,S21)
    integer,intent(in) :: Ndim
    real(kind=rkind),intent(in)    :: T11,T22,T33
    real(kind=rkind),intent(in)    :: T23,T13,T12
    real(kind=rkind),intent(in)    :: E11(Ndim),E22(Ndim),E33(Ndim)
    real(kind=rkind),intent(in)    :: E23(Ndim),E13(Ndim),E12(Ndim)
    real(kind=rkind),intent(in)    :: E32(Ndim),E31(Ndim),E21(Ndim)
    real(kind=rkind),intent(inout) :: S11(Ndim),S22(Ndim),S33(Ndim)
    real(kind=rkind),intent(inout) :: S23(Ndim),S13(Ndim),S12(Ndim)
    real(kind=rkind),intent(inout) :: S32(Ndim),S31(Ndim),S21(Ndim)

    S11=S11              +(T12    )*E21+(T13    )*E31&
           -(T12    )*E12+(T11+T22)*E22+(T23    )*E32&
           -(T13    )*E13+(T23    )*E23+(T11+T33)*E33
    S21=S21+(T12    )*E11-(T11-T22)*E21+(T23    )*E31&
           -(T11+T22)*E12-(T12    )*E22-(T13    )*E32&
           -(T23    )*E13-(T13    )*E23+(T12    )*E33
    S31=S31+(T13    )*E11+(T23    )*E21-(T11-T33)*E31&
           -(T23    )*E12+(T13    )*E22-(T12    )*E32&
           -(T11+T33)*E13-(T12    )*E23-(T13    )*E33
    S12=S12-(T12    )*E11-(T22+T11)*E21-(T23    )*E31&
           +(T11-T22)*E12+(T12    )*E22+(T13    )*E32&
           -(T23    )*E13-(T13    )*E23+(T12    )*E33
    S22=S22+(T22+T11)*E11-(T12    )*E21+(T13    )*E31&
           +(T12    )*E12              +(T23    )*E32&
           +(T13    )*E13-(T23    )*E23+(T22+T33)*E33
    S32=S32+(T23    )*E11-(T13    )*E21-(T12    )*E31&
           +(T13    )*E12+(T23    )*E22-(T22-T33)*E32&
           -(T12    )*E13-(T22+T33)*E23-(T23    )*E33
    S13=S13-(T13    )*E11-(T23    )*E21-(T33+T11)*E31&
           -(T23    )*E12+(T13    )*E22-(T12    )*E32&
           +(T11-T33)*E13+(T12    )*E23+(T13    )*E33
    S23=S23+(T23    )*E11-(T13    )*E21-(T12    )*E31&
           -(T13    )*E12-(T23    )*E22-(T33+T22)*E32&
           +(T12    )*E13+(T22-T33)*E23+(T23    )*E33
    S33=S33+(T33+T11)*E11+(T12    )*E21-(T13    )*E31&
           +(T12    )*E12+(T33+T22)*E22-(T23    )*E32&
           +(T13    )*E13+(T23    )*E23

end subroutine stress_TPK1

end module solve_mod
