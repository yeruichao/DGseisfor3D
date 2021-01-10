!*******************************************************************!
!*  This module obtains semidiscretized wave equation              *!
!*                                                                 *!
!*  Author: Ruichao Ye    Email: Ruichao.ye@gmail.com              *!
!*******************************************************************!

!--------------------------------------------------------------------
module solvePML_mod
!--------------------------------------------------------------------

    use para_mod,     only : def_pOrder,pNp,Nfp,Nfp4,&
                             rpml,convtest,withrupt
    use datatype_mod, only : rkind,matrices,sources,rupture,&
                             tetmesh_geometry,tetmesh_domain,&
                             tetmesh_material,&
                             vector_array,tensor_array,wavefield,&
                             auxilary_array,reset_vector,&
                             reset_wavefield,reset_auxilary_array,&
                             BLK_axpy
    use surface_mod,  only : rupture_bnd
    use solve_mod,    only : strain_stress!,penalty_flux_full

!--------------------------------------------------------------------

    implicit none

!--------------------------------------------------------------------
    contains
!--------------------------------------------------------------------

subroutine bilinear_Elastic3D_full(matrix,mesh,subdomain,material,&
        srcs,rhsW,divW,tmpW,localtime,flag)
! inout:
    type(matrices)         :: matrix
    type(tetmesh_geometry) :: mesh
    type(tetmesh_Domain)   :: subdomain
    type(tetmesh_material) :: material
    type(sources)          :: srcs
    type(wavefield)        :: rhsW
    type(auxilary_array)   :: divW,tmpW
    real(kind=rkind),intent(in) :: localtime
    logical,intent(in)     :: flag(:)
! auxilary:
    integer :: head,fhead,fh,offset,i,j
    integer :: ierr,Ndim

!!!!!!!!!!!!!!!! compute volume terms !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call reset_auxilary_array(divW)

    call DGEMM('n','n',pNp,mesh%Nele*12,pNp,1d0,&
        matrix%Dr,pNp,tmpW%array,pNp,0d0,divW%array,pNp)

    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,1),divW%V%x,rhsW%E%xx)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,1),divW%V%y,rhsW%E%xy)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,1),divW%V%z,rhsW%E%xz)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,4),divW%V%x,rhsW%E%yx)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,4),divW%V%y,rhsW%E%yy)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,4),divW%V%z,rhsW%E%yz)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,7),divW%V%x,rhsW%E%zx)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,7),divW%V%y,rhsW%E%zy)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,7),divW%V%z,rhsW%E%zz)
                                               
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,1),divW%S%xx,rhsW%V%x)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,1),divW%S%xy,rhsW%V%y)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,1),divW%S%xz,rhsW%V%z)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,4),divW%S%yx,rhsW%V%x)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,4),divW%S%yy,rhsW%V%y)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,4),divW%S%yz,rhsW%V%z)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,7),divW%S%zx,rhsW%V%x)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,7),divW%S%zy,rhsW%V%y)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,7),divW%S%zz,rhsW%V%z)

    call DGEMM('n','n',pNp,mesh%Nele*12,pNp,1d0,&
        matrix%Ds,pNp,tmpW%array,pNp,0d0,divW%array,pNp)

    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,2),divW%V%x,rhsW%E%xx)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,2),divW%V%y,rhsW%E%xy)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,2),divW%V%z,rhsW%E%xz)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,5),divW%V%x,rhsW%E%yx)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,5),divW%V%y,rhsW%E%yy)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,5),divW%V%z,rhsW%E%yz)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,8),divW%V%x,rhsW%E%zx)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,8),divW%V%y,rhsW%E%zy)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,8),divW%V%z,rhsW%E%zz)
                                                
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,2),divW%S%xx,rhsW%V%x)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,2),divW%S%xy,rhsW%V%y)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,2),divW%S%xz,rhsW%V%z)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,5),divW%S%yx,rhsW%V%x)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,5),divW%S%yy,rhsW%V%y)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,5),divW%S%yz,rhsW%V%z)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,8),divW%S%zx,rhsW%V%x)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,8),divW%S%zy,rhsW%V%y)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,8),divW%S%zz,rhsW%V%z)

    call DGEMM('n','n',pNp,mesh%Nele*12,pNp,1d0,&
        matrix%Dt,pNp,tmpW%array,pNp,0d0,divW%array,pNp)

    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,3),divW%V%x,rhsW%E%xx)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,3),divW%V%y,rhsW%E%xy)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,3),divW%V%z,rhsW%E%xz)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,6),divW%V%x,rhsW%E%yx)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,6),divW%V%y,rhsW%E%yy)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,6),divW%V%z,rhsW%E%yz)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,9),divW%V%x,rhsW%E%zx)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,9),divW%V%y,rhsW%E%zy)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,9),divW%V%z,rhsW%E%zz)
                                               
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,3),divW%S%xx,rhsW%V%x)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,3),divW%S%xy,rhsW%V%y)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,3),divW%S%xz,rhsW%V%z)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,6),divW%S%yx,rhsW%V%x)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,6),divW%S%yy,rhsW%V%y)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,6),divW%S%yz,rhsW%V%z)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,9),divW%S%zx,rhsW%V%x)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,9),divW%S%zy,rhsW%V%y)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,9),divW%S%zz,rhsW%V%z)

    call DSCAL(mesh%Nele*pNp*3,0.5d0,rhsW%E%yz,1)

!!!!!!!!!!!!!!!! compute flux !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call reset_auxilary_array(divW)

    do i=1,mesh%Nhele
        ! Further than first arrival
        if(localtime.lt.srcs%t_direct(i))cycle
        ! Not computed if indicated 
        if(flag(i))cycle
        head=(i-1)*pNp+1;offset=(i-1)*Nfp*4
        do j=1,4
            fhead=(j-1)*Nfp+1
            fh=fhead+offset
            call penalty_flux_full(Nfp,&
                mesh%fbctype(1:,j,i),material%k_media(i),&
                mesh%nx(j,i),mesh%ny(j,i),mesh%nz(j,i),&
                mesh%Fscale(j,i),mesh%vmapM(fhead:),mesh%vmapP(fh:),&
                tmpW%V%x (head:),tmpW%V%y (head:),tmpW%V%z (head:),&
                tmpW%S%xx(head:),tmpW%S%xy(head:),tmpW%S%xz(head:),&
                tmpW%S%yx(head:),tmpW%S%yy(head:),tmpW%S%yz(head:),&
                tmpW%S%zx(head:),tmpW%S%zy(head:),tmpW%S%zz(head:),&
                tmpW%V%x        ,tmpW%V%y        ,tmpW%V%z        ,&
                tmpW%S%xx       ,tmpW%S%xy       ,tmpW%S%xz       ,&
                tmpW%S%yx       ,tmpW%S%yy       ,tmpW%S%yz       ,&
                tmpW%S%zx       ,tmpW%S%zy       ,tmpW%S%zz       ,&
                divW%fV%x (fh:) ,divW%fV%y (fh:) ,divW%fV%z (fh:) ,&
                divW%fE%xx(fh:) ,divW%fE%yy(fh:) ,divW%fE%zz(fh:) ,&
                divW%fE%yz(fh:) ,divW%fE%xz(fh:) ,divW%fE%xy(fh:) ,&
                mesh)
        enddo
    enddo

!    call DGEMM('n','n',pNp,mesh%Nele*3,4*Nfp,1d0,&
!        matrix%LIFT,pNp,divW%fV%x,4*Nfp,1d0,rhsW%V%x,pNp)
!    call DGEMM('n','n',pNp,mesh%Nele*3,4*Nfp,1d0,&
!        matrix%LIFT,pNp,divW%fE%xx,4*Nfp,1d0,rhsW%E%xx,pNp)
!    call DGEMM('n','n',pNp,mesh%Nele*3,4*Nfp,1d0,&
!        matrix%LIFT,pNp,divW%fE%yz,4*Nfp,1d0,rhsW%E%yz,pNp)
    call DGEMM('n','n',pNp,mesh%Nele*9,4*Nfp,1d0,&
        matrix%LIFT,pNp,divW%array,4*Nfp,1d0,rhsW%array,pNp)

    rhsW%V%x=rhsW%V%x/material%rho
    rhsW%V%y=rhsW%V%y/material%rho
    rhsW%V%z=rhsW%V%z/material%rho

end subroutine bilinear_Elastic3D_full

subroutine bilinear_Elastic3D_trE(matrix,mesh,subdomain,material,&
        srcs,rhsW,divW,tmpW,localtime,flag)
! inout:
    type(matrices)         :: matrix
    type(tetmesh_geometry) :: mesh
    type(tetmesh_Domain)   :: subdomain
    type(tetmesh_material) :: material
    type(sources)          :: srcs
    type(wavefield)        :: rhsW
    type(auxilary_array)   :: divW,tmpW
    real(kind=rkind),intent(in) :: localtime
    logical,intent(in)     :: flag(:)
! auxilary:
    integer :: head,fhead,fh,offset,i,j
    integer :: ierr,Ndim

!!!!!!!!!!!!!!!! compute volume terms !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call reset_auxilary_array(divW)

    call DGEMM('n','n',pNp,mesh%Nele*3,pNp,1d0,&
        matrix%Dr,pNp,tmpW%V%array,pNp,0d0,divW%V%array,pNp)

    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,1),divW%V%x,rhsW%E%xx)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,4),divW%V%y,rhsW%E%yy)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,7),divW%V%z,rhsW%E%zz)

    call DGEMM('n','n',pNp,mesh%Nele*3,pNp,1d0,&
        matrix%Ds,pNp,tmpW%V%array,pNp,0d0,divW%V%array,pNp)

    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,2),divW%V%x,rhsW%E%xx)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,5),divW%V%y,rhsW%E%yy)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,8),divW%V%z,rhsW%E%zz)

    call DGEMM('n','n',pNp,mesh%Nele*3,pNp,1d0,&
        matrix%Dt,pNp,tmpW%V%array,pNp,0d0,divW%V%array,pNp)

    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,3),divW%V%x,rhsW%E%xx)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,6),divW%V%y,rhsW%E%yy)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,9),divW%V%z,rhsW%E%zz)

!!!!!!!!!!!!!!!! compute flux !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call reset_auxilary_array(divW)

    do i=1,mesh%Nhele
        ! Further than first arrival
        if(localtime.lt.srcs%t_direct(i))cycle
        ! Not computed if indicated 
        if(flag(i))cycle
        head=(i-1)*pNp+1;offset=(i-1)*Nfp*4
        do j=1,4
            fhead=(j-1)*Nfp+1
            fh=fhead+offset
            call penalty_flux_full(Nfp,&
                mesh%fbctype(1:,j,i),material%k_media(i),&
                mesh%nx(j,i),mesh%ny(j,i),mesh%nz(j,i),&
                mesh%Fscale(j,i),mesh%vmapM(fhead:),mesh%vmapP(fh:),&
                tmpW%V%x (head:),tmpW%V%y (head:),tmpW%V%z (head:),&
                tmpW%S%xx(head:),tmpW%S%xy(head:),tmpW%S%xz(head:),&
                tmpW%S%yx(head:),tmpW%S%yy(head:),tmpW%S%yz(head:),&
                tmpW%S%zx(head:),tmpW%S%zy(head:),tmpW%S%zz(head:),&
                tmpW%V%x        ,tmpW%V%y        ,tmpW%V%z        ,&
                tmpW%S%xx       ,tmpW%S%xy       ,tmpW%S%xz       ,&
                tmpW%S%yx       ,tmpW%S%yy       ,tmpW%S%yz       ,&
                tmpW%S%zx       ,tmpW%S%zy       ,tmpW%S%zz       ,&
                divW%fV%x (fh:) ,divW%fV%y (fh:) ,divW%fV%z (fh:) ,&
                divW%fE%xx(fh:) ,divW%fE%yy(fh:) ,divW%fE%zz(fh:) ,&
                divW%fE%yz(fh:) ,divW%fE%xz(fh:) ,divW%fE%xy(fh:) ,&
                mesh)
        enddo
    enddo

    call DGEMM('n','n',pNp,mesh%Nele*3,4*Nfp,1d0,&
        matrix%LIFT,pNp,divW%fE%array,4*Nfp,1d0,rhsW%E%array,pNp)

end subroutine bilinear_Elastic3D_TrE

subroutine bilinear_Elastic3D_PML(matrix,mesh,subdomain,material,&
        srcs,rhsW,divW,tmpW,localtime,flag)
! inout:
    type(matrices)         :: matrix
    type(tetmesh_geometry) :: mesh
    type(tetmesh_Domain)   :: subdomain
    type(tetmesh_material) :: material
    type(sources)          :: srcs
    type(wavefield)        :: rhsW
    type(auxilary_array)   :: divW,tmpW
    real(kind=rkind),intent(in) :: localtime
    logical,intent(in)     :: flag(:)
! auxilary:
    integer :: head,fhead,fh,offset,i,j
    integer :: ierr,Ndim

!!!!!!!!!!!!!!!! compute volume terms !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call reset_auxilary_array(divW)

    call DGEMM('n','n',pNp,mesh%Nele*12,pNp,1d0,&
        matrix%Dr,pNp,tmpW%array,pNp,0d0,divW%array,pNp)

    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,1),divW%V%y,rhsW%E%xy)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,1),divW%V%z,rhsW%E%xz)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,4),divW%V%x,rhsW%E%yx)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,4),divW%V%z,rhsW%E%yz)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,7),divW%V%x,rhsW%E%zx)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,7),divW%V%y,rhsW%E%zy)
                                               
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,1),divW%S%xx,rhsW%V%x)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,1),divW%S%xy,rhsW%V%y)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,1),divW%S%xz,rhsW%V%z)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,4),divW%S%yx,rhsW%V%x)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,4),divW%S%yy,rhsW%V%y)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,4),divW%S%yz,rhsW%V%z)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,7),divW%S%zx,rhsW%V%x)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,7),divW%S%zy,rhsW%V%y)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,7),divW%S%zz,rhsW%V%z)

    call DGEMM('n','n',pNp,mesh%Nele*12,pNp,1d0,&
        matrix%Ds,pNp,tmpW%array,pNp,0d0,divW%array,pNp)

    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,2),divW%V%y,rhsW%E%xy)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,2),divW%V%z,rhsW%E%xz)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,5),divW%V%x,rhsW%E%yx)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,5),divW%V%z,rhsW%E%yz)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,8),divW%V%x,rhsW%E%zx)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,8),divW%V%y,rhsW%E%zy)
                                               
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,2),divW%S%xx,rhsW%V%x)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,2),divW%S%xy,rhsW%V%y)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,2),divW%S%xz,rhsW%V%z)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,5),divW%S%yx,rhsW%V%x)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,5),divW%S%yy,rhsW%V%y)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,5),divW%S%yz,rhsW%V%z)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,8),divW%S%zx,rhsW%V%x)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,8),divW%S%zy,rhsW%V%y)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,8),divW%S%zz,rhsW%V%z)

    call DGEMM('n','n',pNp,mesh%Nele*12,pNp,1d0,&
        matrix%Dt,pNp,tmpW%array,pNp,0d0,divW%array,pNp)

    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,3),divW%V%y,rhsW%E%xy)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,3),divW%V%z,rhsW%E%xz)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,6),divW%V%x,rhsW%E%yx)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,6),divW%V%z,rhsW%E%yz)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,9),divW%V%x,rhsW%E%zx)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,9),divW%V%y,rhsW%E%zy)
                                               
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,3),divW%S%xx,rhsW%V%x)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,3),divW%S%xy,rhsW%V%y)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,3),divW%S%xz,rhsW%V%z)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,6),divW%S%yx,rhsW%V%x)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,6),divW%S%yy,rhsW%V%y)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,6),divW%S%yz,rhsW%V%z)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,9),divW%S%zx,rhsW%V%x)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,9),divW%S%zy,rhsW%V%y)
    call BLK_axpy(pNp,mesh%Nhele,mesh%invJac(:,9),divW%S%zz,rhsW%V%z)

    call DSCAL(mesh%Nele*pNp*3,0.5d0,rhsW%E%yz,1)

!!!!!!!!!!!!!!!! compute flux !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call reset_auxilary_array(divW)

    do i=1,mesh%Nhele
        ! Further than first arrival
        if(localtime.lt.srcs%t_direct(i))cycle
        ! Not computed if indicated 
        if(flag(i))cycle
        head=(i-1)*pNp+1;offset=(i-1)*Nfp*4
        do j=1,4
            fhead=(j-1)*Nfp+1
            fh=fhead+offset
            call penalty_flux_full(Nfp,&
                mesh%fbctype(1:,j,i),material%k_media(i),&
                mesh%nx(j,i),mesh%ny(j,i),mesh%nz(j,i),&
                mesh%Fscale(j,i),mesh%vmapM(fhead:),mesh%vmapP(fh:),&
                tmpW%V%x (head:),tmpW%V%y (head:),tmpW%V%z (head:),&
                tmpW%S%xx(head:),tmpW%S%xy(head:),tmpW%S%xz(head:),&
                tmpW%S%yx(head:),tmpW%S%yy(head:),tmpW%S%yz(head:),&
                tmpW%S%zx(head:),tmpW%S%zy(head:),tmpW%S%zz(head:),&
                tmpW%V%x        ,tmpW%V%y        ,tmpW%V%z        ,&
                tmpW%S%xx       ,tmpW%S%xy       ,tmpW%S%xz       ,&
                tmpW%S%yx       ,tmpW%S%yy       ,tmpW%S%yz       ,&
                tmpW%S%zx       ,tmpW%S%zy       ,tmpW%S%zz       ,&
                divW%fV%x (fh:) ,divW%fV%y (fh:) ,divW%fV%z (fh:) ,&
                divW%fE%xx(fh:) ,divW%fE%yy(fh:) ,divW%fE%zz(fh:) ,&
                divW%fE%yz(fh:) ,divW%fE%xz(fh:) ,divW%fE%xy(fh:) ,&
                mesh)
        enddo
    enddo

    call DGEMM('n','n',pNp,mesh%Nele*3,4*Nfp,1d0,&
        matrix%LIFT,pNp,divW%fV%array,4*Nfp,1d0,rhsW%V%array,pNp)
    call DGEMM('n','n',pNp,mesh%Nele*3,4*Nfp,1d0,&
        matrix%LIFT,pNp,divW%fE%yz,4*Nfp,1d0,rhsW%E%yz,pNp)

    rhsW%V%x=rhsW%V%x/material%rho
    rhsW%V%y=rhsW%V%y/material%rho
    rhsW%V%z=rhsW%V%z/material%rho

end subroutine bilinear_Elastic3D_PML

subroutine penalty_flux_full(&
    Nfp,fbctype,k_media,nx,ny,nz,Fscale,mapM,mapP,&
    v1 ,v2 ,v3 ,S11 ,S12 ,S13 ,S21 ,S22 ,S23 ,S31 ,S32 ,S33 ,&
    v1g,v2g,v3g,S11g,S12g,S13g,S21g,S22g,S23g,S31g,S32g,S33g,&
    fv1,fv2,fv3,fE11,fE22,fE33,fE23,fE13,fE12,&
    mesh)
! input
    integer,intent(in) :: Nfp,k_media,fbctype(2)
    real(kind=rkind),intent(in) :: nx,ny,nz,Fscale
    type(tetmesh_geometry) :: mesh
    real(kind=rkind),intent(in) :: &
             v1 (1), v2 (1), v3 (1),S11 (1),S12 (1),S13 (1),&
            S21 (1),S22 (1),S23 (1),S31 (1),S32 (1),S33 (1),&
             v1g(1), v2g(1), v3g(1),S11g(1),S12g(1),S13g(1),&
            S21g(1),S22g(1),S23g(1),S31g(1),S32g(1),S33g(1)
    integer,intent(in) :: mapM(Nfp),mapP(Nfp)
! output
    real(kind=rkind) ::  fv1(Nfp), fv2(Nfp), fv3(Nfp),&
                        fE11(Nfp),fE22(Nfp),fE33(Nfp),&
                        fE23(Nfp),fE13(Nfp),fE12(Nfp)
! auxilary
    real(kind=rkind) :: &
               S11m(Nfp), S12m(Nfp), S13m(Nfp),&
               S21m(Nfp), S22m(Nfp), S23m(Nfp),&
               S31m(Nfp), S32m(Nfp), S33m(Nfp),&
               S11p(Nfp), S12p(Nfp), S13p(Nfp),&
               S21p(Nfp), S22p(Nfp), S23p(Nfp),&
               S31p(Nfp), S32p(Nfp), S33p(Nfp),&
                dv1(Nfp),  dv2(Nfp),  dv3(Nfp),&
               dSn1(Nfp), dSn2(Nfp), dSn3(Nfp),&
               St1m(Nfp), St2m(Nfp), St3m(Nfp),&
               St1p(Nfp), St2p(Nfp), St3p(Nfp),&
               Sn1m(Nfp), Sn2m(Nfp), Sn3m(Nfp),&
               Sn1p(Nfp), Sn2p(Nfp), Sn3p(Nfp),&
               dv1p(Nfp), dv2p(Nfp), dv3p(Nfp),&
              dSn1p(Nfp),dSn2p(Nfp),dSn3p(Nfp),&
               tau1(Nfp), tau2(Nfp), tau3(Nfp),&
               ndSn(Nfp), ndv(Nfp), ndSnp(Nfp), ndvp(Nfp)
    real(kind=rkind),parameter :: alpha=0.5d0,alpha_rupt=0.5d0
    integer :: j,fhead,ftail
    logical :: k_rupt(Nfp)

    if(fbctype(1).eq.0)then 
        dv1 = v1g(mapP) - v1(mapM) 
        dv2 = v2g(mapP) - v2(mapM) 
        dv3 = v3g(mapP) - v3(mapM)
        S11p = S11g(mapP); S11m = S11(mapM)
        S12p = S12g(mapP); S12m = S12(mapM)
        S13p = S13g(mapP); S13m = S13(mapM)
        S21p = S21g(mapP); S21m = S21(mapM)
        S22p = S22g(mapP); S22m = S22(mapM)
        S23p = S23g(mapP); S23m = S23(mapM)
        S31p = S31g(mapP); S31m = S31(mapM)
        S32p = S32g(mapP); S32m = S32(mapM)
        S33p = S33g(mapP); S33m = S33(mapM)
        dSn1 = nx*(S11p-S11m) + ny*(S21p-S21m) + nz*(S31p-S31m) 
        dSn2 = nx*(S12p-S12m) + ny*(S22p-S22m) + nz*(S32p-S32m) 
        dSn3 = nx*(S13p-S13m) + ny*(S23p-S23m) + nz*(S33p-S33m) 
        if(k_media.ge.1)then ! solid media
            if(fbctype(2).eq.0)then ! solid-solid interface
                dSn1p=dSn1+alpha*dv1
                dSn2p=dSn2+alpha*dv2
                dSn3p=dSn3+alpha*dv3
                 dv1p=dv1+alpha*dSn1
                 dv2p=dv2+alpha*dSn2
                 dv3p=dv3+alpha*dSn3
            elseif(fbctype(2).eq.1)then ! solid-fluid interface
                  ndv=dv1*nx+dv2*ny+dv3*nz 
                 dv1p=ndv*nx
                 dv2p=ndv*ny
                 dv3p=ndv*nz
                dSn1p=dSn1+alpha*dv1p
                dSn2p=dSn2+alpha*dv2p
                dSn3p=dSn3+alpha*dv3p
                 dv1p=dv1p+alpha*dSn1
                 dv2p=dv2p+alpha*dSn2
                 dv3p=dv3p+alpha*dSn3
            endif
        else ! fluid media
                 ndSn=dSn1*nx+dSn2*ny+dSn3*nz
                  ndv= dv1*nx+ dv2*ny+ dv3*nz
            if(fbctype(2).eq.2)then ! fluid-solid interface
                ndSnp=ndSn+alpha* ndv
                 ndvp= ndv+alpha*ndSn
            else ! fluid-fluid interface
                ndSnp=ndSn+alpha* ndv
                 ndvp= ndv+alpha*ndSn
            endif
                dSn1p=ndSnp*nx
                dSn2p=ndSnp*ny
                dSn3p=ndSnp*nz
                 dv1p= ndvp*nx
                 dv2p= ndvp*ny
                 dv3p= ndvp*nz
        endif
    elseif(fbctype(1).eq.1)then
        S11m = S11(mapM)
        S12m = S12(mapM)
        S13m = S13(mapM)
        S21m = S21(mapM)
        S22m = S22(mapM)
        S23m = S23(mapM)
        S31m = S31(mapM)
        S32m = S32(mapM)
        S33m = S33(mapM)
        dSn1 = - nx*S11m - ny*S21m - nz*S31m
        dSn2 = - nx*S12m - ny*S22m - nz*S32m
        dSn3 = - nx*S13m - ny*S23m - nz*S33m
        if(fbctype(2).eq.1)then ! Free surface boundary
            dv1=dSn1
            dv2=dSn2
            dv3=dSn3
        elseif(fbctype(2).eq.2)then ! Absrobing boundary
            dv1=dSn1-v1(mapM)
            dv2=dSn2-v2(mapM)
            dv3=dSn3-v3(mapM)
            dSn1=0d0
            dSn2=0d0
            dSn3=0d0
        endif
        if(k_media.ge.1)then ! solid media
                dSn1p=dSn1+alpha*dv1
                dSn2p=dSn2+alpha*dv2
                dSn3p=dSn3+alpha*dv3
                 dv1p=dv1+alpha*dSn1
                 dv2p=dv2+alpha*dSn2
                 dv3p=dv3+alpha*dSn3
        else ! fluid media
                 ndSn=dSn1*nx+dSn2*ny+dSn3*nz
                  ndv= dv1*nx+ dv2*ny+ dv3*nz
                ndSnp=ndSn+alpha* ndv
                 ndvp= ndv+alpha*ndSn
                dSn1p=ndSnp*nx
                dSn2p=ndSnp*ny
                dSn3p=ndSnp*nz
                 dv1p= ndvp*nx
                 dv2p= ndvp*ny
                 dv3p= ndvp*nz
        endif
    elseif(fbctype(1).eq.2)then ! rupture
         dv1 = v1g(mapP) - v1(mapM) 
         dv2 = v2g(mapP) - v2(mapM) 
         dv3 = v3g(mapP) - v3(mapM)
         ndv = dv1*nx+dv2*ny+dv3*nz 
        S11p = S11g(mapP); S11m = S11(mapM)
        S12p = S12g(mapP); S12m = S12(mapM)
        S13p = S13g(mapP); S13m = S13(mapM)
        S21p = S21g(mapP); S21m = S21(mapM)
        S22p = S22g(mapP); S22m = S22(mapM)
        S23p = S23g(mapP); S23m = S23(mapM)
        S31p = S31g(mapP); S31m = S31(mapM)
        S32p = S32g(mapP); S32m = S32(mapM)
        S33p = S33g(mapP); S33m = S33(mapM)
        dSn1 = nx*S11m + ny*S21m + nz*S31m
        dSn2 = nx*S12m + ny*S22m + nz*S32m
        dSn3 = nx*S13m + ny*S23m + nz*S33m
        ndSn = dSn1*nx + dSn2*ny + dSn3*nz
        St1m = dSn1 - ndSn * nx
        St2m = dSn2 - ndSn * ny
        St3m = dSn3 - ndSn * nz
        call rupture_bnd(Nfp,fbctype(2),mesh%rupt,&
            nx,ny,nz,dv1,dv2,dv3,ndv,&
            St1m,St2m,St3m,ndSn,&
            dSn1,dSn2,dSn3,dSn1p,dSn2p,dSn3p)
        dv1p =dv1 + alpha_rupt*dSn1p
        dv2p =dv2 + alpha_rupt*dSn2p
        dv3p =dv3 + alpha_rupt*dSn3p
        dSn1p=dSn1 + alpha_rupt*dv1
        dSn2p=dSn2 + alpha_rupt*dv2
        dSn3p=dSn3 + alpha_rupt*dv3
    elseif(fbctype(1).eq.3)then ! external Neumann surface source
        S11m = S11(mapM)
        S12m = S12(mapM)
        S13m = S13(mapM)
        S21m = S21(mapM)
        S22m = S22(mapM)
        S23m = S23(mapM)
        S31m = S31(mapM)
        S32m = S32(mapM)
        S33m = S33(mapM)
        j=fbctype(2)
        fhead=(j-1)*Nfp+1;ftail=j*Nfp
        dSn1 = - nx*S11m - ny*S21m - nz*S31m
        dSn2 = - nx*S12m - ny*S22m - nz*S32m
        dSn3 = - nx*S13m - ny*S23m - nz*S33m
        dSn1=dSn1+mesh%surf%Sn%x(fhead:ftail)
        dSn2=dSn2+mesh%surf%Sn%y(fhead:ftail)
        dSn3=dSn3+mesh%surf%Sn%z(fhead:ftail)
        dv1=0d0!dSn1
        dv2=0d0!dSn2
        dv3=0d0!dSn3
        dv1p=dv1+alpha*dSn1
        dv2p=dv2+alpha*dSn2
        dv3p=dv3+alpha*dSn3
        dSn1p=dSn1+alpha*dv1
        dSn2p=dSn2+alpha*dv2
        dSn3p=dSn3+alpha*dv3
    endif

    fv1  = Fscale * 0.5d0  *  dSn1p
    fv2  = Fscale * 0.5d0  *  dSn2p
    fv3  = Fscale * 0.5d0  *  dSn3p
    fE11 = Fscale * 0.5d0  *  nx*dv1p
    fE22 = Fscale * 0.5d0  *  ny*dv2p
    fE33 = Fscale * 0.5d0  *  nz*dv3p
    fE23 = Fscale * 0.25d0 *( ny*dv3p + nz*dv2p )
    fE13 = Fscale * 0.25d0 *( nx*dv3p + nz*dv1p )
    fE12 = Fscale * 0.25d0 *( nx*dv2p + ny*dv1p )

end subroutine penalty_flux_full

end module solvePML_mod
