!*******************************************************************!
!*  This module solves semi-discretized wave ODE                   *!
!*                                                                 *!
!*  Author: Ruichao Ye    Email: Ruichao.ye@gmail.com              *!
!*******************************************************************!

!--------------------------------------------------------------------
module ODEsolvePML_mod
!--------------------------------------------------------------------

	USE mpi
    use para_mod,     only : def_pOrder,pNp,Nfp,withrupt,friction,&
                             timescheme,&
                             rk4a,rk4b,rk4c,&
                             CNRKW3_IMa,CNRKW3_IMb,&
                             CNRKW3_EXa,CNRKW3_EXb,CNRKW3_c, &
                             IMEXRKCB3a_IMa,IMEXRKCB3a_EXa,&
                             IMEXRKCB3a_b,IMEXRKCB3a_c, &
                             Euler_IMa,Euler_IMb,Euler_IMc,&
                             Euler_EXa,Euler_EXb,Euler_EXc
    use datatype_mod, only : rkind,matrices,sources,tetmesh_domain,&
                             tetmesh_geometry,tetmesh_material,&
                             PML_geometry,&
                             vector_array,tensor_array,&
                             wavefield,auxilary_array,&
                             ini_tensor,ini_vector,del_vector,&
                             ini_wavefield,ini_auxilary_array,&
                             reset_wavefield,reset_auxilary_array,&
                             del_wavefield,del_auxilary_array,&
                             wavefield_AXPY,wavefield_COPY,&
                             wavefield_SCAL,Array_scale,Array_axypz,&
                             vector_SCAL,vector_AXPY,vector_COPY
    use solve_mod,    only : strain_stress
    use solvePML_mod, only : bilinear_Elastic3D_trE,&
                             bilinear_Elastic3D_PML,&
                             bilinear_Elastic3D_full
    use surface_mod,  only : rupture_update,rupture_law
    use source_mod,   only : body_source,rupt_source,surf_source
    use parallel_mod, only : domain_exchange

    implicit none

!    include 'mpif.h'

    type(wavefield) :: EX_W,IM_W,resW,resWs,resWss
    type(auxilary_array) :: divW,tmpW
    type(vector_array) :: dampA,dampB,dampS,resVsss

!--------------------------------------------------------------------
    contains
!--------------------------------------------------------------------

subroutine init_RKvariable_PML(mesh,pNp,Nfp)
    type(tetmesh_geometry) :: mesh
    integer,intent(in) :: pNp,Nfp
    integer :: Ndim
    type(vector_array),pointer :: damp

    Ndim=mesh%Nele*pNp
    if(trim(timescheme).eq.'EXRK4')then
        call ini_wavefield(Ndim,resW)
        call ini_wavefield(Ndim,resWs)
        call ini_wavefield(Ndim,resWss)
    endif
    call ini_auxilary_array(tmpW,pNp,Nfp,mesh%Nele,.false.)
    call ini_auxilary_array(divW,pNp,Nfp,mesh%Nele,.false.)
    call ini_vector(Ndim,dampA)
    call ini_vector(Ndim,dampB)
    call ini_vector(Ndim,dampS)
    call ini_vector(Ndim,resVsss)
    call vector_COPY(Ndim,mesh%PMLinfo%damp,resVsss)
    dampA%x = resVsss%y + resVsss%z             ! Axx
    dampA%y = resVsss%x + resVsss%z             ! Ayy
    dampA%z = resVsss%x + resVsss%y             ! Azz
    dampB%x = resVsss%y * resVsss%z             ! Bxx
    dampB%y = resVsss%x * resVsss%z             ! Byy
    dampB%z = resVsss%x * resVsss%y             ! Bzz
    dampS%x = resVsss%x + resVsss%y + resVsss%z ! a
    dampS%z = resVsss%x * resVsss%y * resVsss%z ! c
    dampS%y =   dampB%x +   dampB%y +   dampB%z ! b
end subroutine init_RKvariable_PML

subroutine del_RKvariable_PML
    call del_wavefield(resW)
    call del_wavefield(resWs)
    call del_wavefield(resWss)
    call del_auxilary_array(tmpW)
    call del_auxilary_array(divW)
    call del_vector(dampA)
    call del_vector(dampB)
    call del_vector(dampS)
end subroutine del_RKvariable_PML

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine RHS_ElasticPML(matrix,mesh,subdomain,material,PMLinfo,&
        srcs,rhsW,W,Ws,Wss,divW,tmpW,localtime,flag)
! inout:
    type(matrices)         :: matrix
    type(tetmesh_geometry) :: mesh
    type(tetmesh_Domain)   :: subdomain
    type(tetmesh_material) :: material
    type(PML_geometry)     :: PMLinfo
    type(sources)          :: srcs
    type(wavefield)        :: rhsW,W,Ws,Wss
    type(auxilary_array)   :: divW,tmpW
    real(kind=rkind),intent(in) :: localtime
    logical,intent(in)     :: flag(:)
! auxilary
    integer :: Ndim,Ndim3
    integer :: ierr

    Ndim=mesh%Nele*pNp
    Ndim3=Ndim*3

    call reset_wavefield(rhsW)

    ! generate Sij

    call strain_stress(pNp,mesh%Nele,material,W%E,tmpW%S)
    call DCOPY(Ndim3,W%V%array,1,tmpW%V%array,1)
    call domain_exchange(subdomain,subdomain%SRtype,tmpW%array,&
        subdomain%SRtype%send_req,subdomain%SRtype%recv_req,100)
    call MPI_Waitall(subdomain%SRtype%host%N_DD_Conn, &
        subdomain%SRtype%send_req,subdomain%SRtype%sta,ierr)
    call MPI_Waitall(subdomain%SRtype%gost%N_DD_Conn,&
        subdomain%SRtype%recv_req,subdomain%SRtype%sta,ierr)
    call DCOPY(Ndim3,tmpW%S%yz,1,tmpW%S%zy,1)
    
    ! generate Sij*

    call strain_stress(pNp,mesh%Nele,material,Ws%E,divW%S)
    call DCOPY(Ndim3,Ws%V%array,1,divW%V%array,1)
    call domain_exchange(subdomain,PMLinfo%SRtype,divW%array,&
        subdomain%SRtype%send_req,subdomain%SRtype%recv_req,200)
    call MPI_Waitall(subdomain%SRtype%host%N_DD_Conn, &
        subdomain%SRtype%send_req,subdomain%SRtype%sta,ierr)
    call MPI_Waitall(subdomain%SRtype%gost%N_DD_Conn,&
        subdomain%SRtype%recv_req,subdomain%SRtype%sta,ierr)
    call DCOPY(Ndim3,divW%S%yz   ,1,divW%S%zy ,1)
    call DCOPY(Ndim3,divW%V%array,1,Ws%V%array,1)
    
    ! Sij = Sij + Aii Sij*
    call Array_axypz(Ndim,1d0,dampA%x,divW%S%xx,tmpW%S%xx) 
    call Array_axypz(Ndim,1d0,dampA%x,divW%S%xy,tmpW%S%xy) 
    call Array_axypz(Ndim,1d0,dampA%x,divW%S%xz,tmpW%S%xz) 
    call Array_axypz(Ndim,1d0,dampA%y,divW%S%yx,tmpW%S%yx) 
    call Array_axypz(Ndim,1d0,dampA%y,divW%S%yy,tmpW%S%yy) 
    call Array_axypz(Ndim,1d0,dampA%y,divW%S%yz,tmpW%S%yz) 
    call Array_axypz(Ndim,1d0,dampA%z,divW%S%zx,tmpW%S%zx) 
    call Array_axypz(Ndim,1d0,dampA%z,divW%S%zy,tmpW%S%zy) 
    call Array_axypz(Ndim,1d0,dampA%z,divW%S%zz,tmpW%S%zz) 

    ! generate Sij**

    call strain_stress(pNp,mesh%Nele,material,Wss%E,divW%S)
    call DCOPY(Ndim3,Wss%V%array,1,divW%V%array,1)
    call domain_exchange(subdomain,PMLinfo%SRtype,divW%array,&
        subdomain%SRtype%send_req,subdomain%SRtype%recv_req,200)
    call MPI_Waitall(subdomain%SRtype%host%N_DD_Conn, &
        subdomain%SRtype%send_req,subdomain%SRtype%sta,ierr)
    call MPI_Waitall(subdomain%SRtype%gost%N_DD_Conn,&
        subdomain%SRtype%recv_req,subdomain%SRtype%sta,ierr)
    call DCOPY(Ndim3,divW%S%yz,1,divW%S%zy,1)
    call DCOPY(Ndim3,divW%V%array,1,Wss%V%array,1)

    ! Sij = Sij + Bii Sij**
    call Array_axypz(Ndim,1d0,dampB%x,divW%S%xx,tmpW%S%xx)
    call Array_axypz(Ndim,1d0,dampB%x,divW%S%xy,tmpW%S%xy)
    call Array_axypz(Ndim,1d0,dampB%x,divW%S%xz,tmpW%S%xz)
    call Array_axypz(Ndim,1d0,dampB%y,divW%S%yx,tmpW%S%yx)
    call Array_axypz(Ndim,1d0,dampB%y,divW%S%yy,tmpW%S%yy)
    call Array_axypz(Ndim,1d0,dampB%y,divW%S%yz,tmpW%S%yz)
    call Array_axypz(Ndim,1d0,dampB%z,divW%S%zx,tmpW%S%zx)
    call Array_axypz(Ndim,1d0,dampB%z,divW%S%zy,tmpW%S%zy)
    call Array_axypz(Ndim,1d0,dampB%z,divW%S%zz,tmpW%S%zz)

    ! apply surface source
    call surf_source(srcs%N_surfsrc,srcs%surfsrc,mesh%surf,&
        localtime)
    call rupture_update(mesh,subdomain,tmpW%V,tmpW%S,flag)
    call rupt_source(srcs%N_ruptsrc,srcs%ruptsrc,mesh%rupt,&
        mesh%rupt%dtau,0d0,localtime)
    call rupture_law(Nfp,mesh%rupt)

    call bilinear_Elastic3D_trE(matrix,mesh,subdomain,material,&
        srcs,rhsW,divW,tmpW,localtime,flag)

    ! Vj = Vj + damp_j Vj*
    call Array_axypz(Ndim3,1d0,PMLinfo%damp%x,Ws%V%x,tmpW%V%x) 

    call bilinear_Elastic3D_PML(matrix,mesh,subdomain,material,&
        srcs,rhsW,divW,tmpW,localtime,flag)

    ! rhsEij = rhsEij - (damp_i+damp_j) Eij
    call Array_axypz(Ndim3,-1d0,dampA%x,W%E%yz,rhsW%E%yz) 
    ! rhsEij = rhsEij - (damp_i*damp_j) Eij*
    call Array_axypz(Ndim3,-1d0,dampB%x,Ws%E%yz,rhsW%E%yz) 

    ! rhsEjj = rhsEjj - damp_j Ejj
    call Array_axypz(Ndim3,-1d0,PMLinfo%damp%x,W%E%xx,rhsW%E%xx) 

    ! rhsVj = rhsVj - a Vj
    call Array_axypz(Ndim,-1d0,dampS%x,W%V%x,rhsW%V%x) 
    call Array_axypz(Ndim,-1d0,dampS%x,W%V%y,rhsW%V%y) 
    call Array_axypz(Ndim,-1d0,dampS%x,W%V%z,rhsW%V%z) 

    ! rhsVj = rhsVj - b Vj*
    call Array_axypz(Ndim,-1d0,dampS%y,Ws%V%x,rhsW%V%x) 
    call Array_axypz(Ndim,-1d0,dampS%y,Ws%V%y,rhsW%V%y) 
    call Array_axypz(Ndim,-1d0,dampS%y,Ws%V%z,rhsW%V%z) 

    ! rhsVj = rhsVj - c Vj**
    call Array_axypz(Ndim,-1d0,dampS%z,Wss%V%x,rhsW%V%x) 
    call Array_axypz(Ndim,-1d0,dampS%z,Wss%V%y,rhsW%V%y) 
    call Array_axypz(Ndim,-1d0,dampS%z,Wss%V%z,rhsW%V%z) 

    call body_source(srcs%N_bodysrc,mesh%Nhele,srcs%bodysrc,&
        rhsW%V,rhsW%E,localtime)

end subroutine RHS_ElasticPML

!--------------------------------------------------------------------

subroutine UpdateElastic3DPML_EXRK4(matrix,mesh,subdomain,material,&
        srcs,rhsW,W,Ws,Wss,time,dt)
! inout:
    type(matrices)         :: matrix
    type(tetmesh_geometry) :: mesh
    type(tetmesh_Domain)   :: subdomain
    type(tetmesh_material) :: material
    type(sources)          :: srcs
    type(wavefield)        :: rhsW,W,Ws,Wss
    real(kind=rkind)       :: time,dt
! auxilary:
    real(kind=rkind) :: localtime
    integer :: k
    integer :: Ndim

    do k=1,5

        localtime=time+rk4c(k)*dt

        call RHS_ElasticPML(matrix,mesh,subdomain,material,&
            mesh%PMLinfo,srcs,rhsW,W,Ws,Wss,divW,tmpW,&
            localtime,mesh%EXflag)

        Ndim=mesh%Nele*pNp

        !resU = rk4a(rkstage)*resU 
        call wavefield_SCAL(Ndim,rk4a(k),resW)
        !resU = resU + delt*rhsU
        call wavefield_AXPY(Ndim,dt,rhsW,resW)
        !U = U + rk4b(rkstage)*resU
        call wavefield_AXPY(Ndim,rk4b(k),resW,W)

        !resU = rk4a(rkstage)*resU 
        call wavefield_SCAL(Ndim,rk4a(k),resWs)
        !resU = resU + delt*rhsU
        call wavefield_AXPY(Ndim,dt,W,resWs)
        !U = U + rk4b(rkstage)*resU
        call wavefield_AXPY(Ndim,rk4b(k),resWs,Ws)

        !resU = rk4a(rkstage)*resU 
        call wavefield_SCAL(Ndim,rk4a(k),resWss)
        !resU = resU + delt*rhsU
        call wavefield_AXPY(Ndim,dt,Ws,resWss)
        !U = U + rk4b(rkstage)*resU
        call wavefield_AXPY(Ndim,rk4b(k),resWss,Wss)

        if(trim(friction).eq.'aging')then
            mesh%rupt%rsi=rk4a(k)*mesh%rupt%rsi+dt*mesh%rupt%hsi
            mesh%rupt%psi=rk4b(k)*mesh%rupt%rsi+   mesh%rupt%psi
        elseif(trim(friction).eq.'linearSW')then
            mesh%rupt%rUt=rk4a(k)*mesh%rupt%rUt+dt*mesh%rupt%Vt
            mesh%rupt%dUt=rk4b(k)*mesh%rupt%rUt+   mesh%rupt%dUt
        endif

    enddo

end subroutine UpdateElastic3DPML_EXRK4

subroutine UpdateElastic3D_CFSPML_EXRK4(matrix,mesh,subdomain,&
        material,srcs,rhsW,W,Ws,Wss,Vsss,time,dt)
! inout:
    type(matrices)         :: matrix
    type(tetmesh_geometry) :: mesh
    type(tetmesh_Domain)   :: subdomain
    type(tetmesh_material) :: material
    type(sources)          :: srcs
    type(wavefield)        :: rhsW,W,Ws,Wss
    type(vector_array)     :: Vsss
    real(kind=rkind)       :: time,dt
! auxilary:
    real(kind=rkind) :: localtime
    integer :: k
    integer :: Ndim
    real(kind=rkind),parameter :: PMLr=40d0

    do k=1,5

        localtime=time+rk4c(k)*dt

        Ndim=mesh%Nele*pNp

        call RHS_ElasticPML(matrix,mesh,subdomain,material,&
            mesh%PMLinfo,srcs,rhsW,W,Ws,Wss,divW,tmpW,&
            localtime,mesh%EXflag)

!        ! rhsVj += \gamma a Vj*
!        call Array_axypz(Ndim,PMLr,dampS%x,Ws%V%x,rhsW%V%x) 
!        call Array_axypz(Ndim,PMLr,dampS%x,Ws%V%y,rhsW%V%y) 
!        call Array_axypz(Ndim,PMLr,dampS%x,Ws%V%z,rhsW%V%z) 
!        ! rhsVj += \gamma b Vj**
!        call Array_axypz(Ndim,PMLr,dampS%y,Wss%V%x,rhsW%V%x) 
!        call Array_axypz(Ndim,PMLr,dampS%y,Wss%V%y,rhsW%V%y) 
!        call Array_axypz(Ndim,PMLr,dampS%y,Wss%V%z,rhsW%V%z) 
!        ! rhsVj += \gamma c Vj***
!        call Array_axypz(Ndim,PMLr,dampS%z,Vsss%x,rhsW%V%x) 
!        call Array_axypz(Ndim,PMLr,dampS%z,Vsss%y,rhsW%V%y) 
!        call Array_axypz(Ndim,PMLr,dampS%z,Vsss%z,rhsW%V%z) 
!
!        ! rhsEij += \gamma (damp_i+damp_j) Eij*
!        call Array_axypz(Ndim*3,PMLr,dampA%x,Ws%E%yz,rhsW%E%yz) 
!        ! rhsEij += \gamma (damp_i damp_j) Eij**
!        call Array_axypz(Ndim*3,PMLr,dampB%x,Wss%E%yz,rhsW%E%yz) 
!        ! rhsEjj += \gamma damp_j Ejj*
!        call Array_axypz(Ndim*3,PMLr,mesh%PMLinfo%damp%x,&
!            Ws%E%xx,rhsW%E%xx) 

        !resU = rk4a(rkstage)*resU 
        call wavefield_SCAL(Ndim,rk4a(k),resW)
        !resU = resU + delt*rhsU
        call wavefield_AXPY(Ndim,dt,rhsW,resW)
        !U = U + rk4b(rkstage)*resU
        call wavefield_AXPY(Ndim,rk4b(k),resW,W)

        !resU = rk4a(rkstage)*resU 
        call wavefield_SCAL(Ndim,rk4a(k),resWs)
        !resU = resU + delt*rhsU
        call wavefield_AXPY(Ndim,dt,W,resWs)
        call wavefield_AXPY(Ndim,-PMLr*dt,Ws,resWs)
        !U = U + rk4b(rkstage)*resU
        call wavefield_AXPY(Ndim,rk4b(k),resWs,Ws)

        !resU = rk4a(rkstage)*resU 
        call wavefield_SCAL(Ndim,rk4a(k),resWss)
        !resU = resU + delt*rhsU
        call wavefield_AXPY(Ndim,dt,Ws,resWss)
        call wavefield_AXPY(Ndim,-PMLr*dt,Wss,resWss)
        !U = U + rk4b(rkstage)*resU
        call wavefield_AXPY(Ndim,rk4b(k),resWss,Wss)

        !resU = rk4a(rkstage)*resU 
        call vector_SCAL(Ndim,rk4a(k),resVsss)
        !resU = resU + delt*rhsU
        call vector_AXPY(Ndim,dt,Wss%V,resVsss)
        call vector_AXPY(Ndim,-PMLr*dt,Vsss,resVsss)
        !U = U + rk4b(rkstage)*resU
        call vector_AXPY(Ndim,rk4b(k),resVsss,Vsss)

        if(trim(friction).eq.'aging')then
            mesh%rupt%rsi=rk4a(k)*mesh%rupt%rsi+dt*mesh%rupt%hsi
            mesh%rupt%psi=rk4b(k)*mesh%rupt%rsi+   mesh%rupt%psi
        elseif(trim(friction).eq.'linearSW')then
            mesh%rupt%rUt=rk4a(k)*mesh%rupt%rUt+dt*mesh%rupt%Vt
            mesh%rupt%dUt=rk4b(k)*mesh%rupt%rUt+   mesh%rupt%dUt
        endif

    enddo

end subroutine UpdateElastic3D_CFSPML_EXRK4

!--------------------------------------------------------------------

end module ODEsolvePML_mod

