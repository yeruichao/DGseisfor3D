!*******************************************************************!
!*  This module solves semi-discretized wave ODE                   *!
!*                                                                 *!
!*  Author: Ruichao Ye    Email: Ruichao.ye@gmail.com              *!
!*******************************************************************!

!--------------------------------------------------------------------
module ODEsolve_mod
!--------------------------------------------------------------------

	USE mpi
    use para_mod,     only : def_pOrder,pNp,Nfp,withrupt,friction,&
                             timescheme,rk4a,rk4b,rk4c,convtest,&
                             prestress,rupt_gamma,grav
    use datatype_mod, only : rkind,matrices,sources,tetmesh_domain,&
                             tetmesh_geometry,tetmesh_material,&
                             vector_array,tensor_array,&
                             wavefield,auxilary_array,&
                             ini_tensor,reset_tensor,tensor_COPY,&
                             ini_vector,vector_COPY,&
                             ini_wavefield,ini_auxilary_array,&
                             reset_wavefield,reset_auxilary_array,&
                             del_wavefield,del_auxilary_array,&
                             vector_AXPY,vector_COPY,vector_scal,&
                             tensor_AXPY,tensor_COPY,tensor_scal,&
                             wavefield_AXPY,wavefield_COPY,&
                             wavefield_SCAL,wavefield_INFNORM,&
                             Array_axypz
    use solve_mod,    only : bilinear_Elastic3D,strain_stress,&
                             flux_U,grav_T0Erhs
    use source_mod,   only : body_source,surf_source
    use parallel_mod, only : domain_exchange
    use surface_mod,  only : rupture_average
    use ruptsolve_mod,only : rupt_nonslip_RHS,form_rupt_matrices1,&
                             rupt_update_psi_RK4,&
                             rupt_Newton, rupt_on
    use grav_mod,     only : apply_grav_blocMV

    implicit none
    public :: VS_W,T0,resW

!    include 'mpif.h'

    type(wavefield) :: VS_W,EX_W,IM_W,resW
    real(kind=rkind),allocatable :: damp(:)
    type(auxilary_array) :: divW
    type(tensor_array) :: T0
    real(kind=rkind),allocatable :: psi0(:),sig0(:),Vt0(:)

!--------------------------------------------------------------------
    contains
!--------------------------------------------------------------------

subroutine IMEX_domains(mesh)
    type(tetmesh_geometry) :: mesh
    integer :: iele
    allocate(mesh%EXflag(mesh%Nhele))
    allocate(mesh%IMflag(mesh%Nhele))
    mesh%EXflag=.false.
    mesh%IMflag=.true. 
end subroutine IMEX_domains

subroutine init_RKvariable(mesh,pNp,Nfp,symm)
    type(tetmesh_geometry) :: mesh
    integer,intent(in) :: pNp,Nfp
    logical,intent(in) :: symm
    integer :: Ndim,i

    Ndim=mesh%Nele*pNp
    call ini_wavefield(Ndim,VS_W,symm)
    call ini_wavefield(Ndim,resW,symm)
    call reset_wavefield(resW)
!    call ini_wavefield(Ndim,EX_W,symm)
    call ini_wavefield(Ndim,IM_W,symm)
    call ini_auxilary_array(divW,pNp,Nfp,mesh%Nele,symm)
    if(mesh%PMLinfo%pml)then
        allocate(damp(Ndim))
        do i=1,Ndim
            damp(i)=-max(&
                mesh%PMLinfo%damp%x(i),&
                mesh%PMLinfo%damp%y(i),&
                mesh%PMLinfo%damp%z(i))
        enddo
    endif
    allocate(psi0(mesh%rupt%Nface*Nfp))
    allocate(sig0(mesh%rupt%Nface*Nfp))
    allocate(Vt0 (mesh%rupt%Nface*Nfp))
end subroutine init_RKvariable

subroutine del_RKvariable
    call del_wavefield(resW)
    call del_wavefield(VS_W)
    call del_auxilary_array(divW)
end subroutine del_RKvariable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine UpdateElastic3D_EXRK4(&
        matrix,mesh,subdomain,material,srcs,rhsW,W,time,dt)
! inout:
    type(matrices)         :: matrix
    type(tetmesh_geometry) :: mesh
    type(tetmesh_Domain)   :: subdomain
    type(tetmesh_material) :: material
    type(sources)           :: srcs
    type(wavefield)        :: rhsW,W
    real(kind=rkind)       :: time,dt
! auxilary:
    real(kind=rkind) :: localtime
    integer :: k,i,ierr
    integer :: Ndim,Nvar

    Ndim=mesh%Nele*pNp
    if(W%E%symmetric)then
        Nvar=9
    else
        Nvar=12
    endif

    do k=1,5

        localtime=time+rk4c(k)*dt

        call strain_stress(pNp,mesh%Nele,material,W%E,VS_W%S)
        call vector_COPY(Ndim,W%V,VS_W%V)

        if(srcs%N_surfsrc.gt.0.or.convtest)&
        call surf_source(srcs%N_surfsrc,srcs%surfsrc,mesh%surf,&
            localtime)

        call body_source(srcs%N_bodysrc,mesh%Nele*pNp,srcs%bodysrc,&
            VS_W%V,VS_W%S,localtime)

        call bilinear_Elastic3D(matrix,mesh,subdomain,material,&
            srcs,rhsW,divW,VS_W,localtime)

        if(withrupt)&
        call rupt_nonslip_RHS(pNp,Nfp,mesh,mesh%rupt,&
            W%E,W%V,rhsW%E,rhsW%V,dt)

        call domain_exchange(subdomain,subdomain%SRtype,rhsW%array,&
            subdomain%SRtype%send_req,subdomain%SRtype%recv_req,100)
        if(subdomain%SRtype%host%N_DD_Conn.gt.0)&
        call MPI_Waitall(subdomain%SRtype%host%N_DD_Conn, &
            subdomain%SRtype%send_req,subdomain%SRtype%sta,ierr)
        if(subdomain%SRtype%gost%N_DD_Conn.gt.0)&
        call MPI_Waitall(subdomain%SRtype%gost%N_DD_Conn,&
            subdomain%SRtype%recv_req,subdomain%SRtype%sta,ierr)

        if(mesh%PMLinfo%pml)then
            do i=0,Nvar-1
                call Array_axypz(Ndim,1d0,damp,&
                    W%array(i*Ndim+1:),rhsW%array(i*Ndim+1:))
            enddo
        endif

        !resU = rk4a(rkstage)*resU 
        call wavefield_SCAL(Ndim,rk4a(k),resW)
        !resU = resU + delt*rhsU
        call wavefield_AXPY(Ndim,dt,rhsW,resW)
        !U = U + rk4b(rkstage)*resU
        call wavefield_AXPY(Ndim,rk4b(k),resW,W)

        call vector_SCAL(Ndim,rk4a(k),resW%U)
        call vector_AXPY(Ndim,dt,W%V,resW%U)
        call vector_AXPY(Ndim,rk4b(k),resW%U,W%U)

    enddo

end subroutine UpdateElastic3D_EXRK4

!--------------------------------------------------------------------

subroutine UpdateElastic3D_Euler(&
        matrix,mesh,subdomain,material,srcs,rhsW,W,time,dt)
! inout:
    type(matrices)         :: matrix
    type(tetmesh_geometry) :: mesh
    type(tetmesh_Domain)   :: subdomain
    type(tetmesh_material) :: material
    type(sources)           :: srcs
    type(wavefield)        :: rhsW,W
    real(kind=rkind)       :: time,dt,res,rtmp,err,maxerr
    integer :: Ndim,Nhdim,Nvar,ierr,iter,i,j
    logical :: nocheat = .false.

    Ndim=mesh%Nele*pNp
    Nhdim=mesh%Nhele*pNp
    if(W%E%symmetric)then
        Nvar=9
    else
        Nvar=12
    endif

!!  iterative coupling  
!    call wavefield_COPY(Ndim,W,resW)
!    call wavefield_COPY(Ndim,W,IM_W)

!    if(mesh%rupt%Nface.gt.0)then
!        sig0=mesh%rupt%sigma
!        Vt0=mesh%rupt%dVt
!        psi0=mesh%rupt%psi
!    endif

!! ! start iterative coupling 
!do iter=1,4

    call reset_tensor(VS_W%E)
!    if(rupt_gamma.gt.0)then
!        call strain_stress(pNp,mesh%Nele,material,resW%E,VS_W%S)
!        call tensor_COPY(Ndim,W%E,resW%E)
!    endif
    call strain_stress(pNp,mesh%Nele,material,W%E,VS_W%S)
    call vector_COPY(Ndim,W%V,VS_W%V)

    if(srcs%N_surfsrc.gt.0.or.convtest)&
        call surf_source(srcs%N_surfsrc,srcs%surfsrc,mesh%surf,&
            time)

    if(srcs%N_bodysrc.gt.0)&
        call body_source(srcs%N_bodysrc,mesh%Nele*pNp,srcs%bodysrc,&
            VS_W%V,VS_W%S,time)

    call reset_wavefield(rhsW)
!    if(rupt_gamma.gt.0)then
!        call wavefield_AXPY(Ndim,-rupt_gamma,W,rhsW)
!    endif

    call bilinear_Elastic3D(matrix,mesh,subdomain,material,srcs,&
        rhsW,divW,VS_W,time)

    if(grav)then
!        call domain_exchange(subdomain,subdomain%SRtype1,W%U%array,&
!            subdomain%SRtype1%send_req,subdomain%SRtype1%recv_req,400)
!        if(subdomain%SRtype1%host%N_DD_Conn.gt.0)&
!        call MPI_Waitall(subdomain%SRtype1%host%N_DD_Conn, &
!            subdomain%SRtype1%send_req,subdomain%SRtype1%sta,ierr)
!        if(subdomain%SRtype1%gost%N_DD_Conn.gt.0)&
!        call MPI_Waitall(subdomain%SRtype1%gost%N_DD_Conn,&
!            subdomain%SRtype1%recv_req,subdomain%SRtype1%sta,ierr)

        call apply_grav_blocMV(pNp,mesh%Nhele,&
            rhsW%V%x,rhsW%V%y,rhsW%V%z,&
            W%U%x,W%U%y,W%U%z,matrix%G0)

!        call flux_U(pNp,Nfp,matrix,mesh,material%T0,material%rho,W%U,&
!                rhsW%V,dt)

        call grav_T0Erhs(pNp,mesh,matrix,material%rho,material%T0,&
                material%g0,W%E,rhsW%V)
    endif

    if(mesh%PMLinfo%pml)then
        do i=0,Nvar-1
            call Array_axypz(Ndim,1d0,damp,&
                W%array(i*Ndim+1:),rhsW%array(i*Ndim+1:))
        enddo
    endif

    if(convtest)then

        if(withrupt)then
            call rupt_nonslip_RHS(pNp,Nfp,mesh,mesh%rupt,&
                W%E,W%V,rhsW%E,rhsW%V,dt)
            call wavefield_COPY(Ndim,rhsW,resW)
            rupt_on=.false.
            call bilinear_Elastic3D(matrix,mesh,subdomain,material,&
                srcs,rhsW,divW,VS_W,time)
            rupt_on=.true.
            err=max(&
                maxval(abs(rhsW%V%x (1:Nhdim)-resW%V%x (1:Nhdim))),&
                maxval(abs(rhsW%V%y (1:Nhdim)-resW%V%y (1:Nhdim))),&
                maxval(abs(rhsW%V%z (1:Nhdim)-resW%V%z (1:Nhdim))))
            call MPI_REDUCE(err,maxerr,1,MPI_DOUBLE_PRECISION,&
                MPI_MAX,0,MPI_COMM_WORLD,ierr)
            if(subdomain%pid.eq.0.and.maxerr.gt.1d-12)&
                print*,'rhsVerr: ',maxerr
            err=max(&
                maxval(abs(rhsW%E%xx(1:Nhdim)-resW%E%xx(1:Nhdim))),&
                maxval(abs(rhsW%E%yx(1:Nhdim)-resW%E%yx(1:Nhdim))),&
                maxval(abs(rhsW%E%zx(1:Nhdim)-resW%E%zx(1:Nhdim))),&
                maxval(abs(rhsW%E%xy(1:Nhdim)-resW%E%xy(1:Nhdim))),&
                maxval(abs(rhsW%E%yy(1:Nhdim)-resW%E%yy(1:Nhdim))),&
                maxval(abs(rhsW%E%zy(1:Nhdim)-resW%E%zy(1:Nhdim))),&
                maxval(abs(rhsW%E%xz(1:Nhdim)-resW%E%xz(1:Nhdim))),&
                maxval(abs(rhsW%E%yz(1:Nhdim)-resW%E%yz(1:Nhdim))),&
                maxval(abs(rhsW%E%zz(1:Nhdim)-resW%E%zz(1:Nhdim))))
            call MPI_REDUCE(err,maxerr,1,MPI_DOUBLE_PRECISION,&
                MPI_MAX,0,MPI_COMM_WORLD,ierr)
            if(subdomain%pid.eq.0.and.maxerr.gt.1d-12)&
                print*,'rhsEerr: ',maxerr
        endif

        call wavefield_AXPY(Ndim,dt,rhsW,W)

        call domain_exchange(subdomain,subdomain%SRtype,W%array,&
            subdomain%SRtype%send_req,subdomain%SRtype%recv_req,100)
        if(subdomain%SRtype%host%N_DD_Conn.gt.0)then
            call MPI_Waitall(subdomain%SRtype%host%N_DD_Conn, &
            subdomain%SRtype%send_req,subdomain%SRtype%sta,ierr)
        endif
        if(subdomain%SRtype%gost%N_DD_Conn.gt.0)then
            call MPI_Waitall(subdomain%SRtype%gost%N_DD_Conn,&
            subdomain%SRtype%recv_req,subdomain%SRtype%sta,ierr)
        endif

    elseif(withrupt.and.prestress.gt.0)then

!        call tensor_COPY(Ndim,material%T0,T0)
        call reset_tensor(T0)
        if(srcs%N_ruptsrc.gt.0)&
            call body_source(srcs%N_ruptsrc,mesh%Nele*pNp,&
                srcs%ruptsrc,W%V,T0,time)

!!  iterative coupling  
!        call wavefield_COPY(Ndim,resW,W)

        call wavefield_AXPY(Ndim,dt,rhsW,W)

        call domain_exchange(subdomain,subdomain%SRtype,W%array,&
            subdomain%SRtype%send_req,subdomain%SRtype%recv_req,100)
        if(subdomain%SRtype%host%N_DD_Conn.gt.0)&
        call MPI_Waitall(subdomain%SRtype%host%N_DD_Conn, &
            subdomain%SRtype%send_req,subdomain%SRtype%sta,ierr)
        if(subdomain%SRtype%gost%N_DD_Conn.gt.0)&
        call MPI_Waitall(subdomain%SRtype%gost%N_DD_Conn,&
            subdomain%SRtype%recv_req,subdomain%SRtype%sta,ierr)

        call rupt_Newton(pNp,Nfp,mesh,mesh%rupt,dt,W,T0,matrix)

        call domain_exchange(subdomain,subdomain%SRtype,W%array,&
            subdomain%SRtype%send_req,subdomain%SRtype%recv_req,200)
        if(subdomain%SRtype%host%N_DD_Conn.gt.0)&
        call MPI_Waitall(subdomain%SRtype%host%N_DD_Conn, &
            subdomain%SRtype%send_req,subdomain%SRtype%sta,ierr)
        if(subdomain%SRtype%gost%N_DD_Conn.gt.0)&
        call MPI_Waitall(subdomain%SRtype%gost%N_DD_Conn,&
            subdomain%SRtype%recv_req,subdomain%SRtype%sta,ierr)

        call vector_AXPY(Ndim,dt,W%V,W%U)
!
!        call domain_exchange(subdomain,subdomain%SRtype1,W%U%array,&
!            subdomain%SRtype1%send_req,subdomain%SRtype1%recv_req,300)
!        if(subdomain%SRtype1%host%N_DD_Conn.gt.0)&
!        call MPI_Waitall(subdomain%SRtype1%host%N_DD_Conn, &
!            subdomain%SRtype1%send_req,subdomain%SRtype1%sta,ierr)
!        if(subdomain%SRtype1%gost%N_DD_Conn.gt.0)&
!        call MPI_Waitall(subdomain%SRtype1%gost%N_DD_Conn,&
!            subdomain%SRtype1%recv_req,subdomain%SRtype1%sta,ierr)

!        if(rupt_gamma.gt.0)then
!            call tensor_AXPY(Ndim,-1d0,W%E,resW%E)
!            call tensor_SCAL(Ndim,-rupt_gamma/dt,resW%E)
!        endif
    endif

!! check force balance
!    call strain_stress(pNp,mesh%Nele,material,W%E,VS_W%S)
!
!    call rupt_forcebalance_check(pNp,Nfp,mesh%rupt,mesh,W%E,VS_W%S,T0)

!    if(withrupt.and.mesh%rupt%Nface.gt.0)&
!        call rupt_update_psi_RK4(&
!            mesh%rupt%Nhface*Nfp,mesh%rupt,psi0,mesh%rupt%psi,&
!            Vt0,mesh%rupt%dVt,sig0,mesh%rupt%sigma,dt)

!!  iterative coupling  
!    rtmp=maxval(abs(W%array-IM_W%array))
!    call MPI_ALLREDUCE(rtmp,res,1,MPI_DOUBLE_PRECISION,&
!            MPI_MAX,MPI_COMM_WORLD,ierr)
!    if(res.lt.1e-10)exit
!    if(subdomain%pid.eq.0)write(*,'(A1\ )')'*'
!    call wavefield_COPY(Ndim,W,IM_W)
!enddo
!! end iterative coupling 


end subroutine UpdateElastic3D_Euler

!--------------------------------------------------------------------

end module ODEsolve_mod

