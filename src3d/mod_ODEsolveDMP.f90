!*******************************************************************!
!*  This module solves semi-discretized wave ODE                   *!
!*                                                                 *!
!*  Author: Ruichao Ye    Email: Ruichao.ye@gmail.com              *!
!*******************************************************************!

!--------------------------------------------------------------------
module ODEsolveDMP_mod
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
                             vector_array,tensor_array,&
                             wavefield,auxilary_array,&
                             ini_tensor,ini_vector,&
                             ini_wavefield,ini_auxilary_array,&
                             reset_wavefield,reset_auxilary_array,&
                             del_wavefield,del_auxilary_array,&
                             wavefield_AXPY,wavefield_COPY,&
                             wavefield_SCAL,Array_axypz
    use solve_mod,    only : bilinear_Elastic3D,strain_stress
    use ODEsolve_mod, only : RHS_Elastic3D
    use surface_mod,  only : rupture_update,rupture_law
    use source_mod,   only : body_source,rupt_source,surf_source
    use parallel_mod, only : domain_exchange

    implicit none

!    include 'mpif.h'

    type(wavefield) :: EX_W,IM_W,resW,resWs
    type(auxilary_array) :: divW,tmpW
    type(wavefield) :: damp

!--------------------------------------------------------------------
    contains
!--------------------------------------------------------------------

subroutine init_RKvariable_DMP(mesh,pNp,Nfp)
    type(tetmesh_geometry) :: mesh
    integer,intent(in) :: pNp,Nfp
    integer :: Ndim,i

    Ndim=mesh%Nele*pNp
    if(trim(timescheme).eq.'EXRK4')then
        call ini_wavefield(Ndim,resW)
    endif
    call ini_auxilary_array(tmpW,pNp,Nfp,mesh%Nele,.true.)
    call ini_auxilary_array(divW,pNp,Nfp,mesh%Nele,.true.)
    call ini_wavefield(Ndim,damp)
    do i=1,Ndim
        damp%V%x(i)=-max(&
            mesh%PMLinfo%damp%x(i),&
            mesh%PMLinfo%damp%y(i),&
            mesh%PMLinfo%damp%z(i))
    enddo
    call DCOPY(Ndim  ,damp%array(1:),1,damp%array(       1:),1)
    call DCOPY(Ndim*2,damp%array(1:),1,damp%array(2*Ndim+1:),1)
    call DCOPY(Ndim*4,damp%array(1:),1,damp%array(4*Ndim+1:),1)
    call DCOPY(Ndim  ,damp%array(1:),1,damp%array(8*Ndim+1:),1)
end subroutine init_RKvariable_DMP

subroutine del_RKvariable_DMP
    call del_wavefield(resW)
    call del_auxilary_array(tmpW)
    call del_auxilary_array(divW)
    call del_wavefield(damp)
end subroutine del_RKvariable_DMP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine UpdateElastic3D_DMP_EXRK4(matrix,mesh,subdomain,material,&
        srcs,rhsW,W,time,dt)
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
    integer :: k
    integer :: Ndim

    do k=1,5

        localtime=time+rk4c(k)*dt

        call RHS_Elastic3D(matrix,mesh,subdomain,material,srcs,&
            rhsW,W,divW,tmpW,localtime,mesh%EXflag)

        Ndim=mesh%Nele*pNp

        call Array_axypz(Ndim*9,1d0,damp%array,W%array,rhsW%array)

        !resU = rk4a(rkstage)*resU 
        call wavefield_SCAL(Ndim,rk4a(k),resW)
        !resU = resU + delt*rhsU
        call wavefield_AXPY(Ndim,dt,rhsW,resW)
        !U = U + rk4b(rkstage)*resU
        call wavefield_AXPY(Ndim,rk4b(k),resW,W)

        if(trim(friction).eq.'aging')then
            mesh%rupt%rsi=rk4a(k)*mesh%rupt%rsi+dt*mesh%rupt%hsi
            mesh%rupt%psi=rk4b(k)*mesh%rupt%rsi+   mesh%rupt%psi
        elseif(trim(friction).eq.'linearSW')then
            mesh%rupt%rUt=rk4a(k)*mesh%rupt%rUt+dt*mesh%rupt%Vt
            mesh%rupt%dUt=rk4b(k)*mesh%rupt%rUt+   mesh%rupt%dUt
        endif

    enddo

end subroutine UpdateElastic3D_DMP_EXRK4

end module ODEsolveDMP_mod

