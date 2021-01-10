PROGRAM main

!*******************************************************************!
! 3D Discontinuous Galerkin acoustic-wave propagation
! Created by Ruichao Ye, March 7, 2012
! latest version on April 13, 2016
!*******************************************************************!

    use mpi
    use para_mod,       only: startTime,finalTime,dt,src_radius,&
                              Fnamelen,basename,sourcefname, &
                              snapfname,recordfname,&
                              recvfname,ruptrecvfname,&
                              snapintv,recordintv,rec_Ncomp,max_c,&
                              start_ckp,def_pOrder,pNp,Nfp,Nsele, &
                              refer_Cp,refer_Cs,refer_rho,refer_cc,&
                              PMLAlpha,Ndomain,task,job,friction,&
                              convtest,withrupt,verbose,logfid,&
                              para_init,symm,Nprocs,Ddecomp,myrank,&
                              rupt_gamma,grav
    use datatype_mod,   only: rkind,vector_array,wavefield,&
                              tetmesh_geometry,tetmesh_Domain,&
                              tetmesh_material,PML_geometry,&
                              matrices,sources,receivers,rupture,&
                              ini_wavefield,ini_vector,&
                              del_wavefield,del_vector,&
                              reset_wavefield,reset_vector,&
                              vector_SCAL,tensor_COPY,tensor_AXPY,&
                              ini_tensor,reset_tensor
    use Meshfile_mod,   only: readin_meshfile,init_material,&
                              nodes,tet,neigh,glob_Nele,glob_Nnode,&
                              xmax,xmin,ymax,ymin,zmax,zmin
    use geometry_mod,   only: Build_geometry
    use source_mod,     only: init_source,init_receiver,&
                              insert_source,point_initcond,&
                              recording,body_source,&
                              write_rec_to_file,first_arrival_time
    use domain_mod,     only: domain_decomp3D,metis_decomp3D
    use parallel_mod,   only: Domain_map3D,init_rupt_exchange_type,&
                              init_SRtype,check_map,check_blowup,&
                              My_barrier1,My_barrier2
    use surface_mod,    only: rupture_init,rupture_geometry_check,&
                              extsurf_init,rupt_recording,&
                              surface_Jac,check_rupt_map,&
                              rupture_crack_time,rupture_malloc
    use PML_mod,        only: init_PML
    use ODEsolve_mod,   only: init_RKvariable,IMEX_domains,&
                              UpdateElastic3D_EXRK4,&
                              UpdateElastic3D_Euler,&
                              VS_W,T0,resW
    use solve_mod,      only: strain_stress,flux_U_malloc
    use checkpoint_mod, only: subdomain_file,time_subdomain_file
    use conv_mpi_mod,   only: wavetype,conv_body,conv_model,conv_err,&
                              init_convrecv,conv_recording
    use ruptsolve_mod,  only: form_rupt_matrices0,&
                              form_rupt_matrices1,&
                              init_quasi_static,&
                              rupt_conv_err,&
                              rupt_Newton_allocate
    use PREM_mod,       only: prem_gravity 
    use grav_mod,       only: grav_blocMatrix,rupt_grav_blocMatrix
        
!    implicit none
!    include 'mpif.h'
!    For parameter_mod
    character*100 :: argv,parfile
    character*10 :: ctmp,ctmp1
    character (len=Fnamelen+10) :: filename
    character (len=Fnamelen+10) :: srcfname,recfname
    integer :: argc
!    for geometry_mod
    type(tetmesh_geometry)  :: mesh
    type(matrices)          :: matrix
    type(tetmesh_Domain)    :: subdomain
    type(tetmesh_material)  :: material
!    for Source_mod
    type(sources) :: srcs
    type(receivers) :: recvs
    integer :: irec,irec_rupt
!    auxilary
    integer ierr
    integer :: level
    logical alive
    double precision :: start_time0, start_time, finish_time
!    for solve_mod
    type(wavefield) :: rhsW,W,Ws,Wss
    type(vector_array) :: Vsss
    integer :: istep,iele
    real(kind=rkind) :: time
    real(kind=rkind) :: v1max,v2max,v3max,v1min,v2min,v3min
!    for checkpoint_mod
    real(kind=rkind) :: snaptime,recotime
    integer :: isnap
!    auxilary
    logical :: nonblowup=.true.
    integer :: i,j,ishot,NPshot,shot0,start_shot=1
    real(kind=rkind) :: l2err2(9),linfty(9),l2er(9),linf(9)
    integer :: itmp1,itmp2,itmp3
    real(kind=rkind) :: rtmp1,rtmp2,rtmp3
    integer :: Nvarcomp
    type(rupture) :: rupt_tmp

!--------------------------------------------------------------

    call MPI_Init(ierr)
    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr )
    call MPI_Comm_size( MPI_COMM_WORLD, Nprocs, ierr )
    subdomain%pid=myrank

    start_time=MPI_WTIME();start_time0=start_time

!    argc=iargc()
!    call getarg(argc,argv)
    call getarg(1,argv)
    read(argv,'(a)')parfile
    if(myrank.eq.0)then
	    write(*,*)parfile
	endif

    call para_init(parfile)

    subdomain%id=mod(myrank,Ndomain)+1
    NPshot=Nprocs/Ndomain
    shot0=myrank/Ndomain+1

    if(verbose)then
        write(ctmp,'(I6)')myrank
        ctmp=adjustl(ctmp)
        logfid=31900+myrank
        open(logfid,file='subdomain'//trim(ctmp)//'.info',&
            status='replace')
    endif

    if(trim(task).ne.'simu')then
        convtest=.true.
        wavetype=trim(task)
    endif

    if(trim(friction).eq.'aging'.or.&
       trim(friction).eq.'slipping')then
        withrupt=.true.
    else
        withrupt=.false.
    endif

!--------------------------------------------------------------

    if(start_ckp.le.0 .and. start_shot .eq.1)then
    
        call readin_meshfile(basename,glob_Nnode,glob_Nele,&
            subdomain)
        if(myrank.eq.0)then
            write(*,*)'readin mesh done'
            write(*,*)'physical domain: x=',xmin,'~',xmax
            write(*,*)'                 y=',ymin,'~',ymax
            write(*,*)'                 z=',zmin,'~',zmax
        endif
    
        level=log(real(Ndomain))/log(2.0)
!print*,Ddecomp
        call domain_decomp3D(level,Ddecomp)
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        if(myrank.eq.0)write(*,*)'domain decompose done'
    
        call Domain_map3D(glob_Nele,pNp,Ndomain,mesh,subdomain)
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        if(myrank.eq.0)write(*,*)'domain map done'
        deallocate(tet,neigh)
    
        call Build_geometry(glob_Nnode,mesh%Nele,nodes,&
            subdomain%local_tet,subdomain%local_neigh,&
            def_pOrder,pNp,Nfp,Nsele,mesh,matrix)
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        if(myrank.eq.0)write(*,*)'build geometry done'
    
        call init_material(mesh%Nele,pNp,Nfp,mesh,subdomain,material)
        if(convtest)call conv_model(pNp,mesh,material)
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        if(myrank.eq.0)write(*,*)'init material done'
    
        if(withrupt)then
            call rupture_init(basename,mesh%rupt,mesh)
            call init_rupt_exchange_type(mesh,subdomain,&
                mesh%rupt%SRtype)
            call init_quasi_static(pNp,Nfp,mesh,mesh%rupt,&
                material%T0,mesh%rupt%dVt)
            call MPI_Barrier(MPI_COMM_WORLD,ierr)
            if(myrank.eq.0)write(*,*)'Rupture initiation done'
        endif

        call init_PML(mesh,subdomain,mesh%PMLinfo)
        if(myrank.eq.0)write(*,*)'init PML done'
    
        if(subdomain%pid.lt.Ndomain)then
            call subdomain_file(snapfname,mesh,subdomain,material,&
                matrix,Ndomain,'save')
        endif
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        if(myrank.eq.0)write(*,*)'save subdomain infomation done'
    
        deallocate(nodes)
        deallocate(subdomain%local_tet,subdomain%local_neigh)
    
    else
    
        call subdomain_file(snapfname,mesh,subdomain,material,&
            matrix,Ndomain,'load')
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        if(myrank.eq.0)write(*,*)'load subdomain infomation done'
        call rupture_geometry_check(Nfp,mesh%rupt,mesh)
    
    endif

    call ini_tensor(pNp*mesh%Nele,T0,.true.)
    call reset_tensor(T0)
    call ini_vector(pNp*mesh%Nele,material%g0)
    call reset_vector(material%g0)

    if(grav)then
        call prem_gravity(pNp*mesh%Nele,&
            mesh%coord%x,mesh%coord%y,mesh%coord%z,&
            0d0,-6371d0,0d0,material%g0%x,material%g0%y,material%g0%z)
!        print*,maxval(material%g0%y),minval(material%g0%y)
        allocate(matrix%G0(3*pNp,3*pNp,mesh%Nele));matrix%G0=0d0
        call grav_blocMatrix(pNp,mesh,matrix,material%T0,material%g0,&
            material%rho)
        call flux_U_malloc(pNp,Nfp,material%T0,mesh)
    endif


!    print*,withrupt,mesh%rupt%Nface
    if(withrupt .and. mesh%rupt%Nface.gt.0)then
        call rupt_Newton_allocate(pNp,Nfp,mesh%rupt)
!        call surface_Jac(def_pOrder,Nfp,mesh%rupt%Nface,&
!            mesh%rupt%coord,mesh%rupt%Jac,mesh%rupt%invJac)
        if(.not.convtest)&
        call init_quasi_static(pNp,Nfp,mesh,mesh%rupt,&
            material%T0,mesh%rupt%dVt)
        call form_rupt_matrices0(pNp,Nfp,mesh,mesh%rupt,&
            material,matrix)
        call form_rupt_matrices1(pNp,Nfp,mesh,mesh%rupt,T0,matrix,dt)
        if(grav)then
            call rupt_grav_blocMatrix(pNp,Nfp,mesh,mesh%rupt,&
                matrix,material%T0)
        endif
    endif
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    if(withrupt .and. myrank.eq.0)&
        write(*,*)'initiate rupture matrices done'

!    print*,'----------------'
!    print*,mesh%rupt%Rmv0(:,1,1)
!    print*,'----------------'
!    stop

!     open(unit=9341,file='Rm0.dat',access='stream',status='replace')
!     write(9341)mesh%rupt%Rm0
!     close(9341)
!     open(unit=9341,file='Rp0.dat',access='stream',status='replace')
!     write(9341)mesh%rupt%Rp0
!     close(9341)
!     open(unit=9341,file='Rmv0.dat',access='stream',status='replace')
!     write(9341)mesh%rupt%Rmv0
!     close(9341)
!     open(unit=9341,file='Rpv0.dat',access='stream',status='replace')
!     write(9341)mesh%rupt%Rpv0
!     close(9341)
!     open(unit=9341,file='Rm1.dat',access='stream',status='replace')
!     write(9341)mesh%rupt%Rm1
!     close(9341)
!     open(unit=9341,file='Rp1.dat',access='stream',status='replace')
!     write(9341)mesh%rupt%Rp1
!     close(9341)
!     open(unit=9341,file='Rmv1.dat',access='stream',status='replace')
!     write(9341)mesh%rupt%Rmv1
!     close(9341)
!     open(unit=9341,file='Rpv1.dat',access='stream',status='replace')
!     write(9341)mesh%rupt%Rpv1
!     close(9341)

!--------------------------------------------------------------

    if(myrank.eq.0)then
        write(*,*)'Start Time = ',startTime
        write(*,*)'Final Time = ',finalTime
        write(*,*)'Time step  = ',dt
        write(*,*)'Iteration  = ',aint((finalTime - startTime)/dt)
    endif

!    call My_barrier1(myrank)
!    print*,'pid=',myrank
!    print*,material%job
!    call My_barrier2(myrank)

    if(symm)then
        Nvarcomp=9
    else
        Nvarcomp=12
    endif

    call init_SRtype(subdomain%SRtype,pNp,mesh%Nele,Nvarcomp)

    call ini_wavefield(mesh%Nele*pNp,W   ,symm)
    call reset_wavefield(W)
    call ini_wavefield(mesh%Nele*pNp,rhsW,symm)
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    if(myrank.eq.0)write(*,*)'init wavefield done'

!--------------------------------------------------------------

    call DCOPY(mesh%Nele*pNp*3,mesh%coord%array,1,W%V%array,1)
    call DCOPY(mesh%Nele*pNp*3,mesh%coord%array,1,W%E%array,1)
    call DCOPY(mesh%Nele*pNp*3,mesh%coord%array,1,W%E%yz   ,1)
    call DCOPY(mesh%Nele*pNp*3,mesh%coord%array,1,W%E%zy   ,1)
    call check_map(subdomain,subdomain%SRtype,&
        W%array,rhsW%array,mesh%Nele*pNp,mesh%Nhele*pNp,Nvarcomp)
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    if(myrank.eq.0)write(*,*)'check mesh communication done'

    call check_rupt_map(Nfp,subdomain,mesh)
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    if(myrank.eq.0)write(*,*)'check rupture communication done'

    if(mesh%PMLinfo%pml)then
        call vector_SCAL(mesh%Nele*pNp,PMLAlpha,mesh%PMLinfo%damp)
!        do i=1,mesh%Nele*pNp
!            if( mesh%PMLinfo%damp%x(i).le.0d0 .and. &
!                mesh%PMLinfo%damp%y(i).le.0d0 .and. &
!                mesh%PMLinfo%damp%z(i).le.0d0 ) then
!                W%V%x (i)=0d0;W%V%y (i)=0d0;W%V%z (i)=0d0
!                W%E%xx(i)=0d0;W%E%yy(i)=0d0;W%E%zz(i)=0d0
!                W%E%yz(i)=0d0;W%E%xz(i)=0d0;W%E%xy(i)=0d0
!            endif
!        enddo
!        call reset_wavefield(rhsW)
!        call check_map(subdomain,mesh%PMLinfo%SRtype,&
!            W%array,rhsW%array,mesh%Nele*pNp,mesh%Nhele*pNp,Nvarcomp)
    endif
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    if(myrank.eq.0)write(*,*)'check PML communication done'

    call reset_wavefield(W)
    call reset_wavefield(rhsW)
!    if(mesh%PMLinfo%pml)then
!        call ini_wavefield(mesh%Nele*pNp,Ws ,symm)
!        call ini_wavefield(mesh%Nele*pNp,Wss,symm)
!        call ini_vector(mesh%Nele*pNp,Vsss)
!        call reset_wavefield(Ws)
!        call reset_wavefield(Wss)
!        call reset_vector(Vsss)
!    endif

!--------------------------------------------------------------

    if(convtest)then
        call init_convrecv(subdomain%id,startTime,finalTime,snapintv,&
            recordintv,recvs)
    else
        call init_receiver(recvfname,&
            matrix,startTime,finalTime,snapintv,recordintv,&
            mesh%coord,mesh%Nele,mesh%Nhele,&
            mesh%coord,mesh%Nele,mesh%Nhele,&
            recvs,rec_Ncomp,3,subdomain)
        call init_receiver(ruptrecvfname,&
            matrix,startTime,finalTime,snapintv,recordintv,&
            mesh%rupt%coord,mesh%rupt%Nface,mesh%rupt%Nhface,&
            mesh%coord,mesh%Nele,mesh%Nhele,&
            mesh%rupt%recvs,9,2,subdomain)
    endif
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    if(myrank.eq.0)write(*,*)'init receiver done,'

!--------------------------------------------------------------

    allocate(srcs%t_direct(mesh%Nele))
    if(convtest)then
        call extsurf_init(pNp,Nfp,mesh)
        allocate(srcs%surfsrc(0))
        allocate(srcs%bodysrc(0))
        srcs%t_direct=0d0
    else
        srcs%t_direct=FinalTime
        call init_source(sourcefname,srcs%N_bodysrc,srcs%bodysrc,&
            mesh%Nele*pNp,mesh%coord,src_radius,&
            mesh,srcs%t_direct,max_c,StartTime)
        if(myrank.eq.0)then
            write(*,*)'Totally ', srcs%N_bodysrc,' srcs found.'
        endif
    endif
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    if(myrank.eq.0)write(*,*)'init source done'

if(.not.convtest .and. withrupt)then
    if(srcs%N_bodysrc.gt.0)then
        if(mesh%rupt%Nface.gt.0)then
        srcs%N_ruptsrc=srcs%N_bodysrc
        allocate(srcs%ruptsrc(srcs%N_ruptsrc))
        do i=1,srcs%N_bodysrc
            srcs%ruptsrc(i)%stime  =srcs%bodysrc(i)%stime
            srcs%ruptsrc(i)%ftime  =srcs%bodysrc(i)%ftime
            srcs%ruptsrc(i)%x      =srcs%bodysrc(i)%x
            srcs%ruptsrc(i)%y      =srcs%bodysrc(i)%y
            srcs%ruptsrc(i)%z      =srcs%bodysrc(i)%z
            srcs%ruptsrc(i)%q      =srcs%bodysrc(i)%q
            srcs%ruptsrc(i)%lag    =srcs%bodysrc(i)%lag
            srcs%ruptsrc(i)%st_type=srcs%bodysrc(i)%st_type
            srcs%ruptsrc(i)%amp    =srcs%bodysrc(i)%amp
            srcs%ruptsrc(i)%freq   =srcs%bodysrc(i)%freq
!        print*,srcs%ruptsrc(i)%st_type
            call insert_source(&
                srcs%bodysrc(i)%x,&
                srcs%bodysrc(i)%y,&
                srcs%bodysrc(i)%z,&
                mesh%Nele*pNp,&
                mesh%coord,src_radius,&
                srcs%ruptsrc(i)%Nnod,&
                srcs%ruptsrc(i)%src2nod,&
                srcs%ruptsrc(i)%sweight)
        enddo
        endif
        srcs%N_bodysrc=0
    endif
!    srcs%t_direct=0d0
endif

!--------------------------------------------------------------
! Check quasi-static status
!    if(withrupt.and..not.convtest)then
!        call rupt_form_RHS(pNp,Nfp,mesh%rupt,mesh,W%E,W%V,&
!                material%T0,dt)
!        call rupt_Newton_solve(Nfp,mesh%rupt,dt)
!    endif
!--------------------------------------------------------------

    if(start_ckp.le.0)then
        time=startTime;istep=0;isnap=0;irec=0;irec_rupt=0
        snaptime=startTime
        recotime=startTime
        if(convtest)then
            call conv_body(pNp,mesh,W%V,W%E,StartTime)
            if(.not.symm)&
                call DCOPY(mesh%Nele*pNp*3,W%E%yz,1,W%E%zy,1)
        else
            call point_initcond(srcs%N_bodysrc,mesh%Nele,&
                srcs%bodysrc,W%V,W%E)
        endif
    else
        isnap=start_ckp
        call time_subdomain_file(snapfname,mesh,&
            isnap,time,istep,irec,W,Ws,Wss,Vsss,&
            subdomain%id-1,'load')
        irec_rupt=irec
        if(myrank.eq.0)print*,time,istep,irec,irec_rupt
        snaptime=time+snapintv
        recotime=time+recordintv
    endif

if(convtest .and. withrupt .and. mesh%rupt%Nhface.gt.0)then
!    do i=1,mesh%rupt%Nface*Nfp
!        mesh%rupt%Vt%array(i)=W%V%array(i)
!    enddo
    mesh%rupt%Vt%array=0d0
    mesh%rupt%dVt=0d0
endif

    finish_time=MPI_WTIME()
    if(myrank.eq.0)&
        write(*,*)'initiation time: ',finish_time-start_time
    if(myrank.eq.0 .and. withrupt)&
        write(*,*)'rupture gamma: ',rupt_gamma*dt
    start_time=finish_time

    call IMEX_domains(mesh)
    call init_RKvariable(mesh,pNp,Nfp,symm)

!    if(start_ckp.gt.0)then
!        if(srcs%N_ruptsrc.gt.0)&
!            call body_source(srcs%N_ruptsrc,mesh%Nele*pNp,&
!                srcs%ruptsrc,W%V,T0,time)
!        call strain_stress(pNp,mesh%Nele,material,W%E,VS_W%S)
!        call tensor_AXPY(mesh%Nele*pNp,1d0,material%T0,T0)
!        call tensor_AXPY(mesh%Nele*pNp,1d0,VS_W%S,T0)
!        call grav_blocMatrix(pNp,mesh,matrix,T0,material%g0,&
!            material%rho)
!        call rupt_grav_blocMatrix(pNp,Nfp,mesh,mesh%rupt,&
!            matrix,T0)
!    endif

    
    do while(time.le.finalTime)

        if(time .ge. snaptime .and. snapintv .gt. 0)then

            call write_rec_to_file(recordfname,recvs,irec)
            call write_rec_to_file(trim(recordfname)//'rupt',&
                mesh%rupt%recvs,irec_rupt)

            isnap=isnap+1
            snaptime=snaptime+snapintv
            call time_subdomain_file(snapfname,mesh,&
                isnap,time,istep,irec,W,Ws,Wss,Vsss,&
                subdomain%id-1,'save')
            if(myrank.eq.0)then
                write(*,*)
                print*,time,istep,irec,irec_rupt
            endif


            call MPI_REDUCE(maxval(W%V%x),v1max,1,&
                MPI_DOUBLE_PRECISION,&
                MPI_MAX,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(maxval(W%V%y),v2max,1,&
                MPI_DOUBLE_PRECISION,&
                MPI_MAX,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(maxval(W%V%z),v3max,1,&
                MPI_DOUBLE_PRECISION,&
                MPI_MAX,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(minval(W%V%x),v1min,1,&
                MPI_DOUBLE_PRECISION,&
                MPI_MIN,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(minval(W%V%y),v2min,1,&
                MPI_DOUBLE_PRECISION,&
                MPI_MIN,0,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(minval(W%V%z),v3min,1,&
                MPI_DOUBLE_PRECISION,&
                MPI_MIN,0,MPI_COMM_WORLD,ierr)

            finish_time=MPI_WTIME()
            if(myrank.eq.0)then
                write(*,*)
                write(*,*)'Snapshot No.',isnap,' complete.'
                write(*,*)'MPI wall time: ', finish_time-start_time
                write(*,*)'v1: ',v1min,v1max
                write(*,*)'v2: ',v2min,v2max
                write(*,*)'v3: ',v3min,v3max
            endif
            start_time=finish_time
!    if(grav)then
!        print*,maxval(rhsW%V%array),minval(rhsW%V%array),&
!            maxval(resW%V%array),minval(resW%V%array)
!    endif

        endif

        if(recordintv.gt.0 .and. time .ge. recotime)then
            recotime=recotime+recordintv
            if(convtest)then
                call conv_err(mesh%Nhele,mesh%detJ,time,matrix%M3D,&
                    mesh%coord%x,mesh%coord%y,mesh%coord%z,&
                    W%V%x,W%V%y,W%V%z,W%E%xx,W%E%yy,W%E%zz,&
                    W%E%yz,W%E%xz,W%E%xy,l2err2,linfty)
                call MPI_REDUCE(l2err2,l2er,9,MPI_DOUBLE_PRECISION,&
                    MPI_SUM,0,MPI_COMM_WORLD,ierr)
                call MPI_REDUCE(linfty,linf,9,MPI_DOUBLE_PRECISION,&
                    MPI_max,0,MPI_COMM_WORLD,ierr)
                if(subdomain%pid.eq.0)then
                    l2er=sqrt(l2err2)
                    print*,l2er(1:3)
                    call conv_recording(recvs,l2er,linf)
                endif
                if(withrupt)then
                    linf=0d0
                    call strain_stress(pNp,mesh%Nele,material,&
                        W%E,VS_W%S)
                    call rupt_conv_err(pNp,Nfp,&
                        mesh,mesh%rupt,W%E,VS_W%S,time,linfty)
                    call MPI_REDUCE(linfty,linf,4,&
                        MPI_DOUBLE_PRECISION,&
                        MPI_MAX,0,MPI_COMM_WORLD,ierr)
                    if(subdomain%pid.eq.0.and.maxval(linf).gt.1d-14)&
                        print*,linf(1:4)
                endif
            else
                call recording(recvs,W%V,W%E)
                call rupt_recording(mesh%rupt%recvs,mesh%rupt,Nfp)
                if(mesh%rupt%Nhface.gt.0)&
                    call rupture_crack_time(mesh%rupt%Nhface*Nfp,&
                        mesh%rupt%dVt,mesh%rupt%crack_t,time)
                if(nonblowup .and. &
                    check_blowup(subdomain,W%V%array))then
                    if(myrank.eq.0)print*,'wavefield blowup at',istep
                    nonblowup=.false.
                endif
            endif

!        if(grav)then
!            call tensor_AXPY(mesh%Nele*pNp,1d0,material%T0,T0)
!            call tensor_AXPY(mesh%Nele*pNp,1d0,VS_W%S,T0)
!            call grav_blocMatrix(pNp,mesh,matrix,T0,material%g0,&
!                material%rho)
!            call rupt_grav_blocMatrix(pNp,Nfp,mesh,mesh%rupt,&
!                matrix,T0)
!        endif

        endif

!        call UpdateElastic3D_EXRK4(&
        call UpdateElastic3D_Euler(&
            matrix,mesh,subdomain,material,srcs,&
            rhsW,W,time,dt)

        time=time+dt;istep=istep+1

        !if(myrank.eq.0)write(*,'(A1\ )')'.'
        if(myrank.eq.0)write(*,'(A1)')'.'

    enddo

    finish_time=MPI_WTIME()
    if(myrank.eq.0)then
        write(*,*)
        write(*,*)'Total run time: ', finish_time-start_time0
    endif

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call MPI_Finalize(ierr)

end PROGRAM main
