!********************************************************************!
!*  This module converts inputed mesh file to DG mesh structure     *!
!*                                                                  *!
!*  Author: Ruichao Ye    Email: Ruichao.ye@gmail.com               *!
!********************************************************************!

!--------------------------------------------------------------------!
module Meshfile_mod
!--------------------------------------------------------------------!

    use datatype_mod, only : rkind,tetmesh_geometry,tetmesh_material,&
                             tetmesh_Domain,tensor_array,ini_tensor
    use para_mod,     only : Fnamelen,basename,def_pOrder,job,pNp,&
                             refer_Cp,refer_Cs,refer_rho,refer_cc,&
                             PMLHx1,PMLHy1,PMLHz1,&
                             PMLHx2,PMLHy2,PMLHz2,&
                             PMLx1_p,PMLy1_p,PMLz1_p,&
                             PMLx2_p,PMLy2_p,PMLz2_p,&
                             PMLx1_type,PMLy1_type,PMLz1_type,&
                             PMLx2_type,PMLy2_type,PMLz2_type,&
                             pmlgeometry,PMLM,withrupt,prestress,&
                             refer_T0,verbose,logfid,Nprocs

    implicit none
    public :: readin_meshfile
    public :: nodes,tet,neigh,F2T, triface, Fatt
    public :: xmax, xmin, ymax, ymin, zmax, zmin, &
              glob_Nele, glob_Nnode, glob_Nface
    private :: alive,error
!--------------------------------------------------------------------!
    include 'mpif.h'
    real(kind=rkind) :: xmax, xmin, ymax, ymin, zmax, zmin
    integer,allocatable :: tet(:,:),neigh(:,:)
    integer,allocatable :: triface(:,:),F2T(:,:),Fatt(:)
    real(kind=rkind),allocatable :: nodes(:,:)
    logical :: alive
    integer :: glob_Nele,glob_Nnode,glob_Nface
    integer :: error
    real(kind=rkind),parameter :: mat_TOL=1d-4

!--------------------------------------------------------------------!
    contains
!--------------------------------------------------------------------!

subroutine readin_meshfile(basename,Nnode,Nele,subdomain)
    ! Read in mesh files in Tetgen format
    ! input
    character(len=*),intent(in) :: basename
    type(tetmesh_Domain)   :: subdomain
    ! output
    integer :: Nele,Nnode
    ! auxilary
    character(len=Fnamelen) :: filename
    integer :: fid
    integer :: tmp1,tmp2,tmp3,tmp4,i,ierr
    real :: rtmp1,rtmp2,rtmp3
    integer :: mproc1,mproc2,mproc3,mproc4

    ! ............... Read in nodes ............... !
    mproc1=max(0,Nprocs-1)

    fid=1002
    filename=trim(basename)//'.node'
    inquire(file=trim(filename),exist=alive)
    if(.not. alive)then
        write(*,*)'File "',trim(filename),&
            '" does not exist.'
        stop
    endif
    open(fid,file=trim(filename),status='old',action='read',&
        position='REWIND')
    read(fid,*,iostat=error)Nnode,tmp1,tmp2,tmp3
    if(Nnode .lt. 1)then
        write(*,*)'Something wrong with nodes file.'
        stop
    else if(tmp1.ne.3)then
        write(*,*)'The model is not 3-Dimensional.'
        stop
    endif
    allocate(nodes(3,Nnode))

    if(subdomain%pid.eq.mproc1)then
    
        if(tmp3.gt.0)then
            do i=1,Nnode
                read(fid,*,iostat=error)tmp1,rtmp1,rtmp2,rtmp3,tmp2
                nodes(1,i)=rtmp1;nodes(2,i)=rtmp2;nodes(3,i)=rtmp3
            end do
        else
            do i=1,Nnode
                read(fid,*,iostat=error)tmp1,rtmp1,rtmp2,rtmp3
                nodes(1,i)=rtmp1;nodes(2,i)=rtmp2;nodes(3,i)=rtmp3
            enddo
        endif
        write(*,*)'Read in nodes done.'
    endif
    close(fid)

    ! ............... Read in tetrahedras ............... !
    mproc2=max(0,Nprocs-2)

    fid=1003
    filename=trim(basename)//'.ele'
    inquire(file=trim(filename),exist=alive)
    if(.not. alive)then
        write(*,*)'Error 1003: File "',trim(filename),&
            '" does not exist.'
        stop
    endif
    open(fid,file=trim(filename),status='old',action='read',&
        position='REWIND')
    read(fid,*,iostat=error)Nele,tmp1,tmp2
    if(Nele.lt.1)then
        write(*,*)'Something wrong with elements file.'
        stop
    else if(tmp1.ne.4)then
        write(*,*)'First order tetrahedron element required.'
        stop
    endif
    glob_Nele=Nele
    allocate(tet(4,Nele))

    if(subdomain%pid.eq.mproc2)then
    
        if(tmp2.gt.0)then
            do i=1,Nele
                read(fid,*,iostat=error)tmp1,tet(1,i),tet(2,i),&
                    tet(3,i),tet(4,i),tmp3
            enddo
        else
            do i=1,Nele
                read(fid,*,iostat=error)tmp1,tet(1,i),tet(2,i),&
                    tet(3,i),tet(4,i)
            end do
        endif
        write(*,*)'Read in elements done.'

    endif
    close(fid)

    ! ............... Read in neighbours ............... !
    mproc3=max(0,Nprocs-3)

    allocate(neigh(4,Nele))

    if(subdomain%pid.eq.mproc3)then
    
        fid=1004
        filename=trim(basename)//'.neigh'
        inquire(file=trim(filename),exist=alive)
        if(.not. alive)then
            write(*,*)'Error 1005: File "',trim(filename),&
                '" does not exist.'
            stop
        endif
        open(fid,file=trim(filename),status='old',action='read',&
            position='REWIND')
        read(fid,*,iostat=error)tmp1,tmp2
        if(Nele.ne.tmp1)then
            write(*,*)'Something wrong with neigh file.'
            stop
        else if(tmp2.ne.4)then
            write(*,*)'First order tetrahedron element required.'
            stop
        endif
        do i=1,Nele
            read(fid,*,iostat=error)tmp1,neigh(3,i),neigh(4,i),&
                neigh(2,i),neigh(1,i)
        end do
        write(*,*)'Read in neighbours done.'
        close(fid)
    endif

    ! ............... Read in face identifiers ............... !
    mproc4=max(0,Nprocs-4)

    if(withrupt)then
        fid=1005
        filename=trim(basename)//'.face'
        inquire(file=trim(filename),exist=alive)
        if(.not. alive)then
            write(*,*)'Error 1005: File "',trim(filename),&
                '" does not exist.'
            stop
        endif
        open(fid,file=trim(filename),status='old',action='read',&
            position='REWIND')
        read(fid,*,iostat=error)glob_Nface
        allocate(triface(3,glob_Nface),F2T(2,glob_Nface))
        allocate(Fatt(glob_Nface))

        if(subdomain%pid.eq.mproc4)then
            do i=1,glob_Nface
                read(fid,*,iostat=error)&
                    triface(1,i),triface(2,i),triface(3,i),&
                    F2T(1,i),F2T(2,i),Fatt(i)
            end do
            write(*,*)'Read in rupture facets done.'
        endif
        close(fid)
        call MPI_BCAST(triface,3*glob_Nface,MPI_INTEGER,mproc4,&
            MPI_COMM_WORLD,ierr)
        call MPI_BCAST(F2T,2*glob_Nface,MPI_INTEGER,mproc4,&
            MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Fatt,glob_Nface,MPI_INTEGER,mproc4,&
            MPI_COMM_WORLD,ierr)
    endif

    call MPI_BCAST(nodes,3*Nnode,MPI_DOUBLE_PRECISION,mproc1,&
        MPI_COMM_WORLD,ierr)
    call MPI_BCAST(tet,4*glob_Nele,MPI_INTEGER,mproc2,&
        MPI_COMM_WORLD,ierr)
    call MPI_BCAST(neigh,4*glob_Nele,MPI_INTEGER,mproc3,&
        MPI_COMM_WORLD,ierr)

    xmin=minval(nodes(1,:));xmax=maxval(nodes(1,:))
    ymin=minval(nodes(2,:));ymax=maxval(nodes(2,:))
    zmin=minval(nodes(3,:));zmax=maxval(nodes(3,:))

    ! .............. generate rotational PML .................. !

    if(trim(pmlgeometry) .eq. 'rotate')then
        PMLM( 1)=rotation_PML_bndvalue(&
            PMLx1_type,PMLx1_p,Nnode,&
            nodes(1,:),nodes(2,:),nodes(3,:),.true.)
        PMLM( 2)=rotation_PML_bndvalue(&
            PMLx2_type,PMLx2_p,Nnode,&
            nodes(1,:),nodes(2,:),nodes(3,:),.true.)
        PMLM( 3)=rotation_PML_bndvalue(&
            PMLy1_type,PMLy1_p,Nnode,&
            nodes(1,:),nodes(2,:),nodes(3,:),.true.)
        PMLM( 4)=rotation_PML_bndvalue(&
            PMLy2_type,PMLy2_p,Nnode,&
            nodes(1,:),nodes(2,:),nodes(3,:),.true.)
        PMLM( 5)=rotation_PML_bndvalue(&
            PMLz1_type,PMLz1_p,Nnode,&
            nodes(1,:),nodes(2,:),nodes(3,:),.true.)
        PMLM( 6)=rotation_PML_bndvalue(&
            PMLz2_type,PMLz2_p,Nnode,&
            nodes(1,:),nodes(2,:),nodes(3,:),.true.)
        PMLM( 7)=rotation_PML_bndvalue(&
            PMLx1_type,PMLx1_p,Nnode,&
            nodes(1,:),nodes(2,:),nodes(3,:),.false.)
        PMLM( 8)=rotation_PML_bndvalue(&
            PMLx2_type,PMLx2_p,Nnode,&
            nodes(1,:),nodes(2,:),nodes(3,:),.false.)
        PMLM( 9)=rotation_PML_bndvalue(&
            PMLy1_type,PMLy1_p,Nnode,&
            nodes(1,:),nodes(2,:),nodes(3,:),.false.)
        PMLM(10)=rotation_PML_bndvalue(&
            PMLy2_type,PMLy2_p,Nnode,&
            nodes(1,:),nodes(2,:),nodes(3,:),.false.)
        PMLM(11)=rotation_PML_bndvalue(&
            PMLz1_type,PMLz1_p,Nnode,&
            nodes(1,:),nodes(2,:),nodes(3,:),.false.)
        PMLM(12)=rotation_PML_bndvalue(&
            PMLz2_type,PMLz2_p,Nnode,&
            nodes(1,:),nodes(2,:),nodes(3,:),.false.)
        if(verbose .and. logfid .gt. 0)&
            write(logfid,*)'PMLM=',PMLM

    endif

    return

end subroutine readin_meshfile

subroutine init_material(Nele,pNp,Nfp,mesh,subdomain,material)
    integer,intent(in) :: Nele
    integer,intent(in) :: pNp,Nfp
    type(tetmesh_material) :: material
    type(tetmesh_geometry) :: mesh
    type(tetmesh_Domain)   :: subdomain

    character(len=Fnamelen) :: filename
    character(len=6) :: ctmp
    integer :: rec_len, fidcp, fidcs, fidcc, fidrho, fidT0, error
    integer :: i,j,k
    real(kind = rkind),pointer :: Cp(:),Cs(:),rho(:),CC(:,:)
    real :: rtmp=0.0
    inquire(iolength=rec_len)rtmp

    write(ctmp,'(I6)')def_pOrder
    ctmp=adjustl(ctmp)

    allocate(material%rho(pNp*Nele));rho=>material%rho
    allocate(material%Cp (pNp*Nele));Cp =>material%Cp
    allocate(material%Cs (pNp*Nele));Cs =>material%Cs
    if(job==0)then
        allocate(material%C(pNp*Nele,1));CC=>material%C
        allocate(material%k_media(Nele))
    elseif(job==1)then
        allocate(material%C(pNp*Nele,2));CC=>material%C
        allocate(material%k_media(Nele))
    elseif(job==2)then
        allocate(material%C(pNp*Nele,9));CC=>material%C
        allocate(material%k_media(Nele))
    else
        allocate(material%C(pNp*Nele,21));CC=>material%C
        allocate(material%k_media(Nele))
    endif

    fidrho=1005
    filename=trim(basename)//'.rho.p'//trim(ctmp)
    inquire(file=trim(filename),exist=alive)
    if(.not. alive)then
        if(subdomain%pid.eq.0)&
            write(*,*)'File "',trim(filename),&
            '" does not exist, assuming rho=',refer_rho
        rho=refer_rho
        fidrho=-1
    else
        open( unit = fidrho, file = filename, &
            form = "unformatted", status = "old", &
            access = "direct", recl = rec_len )
    endif

    fidcp=-1
    if(job.le.1)then
        fidcp=1006
        filename=trim(basename)//'.cp.p'//trim(ctmp)
        inquire(file=trim(filename),exist=alive)
        if(.not. alive)then
            if(subdomain%pid.eq.0)&
                write(*,*)'File "',trim(filename),&
                '" does not exist, assuming Cp=',refer_Cp
            Cp=refer_Cp
            fidcp=-1
        else
            open( unit = fidcp, file = filename, &
                form = "unformatted", status = "old", &
                access = "direct", recl = rec_len )
        endif
    endif

    fidcs=-1
    if(job.eq.1)then
        fidcs=1007
        filename=trim(basename)//'.cs.p'//trim(ctmp)
        inquire(file=trim(filename),exist=alive)
        if(.not. alive)then
            if(subdomain%pid.eq.0)&
                write(*,*)'File "',trim(filename),&
                '" does not exist, assuming Cs=',refer_Cs
            fidcs=-1
            Cs=refer_Cs
        else
            open( unit = fidcs, file = filename, &
                form = "unformatted", status = "old", &
                access = "direct", recl = rec_len )
        endif
    endif

    fidcc=-1
    if(job.ge.2)then
        fidcc=1008
        filename=trim(basename)//'.cij.p'//trim(ctmp)
        inquire(file=trim(filename),exist=alive)
        if(.not. alive)then
            if(subdomain%pid.eq.0)&
                write(*,*)'File "',trim(filename),&
                '" does not exist, assuming Cij as reference.'
            fidcc=-1
            do i=1,9
                CC(:,i)=refer_cc(i)
            enddo
            if(job.gt.2)then
            do i=10,21
                CC(:,i)=refer_cc(i)
            enddo
            endif
        else
            open( unit = fidcc, file = filename, &
                form = "unformatted", status = "old", &
                access = "direct", recl = rec_len )
        endif
    endif

    call nod_wise_mat(mesh,Nele,pNp,Cp,Cs,rho,CC,&
              fidcp,fidcs,fidcc,fidrho)

    if(verbose .and. logfid.gt.0)write(logfid,*)'read file done'

    if(job.eq.0)then
        CC(:,1)=(Cp**2)*rho
        deallocate(Cp)
        material%k_media=0
    elseif(job.eq.1)then
        CC(:,2)=(Cs**2)*rho*2d0
        CC(:,1)=(Cp**2)*rho-CC(:,2)
        material%k_media=job
        do i=1,Nele
            if(all(Cs((i-1)*pNp+1:i*pNp).le.mat_TOL)) &
                material%k_media(i)=0
        enddo
    endif

    if(job.ge.2)then
        material%k_media=job
        do i=1,Nele
            if(    all(abs(CC((i-1)*pNp+1:i*pNp,4)&
                          ).le.mat_TOL).and.&
                   all(abs(CC((i-1)*pNp+1:i*pNp,5)&
                          ).le.mat_TOL).and.&
                   all(abs(CC((i-1)*pNp+1:i*pNp,6)&
                          ).le.mat_TOL))then
                material%k_media(i)=0
            elseif(all(abs(CC((i-1)*pNp+1:i*pNp,1)&
                          -CC((i-1)*pNp+1:i*pNp,2)&
                          ).le.mat_TOL).and.&
                   all(abs(CC((i-1)*pNp+1:i*pNp,2)&
                          -CC((i-1)*pNp+1:i*pNp,3)&
                          ).le.mat_TOL))then
                material%k_media(i)=1
            endif
        enddo
    endif

    if(job.ge.3)then
        do i=1,Nele
            if(material%k_media(i).ge.2 .and. & 
                all(CC((i-1)*pNp+1:i*pNp,10:21).le.mat_TOL))&
                material%k_media(i)=2
        enddo
    endif

    material%job=maxval(material%k_media)
    if(material%job.lt.min(job,3))then
        CC=>null()
        if(material%job.eq.0)then
            allocate(CC(pNp*Nele,1))
            CC(:,1)=material%C(:,1)
        elseif(material%job.eq.1)then
            allocate(CC(pNp*Nele,2))
            CC(:,1:2)=material%C(:,1:2)
        elseif(material%job.eq.2)then
            allocate(CC(pNp*Nele,9))
            CC(:,1:9)=material%C(:,1:9)
        endif
        deallocate(material%C)
        material%C=>CC
    endif

    do i=1,Nele
        do j=1,4
            k=(i-1)*Nfp*4+(j-1)*Nfp+1
            if(mesh%fbctype(1,j,i) .le. 1)then
                if(mesh%vmapP(k).eq.&
                        mesh%vmapM((j-1)*Nfp+1)+(i-1)*pNp)then
                    mesh%fbctype(1,j,i)=1
                    mesh%fbctype(2,j,i)=1
                elseif(material%k_media(i).ge.1)then ! solid 
                    k=(mesh%vmapP(k)-1)/pNp+1 !neigheleID
                    if(material%k_media(k).eq.0)then ! solid-fluid bnd
                        mesh%fbctype(2,j,i)=1
                    endif
                else ! fluid
                    k=(mesh%vmapP(k)-1)/pNp+1 !neigheleID
                    if(material%k_media(k).ge.1)then ! fluid-solid bnd
                        mesh%fbctype(2,j,i)=2
                    else ! fluid-fluid bnd
                        mesh%fbctype(2,j,i)=3
                    endif
                endif
            endif
        enddo
    enddo
    ! boundary source or ruptures are initiated in mod_source
    if(prestress)then
        call ini_tensor(pNp*Nele,material%T0,.true.)
        fidT0=1010
        filename=trim(basename)//'.T0.p'//trim(ctmp)
        inquire(file=trim(filename),exist=alive)
        if(.not. alive)then
            if(subdomain%pid.eq.0)&
                write(*,*)'File "',trim(filename),&
                '" does not exist, assuming T0 as reference.'
            fidT0=-1
            material%T0%xx=refer_T0(1)
            material%T0%yy=refer_T0(2)
            material%T0%zz=refer_T0(3)
            material%T0%yz=refer_T0(4)
            material%T0%xz=refer_T0(5)
            material%T0%xy=refer_T0(6)
        else
            open( unit = fidT0, file = filename, &
                form = "unformatted", status = "old", &
                access = "direct", recl = rec_len )
            call nod_wise_T0(mesh,Nele,pNp,material%T0,fidT0)
        endif
    endif

end subroutine init_material

subroutine nod_wise_mat(mesh,Nele,pNp,&
            Cp,Cs,rho,CC,&
            fidcp,fidcs,fidcc,fidrho)
    integer,intent(in) :: Nele
    integer,intent(in)      :: pNp 
    integer,intent(in)      :: fidcp,fidcs,fidcc,fidrho
    real(kind=rkind)    :: Cp(:),Cs(:),rho(:),CC(:,:)
    type(tetmesh_geometry)  :: mesh
    real            :: rtmp
    integer,pointer :: globID(:)
    integer :: i,j,inod,offset

    globID=>mesh%globID
    offset=glob_Nele*pNp

    if(fidcp>0)then
        do i=1,Nele
            inod=(globID(i)-1)*pNp
            do j=1,pNp
            read(fidcp,rec=inod+j)rtmp
            Cp((i-1)*pNp+j)=rtmp
            enddo
        enddo
        close(fidcp)
    endif
    if(fidrho>0)then
        do i=1,Nele
            inod=(globID(i)-1)*pNp
            do j=1,pNp
            read(fidrho,rec=inod+j)rtmp
            rho((i-1)*pNp+j)=rtmp
            enddo
        enddo
        close(fidrho)
    endif
    if(job>=1 .and. fidcs>0)then
        do i=1,Nele
            do j=1,pNp
            inod=(globID(i)-1)*pNp+j
            read(fidcs,rec=inod)rtmp
            Cs((i-1)*pNp+j)=rtmp
            enddo
        enddo
        close(fidcs)
    endif
    if(job>=2 .and. fidcc>0)then
        do i=1,Nele
            do j=1,pNp
            inod=(globID(i)-1)*pNp+j
            read(fidcc,rec=inod)rtmp
            CC((i-1)*pNp+j,1)=rtmp ! c11
            read(fidcc,rec=inod+offset)rtmp
            CC((i-1)*pNp+j,2)=rtmp ! c22
            read(fidcc,rec=inod+offset*2)rtmp
            CC((i-1)*pNp+j,3)=rtmp ! c33
            read(fidcc,rec=inod+offset*3)rtmp
            CC((i-1)*pNp+j,4)=rtmp ! c23
            read(fidcc,rec=inod+offset*4)rtmp
            CC((i-1)*pNp+j,5)=rtmp ! c13
            read(fidcc,rec=inod+offset*5)rtmp
            CC((i-1)*pNp+j,6)=rtmp ! c12
            read(fidcc,rec=inod+offset*6)rtmp
            CC((i-1)*pNp+j,7)=rtmp ! c44
            read(fidcc,rec=inod+offset*7)rtmp
            CC((i-1)*pNp+j,8)=rtmp ! c55
            read(fidcc,rec=inod+offset*8)rtmp
            CC((i-1)*pNp+j,9)=rtmp ! c66
            enddo
        enddo
        if(job>=3)then
        do i=1,Nele
            do j=1,pNp
            inod=(globID(i)-1)*pNp+j
            read(fidcc,rec=inod+offset*9)rtmp
            CC((i-1)*pNp+j,10)=rtmp ! c14
            read(fidcc,rec=inod+offset*10)rtmp
            CC((i-1)*pNp+j,11)=rtmp ! c15
            read(fidcc,rec=inod+offset*11)rtmp
            CC((i-1)*pNp+j,12)=rtmp ! c16
            read(fidcc,rec=inod+offset*12)rtmp
            CC((i-1)*pNp+j,13)=rtmp ! c24
            read(fidcc,rec=inod+offset*13)rtmp
            CC((i-1)*pNp+j,14)=rtmp ! c25
            read(fidcc,rec=inod+offset*14)rtmp
            CC((i-1)*pNp+j,15)=rtmp ! c26
            read(fidcc,rec=inod+offset*15)rtmp
            CC((i-1)*pNp+j,16)=rtmp ! c34
            read(fidcc,rec=inod+offset*16)rtmp
            CC((i-1)*pNp+j,17)=rtmp ! c35
            read(fidcc,rec=inod+offset*17)rtmp
            CC((i-1)*pNp+j,18)=rtmp ! c36
            read(fidcc,rec=inod+offset*18)rtmp
            CC((i-1)*pNp+j,19)=rtmp ! c45
            read(fidcc,rec=inod+offset*19)rtmp
            CC((i-1)*pNp+j,20)=rtmp ! c46
            read(fidcc,rec=inod+offset*20)rtmp
            CC((i-1)*pNp+j,21)=rtmp ! c56
            enddo
        enddo
        endif
        close(fidcc)
    endif
end subroutine nod_wise_mat

subroutine nod_wise_T0(mesh,Nele,pNp,T0,fidT0)
    integer,intent(in) :: Nele,pNp,fidT0
    type(tetmesh_geometry) :: mesh
    type(tensor_array) :: T0
    real            :: rtmp
    integer,pointer :: globID(:)
    integer :: i,j,inod,offset

    globID=>mesh%globID
    offset=glob_Nele*pNp
    do i=1,Nele
        do j=1,pNp
            inod=(globID(i)-1)*pNp+j
            read(fidT0,rec=inod)rtmp
            T0%xx((i-1)*pNp+j)=rtmp
            read(fidT0,rec=inod+offset)rtmp
            T0%yy((i-1)*pNp+j)=rtmp
            read(fidT0,rec=inod+offset*2)rtmp
            T0%zz((i-1)*pNp+j)=rtmp
            read(fidT0,rec=inod+offset*3)rtmp
            T0%yz((i-1)*pNp+j)=rtmp
            read(fidT0,rec=inod+offset*4)rtmp
            T0%xz((i-1)*pNp+j)=rtmp
            read(fidT0,rec=inod+offset*5)rtmp
            T0%xy((i-1)*pNp+j)=rtmp
        enddo
    enddo
    close(fidT0)
end subroutine nod_wise_T0

function rotation_PML_bndvalue(PML_type,PML_p,Ndim,x,y,z,flag)
    real(kind=rkind) :: rotation_PML_bndvalue
    integer :: PML_type,Ndim
    real(kind=rkind) :: PML_p(3),x(Ndim),y(Ndim),z(Ndim)
    logical :: flag
    if(PML_type .eq. 0)then
        if(flag)then
        rotation_PML_bndvalue=maxval(&
            x*PML_p(1)+y*PML_p(2)+z*PML_p(3))
        else
        rotation_PML_bndvalue=minval(&
            x*PML_p(1)+y*PML_p(2)+z*PML_p(3))
        endif
    elseif(PML_type .eq. 1)then
        if(flag)then
        rotation_PML_bndvalue=sqrt(maxval(&
            (x-PML_p(1))**2+(y-PML_p(2))**2+(z-PML_p(3))**2))
        else
        rotation_PML_bndvalue=sqrt(minval(&
            (x-PML_p(1))**2+(y-PML_p(2))**2+(z-PML_p(3))**2))
        endif
    elseif(PML_type .eq. -1)then
        if(flag)then
        rotation_PML_bndvalue=sqrt(minval(&
            (x-PML_p(1))**2+(y-PML_p(2))**2+(z-PML_p(3))**2))
        else
        rotation_PML_bndvalue=sqrt(maxval(&
            (x-PML_p(1))**2+(y-PML_p(2))**2+(z-PML_p(3))**2))
        endif
    else
        rotation_PML_bndvalue=0d0
    endif
end function rotation_PML_bndvalue

end module Meshfile_mod

