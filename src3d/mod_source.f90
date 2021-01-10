!********************************************************************!
!*  This module declares and read in parameters variables           *!
!*                                                                  *!
!*  Author: Ruichao Ye    Email: Ruichao.ye@gmail.com               *!
!********************************************************************!

!---------------------------------------------------------------------
module source_mod
!---------------------------------------------------------------------

	USE mpi
    use string_mod,     only : string_conf
    use datatype_mod,   only : rkind,point_source,surface,&
                               tetmesh_geometry,matrices,&
                               vector_array,tensor_array,&
                               receivers,sources,pointer_vec_int,&
                               pointer_vec_real,rupture,&
                               tetmesh_Domain
    use para_mod,       only : def_pOrder,pNp,Nfp,Fnamelen,max_c,&
                               src_radius,startTime,verbose,logfid,&
                               job,convtest,Ndomain
    use Meshfile_mod,   only : xmin,xmax,ymin,ymax,zmin,zmax
    use jacobi_mod,     only : Basis3d,Basis2d
    use conv_mpi_mod,   only : conv_boundary
    
    implicit none
!    include 'mpif.h'
    public :: ricker,series,init_source,insert_source,glob_Nrecv
!---------------------------------------------------------------------
    real(kind=rkind) PI
    parameter(PI=3.1415926535897932384626)
    real(kind=rkind),parameter :: sTOL=1d-32
    real(kind=rkind)       :: coeff
    real(kind=rkind),parameter :: expansion=1.00d0
    integer :: glob_Nrecv
!---------------------------------------------------------------------
    contains
!---------------------------------------------------------------------

function ricker(amp,freq,t0,t)
    real(kind=rkind)::ricker,amp,freq,t0,t
    ricker=amp*(1d0-2d0*((t-t0)*PI*freq)**2)&
        *exp(-((t-t0)*PI*freq)**2)
end function ricker

function rupt_time(amp,t0,lag,t)
    real(kind=rkind)::rupt_time,amp,lag,t0,t,tt,tmp
    real(kind=rkind),parameter :: tol=1d-32
    tt=t-lag
    if(tt.le.0d0)then
        rupt_time=0d0
    elseif(tt.ge.t0)then
        rupt_time=amp
    else
        rupt_time=amp*exp((tt-t0)**2/(tt*(tt-2*t0)))
    endif
! take time derivative
!    tmp=0d0
!    if(tt.gt.0 .and. tt.lt.t0)then
!        tmp=-2d0*t0**2*(tt-t0)
!        tmp=tmp/(tt*(tt-2d0*t0))**2
!    endif
!    rupt_time=rupt_time*tmp
end function rupt_time

function series(amp,seq,seqpp,N,t,t0,dt,driv)
    integer :: N,driv,i,Nt
    real(kind=rkind) :: series,amp,seq(N),seqpp(N),t,t0,dt,delt
    Nt=N
    i=aint((t-t0)/dt)+1
    if(i.le.0)then
        series=0d0
        return
    endif
    delt=(t-t0)-(i-1)*dt
! The time series starting from t=0.0 at seq(3)
    if(driv.eq.0)then
        series=seq(i)+delt*((seq(i+1)-seq(i))/dt &
           -(seqpp(i+1)/6+seqpp(i)/3)*dt &
           + delt*(seqpp(i)/2 &
           + delt*((seqpp(i+1)-seqpp(i))/(6*dt))))
    elseif(driv.eq.1)then
        series=(seq(i+1)-seq(i))/dt &
           - (seqpp(i+1)/6+seqpp(i)/3)*dt &
           + delt*(seqpp(i) &
           + delt*(0.5*(seqpp(i+1)-seqpp(i))/dt))
    elseif(driv.eq.2)then
        series=seqpp(i)+delt*(seqpp(i+1)-seqpp(i))/dt
    else
        series=0d0
    endif
    series=series*amp
    return
end function series

subroutine init_source(sourcefname,src_num,srcs,Ndim,coord,radius,&
        mesh,t_direct,max_c,stime)
! input:
    integer,intent(in) :: Ndim
    character(len=*),intent(in) :: sourcefname
    ! file containing source information
    type(vector_array) :: coord
    real(kind=rkind),intent(in) :: radius,stime,max_c
    type(tetmesh_geometry)  :: mesh
! output:
    integer :: src_num
    ! number of source
    type(point_source),pointer :: srcs(:)
    ! source list
    real(kind=rkind),intent(inout) :: t_direct(1)
! auxilary:
    integer :: glob_src_num
    ! global number of body source
    character(len=10) :: src_array
    ! single / align
    character(len=10) :: src_time
    ! ricker / signal / initcond
    character(len=Fnamelen) :: signfname
    ! source-time sequence file name
    ! file is written is in the way of:
    !   dim(Nseq,NsignComp,glob_src_num)
    integer :: Nseq,NsignComp
    ! source-time sequence sampling number and components number
    integer :: src_numx,src_numy,src_numz
    ! number of sources in a line
    real(kind=rkind) :: srcx,srcy,srcz
    real(kind=rkind) :: srcx0,srcy0,srcz0,dsrcx,dsrcy,dsrcz 
    integer :: Nsrcx,Nsrcy,Nsrcz,Nnod,Nsign_src
    type(pointer_vec_int),allocatable :: src2nod(:) 
    type(pointer_vec_real),allocatable :: sweight(:)
    integer,allocatable :: globloc(:)
    character(len=Fnamelen) :: fname
    integer :: fid,fid1,i,j,k,m,mm,rec_len,error
    real(kind=rkind) :: recdt, rect
    logical :: alive
    inquire(iolength=rec_len)rect

    fid=2001
    src_num=0

    inquire(file=trim(sourcefname),exist=alive)
    if(.not. alive)then
        allocate(srcs(src_num))
        return
    endif

    open(fid,file=trim(sourcefname),status='old')
    call string_conf(fid,1,'src_array',2,src_array)
    call string_conf(fid,1,'src_num',2,glob_src_num)
    allocate(src2nod(glob_src_num))
    allocate(sweight(glob_src_num))
    allocate(globloc(glob_src_num))
    globloc=-1;k=0
    if(trim(src_array).eq.'single')then
        do i=1,glob_src_num
            call string_conf(fid,1,"src_x",i+1,srcx)
            call string_conf(fid,1,'src_y',i+1,srcy)
            call string_conf(fid,1,'src_z',i+1,srcz)
            call insert_source(srcx,srcy,srcz,&
                Ndim,coord,radius,Nnod,src2nod(i)%p,sweight(i)%p)
            if(Nnod.gt.0)then
                k=k+1;globloc(i)=k
            endif
            call first_arrival_time(pNp,stime,max_c,radius,&
                mesh,srcx,srcy,srcz,t_direct)
        enddo
        src_num=k
        allocate(srcs(src_num))
        if(src_num.eq.0)goto 19
        do i=1,glob_src_num
            k=globloc(i)
            if(k.le.0)cycle
            call string_conf(fid,1,'src_x'    ,k+1,srcs(k)%x)
            call string_conf(fid,1,'src_y'    ,k+1,srcs(k)%y)
            call string_conf(fid,1,'src_z'    ,k+1,srcs(k)%z)
            call string_conf(fid,1,'src_v1'   ,k+1,srcs(k)%q(1))
            call string_conf(fid,1,'src_v2'   ,k+1,srcs(k)%q(2))
            call string_conf(fid,1,'src_v3'   ,k+1,srcs(k)%q(3))
            call string_conf(fid,1,'src_e11'  ,k+1,srcs(k)%q(4))
            call string_conf(fid,1,'src_e22'  ,k+1,srcs(k)%q(5))
            call string_conf(fid,1,'src_e33'  ,k+1,srcs(k)%q(6))
            call string_conf(fid,1,'src_e23'  ,k+1,srcs(k)%q(7))
            call string_conf(fid,1,'src_e13'  ,k+1,srcs(k)%q(8))
            call string_conf(fid,1,'src_e12'  ,k+1,srcs(k)%q(9))
            call string_conf(fid,1,'amplitude',k+1,srcs(k)%amp)
            call string_conf(fid,1,'lag'      ,k+1,srcs(k)%lag)
            call string_conf(fid,1,'begin_t'  ,k+1,srcs(k)%stime)
            call string_conf(fid,1,'end_t'    ,k+1,srcs(k)%ftime)
            call string_conf(fid,1,'src_time' ,k+1,src_time)
            if(trim(src_time).eq.'ricker')then
                call string_conf(fid,1,'frequency',k+1,srcs(k)%freq)
                srcs(k)%st_type=0
            elseif(trim(src_time).eq.'series')then
                srcs(k)%st_type=1
            elseif(trim(src_time).eq.'initcond')then
                srcs(k)%st_type=2
            elseif(trim(src_time).eq.'rupture')then
                call string_conf(fid,1,'frequency',k+1,srcs(k)%freq)
                srcs(k)%st_type=3
            endif
            srcs(k)%globID=i
            srcs(k)%Nnod=size(src2nod(i)%p)
            srcs(k)%src2nod=>src2nod(i)%p
            srcs(k)%sweight=>sweight(i)%p
            if(minval(srcs(k)%src2nod).le.0)&
                print*,'src2nod error:',minval(srcs(k)%src2nod)
        enddo
    elseif(trim(src_array).eq.'align')then
        call string_conf(fid,1,'src_x',2,srcx0)
        call string_conf(fid,1,'src_x',3,dsrcx)
        call string_conf(fid,1,'src_x',4,Nsrcx)
        call string_conf(fid,1,'src_y',2,srcy0)
        call string_conf(fid,1,'src_y',3,dsrcy)
        call string_conf(fid,1,'src_y',4,Nsrcy)
        call string_conf(fid,1,'src_z',2,srcz0)
        call string_conf(fid,1,'src_z',3,dsrcz)
        call string_conf(fid,1,'src_z',4,Nsrcz)
        m=0;mm=0
        srcz=srcz0
        do k=1,Nsrcz
            srcy=srcy0
            do j=1,Nsrcy
                srcx=srcx0
                do i=1,Nsrcx
                    m=m+1
                    call insert_source(srcx,srcy,srcz,Ndim,coord,&
                        radius,Nnod,src2nod(m)%p,sweight(m)%p)
                    if(Nnod.gt.0)then
                        mm=mm+1;globloc(m)=mm
                    endif
                    call first_arrival_time(pNp,stime,max_c,radius,&
                        mesh,srcx,srcy,srcz,t_direct)
                    srcx=srcx+dsrcx
                enddo
                srcy=srcy+dsrcy
            enddo
            srcz=srcz+dsrcz
        enddo
        src_num=mm
        allocate(srcs(src_num))
        if(src_num.eq.0)goto 19
        m=0;mm=0
        srcz=srcz0
        do k=1,Nsrcz
            srcy=srcy0
            do j=1,Nsrcy
                srcx=srcx0
                do i=1,Nsrcx
                    m=m+1
                    if(globloc(m).gt.0)then
                        mm=globloc(m)
                        srcs(mm)%x=srcx
                        srcs(mm)%y=srcy
                        srcs(mm)%z=srcz
                        srcs(mm)%globID=m
                        srcs(mm)%Nnod=size(src2nod(m)%p)
                        srcs(mm)%src2nod=>src2nod(m)%p
                        srcs(mm)%sweight=>sweight(m)%p
                    endif
                    srcx=srcx+dsrcx
                enddo
                srcy=srcy+dsrcy
            enddo
            srcz=srcz+dsrcz
        enddo
        do i=1,glob_src_num
            k=globloc(i)
            if(k.le.0)cycle
            call string_conf(fid,1,'src_x'    ,2,srcs(k)%x)
            call string_conf(fid,1,'src_y'    ,2,srcs(k)%y)
            call string_conf(fid,1,'src_z'    ,2,srcs(k)%z)
            call string_conf(fid,1,'src_v1'   ,2,srcs(k)%q(1))
            call string_conf(fid,1,'src_v2'   ,2,srcs(k)%q(2))
            call string_conf(fid,1,'src_v3'   ,2,srcs(k)%q(3))
            call string_conf(fid,1,'src_e11'  ,2,srcs(k)%q(4))
            call string_conf(fid,1,'src_e22'  ,2,srcs(k)%q(5))
            call string_conf(fid,1,'src_e33'  ,2,srcs(k)%q(6))
            call string_conf(fid,1,'src_e23'  ,2,srcs(k)%q(7))
            call string_conf(fid,1,'src_e13'  ,2,srcs(k)%q(8))
            call string_conf(fid,1,'src_e12'  ,2,srcs(k)%q(9))
            call string_conf(fid,1,'amplitude',2,srcs(k)%amp)
            call string_conf(fid,1,'lag'      ,2,srcs(k)%lag)
            call string_conf(fid,1,'begin_t'  ,2,srcs(k)%stime)
            call string_conf(fid,1,'end_t'    ,2,srcs(k)%ftime)
            call string_conf(fid,1,'src_time' ,2,src_time)
            if(trim(src_time).eq.'ricker')then
                call string_conf(fid,1,'frequency',2,srcs(k)%freq)
                srcs(k)%st_type=0
            elseif(trim(src_time).eq.'series')then
                srcs(k)%st_type=1
            elseif(trim(src_time).eq.'initcond')then
                srcs(k)%st_type=2
            elseif(trim(src_time).eq.'rupture')then
                call string_conf(fid,1,'frequency',2,srcs(k)%freq)
                srcs(k)%st_type=3
            endif
            srcs(k)%globID=i
            srcs(k)%Nnod=size(src2nod(i)%p)
            srcs(k)%src2nod=>src2nod(i)%p
            srcs(k)%sweight=>sweight(i)%p
        enddo
    endif

    if(trim(src_time).eq.'series')then
        call string_conf(fid,1,'signal',2,signfname)
        fid1=2002
        fname=trim(signfname)//'_Para.txt'
        open(fid1,file=trim(fname),status='old')
        read(fid1,*)Nsign_src
        if(glob_src_num.ne.Nsign_src)&
            write(*,*)'Number of signal source error'
        read(fid1,*)NsignComp
        read(fid1,*)Nseq
        read(fid1,*)recdt
        rect=recdt*(Nseq-1)
        close(fid1)
        fname=trim(signfname)//'.dat'
        open(fid1,file=trim(fname),status='old',access='direct',&
                form='unformatted',recl=rec_len,action='read')
        do i=1,src_num
            allocate(srcs(i)%signal(Nseq,NsignComp))
            allocate(srcs(i)%ddsignal(Nseq,NsignComp))
            mm=(srcs(i)%globID-1)*NsignComp*Nseq
            do j=1,NsignComp*Nseq
                mm=mm+1
                read(fid1,rec=mm)srcs(i)%signal(j,1)
            enddo
            srcs(i)%recdt=recdt
            srcs(i)%rect=rect
            srcs(i)%Nseq=Nseq
            srcs(i)%NsignComp=NsignComp
            call second_driv(srcs(i)%signal,&
                srcs(i)%ddsignal,&
                Nseq,NsignComp,recdt)
        enddo
        close(fid1)
    endif

19  deallocate(globloc,src2nod,sweight)

end subroutine init_source

subroutine insert_source(srcx,srcy,srcz,Ndim,coord,radius,&
        Nnod,src2nod,sweight)
! input:
    real(kind=rkind),intent(in) :: srcx,srcy,srcz
    ! sources parameters
    real(kind=rkind),intent(in) :: radius
    ! Gaussian radius of source
    type(vector_array) :: coord
! output:
    integer,intent(out) :: Nnod
    integer,pointer,intent(out) :: src2nod(:) 
    real(kind=rkind),pointer,intent(out) :: sweight(:)
! auxilary
    integer :: Ndim,i,j,k
    real(kind=rkind) :: radius2,dist2,coeff

    Nnod=0
    do i=1,Ndim
        dist2=((coord%x(i)-srcx)/src_radius)**2 &
             +((coord%y(i)-srcy)/src_radius)**2 &
             +((coord%z(i)-srcz)/src_radius)**2
        if(dist2 .le. 1d0)Nnod=Nnod+1
    enddo
    if(Nnod.gt.0)then
        allocate(src2nod(Nnod))
        allocate(sweight(Nnod));sweight=0d0
        k=0
        do i=1,Ndim
            dist2=((coord%x(i)-srcx)/src_radius)**2 &
                 +((coord%y(i)-srcy)/src_radius)**2 &
                 +((coord%z(i)-srcz)/src_radius)**2
            if(dist2 .le. 1d0)then
                k=k+1
                src2nod(k)=i
                sweight(k)=exp(dist2/(dist2-1d0-1d-10))
            endif
        enddo
    endif
end subroutine insert_source

subroutine first_arrival_time(pNp,stime,cmax,srad,mesh,sx,sy,sz,&
        t_direct)
! input:
    integer,intent(in) :: pNp
    real(kind=rkind),intent(in) :: stime,cmax,srad,sx,sy,sz
    type(tetmesh_geometry)  :: mesh
! inout:
    real(kind=rkind),intent(inout) :: t_direct(1)
! auxilary::
    integer :: i,j,head,tail
    real(kind=rkind) :: dist2(pNp),dirT,tmp

    do i=1,mesh%Nhele
        head=(i-1)*pNp+1;tail=i*pNp
        dist2=(mesh%coord%x(head:tail)-sx)**2 &
             +(mesh%coord%y(head:tail)-sy)**2 &
             +(mesh%coord%z(head:tail)-sz)**2 
        dirT=(sqrt(minval(dist2))-srad)/max_c+stime
        if(t_direct(i).gt.dirT)t_direct(i)=dirT
    enddo
end subroutine first_arrival_time

subroutine second_driv(y,ypp,n,n1,dt)
    
!suppose the signal is sampled with constant interval

      implicit none
      integer n,n1,i,j
      real(kind = rkind)a(3,n),y(n,n1),ypp(n,n1)
      real(kind = rkind)xmult,dt

      do j=1,n1

        ypp(1,j) = (y(2,j)-y(1,j))/dt
        a(2,1)=dt/3.0D+00
        a(1,2)=dt/6.0D+00
      do i=2,n-1
        ypp(i,j)=(y(i+1,j)-y(i,j))/dt-(y(i,j)-y(i-1,j))/dt
        a(3,i-1)=dt/6.0D+00
        a(2,i)=2.0*dt/3.0D+00
        a(1,i+1)=dt/6.0D+00
      end do
        ypp(n,j)=(y(n,j)-y(n-1,j))/dt
        a(3,n-1)=dt/6.0D+00
        a(2,n)=dt/3.0D+00
      do i=2,n
        xmult=a(3,i-1)/a(2,i-1)
        a(2,i)=a(2,i)-xmult*a(1,i)
        ypp(i,j)= ypp(i,j)-xmult*ypp(i-1,j)
      end do
      ypp(n,j)=ypp(n,j)/a(2,n)
      do i=n-1,1,-1
        ypp(i,j) =(ypp(i,j)-a(1,i+1)*ypp(i+1,j))/a(2,i)
      end do

      enddo

      return
end subroutine second_driv

subroutine body_source(src_num,Ndim,srcs,V,S,time)
! input: 
    integer :: src_num
    integer :: Ndim
    type(point_source) :: srcs(src_num)
    real(kind=rkind) :: time
! inout
    type(vector_array) :: V
    type(tensor_array) :: S
! auxilary:
    integer :: i,j
    real(kind=rkind) :: srct,rtmp
    integer :: nid

    do i=1,src_num
        if(srcs(i)%stime.gt.time .or. &
           srcs(i)%ftime.lt.time)cycle
        if(srcs(i)%Nnod.eq.0)cycle
        if(srcs(i)%st_type.eq.0)then
            srct=ricker(srcs(i)%amp,srcs(i)%freq,&
                        srcs(i)%lag,time)
        elseif(srcs(i)%st_type.eq.1)then
            srct=series(srcs(i)%amp     ,srcs(i)%signal,&
                        srcs(i)%ddsignal,srcs(i)%Nseq,time,&
                        srcs(i)%lag     ,srcs(i)%recdt,0)
        elseif(srcs(i)%st_type.eq.3)then
            srct=rupt_time(srcs(i)%amp,&
                           srcs(i)%freq,srcs(i)%lag,time)
        endif
        do j=1,srcs(i)%Nnod
            nid =srcs(i)%src2nod(j)
            if(nid.gt.Ndim)cycle
            rtmp=srcs(i)%sweight(j)*srct
            V%x (nid)=V%x (nid)+rtmp*srcs(i)%q(1)
            V%y (nid)=V%y (nid)+rtmp*srcs(i)%q(2)
            V%z (nid)=V%z (nid)+rtmp*srcs(i)%q(3)
            S%xx(nid)=S%xx(nid)+rtmp*srcs(i)%q(4)
            S%yy(nid)=S%yy(nid)+rtmp*srcs(i)%q(5)
            S%zz(nid)=S%zz(nid)+rtmp*srcs(i)%q(6)
            S%yz(nid)=S%yz(nid)+rtmp*srcs(i)%q(7)
            S%xz(nid)=S%xz(nid)+rtmp*srcs(i)%q(8)
            S%xy(nid)=S%xy(nid)+rtmp*srcs(i)%q(9)
        enddo
    enddo
    
end subroutine body_source

subroutine surf_source(src_num,srcs,surf,time)
! input: 
    integer,intent(in) :: src_num
    type(point_source) :: srcs(src_num)
    real(kind=rkind) :: time
! inout
    type(surface) :: surf
! auxilary:
    integer :: i,j,k
    real(kind=rkind) :: srct,rtmp
    integer :: nid

if(convtest)then
    do i=1,surf%Nface
        j=(i-1)*Nfp+1
        call conv_boundary(Nfp,time,&
            surf%coord%x(j:),surf%coord%y(j:),surf%coord%z(j:),&
            surf%nx(i),surf%ny(i),surf%nz(i),&
            surf%Sn%x(j:),surf%Sn%y(j:),surf%Sn%z(j:))
    enddo
else
    do i=1,src_num
! at this moment only consider normal direction
! will update in next version
        if(srcs(i)%stime.gt.time .or. &
           srcs(i)%ftime.lt.time)cycle
        if(srcs(i)%Nnod.eq.0)cycle
        if(srcs(i)%st_type.eq.0)then
            srct=ricker(srcs(i)%amp,srcs(i)%freq,&
                        srcs(i)%lag,time)
        elseif(srcs(i)%st_type.eq.3)then
            srct=rupt_time(srcs(i)%amp,&
                           srcs(i)%freq,srcs(i)%lag,time)
        elseif(srcs(i)%st_type.eq.1)then
            srct=series(srcs(i)%amp     ,srcs(i)%signal,&
                        srcs(i)%ddsignal,srcs(i)%Nseq,time,&
                        srcs(i)%lag     ,srcs(i)%recdt,0)
        endif
        do j=1,srcs(i)%Nnod
            rtmp=srcs(i)%sweight(j)*srct
            nid =srcs(i)%src2nod(j)
            k=(nid-1)/Nfp+1
            surf%Sn%x(nid)=rtmp*surf%nx(k)
            surf%Sn%y(nid)=rtmp*surf%ny(k)
            surf%Sn%z(nid)=rtmp*surf%nz(k)
        enddo
    enddo
endif
end subroutine surf_source

subroutine point_initcond(src_num,Nhele,srcs,V,E)
! input: 
    integer :: src_num
    integer :: Nhele
    type(point_source) :: srcs(src_num)
! inout
    type(vector_array) :: V
    type(tensor_array) :: E
    integer :: i,j
    integer :: nid
    if(src_num.le.0)return
    do i=1,src_num
        if(srcs(i)%Nnod.eq.0)cycle
        if(srcs(i)%st_type.eq.2)then
        do j=1,srcs(i)%Nnod
            nid=srcs(i)%src2nod(j)
            V%x(nid) =V%x(nid) +srcs(i)%q(1)
            V%y(nid) =V%y(nid) +srcs(i)%q(2)
            V%z(nid) =V%z(nid) +srcs(i)%q(3)
            E%xx(nid)=E%xx(nid)+srcs(i)%q(4)
            E%yy(nid)=E%yy(nid)+srcs(i)%q(5)
            E%zz(nid)=E%zz(nid)+srcs(i)%q(6)
            E%yz(nid)=E%yz(nid)+srcs(i)%q(7)
            E%xz(nid)=E%xz(nid)+srcs(i)%q(8)
            E%xy(nid)=E%xy(nid)+srcs(i)%q(9)
        enddo
        endif
    enddo
end subroutine point_initcond

subroutine init_receiver(recvfname,&
        matrix,stime,etime,snapintv,recdt,&
        coord,Nele,Nhele,coord3,Nele3,Nhele3,&
        recvs,Ncomp,k_dim,subdomain)
    character(len=*) :: recvfname
    type(vector_array) :: coord,coord3
    integer,intent(in) :: Nele,Nhele,Ncomp
    integer,intent(in) :: Nele3,Nhele3
    type(matrices) :: matrix
    real(kind=rkind) :: stime,etime,snapintv,recdt
    type(receivers) :: recvs
    integer,intent(in) :: k_dim
    type(tetmesh_Domain)   :: subdomain
    ! auxilary
    integer :: Nrecv,Nrecord,i,j,k,l,sk
    integer,allocatable :: recv2ele(:)
    real(kind=rkind),allocatable :: rnodalw(:,:)
    real(kind=rkind) :: rx,ry,rz,r(1),s(1),t(1),rr,ss,tt
    integer :: error
    logical :: alive
    real(kind=rkind),parameter :: tol=1d-2
    real(kind=rkind),pointer :: invV(:,:)
    integer,allocatable :: recv_domain(:),itmp(:)

    recvs%Ncomp=Ncomp
    recvs%Nrecv=0;Nrecv=0
    recvs%Nt=aint((etime-stime)/recdt)
    recvs%Nrecord=0
    Nrecord=ceiling(snapintv/recdt)
    inquire(file=trim(recvfname),exist=alive)
    if(.not. alive)then
        return
    endif
    open(1233,file=trim(recvfname),status='old')
    error=0
    do while(error .eq. 0)
        read(1233,*,iostat=error)rx,ry,rz
        if(error .eq. 0)Nrecv=Nrecv+1
    enddo
    close(1233)
    if(Nrecv .eq. 0)then
        return
    endif
    glob_Nrecv=Nrecv
    if(Nele.le.0 .or. Nhele.le.0)goto 119
    allocate(recv2ele(Nrecv));recv2ele=0
    if(k_dim.eq.3)then
        allocate(rnodalw(pNp,Nrecv))
    elseif(k_dim.eq.2)then
        allocate(rnodalw(Nfp,Nrecv))
    endif
    rnodalw=0d0


    open(1233,file=trim(recvfname),status='old',position='REWIND')
    do i=1,Nrecv
        read(1233,*)rx,ry,rz
        if(k_dim.eq.3)then
            call inside_tet(coord,Nele,Nhele,rx,ry,rz,&
                tol,recv2ele(i),r(1),s(1),t(1))
        elseif(k_dim.eq.2)then
            call inside_tri3D(coord,Nele,Nhele,rx,ry,rz,&
                tol,recv2ele(i),r(1),s(1))
            if(recv2ele(i).gt.0)then
                call inside_tet(coord3,Nele3,Nhele3,rx,ry,rz,&
                    tol,j,rr,ss,tt)
                if(j.le.0)recv2ele(i)=-1
            endif
        endif
        if(recv2ele(i).le.0)cycle
        if(verbose .and. logfid.gt.0)then
            write(logfid,*)&
                'receiver No.',i,&
                ' inserted in this domain.'
        endif
        recvs%Nrecv=recvs%Nrecv+1
        if(k_dim.eq.3)then
            sk=1
            do j=0,def_pOrder
                do k=0,def_pOrder-j
                    do l=0,def_pOrder-j-k
                        call Basis3D(rnodalw(sk:sk,i),r,s,t,1,j,k,l)
                        sk=sk+1
                    enddo
                enddo
            enddo
        elseif(k_dim.eq.2)then
            sk=1
            do j=0,def_pOrder
                do k=0,def_pOrder-j
                call Basis2D(rnodalw(sk:sk,i),r,s,1,j,k)
                        sk=sk+1
                enddo
            enddo
        endif
    enddo
    close(1233)
    if(recvs%Nrecv.le.0)goto 119
    allocate(recvs%recv2ele(recvs%Nrecv))
    allocate(recvs%globID(recvs%Nrecv))
    if(k_dim.eq.3)then
        allocate(recvs%rnodalw(pNp,recvs%Nrecv))
        invV=>matrix%iV3D
    elseif(k_dim.eq.2)then
        allocate(recvs%rnodalw(Nfp,recvs%Nrecv))
        invV=>matrix%iV2D
    endif
    allocate(recvs%rec_buffer(Ncomp,Nrecord+2,recvs%Nrecv))
    sk=0
    do i=1,Nrecv
        if(recv2ele(i).gt.0)then
            sk=sk+1
            recvs%globID(sk)=i
            recvs%recv2ele(sk)=recv2ele(i)
            recvs%rnodalw(:,sk)=matmul(rnodalw(:,i),invV)
        endif
    enddo

    deallocate(recv2ele);deallocate(rnodalw)

119 continue
    if(recvs%Nrecv.gt.0)then
        allocate(recvs%switchoff(recvs%Nrecv))
        recvs%switchoff=.false.
    endif
    allocate(recv_domain(glob_Nrecv))
    allocate(itmp(glob_Nrecv));itmp=Ndomain+1
    do i=1,recvs%Nrecv
        itmp(recvs%globID(i))=subdomain%id
    enddo
!    print*,subdomain%id,'itmp',itmp
    call MPI_ALLREDUCE(itmp, recv_domain, glob_Nrecv, MPI_INTEGER,&
        MPI_MIN, MPI_COMM_WORLD, error)
!    print*,subdomain%id,'domn',recv_domain
    do i=1,recvs%Nrecv
        if(recv_domain(recvs%globID(i)).ne.subdomain%id)then
!            print*,'ckp'
            recvs%switchoff(i)=.true.
        endif
    enddo
    deallocate(recv_domain);deallocate(itmp)
end subroutine init_receiver

subroutine recording(recvs,V,S)
    type(receivers) :: recvs
    type(vector_array) :: V
    type(tensor_array) :: S
    integer :: i,iele,head,tail

    recvs%Nrecord=recvs%Nrecord+1
    if(recvs%Nrecv.le.0)return
    if(recvs%Ncomp.eq.1)then ! hydrophone
        do i=1,recvs%Nrecv
            iele=recvs%recv2ele(i)
            head=(iele-1)*pNp+1;tail=iele*pNp
            recvs%rec_buffer(1,recvs%Nrecord,i)=&
                dot_product(recvs%rnodalw(:,i),S%xx(head:tail))
        enddo
    elseif(recvs%Ncomp.eq.3)then ! geophone
        do i=1,recvs%Nrecv
            iele=recvs%recv2ele(i)
            head=(iele-1)*pNp+1;tail=iele*pNp
            recvs%rec_buffer(1,recvs%Nrecord,i)=&
                dot_product(recvs%rnodalw(:,i),V%x(head:tail))
            recvs%rec_buffer(2,recvs%Nrecord,i)=&
                dot_product(recvs%rnodalw(:,i),V%y(head:tail))
            recvs%rec_buffer(3,recvs%Nrecord,i)=&
                dot_product(recvs%rnodalw(:,i),V%z(head:tail))
        enddo
    elseif(recvs%Ncomp.eq.4)then ! geophone and hydrophone
        do i=1,recvs%Nrecv
            iele=recvs%recv2ele(i)
            head=(iele-1)*pNp+1;tail=iele*pNp
            recvs%rec_buffer(1,recvs%Nrecord,i)=&
                dot_product(recvs%rnodalw(:,i),V%x(head:tail))
            recvs%rec_buffer(2,recvs%Nrecord,i)=&
                dot_product(recvs%rnodalw(:,i),V%y(head:tail))
            recvs%rec_buffer(3,recvs%Nrecord,i)=&
                dot_product(recvs%rnodalw(:,i),V%z(head:tail))
            recvs%rec_buffer(4,recvs%Nrecord,i)=&
                dot_product(recvs%rnodalw(:,i),S%xx(head:tail))
        enddo
    elseif(recvs%Ncomp.eq.9)then ! full components
        do i=1,recvs%Nrecv
            iele=recvs%recv2ele(i)
            head=(iele-1)*pNp+1;tail=iele*pNp
            recvs%rec_buffer(1,recvs%Nrecord,i)=&
                dot_product(recvs%rnodalw(:,i),V%x(head:tail))
            recvs%rec_buffer(2,recvs%Nrecord,i)=&
                dot_product(recvs%rnodalw(:,i),V%y(head:tail))
            recvs%rec_buffer(3,recvs%Nrecord,i)=&
                dot_product(recvs%rnodalw(:,i),V%z(head:tail))
            recvs%rec_buffer(4,recvs%Nrecord,i)=&
                dot_product(recvs%rnodalw(:,i),S%xx(head:tail))
            recvs%rec_buffer(5,recvs%Nrecord,i)=&
                dot_product(recvs%rnodalw(:,i),S%yy(head:tail))
            recvs%rec_buffer(6,recvs%Nrecord,i)=&
                dot_product(recvs%rnodalw(:,i),S%zz(head:tail))
            recvs%rec_buffer(7,recvs%Nrecord,i)=&
                dot_product(recvs%rnodalw(:,i),S%yz(head:tail))
            recvs%rec_buffer(8,recvs%Nrecord,i)=&
                dot_product(recvs%rnodalw(:,i),S%xz(head:tail))
            recvs%rec_buffer(9,recvs%Nrecord,i)=&
                dot_product(recvs%rnodalw(:,i),S%xy(head:tail))
        enddo
    endif
end subroutine recording

subroutine write_rec_to_file(recfname,recvs,irec_offset)
    character(len=*),intent(in) :: recfname
    type(receivers) :: recvs
    integer :: irec_offset
    character*10 :: ctmp
    integer :: Ncomp,irecv,ierr,i,j,k
    integer :: rec_len
    real(kind=rkind) :: rtmp=0d0
    inquire(iolength=rec_len)rtmp
    if(verbose .and. logfid.gt.0)then
        write(logfid,*)'curr_irec=',recvs%Nrecord+irec_offset
        write(logfid,*)'Nrecord=',recvs%Nrecord
        write(logfid,*)'offset=',irec_offset
    endif
    if(recvs%Nrecord.le.0)goto 99
    if(recvs%Nrecv.le.0)goto 98
    Ncomp=recvs%Ncomp
    if(irec_offset.le.0)then
        do irecv=1,recvs%Nrecv
            if(recvs%switchoff(irecv))cycle
            write(ctmp,'(I6)')recvs%globID(irecv)
            ctmp=adjustl(ctmp)
            open(3311,&
                file=trim(recfname)//'_'//trim(ctmp)//'.dat',&
                access='direct',form='unformatted',recl=rec_len,&
                status='replace')
            k=0
            do j=1,recvs%Nrecord
                do i=1,Ncomp
                    k=k+1
                    write(3311,rec=k)recvs%rec_buffer(i,j,irecv)
                enddo
            enddo
            close(3311)
        enddo
    else
        do irecv=1,recvs%Nrecv
            if(recvs%switchoff(irecv))cycle
            write(ctmp,'(I6)')recvs%globID(irecv)
            ctmp=adjustl(ctmp)
            open(3311,&
                file=trim(recfname)//'_'//trim(ctmp)//'.dat',&
                access='direct',form='unformatted',recl=rec_len,&
                status='old')
            k=irec_offset*Ncomp
            do j=1,recvs%Nrecord
                do i=1,Ncomp
                    k=k+1
                    write(3311,rec=k)recvs%rec_buffer(i,j,irecv)
                enddo
            enddo
            close(3311)
        enddo
    endif
98  irec_offset=irec_offset+recvs%Nrecord
99  recvs%Nrecord=0

end subroutine write_rec_to_file

subroutine inside_tet(coord,Nele,Nhele,vx,vy,vz,tol,inele,r,s,t)
    type(vector_array) :: coord
    integer,intent(in) :: Nele,Nhele
    real(kind=rkind),intent(in) :: vx,vy,vz,tol
    integer,intent(out) :: inele
    real(kind=rkind),intent(out) :: r,s,t
    real(kind=rkind) :: x1,x2,x3,x4,y1,y2,y3,y4,&
                        z1,z2,z3,z4,xx,yy,zz,Jtmp
    integer :: head,iele
    real(kind=rkind) :: r0,res_tmp
    real(kind=rkind) :: tol_tmp

    inele=-1
    tol_tmp=-1d1
    do iele=1,Nele
        head=(iele-1)*pNp
        x1=coord%x(head+1)
        x2=coord%x(head+def_pOrder+1)
        x3=coord%x(head+Nfp)
        x4=coord%x(head+pNp)
        y1=coord%y(head+1)
        y2=coord%y(head+def_pOrder+1)
        y3=coord%y(head+Nfp)
        y4=coord%y(head+pNp)
        z1=coord%z(head+1)
        z2=coord%z(head+def_pOrder+1)
        z3=coord%z(head+Nfp)
        z4=coord%z(head+pNp)
        Jtmp=det4(1d0,1d0,1d0,1d0,&
                   x1, x2, x3, x4,&
                   y1, y2, y3, y4,&
                   z1, z2, z3, z4)
        r0=(abs(Jtmp)**(1d0/3d0))*tol
        if(min(x1,x2,x3,x4)-r0 .gt. vx)cycle
        if(max(x1,x2,x3,x4)+r0 .lt. vx)cycle
        if(min(y1,y2,y3,y4)-r0 .gt. vy)cycle
        if(max(y1,y2,y3,y4)+r0 .lt. vy)cycle
        if(min(z1,z2,z3,z4)-r0 .gt. vz)cycle
        if(max(z1,z2,z3,z4)+r0 .lt. vz)cycle
        xx=det4(1d0,1d0,1d0,1d0,&
                 x1, vx, x3, x4,&
                 y1, vy, y3, y4,&
                 z1, vz, z3, z4)/Jtmp
        if(xx.lt.-tol .or. xx.gt.1d0+tol)cycle
        yy=det4(1d0,1d0,1d0,1d0,&
                 x1, x2, vx, x4,&
                 y1, y2, vy, y4,&
                 z1, z2, vz, z4)/Jtmp
        if(yy.lt.-tol .or. yy.gt.1d0+tol)cycle
        zz=det4(1d0,1d0,1d0,1d0,&
                 x1, x2, x3, vx,&
                 y1, y2, y3, vy,&
                 z1, z2, z3, vz)/Jtmp
        if(zz.lt.-tol .or. zz.gt.1d0+tol)cycle
        if(xx+yy+zz.gt.1d0+tol)cycle
        res_tmp=min(xx,yy,zz,1d0-xx-yy-zz,1d0-xx,1d0-yy,1d0-zz)
        if(tol_tmp.lt.res_tmp)then
            if(iele.le.Nhele)then
                inele=iele
            else
                inele=-1
            endif
            tol_tmp=res_tmp
            r=2d0*xx-1d0
            s=2d0*yy-1d0
            t=2d0*zz-1d0
        endif
    enddo

end subroutine inside_tet

subroutine inside_tri3D(coord,Nele,Nhele,vx,vy,vz,tol,inele,r,s)
    type(vector_array) :: coord
    integer,intent(in) :: Nele,Nhele
    real(kind=rkind),intent(in) :: vx,vy,vz,tol
    integer,intent(out) :: inele
    real(kind=rkind),intent(out) :: r,s
    real(kind=rkind) :: tol_tmp
    real(kind=rkind) :: x1,x2,x3,y1,y2,y3,z1,z2,z3
    real(kind=rkind) :: nx,ny,nz,sx,sy,sz,tx,ty,tz
    real(kind=rkind) :: rr,ss,tt
    integer :: i,j

    inele=-1
    tol_tmp=-1d1
    do i=1,Nele
        j=(i-1)*Nfp
        x1=coord%x(j+1)
        x2=coord%x(j+def_pOrder+1)
        x3=coord%x(j+Nfp)
        y1=coord%y(j+1)
        y2=coord%y(j+def_pOrder+1)
        y3=coord%y(j+Nfp)
        z1=coord%z(j+1)
        z2=coord%z(j+def_pOrder+1)
        z3=coord%z(j+Nfp)
        sx=x2-x1;sy=y2-y1;sz=z2-z1
        ss=sx**2+sy**2+sz**2
        tx=x3-x1;ty=y3-y1;tz=z3-z1
        tt=tx**2+ty**2+tz**2
        tt=sqrt(ss)*2d0
        rr=min(sqrt(ss),sqrt(tt))
        if(min(x1,x2,x3)-rr .gt. vx)cycle
        if(max(x1,x2,x3)+rr .lt. vx)cycle
        if(min(y1,y2,y3)-rr .gt. vy)cycle
        if(max(y1,y2,y3)+rr .lt. vy)cycle
        if(min(z1,z2,z3)-rr .gt. vz)cycle
        if(max(z1,z2,z3)+rr .lt. vz)cycle
        nx=sy*tz-sz*ty
        ny=sz*tx-sx*tz
        nz=sx*ty-sy*tx
        rr=sqrt(nx**2+ny**2+nz**2)
        nx=nx/rr;ny=ny/rr;nz=nz/rr
        x1=vx-x1;y1=vy-y1;z1=vz-z1
        if(abs(nx*x1+ny*y1+nz*z1) .gt. tt)cycle
        tt=det3(sx,sy,sz,tx,ty,tz,nx,ny,nz)
        rr =((ty*nz-ny*tz)*x1 &
            +(nx*tz-tx*nz)*y1 &
            +(tx*ny-nx*ty)*z1)/tt
        if(rr .lt. -tol .or. rr .gt. 1d0+tol)cycle
        ss =((ny*sz-sy*nz)*x1 &
            +(sx*nz-nx*sz)*y1 &
            +(nx*sy-sx*ny)*z1)/tt
        if(ss .lt. -tol .or. ss .gt. 1d0+tol)cycle
        if(rr+ss .gt. 1d0+tol)cycle
        tt=min(rr,ss,1d0-rr-ss,1d0-rr,1d0-ss)
        if(tol_tmp .lt. tt)then
            if(i.le.Nhele)then
                inele=i
            else
                inele=-1
            endif
            tol_tmp=tt
            r=2d0*rr-1d0
            s=2d0*ss-1d0
        endif
    enddo
        
end subroutine inside_tri3D

function det4(a11,a12,a13,a14,&
              a21,a22,a23,a24,&
              a31,a32,a33,a34,&
              a41,a42,a43,a44)
    real(kind=rkind) :: det4,&
        a11,a12,a13,a14,a21,a22,a23,a24,&
        a31,a32,a33,a34,a41,a42,a43,a44
    det4 = a11*det3(a22,a23,a24,&
            a32,a33,a34,&
            a42,a43,a44)&
         - a12*det3(a21,a23,a24,&
            a31,a33,a34,&
            a41,a43,a44)&
         + a13*det3(a21,a22,a24,&
            a31,a32,a34,&
            a41,a42,a44)&
         - a14*det3(a21,a22,a23,&
            a31,a32,a33,&
            a41,a42,a43)
    return
end function det4

function det3(a11,a12,a13,&
              a21,a22,a23,&
              a31,a32,a33)
    real(kind=rkind)::det3,a11,a12,a13,a21,a22,a23,a31,a32,a33
    det3 = a11*a22*a33+a12*a23*a31+a21*a32*a13 &
         - a13*a22*a31-a12*a21*a33-a11*a23*a32
    return
end function det3

subroutine invA3( a11, a12, a13,&
                  a21, a22, a23,&
                  a31, a32, a33,&
                 ia11,ia12,ia13,&
                 ia21,ia22,ia23,&
                 ia31,ia32,ia33)
    real(kind=rkind)::detA,a11,a12,a13,a21,a22,a23,a31,a32,a33
    real(kind=rkind)::ia11,ia12,ia13,ia21,ia22,ia23,ia31,ia32,ia33
    detA = det3(a11,a12,a13,a21,a22,a23,a31,a32,a33)
    ia11=(a22*a33-a23*a32)/detA
    ia21=(a23*a31-a21*a33)/detA
    ia31=(a21*a32-a22*a31)/detA
    ia12=(a13*a32-a12*a33)/detA
    ia22=(a11*a33-a13*a31)/detA
    ia32=(a12*a31-a11*a32)/detA
    ia13=(a12*a23-a13*a22)/detA
    ia23=(a13*a21-a11*a23)/detA
    ia33=(a11*a22-a12*a21)/detA
    return
end subroutine invA3

end module source_mod

