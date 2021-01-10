!*******************************************************************!
!*  This module declares and read in parameters variables          *!
!*                                                                 *!
!*  Author: Ruichao Ye    Email: Ruichao.ye@gmail.com              *!
!*******************************************************************!

!--------------------------------------------------------------------
module para_mod
!--------------------------------------------------------------------

    use string_mod, only : string_conf
    use datatype_mod, only : rkind

    implicit none

    public :: para_init
              ! subroutine to readin para files
    public :: startTime, finalTime, dt
              ! simulation time setting
    public :: basename
              ! readin mesh filename
    public :: sourcefname, src_radius
              ! embedding source 
    public :: recvfname,ruptrecvfname
              ! receiver filenames 
    public :: recordfname, recordintv, rec_Ncomp
              ! recording information 
    public :: snapfname, snapintv
              ! snapshot information 
    public :: job, task, friction, prestress
              ! job identifier
    public :: def_pOrder, pNp, Nfp, Nfp4, Nsele
              ! polynomial order
    public :: PMLHx1,PMLHy1,PMLHz1,&
              PMLHx2,PMLHy2,PMLHz2,&
              PMLx1_p,PMLy1_p,PMLz1_p,&
              PMLx2_p,PMLy2_p,PMLz2_p,&
              PMLx1_type,PMLy1_type,PMLz1_type,&
              PMLx2_type,PMLy2_type,PMLz2_type,&
              PMLAlpha,PMLM,pmlgeometry
              ! absorbing boundary 
    public :: refer_Cp, refer_Cs, refer_rho,&
              refer_cc, max_c,refer_T0
              ! reference elastic material and prestress
    public :: start_ckp, Num_ckp
              ! checkpointing and restarting
    public :: Nprocs,Ndomain,Nshot,Ddecomp
              ! domain decomposition
    public :: Fnamelen
              ! default maximal filename length
    public :: verbose,logfid,convtest,withrupt,rpml,symm,grav
              ! flags
    public :: rk4a,rk4b,rk4c
              ! Runge--kutta constants
    public :: rupt_gamma
    public :: myrank

!--------------------------------------------------------------------

! explicit Runge-Kutta [N2]
    real(kind=rkind),parameter :: &
        rk4a(5) = (/0.0D0                 , &
                   -4.178904744998519D-1  , &
                   -1.192151694642677D0   , &
                   -1.697784692471528D0   , &
                   -1.514183444257156D0/) , &
        rk4b(5) = (/1.496590219992291D-1  , &
                    3.792103129996273D-1  , &
                    8.229550293869817D-1  , &
                    6.994504559491221D-1  , &
                    1.530572479681520D-1/), &
        rk4c(5) = (/0.0D0                 , &
                    1.496590219992291D-1  , &
                    3.704009573642048D-1  , &
                    6.222557631344432D-1  , &
                    9.582821306746903D-1/)

!--------------------------------------------------------------------

    integer,parameter :: Fnamelen=100
    character(len=10) :: pmlgeometry,task,friction,Ddecomp
    character(len=Fnamelen) :: basename,sourcefname,snapfname
    character(len=Fnamelen) :: recvfname,ruptrecvfname
    character(len=Fnamelen) :: recordfname,obsdatfname,residufname
    integer :: def_pOrder, pNp, Nfp, Nfp4, Nsele, job
    integer :: start_ckp, Num_ckp
    integer :: Nprocs, Ndomain, Nshot, rec_Ncomp
    real(kind=rkind) :: startTime, finalTime
    real(kind=rkind) :: src_radius
    real(kind=rkind) :: dt,snapintv,recordintv,corrintv
    integer :: recv_intv
    real(kind=rkind) :: PMLHx1,PMLHx2,&
                        PMLHy1,PMLHy2,&
                        PMLHz1,PMLHz2
    real(kind=rkind) :: PMLx1_p(3),PMLx2_p(3),&
                        PMLy1_p(3),PMLy2_p(3),&
                        PMLz1_p(3),PMLz2_p(3)
    integer ::          PMLx1_type,PMLx2_type,&
                        PMLy1_type,PMLy2_type,&
                        PMLz1_type,PMLz2_type
    real(kind=rkind) :: PMLAlpha,PMLM(12)
    real(kind=rkind) :: refer_Cp, refer_Cs, refer_rho, max_c
    real(kind=rkind) :: refer_cc(21) 
    !ortho: c11,c22,c33,c23,c13,c12,c44,c55,c66
    !TTI  : c14,c15,c16,c24,c25,c26,c34,c35,c36,c45,c46,c56
    integer :: prestress
    real(kind=rkind) :: refer_T0(6)
    ! T11,T22,T33,T23,T13,T12
    logical :: verbose=.false.
    logical :: convtest=.false.
    logical :: withrupt=.false.
    integer :: logfid=0
    logical :: rpml=.false.
    logical :: symm=.true.
    logical :: grav=.true.
    character(len=10),parameter :: timescheme='EXRK4'
    real(kind=rkind) :: rupt_gamma=1d1
    integer :: myrank
!--------------------------------------------------------------------
    contains
!--------------------------------------------------------------------

    subroutine para_init(fnm_conf)
        character (len=*) :: fnm_conf
        integer :: fid, N, itmp
        logical :: alive
        fid=1001
        inquire(file=trim(fnm_conf),exist=alive)
        if(.not. alive)then
            write(*,*)'Error 1003: File "',&
                trim(fnm_conf),'" does not exist.'
            stop
        endif

        open(fid,file=trim(fnm_conf),status="old")
        call string_conf(fid,1,'task',2,task)
!        if(trim(task).eq.'simu')then
        call string_conf(fid,1,'friction',2,friction)
!        endif
        call string_conf(fid,1,'Ndomain',2,Ndomain)
        call string_conf(fid,1,'domaindcmp',2,Ddecomp)
        call string_conf(fid,1,'startTime',2,startTime)
        call string_conf(fid,1,'finalTime',2,finalTime)
        call string_conf(fid,1,'dt',2,dt)
        call string_conf(fid,1,'src_radius',2,src_radius)
        call string_conf(fid,1,'basename',2,basename)
        call string_conf(fid,1,'sourcefname',2,sourcefname)
        call string_conf(fid,1,'recvfname',2,recvfname)
        call string_conf(fid,1,'ruptrecvfnm',2,ruptrecvfname)
        call string_conf(fid,1,'snapfname',2,snapfname)
        call string_conf(fid,1,'recordfname',2,recordfname)
        call string_conf(fid,1,'snapintv',2,snapintv)
        call string_conf(fid,1,'recordintv',2,recordintv)
        call string_conf(fid,1,'rec_Ncomp',2,rec_Ncomp)
        call string_conf(fid,1,'def_pOrder',2,def_pOrder)
        call string_conf(fid,1,'max_c',2,max_c)
        call string_conf(fid,1,'PMLgeomty',2,PMLgeometry)
        call string_conf(fid,1,'PMLalpha',2,PMLAlpha)
        if(trim(PMLgeometry).eq.'regula')then
            call string_conf(fid,1,'PMLH',2,PMLHx1)
            call string_conf(fid,1,'PMLH',3,PMLHx2)
            call string_conf(fid,1,'PMLH',4,PMLHy1)
            call string_conf(fid,1,'PMLH',5,PMLHy2)
            call string_conf(fid,1,'PMLH',6,PMLHz1)
            call string_conf(fid,1,'PMLH',7,PMLHz2)
        elseif(trim(PMLgeometry).eq.'rotate')then
!      0: plane, with normal direction as PMLx1_p(1,2,3)
! (+/-)1: sphere, with center at PMLx1_p(1,2,3), 
!               along +/- radius direction
            rpml=.true.
            call string_conf(fid,1,'PMLx1',2, PMLx1_type)
            call string_conf(fid,1,'PMLx1',3, PMLx1_p(1))
            call string_conf(fid,1,'PMLx1',4, PMLx1_p(2))
            call string_conf(fid,1,'PMLx1',5, PMLx1_p(3))
            call string_conf(fid,1,'PMLx1',6,PMLHx1)
            call string_conf(fid,1,'PMLx2',2, PMLx2_type)
            call string_conf(fid,1,'PMLx2',3, PMLx2_p(1))
            call string_conf(fid,1,'PMLx2',4, PMLx2_p(2))
            call string_conf(fid,1,'PMLx2',5, PMLx2_p(3))
            call string_conf(fid,1,'PMLx2',6,PMLHx2)
            call string_conf(fid,1,'PMLy1',2, PMLy1_type)
            call string_conf(fid,1,'PMLy1',3, PMLy1_p(1))
            call string_conf(fid,1,'PMLy1',4, PMLy1_p(2))
            call string_conf(fid,1,'PMLy1',5, PMLy1_p(3))
            call string_conf(fid,1,'PMLy1',6,PMLHy1)
            call string_conf(fid,1,'PMLy2',2, PMLy2_type)
            call string_conf(fid,1,'PMLy2',3, PMLy2_p(1))
            call string_conf(fid,1,'PMLy2',4, PMLy2_p(2))
            call string_conf(fid,1,'PMLy2',5, PMLy2_p(3))
            call string_conf(fid,1,'PMLy2',6,PMLHy2)
            call string_conf(fid,1,'PMLz1',2, PMLz1_type)
            call string_conf(fid,1,'PMLz1',3, PMLz1_p(1))
            call string_conf(fid,1,'PMLz1',4, PMLz1_p(2))
            call string_conf(fid,1,'PMLz1',5, PMLz1_p(3))
            call string_conf(fid,1,'PMLz1',6,PMLHz1)
            call string_conf(fid,1,'PMLz2',2, PMLz2_type)
            call string_conf(fid,1,'PMLz2',3, PMLz2_p(1))
            call string_conf(fid,1,'PMLz2',4, PMLz2_p(2))
            call string_conf(fid,1,'PMLz2',5, PMLz2_p(3))
            call string_conf(fid,1,'PMLz2',6,PMLHz2)
        endif
        call string_conf(fid,1,'job',2,job)
        call string_conf(fid,1,'symm',2,itmp)
        symm=itmp.gt.0
        call string_conf(fid,1,'start_ckp',2,start_ckp)
        call string_conf(fid,1,'refer_Cp',2,refer_Cp)
        call string_conf(fid,1,'refer_rho',2,refer_rho)
        if(job.gt.0)&
            call string_conf(fid,1,'refer_Cs',2,refer_Cs)
        if(job.gt.1)then
            call string_conf(fid,1,'refer_c11',2,refer_cc(1))
            call string_conf(fid,1,'refer_c22',2,refer_cc(2))
            call string_conf(fid,1,'refer_c33',2,refer_cc(3))
            call string_conf(fid,1,'refer_c23',2,refer_cc(4))
            call string_conf(fid,1,'refer_c13',2,refer_cc(5))
            call string_conf(fid,1,'refer_c12',2,refer_cc(6))
            call string_conf(fid,1,'refer_c44',2,refer_cc(7))
            call string_conf(fid,1,'refer_c55',2,refer_cc(8))
            call string_conf(fid,1,'refer_c66',2,refer_cc(9))
        endif
        if(job.gt.2)then
            call string_conf(fid,1,'refer_c14',2,refer_cc(10))
            call string_conf(fid,1,'refer_c15',2,refer_cc(11))
            call string_conf(fid,1,'refer_c16',2,refer_cc(12))
            call string_conf(fid,1,'refer_c24',2,refer_cc(13))
            call string_conf(fid,1,'refer_c25',2,refer_cc(14))
            call string_conf(fid,1,'refer_c26',2,refer_cc(15))
            call string_conf(fid,1,'refer_c34',2,refer_cc(16))
            call string_conf(fid,1,'refer_c35',2,refer_cc(17))
            call string_conf(fid,1,'refer_c36',2,refer_cc(18))
            call string_conf(fid,1,'refer_c45',2,refer_cc(19))
            call string_conf(fid,1,'refer_c46',2,refer_cc(20))
            call string_conf(fid,1,'refer_c56',2,refer_cc(21))
        endif
        call string_conf(fid,1,'prestress',2,prestress)
        if(prestress.gt.0)then
            call string_conf(fid,1,'refer_T11',2,refer_T0(1))
            call string_conf(fid,1,'refer_T22',2,refer_T0(2))
            call string_conf(fid,1,'refer_T33',2,refer_T0(3))
            call string_conf(fid,1,'refer_T23',2,refer_T0(4))
            call string_conf(fid,1,'refer_T13',2,refer_T0(5))
            call string_conf(fid,1,'refer_T12',2,refer_T0(6))
        endif
        call string_conf(fid,1,'rupt_gamma',2,rupt_gamma)

        close(fid)

        rupt_gamma=rupt_gamma/dt

        N=def_pOrder

        pNp=(N+1)*(N+2)*(N+3)/6
        Nfp=(N+1)*(N+2)/2
        Nfp4=Nfp*4
        Nsele=N*(N+1)*(N+2)/6
        if(N >= 2)Nsele=Nsele+(N-1)*N*(N+1)*2/3
        if(N >= 3)Nsele=Nsele+(N-2)*(N-1)*N/6

    end subroutine para_init

!-------------------------------------------------------------------

end module para_mod

