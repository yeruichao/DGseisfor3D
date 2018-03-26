!*******************************************************************!
!*  This module solves nonlinear rupture coupling problem          *!
!*                                                                 *!
!*  Author: Ruichao Ye    Email: Ruichao.ye@gmail.com              *!
!*******************************************************************!

!--------------------------------------------------------------------
module ruptsolve_mod
!--------------------------------------------------------------------

    use para_mod,     only : rk4a,rk4b,rk4c,convtest,rupt_gamma,myrank
    use datatype_mod, only : rkind,matrices,sources,rupture,&
                             tetmesh_geometry,tetmesh_material,&
                             vector_array,tensor_array,&
                             tetmesh_domain,rupture_matrices,&
                             pointer_vec_real,pointer_mat_real,&
                             DiagMM,DiagMMw
    use conv_mpi_mod, only : conv_boundary
    use friction_mod, only : RSfric_F,RSfric_dFdN,RSfric_dFdS,&
                             RSfric_invF,RSfric_G,RSfric_invG,&
                             RSfric_dGdS,RSfric_dFdP,RSfric_dGdP,&
                             RSfric_invG
    use parallel_mod, only : domain_exchange

    implicit none
    include 'mpif.h'

    private :: ntH,ntG,ntS,ntY,ntZ

    private :: iUm1 ,iUm2 ,iUm3 ,iUp1 ,iUp2 ,iUp3 ,&
               iEm11,iEm12,iEm13,iEp11,iEp12,iEp13,&
               iEm21,iEm22,iEm23,iEp21,iEp22,iEp23,&
               iEm31,iEm32,iEm33,iEp31,iEp32,iEp33,&
               iVm1 ,iVm2 ,iVm3 ,iVp1 ,iVp2 ,iVp3 ,&
               iTf1 ,iTf2 ,iTf3 ,iVt1 ,iVt2 ,iVt3 ,&
               iSig

    private :: pUm,pUp,pEm,pEp,pVm,pVp,pTf,pTf1,pTf2,pTf3,pSig,&
               hUm,hUp,hEm,hEp,hVm,hVp,hTf,hTf1,hTf2,hTf3,hSig

    private :: C11,C12,C22,C13,C23,C33,C14,C24,C34,C44,&
               C15,C25,C35,C45,C55,C16,C26,C36,C46,C56,C66,&
               P11,P12,P22,P13,P23,P33,M11,M12,M22,M13,M23,M33,&
               P11P22,P11P33,P22P33,M11M22,M11M33,M22M33,&
               M11P22,M11P33,M22P33,P11M22,P11M33,P22M33,Zero,&
               CC,DD

    public :: rupt_alpha_V, rupt_alpha_E, rupt_on

    real(kind=rkind),allocatable,target :: &
                ntH(:,:),ntG(:),ntS(:),ntS0(:),ntY(:),ntZ(:)

    integer :: iUm1 ,iUm2 ,iUm3 ,iUp1 ,iUp2 ,iUp3 ,&
               iEm11,iEm12,iEm13,iEp11,iEp12,iEp13,&
               iEm21,iEm22,iEm23,iEp21,iEp22,iEp23,&
               iEm31,iEm32,iEm33,iEp31,iEp32,iEp33,&
               iVm1 ,iVm2 ,iVm3 ,iVp1 ,iVp2 ,iVp3 ,&
               iTf1 ,iTf2 ,iTf3 ,iVt1 ,iVt2 ,iVt3 ,&
               iSig

    real(kind=rkind),pointer :: &
                pUm(:),pUp(:),pEm(:),pEp(:),pVm(:),pVp(:),&
                pSig(:),pTf(:),pTf1(:),pTf2(:),pTf3(:),&
                hUm(:),hUp(:),hEm(:),hEp(:),hVm(:),hVp(:),&
                hSig(:),hTf(:),hTf1(:),hTf2(:),hTf3(:),&
                hSig0(:),hT0(:),hT01(:),hT02(:),hT03(:),&
                ntY1(:),ntZ1(:),ntYb(:),ntZb(:),ntZ2(:)

    real(kind=rkind),pointer :: &
               C11(:),&
               C12(:),C22(:),&
               C13(:),C23(:),C33(:),&
               C14(:),C24(:),C34(:),C44(:),&
               C15(:),C25(:),C35(:),C45(:),C55(:),&
               C16(:),C26(:),C36(:),C46(:),C56(:),C66(:),&
               P11(:),&
               P12(:),P22(:),&
               P13(:),P23(:),P33(:),&
               M11(:),&
               M12(:),M22(:),&
               M13(:),M23(:),M33(:),&
               P11P22(:),P11P33(:),P22P33(:),&
               M11M22(:),M11M33(:),M22M33(:),&
               M11P22(:),M11P33(:),M22P33(:),&
               P11M22(:),P11M33(:),P22M33(:),&
               ZERO(:)

    type(pointer_vec_real) :: CC(9,9),DD(9,9),GG(9),&
                pG(3),pS(3),pT(3),pS0(3)
    type(pointer_mat_real) :: pH(3,3),pA(4,3)
!    type(pointer_mat_real),pointer :: pA(:,:)

    real(kind=rkind),parameter :: &
        rupt_alpha_V = 1.0d0, rupt_alpha_E=1.0d0
!        rupt_alpha_V = 1.0d1, rupt_alpha_E=1.0d1

! Remark: the case with rupt_alpha_E>0 needs debugging, 
!         if penalty for E is required

    logical :: rupt_on=.true.
    logical :: debugging=.false.

!--------------------------------------------------------------------
    contains
!--------------------------------------------------------------------

subroutine rupt_Newton_allocate(pNp,Nfp,rupt)
    integer,intent(in) :: pNp,Nfp
    type(rupture) :: rupt
    integer :: i,j

    allocate(ntH(3*Nfp,3*Nfp))
    allocate(ntG(3*Nfp))
    allocate(ntS(3*Nfp))
    allocate(ntS0(3*Nfp))
    allocate(ntY1(30*pNp))
    allocate(ntYb(4*Nfp))
    allocate(ntZ1(30*pNp))
    allocate(ntZ2(4*Nfp))
    allocate(ntZb(4*Nfp))

!if(debugging)then
!    allocate(ntY(7*Nfp))
!endif

    i=0

    iUm1 =i; i=i+pNp; iUm2 =i; i=i+pNp; iUm3 =i; i=i+pNp
    iUp1 =i; i=i+pNp; iUp2 =i; i=i+pNp; iUp3 =i; i=i+pNp

    iEm11=i; i=i+pNp; iEm21=i; i=i+pNp; iEm31=i; i=i+pNp
    iEm12=i; i=i+pNp; iEm22=i; i=i+pNp; iEm32=i; i=i+pNp
    iEm13=i; i=i+pNp; iEm23=i; i=i+pNp; iEm33=i; i=i+pNp
    iEp11=i; i=i+pNp; iEp21=i; i=i+pNp; iEp31=i; i=i+pNp
    iEp12=i; i=i+pNp; iEp22=i; i=i+pNp; iEp32=i; i=i+pNp
    iEp13=i; i=i+pNp; iEp23=i; i=i+pNp; iEp33=i; i=i+pNp

    iVm1 =i; i=i+pNp; iVm2 =i; i=i+pNp; iVm3 =i; i=i+pNp
    iVp1 =i; i=i+pNp; iVp2 =i; i=i+pNp; iVp3 =i; i=i+pNp

    i=0

    iTf1 =i; i=i+Nfp; iTf2 =i; i=i+Nfp; iTf3 =i; i=i+Nfp
    iSig =i; i=i+Nfp
    iVt1 =i; i=i+Nfp; iVt2 =i; i=i+Nfp; iVt3 =i; i=i+Nfp

    if(rupt%Nhface.gt.0)then
        allocate(rupt%M(rupt%Nhface))
        do i=1,rupt%Nhface
            call rupt_matrix_allocate(pNp,Nfp,rupt%M(i))
        enddo
    endif

    pUm  => ntY1(iUm1 +1:iUm3 +pNp)
    pUp  => ntY1(iUp1 +1:iUp3 +pNp)
    pEm  => ntY1(iEm11+1:iEm33+pNp)
    pEp  => ntY1(iEp11+1:iEp33+pNp)
    pVm  => ntY1(iVm1 +1:iVm3 +pNp)
    pVp  => ntY1(iVp1 +1:iVp3 +pNp)
    pTf  => ntYb(iTf1 +1:iTf3 +Nfp)
    pTf1 => ntYb(iTf1 +1:iTf1 +Nfp)
    pTf2 => ntYb(iTf2 +1:iTf2 +Nfp)
    pTf3 => ntYb(iTf3 +1:iTf3 +Nfp)
    pSig => ntYb(iSig +1:iSig +Nfp)

    hUm  => ntZ1(iUm1 +1:iUm3 +pNp)
    hUp  => ntZ1(iUp1 +1:iUp3 +pNp)
    hEm  => ntZ1(iEm11+1:iEm33+pNp)
    hEp  => ntZ1(iEp11+1:iEp33+pNp)
    hVm  => ntZ1(iVm1 +1:iVm3 +pNp)
    hVp  => ntZ1(iVp1 +1:iVp3 +pNp)
    hTf  => ntZb(iTf1 +1:iTf3 +Nfp)
    hTf1 => ntZb(iTf1 +1:iTf1 +Nfp)
    hTf2 => ntZb(iTf2 +1:iTf2 +Nfp)
    hTf3 => ntZb(iTf3 +1:iTf3 +Nfp)
    hSig => ntZb(iSig +1:iSig +Nfp)

    hT0  => ntZ2(iTf1 +1:iTf3 +Nfp)
    hT01 => ntZ2(iTf1 +1:iTf1 +Nfp)
    hT02 => ntZ2(iTf2 +1:iTf2 +Nfp)
    hT03 => ntZ2(iTf3 +1:iTf3 +Nfp)
    hSig0=> ntZ2(iSig +1:iSig +Nfp)

    C11=>ntY1(       1:   pNp);C12=>ntY1(   pNp+1: 2*pNp)
    C22=>ntY1( 2*pNp+1: 3*pNp);C13=>ntY1( 3*pNp+1: 4*pNp)
    C23=>ntY1( 4*pNp+1: 5*pNp);C33=>ntY1( 5*pNp+1: 6*pNp)
    C14=>ntY1( 6*pNp+1: 7*pNp);C24=>ntY1( 7*pNp+1: 8*pNp)
    C34=>ntY1( 8*pNp+1: 9*pNp);C44=>ntY1( 9*pNp+1:10*pNp)
    C15=>ntY1(10*pNp+1:11*pNp);C25=>ntY1(11*pNp+1:12*pNp)
    C35=>ntY1(12*pNp+1:13*pNp);C45=>ntY1(13*pNp+1:14*pNp)
    C55=>ntY1(14*pNp+1:15*pNp);C16=>ntY1(15*pNp+1:16*pNp)
    C26=>ntY1(16*pNp+1:17*pNp);C36=>ntY1(17*pNp+1:18*pNp)
    C46=>ntY1(18*pNp+1:19*pNp);C56=>ntY1(19*pNp+1:20*pNp)
    C66=>ntY1(20*pNp+1:21*pNp)

    P11   =>ntY1(       1:   pNp);P12   =>ntY1(   pNp+1: 2*pNp)
    P22   =>ntY1( 2*pNp+1: 3*pNp);P13   =>ntY1( 3*pNp+1: 4*pNp)
    P23   =>ntY1( 4*pNp+1: 5*pNp);P33   =>ntY1( 5*pNp+1: 6*pNp)
    M11   =>ntY1( 6*pNp+1: 7*pNp);M12   =>ntY1( 7*pNp+1: 8*pNp)
    M22   =>ntY1( 8*pNp+1: 9*pNp);M13   =>ntY1( 9*pNp+1:10*pNp)
    M23   =>ntY1(10*pNp+1:11*pNp);M33   =>ntY1(11*pNp+1:12*pNp)
    P11P22=>ntY1(12*pNp+1:13*pNp);P11P33=>ntY1(13*pNp+1:14*pNp)
    P22P33=>ntY1(14*pNp+1:15*pNp);M11M22=>ntY1(15*pNp+1:16*pNp)
    M11M33=>ntY1(16*pNp+1:17*pNp);M22M33=>ntY1(17*pNp+1:18*pNp)
    M11P22=>ntY1(18*pNp+1:19*pNp);M11P33=>ntY1(19*pNp+1:20*pNp)
    M22P33=>ntY1(20*pNp+1:21*pNp);P11M22=>ntY1(21*pNp+1:22*pNp)
    P11M33=>ntY1(22*pNp+1:23*pNp);P22M33=>ntY1(23*pNp+1:24*pNp)
    ZERO  =>ntY1(24*pNp+1:25*pNp)

    CC(1,1)%p=>C11; CC(1,2)%p=>C16; CC(1,3)%p=>C15
    CC(2,1)%p=>C16; CC(2,2)%p=>C66; CC(2,3)%p=>C56
    CC(3,1)%p=>C15; CC(3,2)%p=>C56; CC(3,3)%p=>C55
    CC(4,1)%p=>C16; CC(4,2)%p=>C66; CC(4,3)%p=>C56
    CC(5,1)%p=>C12; CC(5,2)%p=>C26; CC(5,3)%p=>C25
    CC(6,1)%p=>C14; CC(6,2)%p=>C46; CC(6,3)%p=>C45
    CC(7,1)%p=>C15; CC(7,2)%p=>C56; CC(7,3)%p=>C55
    CC(8,1)%p=>C14; CC(8,2)%p=>C46; CC(8,3)%p=>C45
    CC(9,1)%p=>C13; CC(9,2)%p=>C36; CC(9,3)%p=>C35
    CC(1,4)%p=>C16; CC(1,5)%p=>C12; CC(1,6)%p=>C14
    CC(2,4)%p=>C66; CC(2,5)%p=>C26; CC(2,6)%p=>C46
    CC(3,4)%p=>C56; CC(3,5)%p=>C25; CC(3,6)%p=>C45
    CC(4,4)%p=>C66; CC(4,5)%p=>C26; CC(4,6)%p=>C46
    CC(5,4)%p=>C26; CC(5,5)%p=>C22; CC(5,6)%p=>C24
    CC(6,4)%p=>C46; CC(6,5)%p=>C24; CC(6,6)%p=>C44
    CC(7,4)%p=>C56; CC(7,5)%p=>C25; CC(7,6)%p=>C45
    CC(8,4)%p=>C46; CC(8,5)%p=>C24; CC(8,6)%p=>C44
    CC(9,4)%p=>C36; CC(9,5)%p=>C23; CC(9,6)%p=>C34
    CC(1,7)%p=>C15; CC(1,8)%p=>C14; CC(1,9)%p=>C13
    CC(2,7)%p=>C56; CC(2,8)%p=>C46; CC(2,9)%p=>C36
    CC(3,7)%p=>C55; CC(3,8)%p=>C45; CC(3,9)%p=>C35
    CC(4,7)%p=>C56; CC(4,8)%p=>C46; CC(4,9)%p=>C36
    CC(5,7)%p=>C25; CC(5,8)%p=>C24; CC(5,9)%p=>C23
    CC(6,7)%p=>C45; CC(6,8)%p=>C44; CC(6,9)%p=>C34
    CC(7,7)%p=>C55; CC(7,8)%p=>C45; CC(7,9)%p=>C35
    CC(8,7)%p=>C45; CC(8,8)%p=>C44; CC(8,9)%p=>C34
    CC(9,7)%p=>C35; CC(9,8)%p=>C34; CC(9,9)%p=>C33

    DD(1,1)%p=>ZERO  ;DD(1,2)%p=>P12   ;DD(1,3)%p=>P13   
    DD(2,1)%p=>P12   ;DD(2,2)%p=>M11P22;DD(2,3)%p=>P23   
    DD(3,1)%p=>P13   ;DD(3,2)%p=>P23   ;DD(3,3)%p=>M11P33
    DD(4,1)%p=>M12   ;DD(4,2)%p=>M11M22;DD(4,3)%p=>M23   
    DD(5,1)%p=>P11P22;DD(5,2)%p=>M12   ;DD(5,3)%p=>P13   
    DD(6,1)%p=>P23   ;DD(6,2)%p=>M13   ;DD(6,3)%p=>M12   
    DD(7,1)%p=>M13   ;DD(7,2)%p=>M23   ;DD(7,3)%p=>M11M33
    DD(8,1)%p=>P23   ;DD(8,2)%p=>M13   ;DD(8,3)%p=>M12   
    DD(9,1)%p=>P11P33;DD(9,2)%p=>P12   ;DD(9,3)%p=>M13   
    DD(1,4)%p=>M12   ;DD(1,5)%p=>P11P22;DD(1,6)%p=>P23   
    DD(2,4)%p=>M11M22;DD(2,5)%p=>M12   ;DD(2,6)%p=>M13   
    DD(3,4)%p=>M23   ;DD(3,5)%p=>P13   ;DD(3,6)%p=>M12   
    DD(4,4)%p=>P11M22;DD(4,5)%p=>P12   ;DD(4,6)%p=>P13   
    DD(5,4)%p=>P12   ;DD(5,5)%p=>ZERO  ;DD(5,6)%p=>P23   
    DD(6,4)%p=>P13   ;DD(6,5)%p=>P23   ;DD(6,6)%p=>M22P33
    DD(7,4)%p=>M23   ;DD(7,5)%p=>P13   ;DD(7,6)%p=>M12   
    DD(8,4)%p=>M13   ;DD(8,5)%p=>M23   ;DD(8,6)%p=>M22M33
    DD(9,4)%p=>P12   ;DD(9,5)%p=>P22P33;DD(9,6)%p=>M23   
    DD(1,7)%p=>M13   ;DD(1,8)%p=>P23   ;DD(1,9)%p=>P11P33
    DD(2,7)%p=>M23   ;DD(2,8)%p=>M13   ;DD(2,9)%p=>P12   
    DD(3,7)%p=>M11M33;DD(3,8)%p=>M12   ;DD(3,9)%p=>M13   
    DD(4,7)%p=>M23   ;DD(4,8)%p=>M13   ;DD(4,9)%p=>P12   
    DD(5,7)%p=>P13   ;DD(5,8)%p=>M23   ;DD(5,9)%p=>P22P33 
    DD(6,7)%p=>M12   ;DD(6,8)%p=>M22M33;DD(6,9)%p=>M23   
    DD(7,7)%p=>P11M33;DD(7,8)%p=>P12   ;DD(7,9)%p=>P13   
    DD(8,7)%p=>P12   ;DD(8,8)%p=>P22M33;DD(8,9)%p=>P23   
    DD(9,7)%p=>P13   ;DD(9,8)%p=>P23   ;DD(9,9)%p=>ZERO  

    GG(1)%p=>P11; GG(2)%p=>P12; GG(3)%p=>P13
    GG(4)%p=>P12; GG(5)%p=>P22; GG(6)%p=>P23
    GG(7)%p=>P13; GG(8)%p=>P23; GG(9)%p=>P33

    do i=1,3
        pG(i)%p=>ntG((i-1)*Nfp+1:i*Nfp)
        pS(i)%p=>ntS((i-1)*Nfp+1:i*Nfp)
        pT(i)%p=>pTf((i-1)*Nfp+1:i*Nfp)
        pS0(i)%p=>ntS0((i-1)*Nfp+1:i*Nfp)
        do j=1,3
            pH(j,i)%p=>ntH((j-1)*Nfp+1:j*Nfp,(i-1)*Nfp+1:i*Nfp)
        enddo
    enddo

end subroutine rupt_Newton_allocate

subroutine rupt_matrix_allocate(pNp,Nfp,M)
    integer,intent(in) :: pNp,Nfp
    type(rupture_matrices) :: M
    integer :: i,j

    allocate(M%Dm(pNp,pNp,3));M%Dm=0d0
    allocate(M%Dp(pNp,pNp,3));M%Dp=0d0
    allocate(M%Dn(pNp,pNp,2));M%Dn=0d0
    allocate(M%Jm(pNp,Nfp));M%Jm=0d0
    allocate(M%Jp(pNp,Nfp));M%Jp=0d0
    allocate(M%Lm(Nfp,pNp));M%Lm=0d0
    allocate(M%Lp(Nfp,pNp));M%Lp=0d0
    allocate(M%Cm(pNp,9,9));M%Cm=0d0
    allocate(M%Cp(pNp,9,9));M%Cp=0d0
    allocate(M%irm(pNp));M%irm=0d0
    allocate(M%irp(pNp));M%irp=0d0
    allocate(M%irDDm(pNp*3,pNp*9));M%irDDm=0d0
    allocate(M%irDDp(pNp*3,pNp*9));M%irDDp=0d0
    allocate(M%BVmEm(pNp*3,pNp*9))
    allocate(M%BVmEp(pNp*3,pNp*9))
    allocate(M%BVpEm(pNp*3,pNp*9))
    allocate(M%BVpEp(pNp*3,pNp*9))

    if(convtest)then
        allocate(M%DDm     (pNp*9,pNp*3));M%DDm     =0d0
        allocate(M%DDp     (pNp*9,pNp*3));M%DDp     =0d0
        allocate(M%irJLmm  (pNp*3,pNp*3));M%irJLmm  =0d0
        allocate(M%irJLmp  (pNp*3,pNp*3));M%irJLmp  =0d0
        allocate(M%irJLpm  (pNp*3,pNp*3));M%irJLpm  =0d0
        allocate(M%irJLpp  (pNp*3,pNp*3));M%irJLpp  =0d0
        allocate(M%RJLmm   (pNp*9,pNp*3));M%RJLmm   =0d0
        allocate(M%RJLmp   (pNp*9,pNp*3));M%RJLmp   =0d0
        allocate(M%RJLpm   (pNp*9,pNp*3));M%RJLpm   =0d0
        allocate(M%RJLpp   (pNp*9,pNp*3));M%RJLpp   =0d0
        allocate(M%irDCm   (pNp*3,pNp*9));M%irDCm   =0d0
        allocate(M%irDCp   (pNp*3,pNp*9));M%irDCp   =0d0
        allocate(M%irJLRCmm(pNp*3,pNp*9));M%irJLRCmm=0d0
        allocate(M%irJLRCmp(pNp*3,pNp*9));M%irJLRCmp=0d0
        allocate(M%irJLRCpm(pNp*3,pNp*9));M%irJLRCpm=0d0
        allocate(M%irJLRCpp(pNp*3,pNp*9));M%irJLRCpp=0d0
        allocate(M%irJm    (pNp*3,Nfp*3));M%irJm    =0d0
        allocate(M%irJp    (pNp*3,Nfp*3));M%irJp    =0d0
        allocate(M%RJLRCmm (pNp*9,pNp*9));M%RJLRCmm =0d0
        allocate(M%RJLRCmp (pNp*9,pNp*9));M%RJLRCmp =0d0
        allocate(M%RJLRCpm (pNp*9,pNp*9));M%RJLRCpm =0d0
        allocate(M%RJLRCpp (pNp*9,pNp*9));M%RJLRCpp =0d0
        allocate(M%RJm     (pNp*9,Nfp*3));M%irJm    =0d0
        allocate(M%RJp     (pNp*9,Nfp*3));M%irJp    =0d0
        allocate(M%KLRCm   (Nfp*3,pNp*9));M%KLRCm   =0d0
        allocate(M%KLRCp   (Nfp*3,pNp*9));M%KLRCp   =0d0
    endif
    allocate(M%BLRCm   (Nfp  ,pNp*9));M%BLRCm   =0d0
    allocate(M%BLRCp   (Nfp  ,pNp*9));M%BLRCp   =0d0

    allocate(M%A11(30*pNp,30*pNp));M%A11=0d0
    allocate(M%A12(30*pNp, 7*Nfp));M%A12=0d0
    allocate(M%A21( 4*Nfp,30*pNp));M%A21=0d0
    allocate(M%A22( 4*Nfp, 7*Nfp));M%A22=0d0

    allocate(M%pA11(30*pNp))
    allocate(M%pAb1( 4*Nfp))

    M%AUmUm=>M%A11(iUm1 +1:iUm3 +pNp,iUm1 +1:iUm3 +pNp)
    M%AUpUp=>M%A11(iUp1 +1:iUp3 +pNp,iUp1 +1:iUp3 +pNp)

    M%AUmVm=>M%A11(iUm1 +1:iUm3 +pNp,iVm1 +1:iVm3 +pNp)
    M%AUpVp=>M%A11(iUp1 +1:iUp3 +pNp,iVp1 +1:iVp3 +pNp)

    M%AEmUm=>M%A11(iEm11+1:iEm33+pNp,iUm1 +1:iUm3 +pNp)
    M%AEmUp=>M%A11(iEm11+1:iEm33+pNp,iUp1 +1:iUp3 +pNp)
    M%AEpUm=>M%A11(iEp11+1:iEp33+pNp,iUm1 +1:iUm3 +pNp)
    M%AEpUp=>M%A11(iEp11+1:iEp33+pNp,iUp1 +1:iUp3 +pNp)

    M%AEmEm=>M%A11(iEm11+1:iEm33+pNp,iEm11+1:iEm33+pNp)
    M%AEmEp=>M%A11(iEm11+1:iEm33+pNp,iEp11+1:iEp33+pNp)
    M%AEpEm=>M%A11(iEp11+1:iEp33+pNp,iEm11+1:iEm33+pNp)
    M%AEpEp=>M%A11(iEp11+1:iEp33+pNp,iEp11+1:iEp33+pNp)

    M%AEmVm=>M%A11(iEm11+1:iEm33+pNp,iVm1 +1:iVm3 +pNp)
    M%AEmVp=>M%A11(iEm11+1:iEm33+pNp,iVp1 +1:iVp3 +pNp)
    M%AEpVm=>M%A11(iEp11+1:iEp33+pNp,iVm1 +1:iVm3 +pNp)
    M%AEpVp=>M%A11(iEp11+1:iEp33+pNp,iVp1 +1:iVp3 +pNp)

    M%AEmVt=>M%A12(iEm11+1:iEm33+pNp,iVt1 +1:iVt3 +Nfp)
    M%AEpVt=>M%A12(iEp11+1:iEp33+pNp,iVt1 +1:iVt3 +Nfp)

    M%AVmUm=>M%A11(iVm1 +1:iVm3 +pNp,iUm1 +1:iUm3 +pNp)
    M%AVmUp=>M%A11(iVm1 +1:iVm3 +pNp,iUp1 +1:iUp3 +pNp)
    M%AVpUm=>M%A11(iVp1 +1:iVp3 +pNp,iUm1 +1:iUm3 +pNp)
    M%AVpUp=>M%A11(iVp1 +1:iVp3 +pNp,iUp1 +1:iUp3 +pNp)

    M%AVmEm=>M%A11(iVm1 +1:iVm3 +pNp,iEm11+1:iEm33+pNp)
    M%AVmEp=>M%A11(iVm1 +1:iVm3 +pNp,iEp11+1:iEp33+pNp)
    M%AVpEm=>M%A11(iVp1 +1:iVp3 +pNp,iEm11+1:iEm33+pNp)
    M%AVpEp=>M%A11(iVp1 +1:iVp3 +pNp,iEp11+1:iEp33+pNp)

    M%AVmVm=>M%A11(iVm1 +1:iVm3 +pNp,iVm1 +1:iVm3 +pNp)
    M%AVmVp=>M%A11(iVm1 +1:iVm3 +pNp,iVp1 +1:iVp3 +pNp)
    M%AVpVm=>M%A11(iVp1 +1:iVp3 +pNp,iVm1 +1:iVm3 +pNp)
    M%AVpVp=>M%A11(iVp1 +1:iVp3 +pNp,iVp1 +1:iVp3 +pNp)

    M%AVmVt=>M%A12(iVm1 +1:iVm3 +pNp,iVt1 +1:iVt3 +Nfp)
    M%AVpVt=>M%A12(iVp1 +1:iVp3 +pNp,iVt1 +1:iVt3 +Nfp)

    M%ATfUm=>M%A21(iTf1 +1:iTf3 +Nfp,iUm1 +1:iUm3 +pNp)
    M%ATfUp=>M%A21(iTf1 +1:iTf3 +Nfp,iUp1 +1:iUp3 +pNp)

    M%ATfEm=>M%A21(iTf1 +1:iTf3 +Nfp,iEm11+1:iEm33+pNp)
    M%ATfEp=>M%A21(iTf1 +1:iTf3 +Nfp,iEp11+1:iEp33+pNp)

    M%ATnUm=>M%A21(iSig +1:iSig +Nfp,iUm1 +1:iUm3 +pNp)
    M%ATnUp=>M%A21(iSig +1:iSig +Nfp,iUp1 +1:iUp3 +pNp)

    M%ATnEm=>M%A21(iSig +1:iSig +Nfp,iEm11+1:iEm33+pNp)
    M%ATnEp=>M%A21(iSig +1:iSig +Nfp,iEp11+1:iEp33+pNp)

    M%ATfTf=>M%A22(iTf1 +1:iTf3 +Nfp,iTf1 +1:iTf3 +Nfp)
    M%ATnTn=>M%A22(iSig +1:iSig +Nfp,iSig +1:iSig +Nfp)

    M%Ab1=>M%A22(1:Nfp*4,1:Nfp*4)
    M%Ab2=>M%A22(1:Nfp*4,Nfp*4+1:Nfp*7)
    do i=1,3
        do j=1,4
            M%Abb(j,i)%p=>M%Ab2((j-1)*Nfp+1:j*Nfp,(i-1)*Nfp+1:i*Nfp)
        enddo
    enddo

end subroutine rupt_matrix_allocate

subroutine init_quasi_static(pNp,Nfp,mesh,rupt,T0,Vt0)
    integer,intent(in) :: pNp,Nfp
    type(tetmesh_geometry) :: mesh
    type(rupture) :: rupt
    type(tensor_array) :: T0
    real(kind=rkind) :: Vt0(:)
    integer :: map(Nfp),itet,itri,jtet,jtri,iface,i,j,k
    real(kind=rkind) :: S11,S22,S33,S23,S13,S12,n1,n2,n3,nn,&
        t1,t2,t3,sigma,tauf,v1,v2,v3,Vt
    real(kind=rkind) :: a,b,L,v0,f0
    do iface=1,rupt%Nhface
        itet=rupt%T2E(1,iface);itri=rupt%T2E(3,iface)
        n1=mesh%nx(itri,itet)
        n2=mesh%ny(itri,itet)
        n3=mesh%nz(itri,itet)
!if(.true.)then
!    jtet=rupt%T2E(2,iface);jtri=rupt%T2E(4,iface)
!    n2=0d0
!    nn=sqrt(n1**2+n3**2)
!    n1=n1/nn;n3=n3/nn
!    mesh%nx(itri,itet)=n1
!    mesh%ny(itri,itet)=n2
!    mesh%nz(itri,itet)=n3
!    mesh%nx(jtri,jtet)=n1
!    mesh%ny(jtri,jtet)=n2
!    mesh%nz(jtri,jtet)=n3
!endif
        map(rupt%perm(:,1,iface))=&
            mesh%vmapM((itri-1)*Nfp+1:itri*Nfp)+(itet-1)*pNp
        k=(iface-1)*Nfp
        do j=1,Nfp
            S11=T0%xx(map(j));S22=T0%yy(map(j));S33=T0%zz(map(j))
            S23=T0%yz(map(j));S13=T0%xz(map(j));S12=T0%xy(map(j))
            t1=n1*S11+n2*S12+n3*S13
            t2=n1*S12+n2*S22+n3*S23
            t3=n1*S13+n2*S23+n3*S33
            sigma=-n1*t1-n2*t2-n3*t3
            t1=t1+sigma*n1
            t2=t2+sigma*n2
            t3=t3+sigma*n3
            tauf=sqrt(t1**2+t2**2+t3**2)
            rupt%sigma0(k+j)=sigma
            rupt%tau0%x(k+j)=t1
            rupt%tau0%y(k+j)=t2
            rupt%tau0%z(k+j)=t3
            Vt=Vt0(k+j)
            a =rupt%a (k+j)
            b =rupt%b (k+j)
            L =rupt%L (k+j)
            v0=rupt%v0(k+j)
            f0=rupt%f0(k+j)
            rupt%psi(k+j)=RSfric_invF(tauf,sigma,Vt,a,b,L,v0,f0)
            t1=t1/tauf*Vt
            t2=t2/tauf*Vt
            t3=t3/tauf*Vt
            rupt%Vt0%x(k+j)=t1
            rupt%Vt0%y(k+j)=t2
            rupt%Vt0%z(k+j)=t3
        enddo
    enddo
    rupt%sigma=rupt%sigma0
    rupt%tauf%array=rupt%tau0%array
    rupt%Vt%array=rupt%Vt0%array
end subroutine init_quasi_static

subroutine form_rupt_matrices0(pNp,Nfp,mesh,rupt,mat,matrix)
    integer,intent(in) :: pNp,Nfp
    type(tetmesh_geometry) :: mesh
    type(rupture) :: rupt
    type(tetmesh_material) :: mat
    type(matrices) :: matrix
    real(kind=rkind) :: n(3)
    integer :: iface,itet1,itet2,itri1,itri2,i,j,k
    integer :: mapm(Nfp),mapp(Nfp),vmapm(pNp),vmapp(pNp)
    real(kind=rkind),pointer :: irm(:),irp(:),Jm(:,:),Jp(:,:),&
        Lm(:,:),Lp(:,:),Dm(:,:,:),Dp(:,:,:),Dn(:,:,:),&
        Cm(:,:,:),Cp(:,:,:)
    real(kind=rkind) :: invJ(9)
    real(kind=rkind) :: tmp1(Nfp),tmp2(Nfp),maxerr

    do iface=1,rupt%Nhface
        itet1=rupt%T2E(1,iface);itri1=rupt%T2E(3,iface)
        itet2=rupt%T2E(2,iface);itri2=rupt%T2E(4,iface)

        n(1)=mesh%nx(itri1,itet1)
        n(2)=mesh%ny(itri1,itet1)
        n(3)=mesh%nz(itri1,itet1)

        do i=1,pNp
            vmapm(i)=(itet1-1)*pNp+i
            vmapp(i)=(itet2-1)*pNp+i
        enddo

        irm=>rupt%M(iface)%irm
        irp=>rupt%M(iface)%irp
        Jm=>rupt%M(iface)%Jm
        Jp=>rupt%M(iface)%Jp
        Lm=>rupt%M(iface)%Lm
        Lp=>rupt%M(iface)%Lp
        Dm=>rupt%M(iface)%Dm
        Dp=>rupt%M(iface)%Dp
        Dn=>rupt%M(iface)%Dn
        Cm=>rupt%M(iface)%Cm
        Cp=>rupt%M(iface)%Cp

!---------------------------------------------------------------------!

        irm=1d0/mat%rho(vmapm)
        Jm=0d0
        mapm=rupt%perm(:,1,iface)
        Jm(:,mapm) = mesh%Fscale(itri1,itet1) * &
            matrix%Lift(:,(itri1-1)*Nfp+1:itri1*Nfp)
        Lm=0d0
        mapm(rupt%perm(:,1,iface))=&
            mesh%vmapM((itri1-1)*Nfp+1:itri1*Nfp)
        do i=1,Nfp
            Lm(i,mapm(i))=1d0
        enddo

        tmp1=matmul(Lm,mesh%coord%x((itet1-1)*pNp+1:itet1*pNp))
        tmp2=mesh%coord%x(mapm+(itet1-1)*pNp)
        maxerr=maxval(abs(tmp1-tmp2))
        tmp1=matmul(Lm,mesh%coord%y((itet1-1)*pNp+1:itet1*pNp))
        tmp2=mesh%coord%y(mapm+(itet1-1)*pNp)
        maxerr=maxval(abs(tmp1-tmp2))
        if(maxerr.gt.1d-16)print*,'err=',maxerr

        invJ=-mesh%invJac(itet1,1:9)
        Dm(:,:,1)=invJ(1)*matrix%Dr+invJ(2)*matrix%Ds+invJ(3)*matrix%Dt
        Dm(:,:,2)=invJ(4)*matrix%Dr+invJ(5)*matrix%Ds+invJ(6)*matrix%Dt
        Dm(:,:,3)=invJ(7)*matrix%Dr+invJ(8)*matrix%Ds+invJ(9)*matrix%Dt
        Dn(:,:,1)=Dm(:,:,1)*n(1)+Dm(:,:,2)*n(2)+Dm(:,:,3)*n(3);

        call rupt_asign_mat(pNp,mat,vmapm,&
            C11,C12,C22,C13,C23,C33,C14,C24,C34,C44,&
            C15,C25,C35,C45,C55,C16,C26,C36,C46,C56,C66)

        do j=1,9
            do i=1,9
                Cm(:,i,j)=CC(i,j)%p
            enddo
        enddo

!---------------------------------------------------------------------!

        irp=1d0/mat%rho(vmapp)
        Jp=0d0
        mapp=rupt%perm(:,2,iface)
        Jp(:,mapp) = mesh%Fscale(itri2,itet2) * &
            matrix%Lift(:,(itri2-1)*Nfp+1:itri2*Nfp)
        Lp=0d0
        mapp(rupt%perm(:,2,iface))=&
            mesh%vmapM((itri2-1)*Nfp+1:itri2*Nfp)
        do i=1,Nfp
            Lp(i,mapp(i))=1d0
        enddo

        invJ=-mesh%invJac(itet2,1:9)
        Dp(:,:,1)=invJ(1)*matrix%Dr+invJ(2)*matrix%Ds+invJ(3)*matrix%Dt
        Dp(:,:,2)=invJ(4)*matrix%Dr+invJ(5)*matrix%Ds+invJ(6)*matrix%Dt
        Dp(:,:,3)=invJ(7)*matrix%Dr+invJ(8)*matrix%Ds+invJ(9)*matrix%Dt
        Dn(:,:,2)=Dp(:,:,1)*n(1)+Dp(:,:,2)*n(2)+Dp(:,:,3)*n(3);

        call rupt_asign_mat(pNp,mat,vmapp,&
            C11,C12,C22,C13,C23,C33,C14,C24,C34,C44,&
            C15,C25,C35,C45,C55,C16,C26,C36,C46,C56,C66)

        do j=1,9
            do i=1,9
                Cp(:,i,j)=CC(i,j)%p
            enddo
        enddo

    enddo
end subroutine form_rupt_matrices0

subroutine form_rupt_matrices1(pNp,Nfp,mesh,rupt,T0,matrix,dt)
    integer,intent(in) :: pNp,Nfp
    type(tetmesh_geometry) :: mesh
    type(rupture) :: rupt
    type(tensor_array) :: T0
    type(tetmesh_material) :: mat
    type(matrices) :: matrix
    real(kind=rkind),intent(in) :: dt
    real(kind=rkind) :: n(3),nK(3,3),KR(3,9),BR(9)
    integer :: iface,itet,itri,i,j,k,ii,jj,kk,info
    integer :: map(Nfp),vmap(pNp)
    real(kind=rkind),parameter :: eye(3,3) = &
        (/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/)
    real(kind=rkind),pointer :: Cm(:,:,:),Cp(:,:,:),irm(:),irp(:),&
        Dm(:,:,:),Dp(:,:,:),Dn(:,:,:),&
        Jm(:,:),Jp(:,:),Lm(:,:),Lp(:,:),mpt(:)
    real(kind=rkind) :: &
        JmLm(pNp,pNp),JpLm(pNp,pNp),JmLp(pNp,pNp),JpLp(pNp,pNp),&
        DDm     (pNp*9,pNp*3),DDp     (pNp*9,pNp*3),&
        irJLmm  (pNp*3,pNp*3),irJLmp  (pNp*3,pNp*3),&
        irJLpm  (pNp*3,pNp*3),irJLpp  (pNp*3,pNp*3),&
        irnnJLmm(pNp*3,pNp*3),irnnJLmp(pNp*3,pNp*3),&
        irnnJLpm(pNp*3,pNp*3),irnnJLpp(pNp*3,pNp*3),&
        irJm    (pNp*3,Nfp*3),irJp    (pNp*3,Nfp*3),&
        RJLmm   (pNp*9,pNp*3),RJLmp   (pNp*9,pNp*3),&
        RJLpm   (pNp*9,pNp*3),RJLpp   (pNp*9,pNp*3),&
        RJm     (pNp*9,Nfp*3),RJp     (pNp*9,Nfp*3),&
        irDCm   (pNp*3,pNp*9),irDCp   (pNp*3,pNp*9),&
        irDDm   (pNp*3,pNp*9),irDDp   (pNp*3,pNp*9),&
        irJLRCmm(pNp*3,pNp*9),irJLRCmp(pNp*3,pNp*9),&
        irJLRCpm(pNp*3,pNp*9),irJLRCpp(pNp*3,pNp*9),&
        RJLRCmm (pNp*9,pNp*9),RJLRCmp (pNp*9,pNp*9),&
        RJLRCpm (pNp*9,pNp*9),RJLRCpp (pNp*9,pNp*9),&
        KLRCm   (Nfp*3,pNp*9),KLRCp   (Nfp*3,pNp*9),&
        BLRCm   (Nfp  ,pNp*9),BLRCp   (Nfp  ,pNp*9),&
        RJQmm   (pNp*9,pNp*3),RJQmp   (pNp*9,pNp*3),&
        RJQpm   (pNp*9,pNp*3),RJQpp   (pNp*9,pNp*3),&
        irJQmm  (pNp*3,pNp*3),irJQmp  (pNp*3,pNp*3),&
        irJQpm  (pNp*3,pNp*3),irJQpp  (pNp*3,pNp*3),&
        KQm     (Nfp*3,pNp*3),KQp     (Nfp*3,pNp*3),&
        BQm     (Nfp  ,pNp*3),BQp     (Nfp  ,pNp*3),&
        RCm(pNp,3,9),RCp(pNp,3,9),&
        KRCm(pNp,3,9),KRCp(pNp,3,9),BRCm(pNp,9),BRCp(pNp,9),&
        Dx(Nfp,Nfp),Dy(Nfp,Nfp),Dz(Nfp,Nfp),&
        DS(pNp,pNp),LDS(Nfp,pNp),JLDS(pNp,pNp),&
        tau0(pNp,3),tau0t(pNp,3),&
        Q(Nfp,3,3),tmpv(pNp),tmpf(Nfp)
    real(kind=rkind),target :: tmp(pNp,pNp)
    real(kind=rkind),pointer :: dtmp(:)

    dtmp=>tmp(1:pNp**2:pNp+1,1)

    ZERO=0d0
    do iface=1,rupt%Nhface
        itet=rupt%T2E(1,iface);itri=rupt%T2E(3,iface)

        n(1)=mesh%nx(itri,itet)
        n(2)=mesh%ny(itri,itet)
        n(3)=mesh%nz(itri,itet)
        do i=1,3
            do j=1,3
                nK(i,j)=n(i)*n(j)
            enddo
        enddo
        nK=eye-nK

        do i=1,pNp
            vmap(i)=(itet-1)*pNp+i
        enddo
        map(rupt%perm(:,1,iface))=&
            mesh%vmapM((itri-1)*Nfp+1:itri*Nfp)

        P11=T0%xx(vmap);P22=T0%yy(vmap);P33=T0%zz(vmap)
        P23=T0%yz(vmap);P13=T0%xz(vmap);P12=T0%xy(vmap)
        M11=-P11;M22=-P22;M33=-P33
        M23=-P23;M13=-P13;M12=-P12
        P11P22=P11+P22;P11P33=P11+P33;P22P33=P22+P33
        M11M22=M11+M22;M11M33=M11+M33;M22M33=M22+M33
        M11P22=M11+P22;M11P33=M11+P33;M22P33=M22+P33
        P11M22=P11+M22;P11M33=P11+M33;P22M33=P22+M33

        tau0(:,1)=n(1)*P11+n(2)*P12+n(3)*P13
        tau0(:,2)=n(1)*P12+n(2)*P22+n(3)*P23
        tau0(:,3)=n(1)*P13+n(2)*P23+n(3)*P33

        irm=>rupt%M(iface)%irm
        irp=>rupt%M(iface)%irp
        Jm=>rupt%M(iface)%Jm
        Jp=>rupt%M(iface)%Jp
        Lm=>rupt%M(iface)%Lm
        Lp=>rupt%M(iface)%Lp
        Cm=>rupt%M(iface)%Cm
        Cp=>rupt%M(iface)%Cp
        Dm=>rupt%M(iface)%Dm
        Dp=>rupt%M(iface)%Dp
        Dn=>rupt%M(iface)%Dn

!---------------------------------------------------------------------!

        DDm=0d0;DDp=0d0
        do i=1,3
            do j=1,3
                k=(i-1)*3+j
                DDm((k-1)*pNp+1:k*pNp,(i-1)*pNp+1:i*pNp)=Dm(:,:,j)
                DDp((k-1)*pNp+1:k*pNp,(i-1)*pNp+1:i*pNp)=Dp(:,:,j)
            enddo
        enddo

        JmLm=matmul(Jm,Lm);JmLp=matmul(Jm,Lp)
        JpLm=matmul(Jp,Lm);JpLp=matmul(Jp,Lp)

        RJLmm=0d0;RJLpm=0d0;RJLmp=0d0;RJLpp=0d0;RJm=0d0;RJp=0d0

        do i=1,3
            do j=1,3
                k=(i-1)*3+j
                RJLmm((k-1)*pNp+1:k*pNp,(i-1)*pNp+1:i*pNp)=n(j)*JmLm
                RJLpm((k-1)*pNp+1:k*pNp,(i-1)*pNp+1:i*pNp)=n(j)*JpLm
                RJLmp((k-1)*pNp+1:k*pNp,(i-1)*pNp+1:i*pNp)=n(j)*JmLp
                RJLpp((k-1)*pNp+1:k*pNp,(i-1)*pNp+1:i*pNp)=n(j)*JpLp
                RJm  ((k-1)*pNp+1:k*pNp,(i-1)*Nfp+1:i*Nfp)=n(j)*Jm
                RJp  ((k-1)*pNp+1:k*pNp,(i-1)*Nfp+1:i*Nfp)=n(j)*Jp
            enddo
        enddo

!open(1111,file='mat.dat',status='replace')
!write(1111,*)RJLmm
!write(1111,*)RJLpm
!write(1111,*)RJLmp
!write(1111,*)RJLpp
!stop

        irJLmm=0d0;irJLpm=0d0;irJLmp=0d0;irJLpp=0d0;irJm=0d0;irJp=0d0

        do i=1,3
            irJLmm((i-1)*pNp+1:i*pNp,(i-1)*pNp+1:i*pNp)=JmLm
            irJLpm((i-1)*pNp+1:i*pNp,(i-1)*pNp+1:i*pNp)=JpLm
            irJLmp((i-1)*pNp+1:i*pNp,(i-1)*pNp+1:i*pNp)=JmLp
            irJLpp((i-1)*pNp+1:i*pNp,(i-1)*pNp+1:i*pNp)=JpLp
            irJm  ((i-1)*pNp+1:i*pNp,(i-1)*Nfp+1:i*Nfp)=Jm
            irJp  ((i-1)*pNp+1:i*pNp,(i-1)*Nfp+1:i*Nfp)=Jp
        enddo
        do i=1,3
        do j=1,3
            irnnJLmm((j-1)*pNp+1:j*pNp,(i-1)*pNp+1:i*pNp)=JmLm*n(j)*n(i)
            irnnJLpm((j-1)*pNp+1:j*pNp,(i-1)*pNp+1:i*pNp)=JpLm*n(j)*n(i)
            irnnJLmp((j-1)*pNp+1:j*pNp,(i-1)*pNp+1:i*pNp)=JmLp*n(j)*n(i)
            irnnJLpp((j-1)*pNp+1:j*pNp,(i-1)*pNp+1:i*pNp)=JpLp*n(j)*n(i)
        enddo
        enddo

        call DiagMMw('l',pNp,pNp*9,irm,irJLmm,pNp)
        call DiagMMw('l',pNp,pNp*9,irp,irJLpm,pNp)
        call DiagMMw('l',pNp,pNp*9,irm,irJLmp,pNp)
        call DiagMMw('l',pNp,pNp*9,irp,irJLpp,pNp)
        call DiagMMw('l',pNp,Nfp*9,irm,irJm,pNp)
        call DiagMMw('l',pNp,Nfp*9,irp,irJp,pNp)
        call DiagMMw('l',pNp,pNp*9,irm,irnnJLmm,pNp)
        call DiagMMw('l',pNp,pNp*9,irp,irnnJLpm,pNp)
        call DiagMMw('l',pNp,pNp*9,irm,irnnJLmp,pNp)
        call DiagMMw('l',pNp,pNp*9,irp,irnnJLpp,pNp)

!---------------------------------------------------------------------!

        irDCm=0d0
        do j=1,9
            do i=1,3
                ii=(i-1)*pNp+1;jj=(j-1)*pNp+1
                do k=1,3
                    tmpv=Cm(:,(i-1)*3+k,j)+0.5d0*DD((i-1)*3+k,j)%p
                    call DiagMM('r',pNp,pNp,1d0,tmpv,Dm(:,:,k),&
                        irDCm(ii:,jj),pNp,pNp*3)
                enddo
            enddo
        enddo

        irDCp=0d0
        do j=1,9
            do i=1,3
                ii=(i-1)*pNp+1;jj=(j-1)*pNp+1
                do k=1,3
                    tmpv=Cp(:,(i-1)*3+k,j)+0.5d0*DD((i-1)*3+k,j)%p
                    call DiagMM('r',pNp,pNp,1d0,tmpv,Dp(:,:,k),&
                        irDCp(ii:,jj),pNp,pNp*3)
                enddo
            enddo
        enddo

        call DiagMMw('l',pNp,pNp*27,irm,irDCm,pNp)
        call DiagMMw('l',pNp,pNp*27,irp,irDCp,pNp)

        irDDm=0d0;irDDp=0d0
        do i=1,3
            do j=1,3
                k=(i-1)*3+j
                irDDm((i-1)*pNp+1:i*pNp,(k-1)*pNp+1:k*pNp)=Dm(:,:,j)
                irDDp((i-1)*pNp+1:i*pNp,(k-1)*pNp+1:k*pNp)=Dp(:,:,j)
            enddo
        enddo

        call DiagMMw('l',pNp,pNp*27,irm,irDDm,pNp)
        call DiagMMw('l',pNp,pNp*27,irp,irDDp,pNp)

!---------------------------------------------------------------------!

        RCm=0d0;RCp=0d0

        do j=1,9
            do i=1,3
                do k=1,3
                    tmpv=Cm(:,(i-1)*3+k,j)+0.5d0*DD((i-1)*3+k,j)%p
                    RCm(:,i,j)=RCm(:,i,j)+n(k)*tmpv
                    tmpv=Cp(:,(i-1)*3+k,j)+0.5d0*DD((i-1)*3+k,j)%p
                    RCp(:,i,j)=RCp(:,i,j)+n(k)*tmpv
                enddo
            enddo
        enddo

!open(1111,file='mat.dat',status='replace')
!write(1111,*)RCm
!write(1111,*)RCp
!stop

        KRCm=0d0;KRCp=0d0
        do i=1,3
            do j=1,9
                do k=1,3
                    KRCm(:,i,j)=KRCm(:,i,j)+nK(i,k)*RCm(:,k,j)
                    KRCp(:,i,j)=KRCp(:,i,j)+nK(i,k)*RCp(:,k,j)
                enddo
            enddo
        enddo
        BRCm=0d0;BRCp=0d0
        do j=1,9
            do k=1,3
                BRCm(:,j)=BRCm(:,j)+n(k)*RCm(:,k,j)
                BRCp(:,j)=BRCp(:,j)+n(k)*RCp(:,k,j)
            enddo
        enddo

!---------------------------------------------------------------------!

        JmLm=matmul(Jm,Lm);JmLp=matmul(Jm,Lp)
        JpLm=matmul(Jp,Lm);JpLp=matmul(Jp,Lp)

        irJLRCmm=0d0;irJLRCpm=0d0;irJLRCmp=0d0;irJLRCpp=0d0
        KLRCm=0d0;BLRCm=0d0;KLRCp=0d0;BLRCp=0d0
        do j=1,9
            jj=(j-1)*pNp+1
            do i=1,3
                ii=(i-1)*pNp+1
                call DiagMM('r',pNp,pNp,1d0,RCm(:,i,j),JmLm,&
                    irJLRCmm(ii:,jj),pNp,pNp*3)
                call DiagMM('r',pNp,pNp,1d0,RCm(:,i,j),JpLm,&
                    irJLRCpm(ii:,jj),pNp,pNp*3)
                call DiagMM('r',pNp,pNp,1d0,RCp(:,i,j),JmLp,&
                    irJLRCmp(ii:,jj),pNp,pNp*3)
                call DiagMM('r',pNp,pNp,1d0,RCp(:,i,j),JpLp,&
                    irJLRCpp(ii:,jj),pNp,pNp*3)
                ii=(i-1)*Nfp+1
                call DiagMM('r',pNp,Nfp,1d0,KRCm(:,i,j),Lm,&
                    KLRCm(ii:,jj),Nfp,Nfp*3)
                call DiagMM('r',pNp,Nfp,1d0,KRCp(:,i,j),Lp,&
                    KLRCp(ii:,jj),Nfp,Nfp*3)
            enddo
            call DiagMM('r',pNp,Nfp,1d0,BRCm(:,j),Lm,&
                BLRCm(1:,jj),Nfp,Nfp)
            call DiagMM('r',pNp,Nfp,1d0,BRCp(:,j),Lp,&
                BLRCp(1:,jj),Nfp,Nfp)
        enddo

        RJLRCmm=0d0;RJLRCpm=0d0;RJLRCmp=0d0;RJLRCpp=0d0
        do i=1,3
            do j=1,3
                k=(i-1)*3+j
                RJLRCmm((k-1)*pNp+1:k*pNp,:)=&
                    n(j)*irJLRCmm((i-1)*pNp+1:i*pNp,:)
                RJLRCpm((k-1)*pNp+1:k*pNp,:)=&
                    n(j)*irJLRCpm((i-1)*pNp+1:i*pNp,:)
                RJLRCmp((k-1)*pNp+1:k*pNp,:)=&
                    n(j)*irJLRCmp((i-1)*pNp+1:i*pNp,:)
                RJLRCpp((k-1)*pNp+1:k*pNp,:)=&
                    n(j)*irJLRCpp((i-1)*pNp+1:i*pNp,:)
            enddo
        enddo

        call DiagMMw('l',pNp,pNp*27,irm,irJLRCmm,pNp)
        call DiagMMw('l',pNp,pNp*27,irm,irJLRCmp,pNp)
        call DiagMMw('l',pNp,pNp*27,irp,irJLRCpm,pNp)
        call DiagMMw('l',pNp,pNp*27,irp,irJLRCpp,pNp)

!!---------------------------------------------------------------------!
! form Q

        irJQmm=0d0;irJQpm=0d0;BQm=0d0;KQm=0d0;
        irJQmp=0d0;irJQpp=0d0;BQp=0d0;KQp=0d0;

        tmpv=n(1)*tau0(:,1)+n(2)*tau0(:,2)+n(3)*tau0(:,3)
        tau0t(:,1)=tau0(:,1)-tmpv*n(1)
        tau0t(:,2)=tau0(:,2)-tmpv*n(2)
        tau0t(:,3)=tau0(:,3)-tmpv*n(3)

        do j=1,3
            jj=(j-1)*pNp+1
            DS=Dm(:,:,j)-Dn(:,:,1)*n(j)
            LDS=matmul(Lm,DS)
            call DiagMM('r',pNp,Nfp,1d0,tmpv,LDS,&
                BQm(1:,jj),Nfp,Nfp)
            do i=1,3
                ii=(i-1)*Nfp+1
                call DiagMM('r',pNp,Nfp,1d0,tau0t(:,i),LDS,&
                    KQm(ii:,jj),Nfp,Nfp*3)
            enddo
            JLDS=matmul(JmLm,DS)
            do i=1,3
                ii=(i-1)*pNp+1
                call DiagMM('r',pNp,pNp,1d0,tau0(:,i),JLDS,&
                    irJQmm(ii:,jj),pNp,pNp*3)
            enddo
            JLDS=matmul(JpLm,DS)
            do i=1,3
                ii=(i-1)*pNp+1
                call DiagMM('r',pNp,pNp,1d0,tau0(:,i),JLDS,&
                    irJQpm(ii:,jj),pNp,pNp*3)
            enddo
        enddo
        call DiagMMw('l',pNp,pNp*3,irm,irJQmm,pNp)
        call DiagMMw('l',pNp,pNp*3,irp,irJQpm,pNp)

        do j=1,3
            jj=(j-1)*pNp+1
            DS=Dp(:,:,j)-Dn(:,:,2)*n(j)
            LDS=matmul(Lp,DS)
            call DiagMM('r',pNp,Nfp,1d0,tmpv,LDS,&
                BQp(1:,jj),Nfp,Nfp)
            do i=1,3
                ii=(i-1)*Nfp+1
                call DiagMM('r',pNp,Nfp,1d0,tau0t(:,i),LDS,&
                    KQp(ii:,jj),Nfp,Nfp*3)
            enddo
            JLDS=matmul(JmLp,DS)
            do i=1,3
                ii=(i-1)*pNp+1
                call DiagMM('r',pNp,pNp,1d0,tau0(:,i),JLDS,&
                    irJQmp(ii:,jj),pNp,pNp*3)
            enddo
            JLDS=matmul(JpLp,DS)
            do i=1,3
                ii=(i-1)*pNp+1
                call DiagMM('r',pNp,pNp,1d0,tau0(:,i),JLDS,&
                    irJQpp(ii:,jj),pNp,pNp*3)
            enddo
        enddo
        call DiagMMw('l',pNp,pNp*3,irm,irJQmp,pNp)
        call DiagMMw('l',pNp,pNp*3,irp,irJQpp,pNp)

!---------------------------------------------------------------------!

        rupt%M(iface)%A11=0d0;rupt%M(iface)%A12=0d0
        rupt%M(iface)%A21=0d0;rupt%M(iface)%A22=0d0

        do i=1,3*pNp
            rupt%M(iface)%AUmVm(i,i)=-dt
            rupt%M(iface)%AUpVp(i,i)=-dt
        enddo

        rupt%M(iface)%AEmVm=dt*DDm-dt/2d0*RJLmm
        rupt%M(iface)%AEmVp=      -dt/2d0*RJLmp
        rupt%M(iface)%AEmVt=       dt/2d0*RJm
        rupt%M(iface)%AEpVm=       dt/2d0*RJLpm
        rupt%M(iface)%AEpVp=dt*DDp+dt/2d0*RJLpp
        rupt%M(iface)%AEpVt=       dt/2d0*RJp

        rupt%M(iface)%AVmUm= rupt_alpha_V*dt*irnnJLmm
        rupt%M(iface)%AVmUp=-rupt_alpha_V*dt*irnnJLmp
        rupt%M(iface)%AVpUm=-rupt_alpha_V*dt*irnnJLpm
        rupt%M(iface)%AVpUp= rupt_alpha_V*dt*irnnJLpp

        rupt%M(iface)%AVmVm= rupt_alpha_V*dt*irJLmm
        rupt%M(iface)%AVmVp=-rupt_alpha_V*dt*irJLmp
        rupt%M(iface)%AVpVm=-rupt_alpha_V*dt*irJLpm
        rupt%M(iface)%AVpVp= rupt_alpha_V*dt*irJLpp
        rupt%M(iface)%AVmVt= rupt_alpha_V*dt*irJm
        rupt%M(iface)%AVpVt=-rupt_alpha_V*dt*irJp

!        rupt%M(iface)%AVmVm= rupt%M(iface)%AVmVm &
!                           + rupt_gamma*dt*(irJLmm-irnnJLmm)
!        rupt%M(iface)%AVpVp= rupt%M(iface)%AVpVp &
!                           + rupt_gamma*dt*(irJLpp-irnnJLpp)

        rupt%M(iface)%AVmEm=(dt)*(irDCm-0.5d0*irJLRCmm+rupt_gamma*irDDm)
        rupt%M(iface)%AVmEp=(dt)*(     -0.5d0*irJLRCmp                 )
        rupt%M(iface)%AVpEm=(dt)*(      0.5d0*irJLRCpm                 )
        rupt%M(iface)%AVpEp=(dt)*(irDCp+0.5d0*irJLRCpp+rupt_gamma*irDDp)

!        rupt%M(iface)%BVmEm=    rupt_gamma *(irDCm-0.5d0*irJLRCmm)
!        rupt%M(iface)%BVmEp=    rupt_gamma *(     -0.5d0*irJLRCmp)
!        rupt%M(iface)%BVpEm=    rupt_gamma *(      0.5d0*irJLRCpm)
!        rupt%M(iface)%BVpEp=    rupt_gamma *(irDCp+0.5d0*irJLRCpp)

        do i=1,pNp
            rupt%M(iface)%A11(iVm1+i,iVm1+i)=&
            rupt%M(iface)%A11(iVm1+i,iVm1+i)+&
            rupt%M(iface)%irm(i)*rupt_gamma*dt
            rupt%M(iface)%A11(iVm2+i,iVm2+i)=&
            rupt%M(iface)%A11(iVm2+i,iVm2+i)+&
            rupt%M(iface)%irm(i)*rupt_gamma*dt
            rupt%M(iface)%A11(iVm3+i,iVm3+i)=&
            rupt%M(iface)%A11(iVm3+i,iVm3+i)+&
            rupt%M(iface)%irm(i)*rupt_gamma*dt
            rupt%M(iface)%A11(iVp1+i,iVp1+i)=&
            rupt%M(iface)%A11(iVp1+i,iVp1+i)+&
            rupt%M(iface)%irp(i)*rupt_gamma*dt
            rupt%M(iface)%A11(iVp2+i,iVp2+i)=&
            rupt%M(iface)%A11(iVp2+i,iVp2+i)+&
            rupt%M(iface)%irp(i)*rupt_gamma*dt
            rupt%M(iface)%A11(iVp3+i,iVp3+i)=&
            rupt%M(iface)%A11(iVp3+i,iVp3+i)+&
            rupt%M(iface)%irp(i)*rupt_gamma*dt
        enddo

        rupt%M(iface)%AEmEm= rupt_alpha_E*dt*RJLRCmm
        rupt%M(iface)%AEmEp=-rupt_alpha_E*dt*RJLRCmp
        rupt%M(iface)%AEpEm=-rupt_alpha_E*dt*RJLRCpm
        rupt%M(iface)%AEpEp= rupt_alpha_E*dt*RJLRCpp

        rupt%M(iface)%ATfEm=-0.5d0*KLRCm
        rupt%M(iface)%ATfEp=-0.5d0*KLRCp
! iter sigma
        rupt%M(iface)%ATnEm= 0.5d0*BLRCm
        rupt%M(iface)%ATnEp= 0.5d0*BLRCp

!!---------------------------------------------------------------------!
!
!        rupt%M(iface)%AEmUm=-rupt_alpha_E*dt*RJQm
!        rupt%M(iface)%AEmUp= rupt_alpha_E*dt*RJQm
!        rupt%M(iface)%AEpUm= rupt_alpha_E*dt*RJQp
!        rupt%M(iface)%AEpUp=-rupt_alpha_E*dt*RJQp

!        rupt%M(iface)%AVmUm=rupt%M(iface)%AVmUm-0.5d0*dt*irJQmm
!        rupt%M(iface)%AVmUp=rupt%M(iface)%AVmUp+0.5d0*dt*irJQmp
!        rupt%M(iface)%AVpUm=rupt%M(iface)%AVpUm-0.5d0*dt*irJQpm
!        rupt%M(iface)%AVpUp=rupt%M(iface)%AVpUp+0.5d0*dt*irJQpp
!
!        rupt%M(iface)%ATfUm= 0.5d0*KQm
!        rupt%M(iface)%ATfUp= 0.5d0*KQp
!        rupt%M(iface)%ATnUm=-0.5d0*BQm
!        rupt%M(iface)%ATnUp=-0.5d0*BQp

!---------------------------------------------------------------------!

        rupt%M(iface)%irDDm = irDDm 
        rupt%M(iface)%irDDp = irDDp 
        rupt%M(iface)%BLRCm = BLRCm
        rupt%M(iface)%BLRCp = BLRCp
        if(convtest)then
            rupt%M(iface)%DDm      = DDm
            rupt%M(iface)%DDp      = DDp
            rupt%M(iface)%RJLmm    = RJLmm
            rupt%M(iface)%RJLpm    = RJLpm
            rupt%M(iface)%RJLmp    = RJLmp
            rupt%M(iface)%RJLpp    = RJLpp
            rupt%M(iface)%irJLmm   = irJLmm
            rupt%M(iface)%irJLpm   = irJLpm
            rupt%M(iface)%irJLmp   = irJLmp
            rupt%M(iface)%irJLpp   = irJLpp
            rupt%M(iface)%irDCm    = irDCm    
            rupt%M(iface)%irDCp    = irDCp    
            rupt%M(iface)%irJLRCmm = irJLRCmm 
            rupt%M(iface)%irJLRCmp = irJLRCmp 
            rupt%M(iface)%irJLRCpm = irJLRCpm 
            rupt%M(iface)%irJLRCpp = irJLRCpp 
            rupt%M(iface)%irJm     = irJm 
            rupt%M(iface)%irJp     = irJp 
            rupt%M(iface)%RJLRCmm  = RJLRCmm
            rupt%M(iface)%RJLRCmp  = RJLRCmp
            rupt%M(iface)%RJLRCpm  = RJLRCpm
            rupt%M(iface)%RJLRCpp  = RJLRCpp
            rupt%M(iface)%RJm      = RJm 
            rupt%M(iface)%RJp      = RJp 
            rupt%M(iface)%KLRCm    = KLRCm
            rupt%M(iface)%KLRCp    = KLRCp
            cycle
        endif

        do i=1,iVp3+pNp
            rupt%M(iface)%A11(i,i)=rupt%M(iface)%A11(i,i)+1d0
        enddo

        do i=1,iSig+Nfp
            rupt%M(iface)%A22(i,i)=1d0
        enddo

        do i=iEm11+1,iEp33+pNp
            rupt%M(iface)%A11(i,i)=rupt%M(iface)%A11(i,i)+rupt_gamma*dt
        enddo
!!---------------------------------------------------------------------!

!if(not(debugging))then
        call dgetrf(30*pNp,30*pNp,rupt%M(iface)%A11,&
            30*pNp,rupt%M(iface)%pA11,info)
        call dgetrs('N',30*pNp,7*Nfp,rupt%M(iface)%A11,&
            30*pNp,rupt%M(iface)%pA11,rupt%M(iface)%A12,&
            30*pNp,info)
        rupt%M(iface)%A22=rupt%M(iface)%A22&
            -matmul(rupt%M(iface)%A21,rupt%M(iface)%A12)

        call dgetrf(4*Nfp,4*Nfp,rupt%M(iface)%Ab1,4*Nfp,&
            rupt%M(iface)%pAb1,info)
        call dgetrs('N',4*Nfp,3*Nfp,rupt%M(iface)%Ab1,4*Nfp,&
            rupt%M(iface)%pAb1,rupt%M(iface)%Ab2,4*Nfp,info)
!endif

!if(any(rupt%M(iface)%Ab2.ne.rupt%M(iface)%Ab2))then
!print*,rupt%M(iface)%Ab1
!print*,'ckp2: Ab2 nan'
!print*,rupt%M(iface)%Ab2
!stop
!endif
        rupt%M(iface)%Ab2=-rupt%M(iface)%Ab2

    enddo
end subroutine form_rupt_matrices1

subroutine rupt_asign_mat(Np,mat,map,&
        C11,C12,C22,C13,C23,C33,C14,C24,C34,C44,&
        C15,C25,C35,C45,C55,C16,C26,C36,C46,C56,C66)
    integer,intent(in) :: Np,map(Np)
    type(tetmesh_material) :: mat
    real(kind=rkind),intent(out) :: &
        C11(Np),&
        C12(Np),C22(Np),&
        C13(Np),C23(Np),C33(Np),&
        C14(Np),C24(Np),C34(Np),C44(Np),&
        C15(Np),C25(Np),C35(Np),C45(Np),C55(Np),&
        C16(Np),C26(Np),C36(Np),C46(Np),C56(Np),C66(Np)
    integer :: job
    job=mat%job
    if(job.eq.0)then
        C11=mat%C(map,1);C22=C11;C33=C11
        C23=C11;C13=C11;C12=C11
        C14=0d0;C24=0d0;C34=0d0;C44=0d0;C15=0d0
        C25=0d0;C35=0d0;C45=0d0;C55=0d0;C16=0d0
        C26=0d0;C36=0d0;C46=0d0;C56=0d0;C66=0d0
    elseif(job.eq.1)then
        C23=mat%C(map,1);C44=mat%C(map,2);
        C11=C23+C44;C22=C11;C33=C11;C44=C44/2d0
        C13=C23;C12=C23;C55=C44;C66=C44
        C14=0d0;C24=0d0;C34=0d0;C15=0d0
        C25=0d0;C35=0d0;C45=0d0;C16=0d0
        C26=0d0;C36=0d0;C46=0d0;C56=0d0
    elseif(job.eq.2)then
        C11=mat%C(map,1);C22=mat%C(map,2);C33=mat%C(map,3)
        C44=mat%C(map,4);C55=mat%C(map,5);C66=mat%C(map,6)
        C23=mat%C(map,7);C13=mat%C(map,8);C12=mat%C(map,9)
        C14=0d0;C24=0d0;C34=0d0;C15=0d0
        C25=0d0;C35=0d0;C45=0d0;C16=0d0
        C26=0d0;C36=0d0;C46=0d0;C56=0d0
    else
        C11=mat%C(map,1);C22=mat%C(map,2);C33=mat%C(map,3)
        C44=mat%C(map,4);C55=mat%C(map,5);C66=mat%C(map,6)
        C23=mat%C(map,7);C13=mat%C(map,8);C12=mat%C(map,9)
        C14=mat%C(map,10);C15=mat%C(map,11);C16=mat%C(map,12)
        C24=mat%C(map,13);C25=mat%C(map,14);C26=mat%C(map,15)
        C34=mat%C(map,16);C35=mat%C(map,17);C36=mat%C(map,18)
        C45=mat%C(map,19);C46=mat%C(map,20);C56=mat%C(map,21)
    endif
end subroutine rupt_asign_mat

subroutine rupt_nonslip_RHS(pNp,Nfp,mesh,rupt,E,V,rhsE,rhsV,dt)
    integer,intent(in) :: pNp,Nfp
    type(tetmesh_geometry) :: mesh
    type(rupture) :: rupt
    type(tensor_array) :: E,rhsE
    type(vector_array) :: V,rhsV
    integer :: mapm(Nfp),mapp(Nfp),vmapm(pNp),vmapp(pNp),map(Nfp)
    real(kind=rkind) :: n1,n2,n3,dt
    integer :: iface,itet1,itri1,itet2,itri2,i

    ntY1=0d0;ntYb=0d0

    do iface=1,rupt%Nhface
        itet1=rupt%T2E(1,iface)
        itri1=rupt%T2E(3,iface)
        mapm(rupt%perm(:,1,iface))=&
            mesh%vmapM((itri1-1)*Nfp+1:itri1*Nfp)+(itet1-1)*pNp
        itet2=rupt%T2E(2,iface)
        itri2=rupt%T2E(4,iface)
        mapp(rupt%perm(:,2,iface))=&
            mesh%vmapM((itri2-1)*Nfp+1:itri2*Nfp)+(itet2-1)*pNp
        do i=1,pNp
            vmapm(i)=(itet1-1)*pNp+i
            vmapp(i)=(itet2-1)*pNp+i
        enddo
        do i=1,Nfp
            map(i)=(iface-1)*Nfp+i
        enddo
        n1=mesh%nx(itri1,itet1)
        n2=mesh%ny(itri1,itet1)
        n3=mesh%nz(itri1,itet1)

        call rupt_map_asign_E(pNp,E,pEm,vmapm,0d0,1d0,'s')
        call rupt_map_asign_E(pNp,E,pEp,vmapp,0d0,1d0,'s')
        call rupt_map_asign_V(pNp,V,pVm,vmapm,0d0,1d0,'s')
        call rupt_map_asign_V(pNp,V,pVp,vmapp,0d0,1d0,'s')
        call rupt_map_asign_V(Nfp,rupt%Vt,ntS,map,0d0,1d0,'s')

        if(.true.)then
            hVm = &
                  matmul(rupt%M(iface)%irDCm,pEm)&
                - matmul(rupt%M(iface)%irJLRCmm,pEm)*0.5d0 &
                - matmul(rupt%M(iface)%irJLRCmp,pEp)*0.5d0 &
                + matmul(rupt%M(iface)%irJLmm,pVm)*rupt_alpha_V &
                - matmul(rupt%M(iface)%irJLmp,pVp)*rupt_alpha_V &
                + matmul(rupt%M(iface)%irJm  ,ntS)*rupt_alpha_V 

            hVp = &
                  matmul(rupt%M(iface)%irDCp,pEp)&
                + matmul(rupt%M(iface)%irJLRCpm,pEm)*0.5d0 &
                + matmul(rupt%M(iface)%irJLRCpp,pEp)*0.5d0 &
                - matmul(rupt%M(iface)%irJLpm,pVm)*rupt_alpha_V &
                + matmul(rupt%M(iface)%irJLpp,pVp)*rupt_alpha_V &
                - matmul(rupt%M(iface)%irJp  ,ntS)*rupt_alpha_V 

            hEm = &
                  matmul(rupt%M(iface)%DDm,pVm)&
                - matmul(rupt%M(iface)%RJLmm,pVm)*0.5d0 &
                - matmul(rupt%M(iface)%RJLmp,pVp)*0.5d0 &
                + matmul(rupt%M(iface)%RJm  ,ntS)*0.5d0 &
                + matmul(rupt%M(iface)%RJLRCmm,pEm)*rupt_alpha_E &
                - matmul(rupt%M(iface)%RJLRCmp,pEp)*rupt_alpha_E 

            hEp = &
                  matmul(rupt%M(iface)%DDp,pVp)&
                + matmul(rupt%M(iface)%RJLpm,pVm)*0.5d0 &
                + matmul(rupt%M(iface)%RJLpp,pVp)*0.5d0 &
                + matmul(rupt%M(iface)%RJp  ,ntS)*0.5d0 &
                - matmul(rupt%M(iface)%RJLRCpm,pEm)*rupt_alpha_E &
                + matmul(rupt%M(iface)%RJLRCpp,pEp)*rupt_alpha_E 
        else
            ntZ1=matmul(rupt%M(iface)%A11,ntY1)
            ntZ1=ntZ1+matmul(rupt%M(iface)%A12(:,Nfp*4+1:Nfp*7),ntS)
            ntZ1=ntZ1/dt
        endif

        call rupt_map_asign_V(pNp,rhsV,hVm,vmapm,1d0,-1d0,'v')
        call rupt_map_asign_V(pNp,rhsV,hVp,vmapp,1d0,-1d0,'v')
        call rupt_map_asign_E(pNp,rhsE,hEm,vmapm,1d0,-1d0,'v')
        call rupt_map_asign_E(pNp,rhsE,hEp,vmapp,1d0,-1d0,'v')

    enddo

end subroutine rupt_nonslip_RHS

subroutine rupt_conv_err(pNp,Nfp,mesh,rupt,E,S,time,rerr)
    real(kind=rkind) :: rerr(4),time
    integer,intent(in) :: pNp,Nfp
    type(tetmesh_geometry) :: mesh
    type(rupture) :: rupt
    type(tensor_array) :: E,S
    integer :: mapm(Nfp),mapp(Nfp),vmapm(pNp),vmapp(pNp)
    integer :: iface,itet1,itri1,itet2,itri2,i,j
    real(kind=rkind) :: n1,n2,n3,err
    rerr=0d0
    do iface=1,rupt%Nhface
        itet1=rupt%T2E(1,iface)
        itri1=rupt%T2E(3,iface)
        itet2=rupt%T2E(2,iface)
        itri2=rupt%T2E(4,iface)
        n1=mesh%nx(itri1,itet1)
        n2=mesh%ny(itri1,itet1)
        n3=mesh%nz(itri1,itet1)
        mapm(rupt%perm(:,1,iface))=&
            mesh%vmapM((itri1-1)*Nfp+1:itri1*Nfp)+(itet1-1)*pNp
        mapp(rupt%perm(:,2,iface))=&
            mesh%vmapM((itri2-1)*Nfp+1:itri2*Nfp)+(itet2-1)*pNp
        do i=1,pNp
            vmapm(i)=(itet1-1)*pNp+i
            vmapp(i)=(itet2-1)*pNp+i
        enddo
        call rupt_map_asign_E(pNp,E,pEm,vmapm,0d0,1d0,'s')
        call rupt_map_asign_E(pNp,E,pEp,vmapp,0d0,1d0,'s')

        if(.true.)then
            pSig= 0.5d0*matmul(rupt%M(iface)%BLRCm,pEm)&
                 +0.5d0*matmul(rupt%M(iface)%BLRCp,pEp)
            pTf =-0.5d0*matmul(rupt%M(iface)%KLRCm,pEm)&
                 -0.5d0*matmul(rupt%M(iface)%KLRCp,pEp)
        else
            ntYb=matmul(rupt%M(iface)%A21,ntY1)
        endif

        if(.true.)then
            hTf1=      n1*S%xx(mapm)+n2*S%yx(mapm)+n3*S%zx(mapm)
            hTf2=      n1*S%xy(mapm)+n2*S%yy(mapm)+n3*S%zy(mapm)
            hTf3=      n1*S%xz(mapm)+n2*S%yz(mapm)+n3*S%zz(mapm)
            hTf1=(hTf1+n1*S%xx(mapp)+n2*S%yx(mapp)+n3*S%zx(mapp))*0.5d0
            hTf2=(hTf2+n1*S%xy(mapp)+n2*S%yy(mapp)+n3*S%zy(mapp))*0.5d0
            hTf3=(hTf3+n1*S%xz(mapp)+n2*S%yz(mapp)+n3*S%zz(mapp))*0.5d0
        else
            j=(iface-1)*Nfp+1
            call conv_boundary(Nfp,time,&
                rupt%coord%x(j:),rupt%coord%y(j:),rupt%coord%z(j:),&
                n1,n2,n3,hTf1,hTf2,hTf3)
        endif

        hSig=-n1*hTf1-n2*hTf2-n3*hTf3
        hTf1=hTf1+hSig*n1
        hTf2=hTf2+hSig*n2
        hTf3=hTf3+hSig*n3

        rerr(1)=max(rerr(1),maxval(abs(pSig+hSig))) 
        rerr(2)=max(rerr(2),maxval(abs(pTf1+hTf1)))
        rerr(3)=max(rerr(3),maxval(abs(pTf2+hTf2)))
        rerr(4)=max(rerr(4),maxval(abs(pTf3+hTf3)))
    enddo
!    stop
end subroutine rupt_conv_err

subroutine rupt_Newton(pNp,Nfp,mesh,rupt,dt,V,E,dT0)
    integer,intent(in) :: pNp,Nfp
    type(tetmesh_geometry) :: mesh
    type(rupture) :: rupt
    real(kind=rkind),intent(in) :: dt
    type(vector_array) :: V
    type(tensor_array) :: E,dT0
    integer :: iface,i,j,k,l,m,n,itet1,itet2,itri1,itri2,head
    integer :: map(Nfp),mapm(Nfp),mapp(Nfp),&
        vmap(pNp),vmapm(pNp),vmapp(pNp)
    real(kind=rkind) :: nn(3),nt(3,3),a0,b0,L0,v0,f0,Vt,sig,psi,aVt
    real(kind=rkind) :: dV(3*pNp),dE(9*pNp),&
        mgS(Nfp),Psi0(Nfp),pPsi(Nfp),&
        F(Nfp),dFdS(Nfp),dFdN(Nfp),dFdP(Nfp),&
        dGdS(Nfp),dGdP(Nfp),dPdS(Nfp),&
        rFl(Nfp,3),rFS(Nfp),rSl(Nfp),eye(Nfp,Nfp),G0(Nfp),&
        WW(Nfp*3,Nfp*3),Xi(Nfp*3),&
        Rm(Nfp*3,Nfp*9),Rp(Nfp*3,Nfp*9)
    real(kind=rkind),target :: tmp1(Nfp,Nfp),tmp2(Nfp,Nfp)
    real(kind=rkind),pointer :: dtmp1(:),dtmp2(:)
    real(kind=rkind),parameter :: delt(3,3)=&
        (/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/)
    integer :: iter,ipiv(3*Nfp),info,Niter
    integer,parameter :: itermax=60
    real(kind=rkind) :: res(itermax),res1(itermax),&
        smp1(itermax),smp2(itermax),smp3(itermax),&
        smp4(itermax),err,resmax,detH,&
        tmpH(3*Nfp,3*Nfp),tmpG(3*Nfp),mxJ(3*pNp,3*Nfp)
    logical :: flag
    real(kind=rkind) :: ntJac(Nfp,Nfp,3)
    real(kind=rkind) :: dUtmp(Nfp*3),dVtmp(Nfp*3)
    real(kind=rkind) :: hE0(pNp*9)
    real(kind=rkind) :: Etest(Nfp*9),Vtest(Nfp*3),dVtest(pNp*3),&
        dVtest1(Nfp*3)
    real :: xtmp,ytmp,ztmp

    Niter=0
    dtmp1=>tmp1(1:Nfp**2:Nfp+1,1)
    dtmp2=>tmp2(1:Nfp**2:Nfp+1,1)

    rupt%tau%array=-rupt%tau0%array

    do iface=1,rupt%Nhface
    ! face loop
        itet1=rupt%T2E(1,iface);itri1=rupt%T2E(3,iface)
        itet2=rupt%T2E(2,iface);itri2=rupt%T2E(4,iface)
        nn(1)=mesh%nx(itri1,itet1)
        nn(2)=mesh%ny(itri1,itet1)
        nn(3)=mesh%nz(itri1,itet1)
        do i=1,3
            do j=1,3
                nt(i,j)=delt(i,j)-nn(i)*nn(j)
            enddo
        enddo
        head=(iface-1)*Nfp
        do i=1,Nfp
            map(i)=head+i
        enddo
        do i=1,pNp
            vmap(i)=(iface-1)*pNp+i
        enddo
        mapm(rupt%perm(:,1,iface))=&
            mesh%vmapM((itri1-1)*Nfp+1:itri1*Nfp)+(itet1-1)*pNp
        mapp(rupt%perm(:,2,iface))=&
            mesh%vmapM((itri2-1)*Nfp+1:itri2*Nfp)+(itet2-1)*pNp
        do i=1,pNp
            vmapm(i)=(itet1-1)*pNp+i
            vmapp(i)=(itet2-1)*pNp+i
        enddo
!        pA=>rupt%M(iface)%Abb
        do i=1,3
            do j=1,4
            pA(j,i)%p=>rupt%M(iface)%Ab2&
                ((j-1)*Nfp+1:j*Nfp,(i-1)*Nfp+1:i*Nfp)
            enddo
!print*,pA(4,i)%p
        enddo
!stop
        Psi0=rupt%psi(map)

    ! form Z
        call rupt_map_asign_V(pNp,rupt%Um,hUm,vmap,0d0,1d0,'s')
        call rupt_map_asign_V(pNp,rupt%Up,hUp,vmap,0d0,1d0,'s')
        call rupt_map_asign_E(pNp,E,hEm,vmapm,0d0,1d0,'s')
        call rupt_map_asign_E(pNp,E,hEp,vmapp,0d0,1d0,'s')
        call rupt_map_asign_V(pNp,V,hVm,vmapm,0d0,1d0,'s')
        call rupt_map_asign_V(pNp,V,hVp,vmapp,0d0,1d0,'s')

!    if(rupt_gamma.gt.0d0)then
!        call rupt_map_asign_E(pNp,rupt%Em,hE0,vmap,0d0,1d0,'s')
!        hVm=hVm+matmul(rupt%M(iface)%irDDm,hE0)*rupt_gamma
!!        hVm=hVm+matmul(rupt%M(iface)%BVmEm,hE0)
!!        hVp=hVp+matmul(rupt%M(iface)%BVpEm,hE0)
!        call rupt_map_asign_E(pNp,rupt%Ep,hE0,vmap,0d0,1d0,'s')
!        hVp=hVp+matmul(rupt%M(iface)%irDDp,hE0)*rupt_gamma
!!        hVm=hVm+matmul(rupt%M(iface)%BVmEp,hE0)
!!        hVp=hVp+matmul(rupt%M(iface)%BVpEp,hE0)
!    endif

        call rupt_map_asign_E(Nfp,dT0,ntY1,mapm,0d0,1d0,'s')
        hTf =       nn(1)*ntY1(      1:3*Nfp)
        hTf = hTf + nn(2)*ntY1(3*Nfp+1:6*Nfp)
        hTf = hTf + nn(3)*ntY1(6*Nfp+1:9*Nfp)
        hSig=-hTf1*nn(1)-hTf2*nn(2)-hTf3*nn(3)
        hTf1 = hTf1 + hSig*nn(1)
        hTf2 = hTf2 + hSig*nn(2)
        hTf3 = hTf3 + hSig*nn(3)

!! iter sigma
!        hSig0=hSig
!        hSig=rupt%sigma(map)

    ! form Zbar
        ntZ2=ntZb
        ntY1=ntZ1
        call dgetrs('N',30*pNp,1,rupt%M(iface)%A11,&
            30*pNp,rupt%M(iface)%pA11,ntY1,30*pNp,info)
        ntZb=ntZb-matmul(rupt%M(iface)%A21,ntY1)

    ! form initial guess of S
        call rupt_map_asign_V(Nfp,rupt%Vt,ntS,map,0d0,1d0,'s')
        call rupt_map_asign_V(Nfp,rupt%Vt0,ntS0,map,0d0,1d0,'s')

!        mgS=0d0
!        do m=1,3
!            mgS=mgS+pS(m)%p*nn(m)
!        enddo
!        do m=1,3
!            pS(m)%p=pS(m)%p-mgS*nn(m)
!        enddo

!        mgS=sqrt(pS(1)%p**2+pS(2)%p**2+pS(3)%p**2)

    ! start Newton iteration
        res=0d0;resmax=0d0;res1=0d0
        do iter=1,itermax

        ! compute Ybar
            ntYb=ntZb

            call dgetrs('N',4*Nfp,1,rupt%M(iface)%Ab1,4*Nfp,&
                rupt%M(iface)%pAb1,ntYb,4*Nfp,info)
            ntYb=ntYb+matmul(rupt%M(iface)%Ab2,ntS-ntS0)

            mgS=sqrt(pS(1)%p**2+pS(2)%p**2+pS(3)%p**2)

            call rupt_map_asign_V(Nfp,rupt%tau0,pTf,map,1d0,1d0,'s')
            pSig=pSig+rupt%sigma0(map)

        ! compute friction derivatives
            do n=1,Nfp
                a0  = rupt%a (head+n)
                b0  = rupt%b (head+n)
                L0  = rupt%L (head+n)
                v0  = rupt%v0(head+n)
                f0  = rupt%f0(head+n)
                Vt  =  mgS(n)
                sig = pSig(n)
                psi = Psi0(n)
                psi     = RSfric_invG(sig,Vt,psi,a0,b0,L0,v0,f0,dt)
                F   (n) = RSfric_F   (sig,Vt,psi,a0,b0,L0,v0,f0)
                dFdN(n) = RSfric_dFdN(sig,Vt,psi,a0,b0,L0,v0,f0)
                dFdS(n) = RSfric_dFdS(sig,Vt,psi,a0,b0,L0,v0,f0)
                dFdP(n) = RSfric_dFdP(sig,Vt,psi,a0,b0,L0,v0,f0)
                dGdS(n) = RSfric_dGdS(sig,Vt,psi,a0,b0,L0,v0,f0)
                dGdP(n) = RSfric_dGdP(sig,Vt,psi,a0,b0,L0,v0,f0)
                pPsi(n) = psi
            enddo

            ! compute residual
            do l=1,3
                rFl(:,l)=F*pS(l)%p/mgS-pT(l)%p
            enddo
            res(iter)=res(iter)+sum(rFl**2)
            res(iter)=sqrt(res(iter))
            if(res(iter).lt.1d-10)exit

            dPdS=-dt*dGdS/(dt*dGdP+1d0)
            rFS=(dFdS+dFdP*dPdS)*mgS-F
!            rFS=(dFdS)*mgS-F

            ! form Gradient and Hessian
            ntG=0d0;ntH=0d0
            do l=1,3
                tmp2=0d0
                dtmp2=pS(l)%p/mgS*dFdN
                do j=1,3
                    tmp1=0d0
                    dtmp1=delt(l,j)*F/mgS+&
                        pS(l)%p*pS(j)%p*rFS/(mgS**3)
                    tmp1=tmp1+matmul(tmp2,pA(4,j)%p)-pA(l,j)%p
                    ntJac(:,:,j)=tmp1
                    pG(j)%p=pG(j)%p+matmul(transpose(tmp1),rFl(:,l))
                enddo
                do j=1,3
                    do i=1,3
                        pH(i,j)%p=pH(i,j)%p+matmul( &
                            transpose(ntJac(:,:,j)),ntJac(:,:,i) )
                    enddo
                enddo
            enddo

!            smp1(iter)=minval(f)
!!            smp2(iter)=minval(mgS)
!            smp2(iter)=minval(pPsi)
!            smp3(iter)=maxval(mgS)
!            res(iter)=sqrt(res(iter))
!            res1(iter)=sum(ntG**2)
!            resmax=max(resmax,res1(iter))
            if(res(iter).lt.1d-10)exit
            call dgesv(Nfp*3,1,ntH,Nfp*3,ipiv,ntG,Nfp*3,info)
            ntS=ntS-ntG
!            smp4(iter)=minval(Psig)
        enddo

        if(iter.ge.itermax)then
            xtmp=sum(rupt%coord%x(head+1:head+Nfp))/Nfp
            ytmp=sum(rupt%coord%y(head+1:head+Nfp))/Nfp
            ztmp=sum(rupt%coord%z(head+1:head+Nfp))/Nfp
            print*,''
            print*,'(',xtmp,',',ytmp,',',ztmp,')'
            print*,'obj'
            print*,res
!            print*,'grad'
!            print*,res1
!            print*,'fric_min'
!            print*,smp1
!!            print*,'slip_min'
!            print*,'psi_min'
!            print*,smp2
!            print*,'slip_max'
!            print*,smp3
!            print*,'sigma_min'
!            print*,smp4
            stop
        endif
        Niter=max(Niter,iter)

    ! solve Y1
        ntY1=ntZ1
!        ntY1=ntY1-matmul(rupt%M(iface)%A12,ntYb)
        ntY1=ntY1-matmul(rupt%M(iface)%A12(:,Nfp*4+1:Nfp*7),ntS-ntS0)
        call dgetrs('N',30*pNp,1,rupt%M(iface)%A11,&
            30*pNp,rupt%M(iface)%pA11,ntY1,30*pNp,info)

!if(debugging)then
!if(nn(1).lt.-1e-3.and.myrank.lt.8)then
!    print*,'ntY1'
!    print*,ntY1
!    print*,'ntYb'
!    print*,ntYb
!    print*,'ntS'
!    print*,ntS-ntS0
!    stop
!endif
!endif
    ! update wavefield
        call rupt_map_asign_V(pNp,rupt%Um,pUm,vmap,0d0,1d0,'v')
        call rupt_map_asign_V(pNp,rupt%Up,pUp,vmap,0d0,1d0,'v')
        call rupt_map_asign_E(pNp,E,pEm,vmapm,0d0,1d0,'v')
        call rupt_map_asign_E(pNp,E,pEp,vmapp,0d0,1d0,'v')
        call rupt_map_asign_V(pNp,V,pVm,vmapm,0d0,1d0,'v')
        call rupt_map_asign_V(pNp,V,pVp,vmapp,0d0,1d0,'v')

        call rupt_map_asign_E(pNp,rupt%Em,pEm,vmap,0d0,1d0,'v')
        call rupt_map_asign_E(pNp,rupt%Ep,pEp,vmap,0d0,1d0,'v')

! iter sigma
        rupt%sigma(map) = pSig
!        rupt%sigma(map)=hSig0&
!            -0.5d0*matmul(rupt%M(iface)%BLRCm,pEm)&
!            -0.5d0*matmul(rupt%M(iface)%BLRCp,pEp)

        call rupt_map_asign_V(Nfp,rupt%tauf,pTf,map,0d0,1d0,'v')
        call rupt_map_asign_V(Nfp,rupt%Vt,ntS,map,0d0,1d0,'v')

        rupt%dVt(map)=mgS
        rupt%f(map)=F/pSig
        rupt%psi(map)=pPsi

!        call rupt_map_asign_V(Nfp,rupt%tau,pTf,map,0d0,1d0,'v')

    enddo

    ! end face loop
!if(debugging)then
!call MPI_REDUCE(Niter,iter,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,n)
!if(myrank.eq.0)print*,iter
!endif

end subroutine rupt_Newton

subroutine rupt_map_asign_V(Nfp,V,fV,map,a,b,dist)
    integer,intent(in) :: Nfp,map(Nfp)
    real(kind=rkind),intent(in) :: a,b
    type(vector_array) :: V
    real(kind=rkind) :: fV(Nfp*3)
    character :: dist
    if(dist.eq.'S'.or.dist.eq.'s')then
        if(a.eq.0d0)then
            fV(      1:  Nfp)=b*V%x(map)
            fV(  Nfp+1:2*Nfp)=b*V%y(map)
            fV(2*Nfp+1:3*Nfp)=b*V%z(map)
        else
            fV(      1:  Nfp)=a*fV(      1:  Nfp)+b*V%x(map)
            fV(  Nfp+1:2*Nfp)=a*fV(  Nfp+1:2*Nfp)+b*V%y(map)
            fV(2*Nfp+1:3*Nfp)=a*fV(2*Nfp+1:3*Nfp)+b*V%z(map)
        endif
    else
        if(a.eq.0d0)then
            V%x(map)=b*fV(      1:  Nfp)
            V%y(map)=b*fV(  Nfp+1:2*Nfp)
            V%z(map)=b*fV(2*Nfp+1:3*Nfp)
        else
            V%x(map)=a*V%x(map)+b*fV(      1:  Nfp)
            V%y(map)=a*V%y(map)+b*fV(  Nfp+1:2*Nfp)
            V%z(map)=a*V%z(map)+b*fV(2*Nfp+1:3*Nfp)
        endif
    endif
end subroutine rupt_map_asign_V

subroutine rupt_map_asign_E(Nfp,E,fE,map,a,b,dist)
    integer,intent(in) :: Nfp,map(Nfp)
    real(kind=rkind),intent(in) :: a,b
    type(tensor_array) :: E
    real(kind=rkind) :: fE(Nfp*9)
    character :: dist
    if(dist.eq.'S'.or.dist.eq.'s')then
        if(a.eq.0d0)then
            fE(      1:  Nfp)=b*E%xx(map)
            fE(  Nfp+1:2*Nfp)=b*E%yx(map)
            fE(2*Nfp+1:3*Nfp)=b*E%zx(map)
            fE(3*Nfp+1:4*Nfp)=b*E%xy(map)
            fE(4*Nfp+1:5*Nfp)=b*E%yy(map)
            fE(5*Nfp+1:6*Nfp)=b*E%zy(map)
            fE(6*Nfp+1:7*Nfp)=b*E%xz(map)
            fE(7*Nfp+1:8*Nfp)=b*E%yz(map)
            fE(8*Nfp+1:9*Nfp)=b*E%zz(map)
        else
            fE(      1:  Nfp)=a*fE(      1:  Nfp)+b*E%xx(map)
            fE(  Nfp+1:2*Nfp)=a*fE(  Nfp+1:2*Nfp)+b*E%yx(map)
            fE(2*Nfp+1:3*Nfp)=a*fE(2*Nfp+1:3*Nfp)+b*E%zx(map)
            fE(3*Nfp+1:4*Nfp)=a*fE(3*Nfp+1:4*Nfp)+b*E%xy(map)
            fE(4*Nfp+1:5*Nfp)=a*fE(4*Nfp+1:5*Nfp)+b*E%yy(map)
            fE(5*Nfp+1:6*Nfp)=a*fE(5*Nfp+1:6*Nfp)+b*E%zy(map)
            fE(6*Nfp+1:7*Nfp)=a*fE(6*Nfp+1:7*Nfp)+b*E%xz(map)
            fE(7*Nfp+1:8*Nfp)=a*fE(7*Nfp+1:8*Nfp)+b*E%yz(map)
            fE(8*Nfp+1:9*Nfp)=a*fE(8*Nfp+1:9*Nfp)+b*E%zz(map)
        endif
    else
        if(a.eq.0d0)then
            E%xx(map)=b*fE(      1:  Nfp)
            E%yx(map)=b*fE(  Nfp+1:2*Nfp)
            E%zx(map)=b*fE(2*Nfp+1:3*Nfp)
            E%xy(map)=b*fE(3*Nfp+1:4*Nfp)
            E%yy(map)=b*fE(4*Nfp+1:5*Nfp)
            E%zy(map)=b*fE(5*Nfp+1:6*Nfp)
            E%xz(map)=b*fE(6*Nfp+1:7*Nfp)
            E%yz(map)=b*fE(7*Nfp+1:8*Nfp)
            E%zz(map)=b*fE(8*Nfp+1:9*Nfp)
        else
            E%xx(map)=a*E%xx(map)+b*fE(      1:  Nfp)
            E%yx(map)=a*E%yx(map)+b*fE(  Nfp+1:2*Nfp)
            E%zx(map)=a*E%zx(map)+b*fE(2*Nfp+1:3*Nfp)
            E%xy(map)=a*E%xy(map)+b*fE(3*Nfp+1:4*Nfp)
            E%yy(map)=a*E%yy(map)+b*fE(4*Nfp+1:5*Nfp)
            E%zy(map)=a*E%zy(map)+b*fE(5*Nfp+1:6*Nfp)
            E%xz(map)=a*E%xz(map)+b*fE(6*Nfp+1:7*Nfp)
            E%yz(map)=a*E%yz(map)+b*fE(7*Nfp+1:8*Nfp)
            E%zz(map)=a*E%zz(map)+b*fE(8*Nfp+1:9*Nfp)
        endif
    endif
end subroutine rupt_map_asign_E

subroutine rupt_update_psi_RK4(Ndim,rupt,psi0,psi1,V0,V1,sig0,sig1,dt)
    integer,intent(in) :: Ndim
    type(rupture) :: rupt
    real(kind=rkind) :: dt,psi0(Ndim),psi1(Ndim),&
        sig0(Ndim),sig1(Ndim),V0(Ndim),V1(Ndim)
    real(kind=rkind) :: a,b,L,s0,f0
    real(kind=rkind) :: psi,rsi,hsi,dv,v,vt,&
                        sig,sigt,dsig
    integer :: i,j,k
    do i=1,Ndim
        a =rupt%a (i)
        b =rupt%b (i)
        L =rupt%L (i)
        s0=rupt%V0(i)
        f0=rupt%f0(i)
        v=V0(i);dv=V1(i)-v
        psi=psi0(i)
        sig=sig0(i);dsig=sig1(i)-sig
        rsi=0d0
        do k=1,5
            vt=v+rk4c(k)*dv
            sigt=sig+rk4c(k)*dsig
            hsi=-RSfric_G(sigt,vt,psi,a,b,L,s0,f0)
            rsi=rk4a(k)*rsi+dt*hsi
            psi=rk4b(k)*rsi+psi
        enddo
        psi1(i)=psi
    enddo
end subroutine rupt_update_psi_RK4

end module ruptsolve_mod

