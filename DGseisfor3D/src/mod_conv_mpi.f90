!*******************************************************************!
!*  This module aims for convergence test                          *!
!*                                                                 *!
!*  Author: Ruichao Ye    Email: Ruichao.ye@gmail.com              *!
!*******************************************************************!

!--------------------------------------------------------------------
!--------------------------------------------------------------------
module conv_mpi_mod
!--------------------------------------------------------------------
    use para_mod,     only : refer_Cp, refer_Cs, refer_rho, pNp, Nfp
    use datatype_mod, only : rkind,receivers,surface,&
                             tetmesh_geometry,tetmesh_material,&
                             vector_array,tensor_array
    use Meshfile_mod, only : zmax,zmin
    implicit none
    public :: conv_boundary,conv_solution,conv_material,conv_err
!--------------------------------------------------------------------
    double precision,parameter :: PI=3.1415926535897932384626d0
    character(len=10) :: wavetype
!--------------------------------------------------------------------
        contains
!--------------------------------------------------------------------

subroutine init_convrecv(domainid,stime,etime,snapintv,recdt,recvs)
    integer :: domainid
    type(receivers) :: recvs
    real(kind=rkind) :: stime,etime,snapintv,recdt
    integer :: i,Nrecord
    if(domainid.ne.1)then
        recvs%Nrecv=0;return
    endif
    recvs%Nrecv=2
    recvs%Nrecord=0
    recvs%Ncomp=3
    Nrecord=ceiling(snapintv/recdt)
    allocate(recvs%globID(recvs%Nrecv))
    allocate(recvs%rec_buffer(3,Nrecord+2,recvs%Nrecv))
    do i=1,recvs%Nrecv
        recvs%globID(i)=i
    enddo
    allocate(recvs%switchoff(recvs%Nrecv))
    recvs%switchoff=.false.
end subroutine init_convrecv

subroutine conv_recording(recvs,l2er,linf)
    type(receivers) :: recvs
    real(kind=rkind) :: l2er(9),linf(9)
    if(recvs%Nrecv.lt.2)return
    recvs%Nrecord=recvs%Nrecord+1
    recvs%rec_buffer(1:3,recvs%Nrecord,1)=l2er(1:3)
    recvs%rec_buffer(1:3,recvs%Nrecord,2)=linf(1:3)
end subroutine conv_recording

subroutine conv_surf(Nfp,surf,time)
    integer :: Nfp
    real(kind=rkind) :: time
    type(surface) :: surf
    integer :: head,i
    do i=1,surf%Nface
        head=(i-1)*Nfp+1
        call conv_boundary(Nfp,time,&
            surf%coord%x(head:),&
            surf%coord%y(head:),&
            surf%coord%z(head:),&
            surf%nx(i),surf%ny(i),surf%nz(i),&
            surf%Sn%x(head:),surf%Sn%y(head:),surf%Sn%z(head:))
    enddo
end subroutine conv_surf

subroutine conv_body(pNp,mesh,V,E,time)
    integer :: pNp
    real(kind=rkind) :: time
    type(tetmesh_geometry) :: mesh
    type(vector_array) :: V
    type(tensor_array) :: E
    integer :: head,i
    do i=1,mesh%Nele
        head=(i-1)*pNp+1
        call conv_solution(pNp,time,&
            mesh%coord%x(head:),&
            mesh%coord%y(head:),&
            mesh%coord%z(head:),&
             V%x(head:), V%y(head:), V%z(head:),&
            E%xx(head:),E%yy(head:),E%zz(head:),&
            E%yz(head:),E%xz(head:),E%xy(head:))
    enddo
end subroutine conv_body

subroutine conv_model(pNp,mesh,mat)
    integer :: pNp
    type(tetmesh_geometry) :: mesh
    type(tetmesh_material) :: mat
    integer :: head,i,j,k
    real(kind=rkind),pointer :: lambda(:),muX2(:)
    lambda=>mat%C(:,1);muX2=>mat%C(:,2)
    do i=1,mesh%Nele
        head=(i-1)*pNp+1
        call conv_material(pNp,&
            mesh%coord%x(head:),&
            mesh%coord%y(head:),&
            mesh%coord%z(head:),&
            lambda(head:),muX2(head:),mat%rho(head:))
        if(mat%C(head,2).le.1e-6)then
            mat%k_media(i)=0
        else
            mat%k_media(i)=1
        endif
    enddo
    do i=1,mesh%Nele
        do j=1,4
            if(mesh%fbctype(1,j,i).eq.0)then
                k=(i-1)*Nfp*4+(j-1)*Nfp+1
                k=(mesh%vmapP(k)-1)/pNp+1 !neigheleID
                if(mat%k_media(i).eq.mat%k_media(k))then
                    mesh%fbctype(2,j,i)=0
                elseif(mat%k_media(i).gt.0 .and. &
                       mat%k_media(k).eq.0)then
                    mesh%fbctype(2,j,i)=1
                elseif(mat%k_media(i).eq.0 .and. &
                       mat%k_media(k).gt.0)then
                    mesh%fbctype(2,j,i)=2
                endif
            endif
        enddo
    enddo
end subroutine conv_model

subroutine conv_boundary(Nfp,time,xm,ym,zm,nx,ny,nz,SBn1,SBn2,SBn3)
    integer,intent(in) :: Nfp
    real(kind=rkind),intent(in) :: xm(Nfp),ym(Nfp),zm(Nfp),time
    real(kind=rkind),intent(in) :: nx,ny,nz
    real(kind=rkind),intent(out):: SBn1(Nfp),SBn2(Nfp),SBn3(Nfp)
    real(kind=rkind) :: S11B(Nfp),S22B(Nfp),S33B(Nfp),&
                        S23B(Nfp),S13B(Nfp),S12B(Nfp)
    real(kind=rkind) :: E11B(Nfp),E22B(Nfp),E33B(Nfp),&
                        E23B(Nfp),E13B(Nfp),E12B(Nfp),&
                         V1B(Nfp), V2B(Nfp), V3B(Nfp)
    real(kind=rkind) :: lambdaB(1),muB(1),Cp(1),Cs(1),rho(1)

    if(trim(wavetype).eq.'Plane')then
! for plane wave
        call planewave(&
                V1B,V2B,V3B,E11B,E22B,E33B,E23B,E13B,E12B,&
                xm,ym,zm,Nfp,time)
!        lambdaB=(refer_Cp**2-2d0*refer_Cs**2)*refer_rho
!        muB=refer_Cs**2*refer_rho
        lambdaB=2d0;muB=1d0
    elseif(trim(wavetype).eq.'Rayleigh')then
! for Rayleigh wave
        call Rayleigh(&
                V1B,V2B,V3B,E11B,E22B,E33B,E23B,E13B,E12B,&
                xm,ym,zm,Nfp,time)
        lambdaB=2d0;muB=1d0
    elseif(trim(wavetype).eq.'Scholte')then
! for Scholte wave
      if(sum(zm).gt.0)then
        call Scholtemat(cp,cs,rho,1,0) 
        call Scholte(&
                V1B,V2B,V3B,E11B,E22B,E33B,E23B,E13B,E12B,&
                xm,ym,zm,Nfp,time,1)
      else
        call Scholtemat(cp,cs,rho,1,1) 
        call Scholte(&
                V1B,V2B,V3B,E11B,E22B,E33B,E23B,E13B,E12B,&
                xm,ym,zm,Nfp,time,0)
      endif
        muB=cs**2*rho
        lambdaB=cp**2*rho-2d0*muB
    else
! for Stoneley wave
      if(sum(zm).gt.0)then
        call Stoneleymat(cp,cs,rho,1,0) 
        call Stoneley(&
                V1B,V2B,V3B,E11B,E22B,E33B,E23B,E13B,E12B,&
                xm,ym,zm,Nfp,time,1)
      else
        call Stoneleymat(cp,cs,rho,1,1) 
        call Stoneley(&
                V1B,V2B,V3B,E11B,E22B,E33B,E23B,E13B,E12B,&
                xm,ym,zm,Nfp,time,0)
      endif
        muB=cs**2*rho
        lambdaB=cp**2*rho-2d0*muB
    endif
    S33B=lambdaB(1)*(E11B+E22B+E33B)
    S11B=S33B+2d0*muB(1)*E11B
    S22B=S33B+2d0*muB(1)*E22B
    S33B=S33B+2d0*muB(1)*E33B
    S23B=2d0*muB(1)*E23B
    S13B=2d0*muB(1)*E13B
    S12B=2d0*muB(1)*E12B
    SBn1=S11B*nx+S12B*ny+S13B*nz
    SBn2=S12B*nx+S22B*ny+S23B*nz
    SBn3=S13B*nx+S23B*ny+S33B*nz
end subroutine conv_boundary

subroutine conv_solution(pNp,time,xm,ym,zm,exV1,exV2,exV3,&
                exE11,exE22,exE33,exE23,exE13,exE12)
    integer,intent(in) :: pNp
    real(kind=rkind),intent(in) :: xm(pNp),ym(pNp),zm(pNp),time
    real(kind=rkind),intent(out):: exE11(pNp),exE22(pNp),exE33(pNp),&
                                   exE23(pNp),exE13(pNp),exE12(pNp),&
                                    exV1(pNp), exV2(pNp), exV3(pNp)

    if(trim(wavetype).eq.'Plane')then
! for plane wave
        call planewave(&
                exV1,exV2,exV3,exE11,exE22,exE33,exE23,exE13,exE12,&
                xm,ym,zm,pNp,time)
    elseif(trim(wavetype).eq.'Rayleigh')then
! for Rayleigh wave
        call Rayleigh(&
                exV1,exV2,exV3,exE11,exE22,exE33,exE23,exE13,exE12,&
                xm,ym,zm,pNp,time)
    elseif(trim(wavetype).eq.'Scholte')then
! for Scholte wave
      if(sum(zm).gt.0)then
        call Scholte(&
                exV1,exV2,exV3,exE11,exE22,exE33,exE23,exE13,exE12,&
                xm,ym,zm,pNp,time,1)
      else
        call Scholte(&
                exV1,exV2,exV3,exE11,exE22,exE33,exE23,exE13,exE12,&
                xm,ym,zm,pNp,time,0)
      endif
    else
! for Stoneley wave
      if(sum(zm).gt.0)then
        call Stoneley(&
                exV1,exV2,exV3,exE11,exE22,exE33,exE23,exE13,exE12,&
                xm,ym,zm,pNp,time,1)
      else
        call Stoneley(&
                exV1,exV2,exV3,exE11,exE22,exE33,exE23,exE13,exE12,&
                xm,ym,zm,pNp,time,0)
      endif
    endif
end subroutine conv_solution

subroutine conv_material(pNp,xm,ym,zm,lambda,muX2,rho)
    integer,intent(in) :: pNp
    real(kind=rkind),intent(in) :: xm(pNp),ym(pNp),zm(pNp)
    real(kind=rkind) :: lambda(pNp),muX2(pNp),rho(pNp)
    real(kind=rkind) :: Cp(pNp),Cs(pNp)

    if(trim(wavetype).eq.'Plane')then
! for plane wave
        lambda=(refer_Cp**2-2d0*refer_Cs**2)*refer_rho
        muX2=refer_Cs**2*refer_rho*2d0
        rho=refer_rho
    elseif(trim(wavetype).eq.'Rayleigh')then
! for Rayleigh wave
        lambda=2d0;muX2=2d0;rho=1d0
    elseif(trim(wavetype).eq.'Scholte')then
! for Scholte wave
      if(sum(zm).gt.0)then
        call Scholtemat(Cp,Cs,rho,pNp,0) 
      else
        call Scholtemat(Cp,Cs,rho,pNp,1) 
      endif
        muX2=2d0*Cs**2*rho
        lambda=Cp**2*rho-muX2
    else
! for Stoneley wave
      if(sum(zm).gt.0)then
        call Stoneleymat(Cp,Cs,rho,pNp,0) 
      else
        call Stoneleymat(Cp,Cs,rho,pNp,1) 
      endif
        muX2=2d0*Cs**2*rho
        lambda=Cp**2*rho-muX2
    endif
end subroutine conv_material

subroutine conv_err(Nele,detJ,time,Mass,x,y,z,&
        v1,v2,v3,e11,e22,e33,e23,e13,e12,&
        l2err2,linfty)
    integer :: Nele
    real(kind=rkind) :: detJ(Nele),time,Mass(pNp,pNp)
    real(kind=rkind) ::   x(pNp*Nele),  y(pNp*Nele),  z(pNp*Nele)
    real(kind=rkind) ::  V1(pNp*Nele), V2(pNp*Nele), V3(pNp*Nele)
    real(kind=rkind) :: E11(pNp*Nele),E22(pNp*Nele),E33(pNp*Nele)
    real(kind=rkind) :: E23(pNp*Nele),E13(pNp*Nele),E12(pNp*Nele)
    real(kind=rkind) ::  exV1(pNp), exV2(pNp), exV3(pNp)
    real(kind=rkind) :: exE11(pNp),exE22(pNp),exE33(pNp)
    real(kind=rkind) :: exE23(pNp),exE13(pNp),exE12(pNp)
    real(kind=rkind) :: l2err2(9),linfty(9)
    integer :: j,head,tail

    l2err2=0d0
    linfty=0d0
    do j=1,Nele
        head=(j-1)*pNp+1;tail=j*pNp
        call conv_solution(pNp,time,x(head:),y(head:),z(head:),&
            exV1,exV2,exV3,exE11,exE22,exE33,exE23,exE13,exE12)
        exV1 = V1(head:tail)- exV1
        exV2 = V2(head:tail)- exV2
        exV3 = V3(head:tail)- exV3
        exE11=E11(head:tail)-exE11
        exE22=E22(head:tail)-exE22
        exE33=E33(head:tail)-exE33
        exE23=E23(head:tail)-exE23
        exE13=E13(head:tail)-exE13
        exE12=E12(head:tail)-exE12
        l2err2(1)=l2err2(1)+dot_product( exV1,matmul(Mass, exV1))*detJ(j)
        l2err2(2)=l2err2(2)+dot_product( exV2,matmul(Mass, exV2))*detJ(j)
        l2err2(3)=l2err2(3)+dot_product( exV3,matmul(Mass, exV3))*detJ(j)
        l2err2(4)=l2err2(4)+dot_product(exE11,matmul(Mass,exE11))*detJ(j)
        l2err2(5)=l2err2(5)+dot_product(exE22,matmul(Mass,exE22))*detJ(j)
        l2err2(6)=l2err2(6)+dot_product(exE33,matmul(Mass,exE33))*detJ(j)
        l2err2(7)=l2err2(7)+dot_product(exE23,matmul(Mass,exE23))*detJ(j)
        l2err2(8)=l2err2(8)+dot_product(exE13,matmul(Mass,exE13))*detJ(j)
        l2err2(9)=l2err2(9)+dot_product(exE12,matmul(Mass,exE12))*detJ(j)
        linfty(1)=max(linfty(1),maxval(abs( exV1)))
        linfty(2)=max(linfty(2),maxval(abs( exV2)))
        linfty(3)=max(linfty(3),maxval(abs( exV3)))
        linfty(4)=max(linfty(4),maxval(abs(exE11)))
        linfty(5)=max(linfty(5),maxval(abs(exE22)))
        linfty(6)=max(linfty(6),maxval(abs(exE33)))
        linfty(7)=max(linfty(7),maxval(abs(exE23)))
        linfty(8)=max(linfty(8),maxval(abs(exE13)))
        linfty(9)=max(linfty(9),maxval(abs(exE12)))
    enddo

end subroutine conv_err

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!  Analytical solutions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine planewave(V1,V2,V3,E11,E22,E33,E23,E13,E12,x,y,z,&
        Ndof,time)
    integer :: Ndof,i
    real(kind=rkind) ::   x(Ndof),  y(Ndof),  z(Ndof)
    real(kind=rkind) ::  V1(Ndof), V2(Ndof), V3(Ndof)
    real(kind=rkind) :: E11(Ndof),E22(Ndof),E33(Ndof)
    real(kind=rkind) :: E23(Ndof),E13(Ndof),E12(Ndof)
    real(kind=rkind) :: time,cp,cs
    real(kind=rkind) :: sinocptkx,sinocstkx
    real(kind=rkind) :: khat1,khat2,khat3,lhat1,lhat2,lhat3
    real(kind=rkind),parameter :: knorm=0.5d0*pi

    khat1 =  0.0d0;khat2 =  0.0d0;khat3 =  1.0d0
    lhat1 =  0.0d0;lhat2 = -1.0d0;lhat3 =  0.0d0
!    cp=refer_Cp;cs=refer_Cs
    cp=2d0;cs=1d0

   do i=1,Ndof         
       sinocptkx = cos( knorm*cp*time+&
               knorm*(khat1*x(i)+khat2*y(i)+khat3*z(i)));
       sinocstkx = cos(-knorm*cs*time+&
               knorm*(khat1*x(i)+khat2*y(i)+khat3*z(i)));
       V1(i)=khat1*cp*sinocptkx-lhat1*cs*sinocstkx;
       V2(i)=khat2*cp*sinocptkx-lhat2*cs*sinocstkx;
       V3(i)=khat3*cp*sinocptkx-lhat3*cs*sinocstkx;
       E11(i)=khat1**2*sinocptkx+khat1*lhat1*sinocstkx;
       E22(i)=khat2**2*sinocptkx+khat2*lhat2*sinocstkx;
       E33(i)=khat3**2*sinocptkx+khat3*lhat3*sinocstkx;
       E12(i)=khat1*khat2*sinocptkx&
               +0.5d0*(khat1*lhat2+lhat1*khat2)*sinocstkx;
       E23(i)=khat2*khat3*sinocptkx&
               +0.5d0*(khat2*lhat3+lhat2*khat3)*sinocstkx;
       E13(i)=khat1*khat3*sinocptkx&
               +0.5d0*(khat1*lhat3+lhat1*khat3)*sinocstkx;
    enddo
end subroutine planewave

subroutine Rayleigh(V1,V2,V3,E11,E22,E33,E23,E13,E12,x,y,z,Ndof,time)
    integer :: Ndof,i
    real(kind=rkind) :: x(:),y(:),z(:),V1(:),V2(:),V3(:)
    real(kind=rkind) :: E11(:),E22(:),E33(:),E23(:),E13(:),E12(:)
    real(kind=rkind) :: time
    real(kind=rkind),parameter:: cp=2.0d0,cs=1.0d0,rho=1.0d0,&
            c=9.325259059311549d-01,k=1.0d0*PI,omega=k*c,&
            bl=sqrt(1d0-c**2/cp**2),bs=sqrt(1d0-c**2/cs**2)
    real(kind=rkind) :: ecosl,esinl,ecoss,esins
    real(kind=rkind),parameter :: A=1.0d0,&
            B=-2.0d0*cs**2/(2*cs**2-c**2)*bl*A 
    real(kind=rkind) :: d2phidxdt,d2phidzdt,d2phidxdx,&
                        d2phidzdz,d2phidxdz,d2psidxdt,&
                        d2psidzdt,d2psidxdx,d2psidzdz,&
                        d2psidxdz

    do i=1,Ndof

        ecosl = exp(-k*bl*(z(i)-zmin))*cos(-k*x(i)+omega*time)
        esinl = exp(-k*bl*(z(i)-zmin))*sin(-k*x(i)+omega*time)
        ecoss = exp(-k*bs*(z(i)-zmin))*cos(-k*x(i)+omega*time)
        esins = exp(-k*bs*(z(i)-zmin))*sin(-k*x(i)+omega*time)

        d2phidxdt = k*omega*   ( A*ecosl)
        d2phidzdt = k*omega*bl*( A*esinl)
        d2phidxdx = k*k*       (-A*ecosl)
        d2phidzdz = k*k*bl*bl* ( A*ecosl)
        d2phidxdz = k*k*bl*    (-A*esinl)
        d2psidxdt = k*omega*   ( B*ecoss)
        d2psidzdt = k*omega*bs*( B*esins)
        d2psidxdx = k*k*       (-B*ecoss)
        d2psidzdz = k*k*bs*bs* ( B*ecoss)
        d2psidxdz = k*k*bs*    (-B*esins)

        V1(i)  = d2phidxdt - d2psidzdt
        V2(i)  = 0.0d0
        V3(i)  = d2phidzdt + d2psidxdt
        E11(i) = d2phidxdx - d2psidxdz
        E22(i) = 0.0d0
        E33(i) = d2phidzdz + d2psidxdz
        E12(i) = 0.0d0
        E23(i) = 0.0d0
        E13(i) = d2phidxdz + 0.5d0 * (d2psidxdx - d2psidzdz)
    enddo

end subroutine Rayleigh

subroutine Scholtemat(cp,cs,rho,Ndof,k_media)
    integer :: Ndof,k_media
    real(kind=rkind) :: cp(Ndof),cs(Ndof),rho(Ndof)
    real(kind=rkind),parameter:: &
        lambda1=1.11d0,mu1=0.0d0,rho1=1.32d0,&
        lambda2=1.2d0,mu2=1.3d0,rho2=1.1d0,&
        cp1=sqrt((lambda1+2*mu1)/rho1),&
        cs1=sqrt(mu1/rho1),&
        cp2=sqrt((lambda2+2*mu2)/rho2),&
        cs2=sqrt(mu2/rho2)

    if(k_media.eq.0)then
            cp(1:)=cp1; cs(1:)=cs1; rho(1:)=rho1
    elseif(k_media.eq.1)then
            cp(1:)=cp2; cs(1:)=cs2; rho(1:)=rho2
    endif
end subroutine Scholtemat

subroutine Scholte(V1,V2,V3,E11,E22,E33,E23,E13,E12,x,y,z,&
        Ndof,time,k_media)
    integer :: Ndof,i,k_media
    real(kind=rkind) :: x(:),y(:),z(:),V1(:),V2(:),V3(:)
    real(kind=rkind) :: E11(:),E22(:),E33(:),E23(:),E13(:),E12(:)
    real(kind=rkind) :: time
    real(kind=rkind),parameter :: &
    lambda1=1.11d0,mu1=0.0d0,rho1=1.32d0,&
    lambda2=1.2d0,mu2=1.3d0,rho2=1.1d0,&
    cp1=sqrt((lambda1+2*mu1)/rho1),&
    cs1=sqrt(mu1/rho1),&
    cp2=sqrt((lambda2+2*mu2)/rho2),&
    cs2=sqrt(mu2/rho2),&
    c=0.712130427093662d0,k=0.25d0*PI,omega=k*c,&
    b1=sqrt(1d0-c**2/cp1**2),bl2=sqrt(1d0-c**2/cp2**2),&
    bs2=sqrt(1d0-c**2/cs2**2)
    real(kind=rkind) :: d2phidxdt,d2phidzdt,d2phidxdx,&
                        d2phidzdz,d2phidxdz,d2psidxdt,&
                        d2psidzdt,d2psidxdx,d2psidzdz,&
                        d2psidxdz
    real(kind=rkind) :: reA1,imA1

    do i=1,Ndof
        if(k_media.le.0)then
            d2phidxdt = k*omega*(1.0d0+bs2**2)/(2.0d0*bl2) &
                    *exp(k*bl2*z(i))*sin(k*x(i)-omega*time)
            d2phidzdt = -k*omega*0.5d0*(1.0d0+bs2**2) &
                    *exp(k*bl2*z(i))*cos(k*x(i)-omega*time)
            d2phidxdx = -k**2*(1.0d0+bs2*bs2)/(2.0d0*bl2) &
                    *exp(k*bl2*z(i))*sin(k*x(i)-omega*time)
            d2phidzdz = k**2*0.5d0*bl2*(1.0d0+bs2**2) &
                    *exp(k*bl2*z(i))*sin(k*x(i)-omega*time)
            d2phidxdz = k**2*0.5d0*(1.0d0+bs2**2) &
                    *exp(k*bl2*z(i))*cos(k*x(i)-omega*time)
            d2psidxdt = omega*k*exp(k*bs2*z(i)) &
                    *cos(k*x(i)-omega*time)
            d2psidzdt = omega*k*bs2*exp(k*bs2*z(i)) &
                    *sin(k*x(i)-omega*time)
            d2psidxdx = -k**2*exp(k*bs2*z(i)) &
                    *cos(k*x(i)-omega*time)
            d2psidzdz = k**2*bs2**2*exp(k*bs2*z(i)) &
                    *cos(k*x(i)-omega*time)
            d2psidxdz = -k**2*bs2*exp(k*bs2*z(i)) &
                    *sin(k*x(i)-omega*time)     
        else
            reA1 = 0.0d0
            imA1 = 0.5d0 * (bs2 * bs2 - 1.0d0) / b1
            d2psidxdt = 0.0d0
            d2psidzdt = 0.0d0
            d2psidxdx = 0.0d0
            d2psidzdz = 0.0d0
            d2psidxdz = 0.0d0
            d2phidxdt = k*omega*reA1*exp(-k*b1*z(i)) &
                    *cos(k*x(i)-omega*time) &
                    - k*omega*imA1*exp(-k*b1*z(i)) &
                    *sin(k*x(i)-omega*time)
            d2phidzdt = -k*omega*b1*imA1*exp(-k*b1*z(i)) &
                    *cos(k*x(i)-omega*time) &
                    + k*omega*b1*reA1*exp(-k*b1*z(i)) &
                    *sin(k*x(i)-omega*time)
            d2phidxdx = -k**2*reA1*exp(-k*b1*z(i)) &
                    *cos(k*x(i)-omega*time) &
                    + k**2*imA1*exp(-k*b1*z(i)) &
                    *sin(k*x(i)-omega*time)
            d2phidzdz = k**2*b1**2*reA1*exp(-k*b1*z(i)) &
                    *cos(k*x(i)-omega*time) &
                    - k**2*b1**2*imA1*exp(-k*b1*z(i)) &
                    *sin(k*x(i)-omega*time)
            d2phidxdz =  k**2*b1*reA1*exp(-k*b1*z(i)) &
                    *sin(k*x(i)-omega*time) &
                    + k**2*b1*imA1*exp(-k*b1*z(i)) &
                    *cos(k*x(i)-omega*time)
        endif
        V1(i)  = d2phidxdt - d2psidzdt
        V2(i)  = 0.0d0
        V3(i)  = d2phidzdt + d2psidxdt
        E11(i) = d2phidxdx - d2psidxdz
        E22(i) = 0.0d0
        E33(i) = d2phidzdz + d2psidxdz
        E12(i) = 0.0d0
        E23(i) = 0.0d0
        E13(i) = d2phidxdz + 0.5d0 * (d2psidxdx - d2psidzdz)
    enddo
end subroutine Scholte

subroutine Stoneleymat(cp,cs,rho,Ndof,k_media)
    integer :: Ndof,k_media
    real(kind=rkind) :: cp(Ndof),cs(Ndof),rho(Ndof)
    real(kind=rkind),parameter::&
        lambda1=3.0d0,mu1=3.0d0,rho1=4.0d0,&
        lambda2=1.2d0,mu2=1.2d0,rho2=1.2d0,&
        cp1=sqrt((lambda1+2d0*mu1)/rho1),&
        cs1=sqrt(mu1/rho1),&
        cp2=sqrt((lambda2+2d0*mu2)/rho2),&
        cs2=sqrt(mu2/rho2)

    if(k_media.eq.0)then
            cp(1:)=cp1; cs(1:)=cs1; rho(1:)=rho1
    elseif(k_media.eq.1)then
            cp(1:)=cp2; cs(1:)=cs2; rho(1:)=rho2
    endif
end subroutine Stoneleymat

subroutine Stoneley(V1,V2,V3,E11,E22,E33,E23,E13,E12,x,y,z,&
        Ndof,time,k_media)
    integer :: Ndof,k_media,i
    real(kind=rkind) :: x(:),y(:),z(:),V1(:),V2(:),V3(:)
    real(kind=rkind) :: E11(:),E22(:),E33(:),E23(:),E13(:),E12(:)
    real(kind=rkind) :: time,bl1,bs1,bl2,bs2
    real(kind=rkind) :: cl1,cs1,cl2,cs2
    real(kind=rkind),parameter:: &
        lambda1=3.0d0,mu1=3.0d0,rho1=4.0d0,&
        lambda2=1.2d0,mu2=1.2d0,rho2=1.2d0,&
        c=8.638082912540790d-01,&
        k=0.25d0*PI,omega=k*c
    real(kind=rkind) :: d2phidxdt,d2phidzdt,d2phidxdx,&
                        d2phidzdz,d2phidxdz,d2psidxdt,&
                        d2psidzdt,d2psidxdx,d2psidzdz,&
                        d2psidxdz,&
                        ecosl,esinl,ecoss,esins
    real(kind=rkind) :: reA1,imA1,reA2,imA2,reB1,imB1,reB2,imB2

    cl1 = sqrt ((lambda1 + 2d0 * mu1) / rho1)
    cs1 = sqrt (mu1 / rho1)
    cl2 = sqrt ((lambda2 + 2d0 * mu2) / rho2)
    cs2 = sqrt (mu2 / rho2)
    bl1 = sqrt (1d0 - c * c / (cl1 * cl1))
    bs1 = sqrt (1d0 - c * c / (cs1 * cs1))
    bl2 = sqrt (1d0 - c * c / (cl2 * cl2))
    bs2 = sqrt (1d0 - c * c / (cs2 * cs2))

    reB2= 2d0*(mu1*bl1+mu2*bl2)*(bl1*bs1-1d0)/(bl1+bl2) &
         -2d0*mu1*bl1*bs1+mu1*(1d0+bs1*bs1)
    imB2=0d0
    
    reB1=-2d0*(mu1*bl1+mu2*bl2)*(1d0+bl1*bs2)/(bl1+bl2) &
         +2d0*mu1*bl1*bs2+mu2*(1d0+bs2*bs2)
    imB1=0d0
    
    reA2 = (bl1 * bs1 - 1d0) / (bl1 + bl2) * imB1 &
      + (1d0 + bl1 * bs2) / (bl1 + bl2) * imB2 
    imA2 = -(bl1 * bs1 - 1d0) / (bl1 + bl2) * reB1 &
      - (1d0 + bl1 * bs2) / (bl1 + bl2) * reB2
    
    reA1 = reA2 - bs1 * imB1 - bs2 * imB2
    imA1 = imA2 + bs1 * reB1 + bs2 * reB2


    do i=1,Ndof
        if(k_media.le.0)then
            ecosl = exp(k*bl2*z(i))*cos(k*x(i)-omega*time)
            esinl = exp(k*bl2*z(i))*sin(k*x(i)-omega*time)
            ecoss = exp(k*bs2*z(i))*cos(k*x(i)-omega*time)
            esins = exp(k*bs2*z(i))*sin(k*x(i)-omega*time)
            d2phidxdt = k*omega*(reA2*ecosl-imA2* esinl)
            d2phidzdt = k*omega*bl2*(imA2*ecosl+reA2*esinl)
            d2phidxdx = k*k*(-reA2*ecosl+imA2*esinl)
            d2phidzdz = k*k*bl2*bl2*(reA2*ecosl-imA2*esinl)
            d2phidxdz = k*k*bl2*(-imA2*ecosl-reA2*esinl)
            d2psidxdt = k*omega*(reB2*ecoss-imB2*esins)
            d2psidzdt = k*omega*bs2*(imB2*ecoss+reB2*esins)
            d2psidxdx = k*k*(-reB2*ecoss+imB2*esins)
            d2psidzdz = k*k*bs2*bs2*(reB2*ecoss-imB2*esins)
            d2psidxdz = k*k*bs2*(-imB2*ecoss-reB2*esins)
        else
            ecosl = exp(-k*bl1*z(i))*cos(k*x(i)-omega*time)
            esinl = exp(-k*bl1*z(i))*sin(k*x(i)-omega*time)
            ecoss = exp(-k*bs1*z(i))*cos(k*x(i)-omega*time)
            esins = exp(-k*bs1*z(i))*sin(k*x(i)-omega*time)
            d2phidxdt = k*omega*(reA1*ecosl-imA1*esinl)
            d2phidzdt = k*omega*bl1*(-imA1*ecosl-reA1*esinl)
            d2phidxdx = k*k*(-reA1*ecosl+imA1*esinl)
            d2phidzdz = k*k*bl1*bl1*(reA1*ecosl-imA1*esinl)
            d2phidxdz = k*k*bl1*(imA1*ecosl+reA1*esinl)
            d2psidxdt = k*omega*(reB1*ecoss-imB1*esins)
            d2psidzdt = k*omega*bs1*(-imB1*ecoss-reB1*esins)
            d2psidxdx = k*k*(-reB1*ecoss+imB1*esins)
            d2psidzdz = k*k*bs1*bs1*(reB1*ecoss-imB1*esins)
            d2psidxdz = k*k*bs1*(imB1*ecoss+reB1*esins)
        endif
        V1(i)  = d2phidxdt - d2psidzdt
        V2(i)  = 0.0d0
        V3(i)  = d2phidzdt + d2psidxdt
        E11(i) = d2phidxdx - d2psidxdz
        E22(i) = 0.0d0
        E33(i) = d2phidzdz + d2psidxdz
        E12(i) = 0.0d0
        E23(i) = 0.0d0
        E13(i) = d2phidxdz + 0.5d0 * (d2psidxdx - d2psidzdz)
    enddo
end subroutine Stoneley

end module conv_mpi_mod

