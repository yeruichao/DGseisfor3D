!*******************************************************************!
!*  This module initiates convolutional perfect matching layer.    *!
!*                                                                 *!
!*  Author: Ruichao Ye    Email: Ruichao.ye@gmail.com              *!
!*******************************************************************!

!--------------------------------------------------------------------
module PML_mod
!--------------------------------------------------------------------

use datatype_mod, only : rkind,tetmesh_geometry,tetmesh_Domain,&
                         send_recv_pack,PML_geometry,&
                         vector_array,ini_vector,reset_vector,&
                         del_vector,vector_COPY,vector_SCAL
use para_mod,     only : pNp,Nfp,job,PMLAlpha,PMLM,&
                         PMLHx1,PMLHy1,PMLHz1,&
                         PMLHx2,PMLHy2,PMLHz2,&
                         PMLx1_p,PMLy1_p,PMLz1_p,&
                         PMLx2_p,PMLy2_p,PMLz2_p,&
                         PMLx1_type,PMLy1_type,PMLz1_type,&
                         PMLx2_type,PMLy2_type,PMLz2_type,&
                         rpml,pmlgeometry,Fnamelen,&
                         logfid,verbose
use parallel_mod, only : sample_SRtype_info,init_SRtype

    implicit none

    real(kind=rkind),parameter :: PI=3.1415926535897932384626d0*2d0
    logical,parameter :: mpml=.true.
    real(kind=rkind),parameter :: mpmlq=0.3d0

!--------------------------------------------------------------------
        contains
!--------------------------------------------------------------------

subroutine init_PML(mesh,subdomain,PMLinfo)
    type(tetmesh_geometry) :: mesh
    type(tetmesh_Domain) :: subdomain
    type(PML_geometry) :: PMLinfo
    real(kind=rkind),pointer :: dtmp(:)
    integer :: i,j,k,head,Ndim
    logical,allocatable :: flag(:)
    logical :: ipml,itmp
    real(kind=rkind) :: rtmp
    real(kind=rkind),allocatable :: dampnx(:,:),dampny(:,:),dampnz(:,:)

    PMLinfo%Nele=0
    Ndim=mesh%Nele*pNp
    call ini_vector(Ndim,PMLinfo%damp)
    call reset_vector(PMLinfo%damp)
    allocate(flag(mesh%Nele))
    flag=.false.

    if(trim(PMLgeometry) .eq. 'regula')then
        call vector_COPY(Ndim,mesh%coord,PMLinfo%damp)
        ipml=.false.
        do i=1,mesh%Nele
            head=(i-1)*pNp+1
    
            call PML_damp(pNp,mesh%xmin,mesh%xmax,PMLHx1,PMLHx2,&
                PMLinfo%damp%x(head:),itmp,&
                mesh%nx(:,i),mesh%fbctype(:,:,i))
            flag(i)=flag(i).or.itmp
    
            call PML_damp(pNp,mesh%ymin,mesh%ymax,PMLHy1,PMLHy2,&
                PMLinfo%damp%y(head:),itmp,&
                mesh%ny(:,i),mesh%fbctype(:,:,i))
            flag(i)=flag(i).or.itmp
    
            call PML_damp(pNp,mesh%zmin,mesh%zmax,PMLHz1,PMLHz2,&
                PMLinfo%damp%z(head:),itmp,&
                mesh%nz(:,i),mesh%fbctype(:,:,i))
            flag(i)=flag(i).or.itmp

            if(flag(i).and.mpml)then
                call MPML_damp(pNp,PMLinfo%damp%x(head:),&
                                   PMLinfo%damp%y(head:),&
                                   PMLinfo%damp%z(head:))
                call MPML_damp(pNp,PMLinfo%damp%y(head:),&
                                   PMLinfo%damp%x(head:),&
                                   PMLinfo%damp%z(head:))
                call MPML_damp(pNp,PMLinfo%damp%z(head:),&
                                   PMLinfo%damp%x(head:),&
                                   PMLinfo%damp%y(head:))
            endif

            ipml=ipml.or.flag(i)
        enddo
        PMLinfo%pml=ipml
    elseif(trim(PMLgeometry) .eq. 'rotate')then
        ipml=.false.
        allocate(dampnx(3,mesh%Nele))
        allocate(dampny(3,mesh%Nele))
        allocate(dampnz(3,mesh%Nele))
        do i=1,mesh%Nele
            head=(i-1)*pNp+1

            call rPML_damp(pNp,PMLx1_type,PMLM(1),PMLx1_p,PMLHx1,&
                mesh%coord%x(head:),&
                mesh%coord%y(head:),&
                mesh%coord%z(head:),&
                PMLinfo%damp%x(head:),dampnx(:,i),itmp,&
                mesh%nx(:,i),mesh%ny(:,i),mesh%nz(:,i),&
                mesh%fbctype(:,:,i))
            flag(i)=flag(i).or.itmp
    
            call rPML_damp(pNp,PMLx2_type,PMLM(2),PMLx2_p,PMLHx2,&
                mesh%coord%x(head:),&
                mesh%coord%y(head:),&
                mesh%coord%z(head:),&
                PMLinfo%damp%x(head:),dampnx(:,i),itmp,&
                mesh%nx(:,i),mesh%ny(:,i),mesh%nz(:,i),&
                mesh%fbctype(:,:,i))
            flag(i)=flag(i).or.itmp
    
            call rPML_damp(pNp,PMLy1_type,PMLM(3),PMLy1_p,PMLHy1,&
                mesh%coord%x(head:),&
                mesh%coord%y(head:),&
                mesh%coord%z(head:),&
                PMLinfo%damp%y(head:),dampny(:,i),itmp,&
                mesh%nx(:,i),mesh%ny(:,i),mesh%nz(:,i),&
                mesh%fbctype(:,:,i))
            flag(i)=flag(i).or.itmp
    
            call rPML_damp(pNp,PMLy2_type,PMLM(4),PMLy2_p,PMLHy2,&
                mesh%coord%x(head:),&
                mesh%coord%y(head:),&
                mesh%coord%z(head:),&
                PMLinfo%damp%y(head:),dampny(:,i),itmp,&
                mesh%nx(:,i),mesh%ny(:,i),mesh%nz(:,i),&
                mesh%fbctype(:,:,i))
            flag(i)=flag(i).or.itmp
    
            call rPML_damp(pNp,PMLz1_type,PMLM(5),PMLz1_p,PMLHz1,&
                mesh%coord%x(head:),&
                mesh%coord%y(head:),&
                mesh%coord%z(head:),&
                PMLinfo%damp%z(head:),dampnz(:,i),itmp,&
                mesh%nx(:,i),mesh%ny(:,i),mesh%nz(:,i),&
                mesh%fbctype(:,:,i))
            flag(i)=flag(i).or.itmp
    
            call rPML_damp(pNp,PMLz2_type,PMLM(6),PMLz2_p,PMLHz2,&
                mesh%coord%x(head:),&
                mesh%coord%y(head:),&
                mesh%coord%z(head:),&
                PMLinfo%damp%z(head:),dampnz(:,i),itmp,&
                mesh%nx(:,i),mesh%ny(:,i),mesh%nz(:,i),&
                mesh%fbctype(:,:,i))
            flag(i)=flag(i).or.itmp
    
            if(flag(i).and.mpml)then
                call MPML_damp(pNp,PMLinfo%damp%x(head:),&
                                   PMLinfo%damp%y(head:),&
                                   PMLinfo%damp%z(head:))
                call MPML_damp(pNp,PMLinfo%damp%y(head:),&
                                   PMLinfo%damp%x(head:),&
                                   PMLinfo%damp%z(head:))
                call MPML_damp(pNp,PMLinfo%damp%z(head:),&
                                   PMLinfo%damp%x(head:),&
                                   PMLinfo%damp%y(head:))
            endif

            ipml=ipml.or.flag(i)
            if(flag(i))PMLinfo%Nele=PMLinfo%Nele+1
        enddo
        allocate(PMLinfo%Jac(3,3,PMLinfo%Nele))
        allocate(PMLinfo%ipml(PMLinfo%Nele))
        k=0
        do i=1,mesh%Nele
            if(flag(i))then
                k=k+1
                PMLinfo%ipml(k)=i
                call rPML_axis(dampnx(:,i),dampny(:,i),dampnz(:,i))
                PMLinfo%Jac(:,1,k)=dampnx(:,i)
                PMLinfo%Jac(:,2,k)=dampny(:,i)
                PMLinfo%Jac(:,3,k)=dampnz(:,i)
            endif
        enddo
        deallocate(dampnx,dampny,dampnz)
    endif

    call sample_SRtype_info(flag,subdomain%SRtype%host,&
        PMLinfo%SRtype%host)
    call sample_SRtype_info(flag,subdomain%SRtype%gost,&
        PMLinfo%SRtype%gost)
    call init_SRtype(PMLinfo%SRtype,pNp,mesh%Nele,9)

    deallocate(flag)

end subroutine init_PML

subroutine PML_damp(pNp,Hmin,Hmax,PMLH1,PMLH2,damp,ipml,n,fbctype)
    integer,intent(in) :: pNp
    real(kind=rkind),intent(in) :: Hmin,Hmax,PMLH1,PMLH2,n(4)
    real(kind=rkind) :: damp(pNp)
    integer :: fbctype(2,4)
    logical :: ipml
    integer :: j

    ipml=.true.
    if(any(damp.lt.Hmin+PMLH1))then
        damp=(Hmin+PMLH1-damp)/PMLH1
    elseif(any(damp.gt.Hmax-PMLH2))then
        damp=(damp-Hmax+PMLH2)/PMLH2
    else
        ipml=.false.
        damp=0d0
    endif
    if(ipml)then
        do j=1,pNp
            damp(j)=max(damp(j),0d0)
        enddo
        do j=1,4
            if(fbctype(1,j).eq.1 .and. abs(n(j)).gt.0.9)then
                fbctype(2,j)=2
            endif
        enddo
    endif
end subroutine PML_damp

subroutine rPML_damp(pNp,PML_type,PMLM,PMLn,PMLH,x,y,z,&
        damp,dampn,ipml,nx,ny,nz,fbctype)
    integer,intent(in) :: pNp,PML_type
    real(kind=rkind),intent(in) :: PMLM,PMLn(3),PMLH,&
        x(pNp),y(pNp),z(pNp),nx(4),ny(4),nz(4)
    logical,intent(inout) :: ipml
    real(kind=rkind) :: damp(pNp),dampn(3),dtmp(pNp)
    integer :: fbctype(2,4)
    integer :: j
    real(kind=rkind) :: n1,n2,n3,rtmp

    ipml=.false.
    if(PML_type .eq. 0)then ! plane
        dtmp=x*PMLn(1)+y*PMLn(2)+z*PMLn(3)
        if(any(dtmp.gt.PMLM-PMLH))then
            ipml=.true.
            do j=1,4
                if(fbctype(1,j) .eq. 1 .and. &
                  abs(nx(j)*PMLn(1)+ny(j)*PMLn(2)+nz(j)*PMLn(3))&
                  .gt. 0.9d0)then
                    fbctype(2,j)=2
                endif
            enddo
            dtmp=(dtmp-PMLM+PMLH)/PMLH
            do j=1,pNp
                damp(j)=max(damp(j),dtmp(j))
            enddo
            rtmp=sqrt(PMLn(1)**2+PMLn(2)**2+PMLn(3)**2)
            dampn=PMLn/rtmp
        endif
    elseif(abs(PML_type) .eq. 1)then ! sphere
        dtmp=sqrt((x-PMLn(1))**2+(y-PMLn(2))**2+(z-PMLn(3))**2)
        if(any(PML_type*(dtmp-PMLM) + PMLH .gt. 0d0))then
            ipml=.true.
            n1=sum(x)/pNp-PMLn(1)
            n2=sum(y)/pNp-PMLn(2)
            n3=sum(z)/pNp-PMLn(3)
            rtmp=sqrt(n1**2+n2**2+n3**2)
            n1=n1/rtmp;n2=n2/rtmp;n3=n3/rtmp
            do j=1,4
                if(fbctype(1,j).eq.1 .and. &
                  abs(n1*nx(j)+n2*ny(j)+n3*nz(j)) .gt. 0.9d0 )then
                    fbctype(2,j)=2
                endif
            enddo
            dtmp=(PML_type*(dtmp-PMLM)+PMLH)/PMLH
            do j=1,pNp
                damp(j)=max(damp(j),dtmp(j))
            enddo
            dampn(1)=n1
            dampn(2)=n2
            dampn(3)=n3
        endif
    endif
end subroutine rPML_damp

subroutine MPML_damp(pNp,damp0,damp1,damp2)
    integer,intent(in) :: pNp
    real(kind=rkind) :: damp0(pNp),damp1(pNp),damp2(pNp)
    integer :: i
    do i=1,pNp
        damp1(i)=max(damp1(i),damp0(i)*mpmlq)
        damp2(i)=max(damp2(i),damp0(i)*mpmlq)
    enddo
end subroutine MPML_damp

subroutine rPML_axis(r,s,t)
    real(kind=rkind) :: r(3),s(3),t(3)
    real(kind=rkind) :: rr,ss,tt,n(3)
    real(kind=rkind),parameter :: tol=1d-1
    integer :: i

    rr=r(1)**2+r(2)**2+r(3)**2
    ss=s(1)**2+s(2)**2+s(3)**2
    tt=t(1)**2+t(2)**2+t(3)**2
    n=0d0
    if(rr.gt.tol .and. ss.le.tol .and. tt.le.tol)then
        i=minloc(abs(r),1)
        n(i)=1d0
        call cross_product(r,n,t)
        tt=sqrt(t(1)**2+t(2)**2+t(3)**2)
        t=t/tt
        call cross_product(t,r,s)
    elseif(ss.gt.tol .and. rr.le.tol .and. tt.le.tol)then
        i=minloc(abs(s),1)
        n(i)=1d0
        call cross_product(s,n,r)
        rr=sqrt(r(1)**2+r(2)**2+r(3)**2)
        r=r/rr
        call cross_product(r,s,t)
    elseif(tt.gt.tol .and. rr.le.tol .and. ss.le.tol)then
        i=minloc(abs(t),1)
        n(i)=1d0
        call cross_product(t,n,s)
        ss=sqrt(s(1)**2+s(2)**2+s(3)**2)
        s=s/ss
        call cross_product(s,t,r)
    elseif(rr.gt.tol .and. ss.gt.tol .and. tt.le.tol)then
        call cross_product(r,s,t)
        tt=sqrt(t(1)**2+t(2)**2+t(3)**2)
        t=t/tt
    elseif(ss.gt.tol .and. tt.gt.tol .and. rr.le.tol)then
        call cross_product(s,t,r)
        rr=sqrt(r(1)**2+r(2)**2+r(3)**2)
        r=r/rr
    elseif(tt.gt.tol .and. rr.gt.tol .and. ss.le.tol)then
        call cross_product(t,r,s)
        ss=sqrt(s(1)**2+s(2)**2+s(3)**2)
        s=s/ss
    endif
end subroutine rPML_axis

subroutine cross_product(r,s,t)
    real(kind=rkind) :: r(3),s(3),t(3)
    t(1)=r(2)*s(3)-r(3)*s(2)
    t(2)=r(3)*s(1)-r(1)*s(3)
    t(3)=r(1)*s(2)-r(2)*s(1)
end subroutine cross_product

end module PML_mod


