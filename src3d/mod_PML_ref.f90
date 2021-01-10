!*******************************************************************!
!*  This module initiates convolutional perfect matching layer.    *!
!*                                                                 *!
!*  Author: Ruichao Ye    Email: Ruichao.ye@gmail.com              *!
!*******************************************************************!

!--------------------------------------------------------------------
module PML_mod
!--------------------------------------------------------------------

use datatype_mod, only : rkind,tetmesh_geometry,&
                         PML_geometry,PML_wavefield,&
                         vector_wavefield,tensor_wavefield
use para_mod,     only : pNp,Nfp,job,PMLAlpha,PMLM,&
                         PMLHx1,PMLHy1,PMLHz1,&
                         PMLHx2,PMLHy2,PMLHz2,&
                         PMLx1_p,PMLy1_p,PMLz1_p,&
                         PMLx2_p,PMLy2_p,PMLz2_p,&
                         PMLx1_type,PMLy1_type,PMLz1_type,&
                         PMLx2_type,PMLy2_type,PMLz2_type,&
                         mpml,rpml,pmlgeometry,Fnamelen,Bproj,&
                         logfid,verbose

    implicit none

    real(kind=rkind),parameter::PI=3.1415926535897932384626d0*2d0
    real(kind=rkind) :: Coeff

!--------------------------------------------------------------------
        contains
!--------------------------------------------------------------------

subroutine init_PML_para(PMLinfo,pNp)
    type(PML_geometry) :: PMLinfo
    integer :: i,j,pNp
    real(kind=rkind) :: rtmp

    Coeff=PMLAlpha*5d2

    do i=1,PMLinfo%Nxb
        do j=1,pNp
            rtmp=PMLinfo%dampx(j,i)
            PMLinfo%alphx(j,i)=(1.0d0-cos(rtmp*PI))/2d0
            PMLinfo%dampx(j,i)=Coeff*(1.0d0-cos(rtmp*PI))/2d0
        enddo
    enddo

    do i=1,PMLinfo%Nyb
        do j=1,pNp
            rtmp=PMLinfo%dampy(j,i)
            PMLinfo%alphy(j,i)=(1.0d0-cos(rtmp*PI))/2d0
            PMLinfo%dampy(j,i)=Coeff*(1.0d0-cos(rtmp*PI))/2d0
        enddo
    enddo

    do i=1,PMLinfo%Nzb
        do j=1,pNp
            rtmp=PMLinfo%dampz(j,i)
            PMLinfo%alphz(j,i)=(1.0d0-cos(rtmp*PI))/2d0
            PMLinfo%dampz(j,i)=Coeff*(1.0d0-cos(rtmp*PI))/2d0
        enddo
    enddo

end subroutine init_PML_para

subroutine init_PML(mesh,PMLinfo)
    type(tetmesh_geometry) :: mesh
    type(PML_geometry) :: PMLinfo
    real(kind=rkind) :: dtmp(pNp),ntmp(3)
    integer :: i,j,head
    logical :: flag,ipml
    real(kind=rkind) :: rtmp
    integer :: nei,fhead,ftail,fhead1,ftail1
    integer :: cx1,cx2,cy1,cy2,cz1,cz2

    allocate(PMLinfo%Xb(mesh%Nele))
    allocate(PMLinfo%Yb(mesh%Nele))
    allocate(PMLinfo%Zb(mesh%Nele))
    PMLinfo%Xb=0; PMLinfo%Yb=0; PMLinfo%Zb=0
    PMLinfo%Nxb=0;PMLinfo%Nyb=0;PMLinfo%Nzb=0

if(trim(PMLgeometry) .eq. 'regula')then

    do i=1,mesh%Nele
        if(any(mesh%x((i-1)*pNp+1:i*pNp).lt.mesh%xmin+PMLHx1).or.&
           any(mesh%x((i-1)*pNp+1:i*pNp).gt.mesh%xmax-PMLHx2))then
                PMLinfo%Nxb=PMLinfo%Nxb+1;PMLinfo%Xb(i)=PMLinfo%Nxb
            do j=1,4
                if(mesh%fbctype(1,j,i).ge.1 .and.&
                        mesh%nx(j,i).gt.0.9d0)then
                    mesh%fbctype(1,j,i)=1
                    mesh%fbctype(2,j,i)=2
                endif
            enddo
        endif
        if(any(mesh%y((i-1)*pNp+1:i*pNp).lt.mesh%ymin+PMLHy1).or.&
           any(mesh%y((i-1)*pNp+1:i*pNp).gt.mesh%ymax-PMLHy2))then
                PMLinfo%Nyb=PMLinfo%Nyb+1;PMLinfo%Yb(i)=PMLinfo%Nyb
            do j=1,4
                if(mesh%fbctype(1,j,i).ge.1 .and.&
                        mesh%ny(j,i).gt.0.9d0)then
                    mesh%fbctype(1,j,i)=1
                    mesh%fbctype(2,j,i)=2
                endif
            enddo
        endif
        if(any(mesh%z((i-1)*pNp+1:i*pNp).lt.mesh%zmin+PMLHz1).or.&
           any(mesh%z((i-1)*pNp+1:i*pNp).gt.mesh%zmax-PMLHz2))then
                PMLinfo%Nzb=PMLinfo%Nzb+1;PMLinfo%Zb(i)=PMLinfo%Nzb
            do j=1,4
                if(mesh%fbctype(1,j,i).ge.1 .and.&
                        mesh%nz(j,i).gt.0.9d0)then
                    mesh%fbctype(1,j,i)=1
                    mesh%fbctype(2,j,i)=2
                endif
            enddo
        endif
    enddo
    if(PMLinfo%Nxb.gt.0)then
        allocate(PMLinfo%alphx(pNp,PMLinfo%Nxb))
        allocate(PMLinfo%dampx(pNp,PMLinfo%Nxb))
    endif
    if(PMLinfo%Nyb.gt.0)then
        allocate(PMLinfo%alphy(pNp,PMLinfo%Nyb))
        allocate(PMLinfo%dampy(pNp,PMLinfo%Nyb))
    endif
    if(PMLinfo%Nzb.gt.0)then
        allocate(PMLinfo%alphz(pNp,PMLinfo%Nzb))
        allocate(PMLinfo%dampz(pNp,PMLinfo%Nzb))
    endif

    do i=1,mesh%Nele
        if(PMLinfo%Xb(i).gt.0)then
            dtmp=mesh%x((i-1)*pNp+1:i*pNp)
            if(any(dtmp.lt.mesh%xmin+PMLHx1))then
                    dtmp=(mesh%xmin+PMLHx1-dtmp)/PMLHx1
            else
                    dtmp=(dtmp-mesh%xmax+PMLHx2)/PMLHx2
            endif
            do j=1,pNp
                    if(dtmp(j) .lt. 0.0d0)dtmp(j)=0.0d0
            enddo
            PMLinfo%dampx(:,PMLinfo%Xb(i)) = dtmp
        endif
        if(PMLinfo%Yb(i).gt.0)then
            dtmp=mesh%y((i-1)*pNp+1:i*pNp)
            if(any(dtmp.lt.mesh%ymin+PMLHy1))then
                    dtmp=(mesh%ymin+PMLHy1-dtmp)/PMLHy1
            else
                    dtmp=(dtmp-mesh%ymax+PMLHy2)/PMLHy2
            endif
            do j=1,pNp
                    if(dtmp(j) .lt. 0.0d0)dtmp(j)=0.0d0
            enddo
            PMLinfo%dampy(:,PMLinfo%Yb(i)) = dtmp
        endif
        if(PMLinfo%Zb(i).gt.0)then
            dtmp=mesh%z((i-1)*pNp+1:i*pNp)
            if(any(dtmp.lt.mesh%zmin+PMLHz1))then
                    dtmp=(mesh%zmin+PMLHz1-dtmp)/PMLHz1
            else
                    dtmp=(dtmp-mesh%zmax+PMLHz2)/PMLHz2
            endif
            do j=1,pNp
                    if(dtmp(j) .lt. 0.0d0)dtmp(j)=0.0d0
            enddo
            PMLinfo%dampz(:,PMLinfo%Zb(i)) = dtmp
        endif
    enddo

    if(verbose .and. logfid.gt.0)then
        write(logfid,*)'PML_Nxb=',PMLinfo%Nxb
        write(logfid,*)'PML_Nyb=',PMLinfo%Nyb
        write(logfid,*)'PML_Nzb=',PMLinfo%Nzb
    endif

elseif(trim(PMLgeometry) .eq. 'rotate')then

    cx1=0;cy1=0;cx2=0;cy2=0;cz1=0;cz2=0
    do i=1,mesh%Nele
        head=(i-1)*pNp+1
        flag=.true.
        call rotation_PML(PMLx1_type,pNp,PMLx1_p,PMLHx1,PMLM(1),&
            mesh%x(head:),mesh%y(head:),mesh%z(head:),&
            mesh%nx(1:,i),mesh%ny(1:,i),mesh%nz(1:,i),flag,&
            mesh%fbctype(1:,1:,i),dtmp,ntmp,ipml)
        if(ipml)then
            cx1=cx1+1
            PMLinfo%Nxb=PMLinfo%Nxb+1;PMLinfo%Xb(i)=PMLinfo%Nxb
        endif
        flag=.true.
        call rotation_PML(PMLx2_type,pNp,PMLx2_p,PMLHx2,PMLM(2),&
            mesh%x(head:),mesh%y(head:),mesh%z(head:),&
            mesh%nx(1:,i),mesh%ny(1:,i),mesh%nz(1:,i),flag,&
            mesh%fbctype(1:,1:,i),dtmp,ntmp,ipml)
        if(ipml)then
            cx2=cx2+1
            PMLinfo%Nxb=PMLinfo%Nxb+1;PMLinfo%Xb(i)=PMLinfo%Nxb
        endif
        flag=.true.
        call rotation_PML(PMLy1_type,pNp,PMLy1_p,PMLHy1,PMLM(3),&
            mesh%x(head:),mesh%y(head:),mesh%z(head:),&
            mesh%nx(1:,i),mesh%ny(1:,i),mesh%nz(1:,i),flag,&
            mesh%fbctype(1:,1:,i),dtmp,ntmp,ipml)
        if(ipml)then
            cy1=cy1+1
            PMLinfo%Nyb=PMLinfo%Nyb+1;PMLinfo%Yb(i)=PMLinfo%Nyb
        endif
        flag=.true.
        call rotation_PML(PMLy2_type,pNp,PMLy2_p,PMLHy2,PMLM(4),&
            mesh%x(head:),mesh%y(head:),mesh%z(head:),&
            mesh%nx(1:,i),mesh%ny(1:,i),mesh%nz(1:,i),flag,&
            mesh%fbctype(1:,1:,i),dtmp,ntmp,ipml)
        if(ipml)then
            cy2=cy2+1
            PMLinfo%Nyb=PMLinfo%Nyb+1;PMLinfo%Yb(i)=PMLinfo%Nyb
        endif
        flag=.true.
        call rotation_PML(PMLz1_type,pNp,PMLz1_p,PMLHz1,PMLM(5),&
            mesh%x(head:),mesh%y(head:),mesh%z(head:),&
            mesh%nx(1:,i),mesh%ny(1:,i),mesh%nz(1:,i),flag,&
            mesh%fbctype(1:,1:,i),dtmp,ntmp,ipml)
        if(ipml)then
            cz1=cz1+1
            PMLinfo%Nzb=PMLinfo%Nzb+1;PMLinfo%Zb(i)=PMLinfo%Nzb
        endif
        flag=.true.
        call rotation_PML(PMLz2_type,pNp,PMLz2_p,PMLHz2,PMLM(6),&
            mesh%x(head:),mesh%y(head:),mesh%z(head:),&
            mesh%nx(1:,i),mesh%ny(1:,i),mesh%nz(1:,i),flag,&
            mesh%fbctype(1:,1:,i),dtmp,ntmp,ipml)
        if(ipml)then
            cz2=cz2+1
            PMLinfo%Nzb=PMLinfo%Nzb+1;PMLinfo%Zb(i)=PMLinfo%Nzb
        endif
    enddo
    if(PMLinfo%Nxb.gt.0)then
        allocate(PMLinfo%alphx(pNp,PMLinfo%Nxb))
        allocate(PMLinfo%dampx(pNp,PMLinfo%Nxb))
        allocate(PMLinfo%xn(3,PMLinfo%Nxb))
    endif
    if(PMLinfo%Nyb.gt.0)then
        allocate(PMLinfo%alphy(pNp,PMLinfo%Nyb))
        allocate(PMLinfo%dampy(pNp,PMLinfo%Nyb))
        allocate(PMLinfo%yn(3,PMLinfo%Nyb))
    endif
    if(PMLinfo%Nzb.gt.0)then
        allocate(PMLinfo%alphz(pNp,PMLinfo%Nzb))
        allocate(PMLinfo%dampz(pNp,PMLinfo%Nzb))
        allocate(PMLinfo%zn(3,PMLinfo%Nzb))
    endif
    do i=1,mesh%Nele
        head=(i-1)*pNp+1
        if(PMLinfo%Xb(i).gt.0)then
            flag=.false.
            call rotation_PML(PMLx1_type,pNp,PMLx1_p,PMLHx1,PMLM(1),&
                mesh%x(head:),mesh%y(head:),mesh%z(head:),&
                mesh%nx(1:,i),mesh%ny(1:,i),mesh%nz(1:,i),flag,&
                mesh%fbctype(1:,1:,i),dtmp,ntmp,ipml)
            rtmp=sum(dtmp)
            if(ipml)then
                PMLinfo%dampx(:,PMLinfo%Xb(i)) = dtmp
                PMLinfo%xn(:,PMLinfo%Xb(i)) = ntmp
            endif
            flag=.false.
            call rotation_PML(PMLx2_type,pNp,PMLx2_p,PMLHx2,PMLM(2),&
                mesh%x(head:),mesh%y(head:),mesh%z(head:),&
                mesh%nx(1:,i),mesh%ny(1:,i),mesh%nz(1:,i),flag,&
                mesh%fbctype(1:,1:,i),dtmp,ntmp,ipml)
            if(ipml .and. sum(dtmp).gt.rtmp)then
                PMLinfo%dampx(:,PMLinfo%Xb(i)) = dtmp
                PMLinfo%xn(:,PMLinfo%Xb(i)) = ntmp
            endif
        endif
        if(PMLinfo%Yb(i).gt.0)then
            flag=.false.
            call rotation_PML(PMLy1_type,pNp,PMLy1_p,PMLHy1,PMLM(3),&
                mesh%x(head:),mesh%y(head:),mesh%z(head:),&
                mesh%nx(1:,i),mesh%ny(1:,i),mesh%nz(1:,i),flag,&
                mesh%fbctype(1:,1:,i),dtmp,ntmp,ipml)
            rtmp=sum(dtmp)
            if(ipml)then
                PMLinfo%dampy(:,PMLinfo%Yb(i)) = dtmp
                PMLinfo%yn(:,PMLinfo%Yb(i)) = ntmp
            endif
            flag=.false.
            call rotation_PML(PMLy2_type,pNp,PMLy2_p,PMLHy2,PMLM(4),&
                mesh%x(head:),mesh%y(head:),mesh%z(head:),&
                mesh%nx(1:,i),mesh%ny(1:,i),mesh%nz(1:,i),flag,&
                mesh%fbctype(1:,1:,i),dtmp,ntmp,ipml)
            if(ipml .and. sum(dtmp).gt.rtmp)then
                PMLinfo%dampy(:,PMLinfo%Yb(i)) = dtmp
                PMLinfo%yn(:,PMLinfo%Yb(i)) = ntmp
            endif
        endif
        if(PMLinfo%Zb(i).gt.0)then
            flag=.false.
            call rotation_PML(PMLz1_type,pNp,PMLz1_p,PMLHz1,PMLM(5),&
                mesh%x(head:),mesh%y(head:),mesh%z(head:),&
                mesh%nx(1:,i),mesh%ny(1:,i),mesh%nz(1:,i),flag,&
                mesh%fbctype(1:,1:,i),dtmp,ntmp,ipml)
            rtmp=sum(dtmp)
            if(ipml)then
                PMLinfo%dampz(:,PMLinfo%Zb(i)) = dtmp
                PMLinfo%zn(:,PMLinfo%Zb(i)) = ntmp
            endif
            flag=.false.
            call rotation_PML(PMLz2_type,pNp,PMLz2_p,PMLHz2,PMLM(6),&
                mesh%x(head:),mesh%y(head:),mesh%z(head:),&
                mesh%nx(1:,i),mesh%ny(1:,i),mesh%nz(1:,i),flag,&
                mesh%fbctype(1:,1:,i),dtmp,ntmp,ipml)
            if(ipml .and. sum(dtmp).gt.rtmp)then
                PMLinfo%dampz(:,PMLinfo%Zb(i)) = dtmp
                PMLinfo%zn(:,PMLinfo%Zb(i)) = ntmp
            endif
        endif
    enddo

    if(verbose .and. logfid.gt.0)then
        write(logfid,*)'PMLM=',PMLM
        write(logfid,*)'PMLx1p=',PMLx1_p
        write(logfid,*)'PMLx2p=',PMLx2_p
        write(logfid,*)'PMLy1p=',PMLy1_p
        write(logfid,*)'PMLy2p=',PMLy2_p
        write(logfid,*)'PMLz1p=',PMLz1_p
        write(logfid,*)'PMLz2p=',PMLz2_p
        write(logfid,*)'PML_Nxb=',PMLinfo%Nxb
        write(logfid,*)'PML_Nyb=',PMLinfo%Nyb
        write(logfid,*)'PML_Nzb=',PMLinfo%Nzb
        write(logfid,*)cx1,cx2,cy1,cy2,cz1,cz2
    endif

endif

    if(PMLinfo%Nxb.gt.0)then
        allocate(PMLinfo%XvmapP(PMLinfo%Nxb*Nfp*4))
        do i=1,mesh%Nele
            if(PMLinfo%Xb(i).gt.0)then
                do j=1,4
                    fhead=(i-1)*Nfp*4+(j-1)*Nfp+1
                    ftail=(i-1)*Nfp*4+j*Nfp
                    fhead1=(PMLinfo%Xb(i)-1)*Nfp*4+(j-1)*Nfp+1
                    ftail1=(PMLinfo%Xb(i)-1)*Nfp*4+j*Nfp
                    nei=(mesh%vmapP(fhead)-1)/pNp
                    if(nei.eq.i.or.nei.gt.mesh%Nele.or.nei.le.0)then
                        head=(PMLinfo%Xb(i)-1)*pNp
                        PMLinfo%XvmapP(fhead1:ftail1)=&
                            mesh%vmapP(fhead:ftail)-i*pNp + head
                    elseif(PMLinfo%Xb(nei) .le. 0)then
                        head=(PMLinfo%Xb(i)-1)*pNp
                        PMLinfo%XvmapP(fhead1:ftail1)=&
                            mesh%vmapP(fhead:ftail)-i*pNp + head
                    else
                        head=(PMLinfo%Xb(nei)-1)*pNp
                        PMLinfo%XvmapP(fhead1:ftail1)=&
                            mesh%vmapP(fhead:ftail)-nei*pNp + head
                    endif
                enddo
            endif
        enddo
    endif
    if(PMLinfo%Nyb.gt.0)then
        allocate(PMLinfo%YvmapP(PMLinfo%Nyb*Nfp*4))
        do i=1,mesh%Nele
            if(PMLinfo%Yb(i).gt.0)then
                do j=1,4
                    fhead=(i-1)*Nfp*4+(j-1)*Nfp+1
                    ftail=(i-1)*Nfp*4+j*Nfp
                    fhead1=(PMLinfo%Yb(i)-1)*Nfp*4+(j-1)*Nfp+1
                    ftail1=(PMLinfo%Yb(i)-1)*Nfp*4+j*Nfp
                    nei=(mesh%vmapP(fhead)-1)/pNp
                    if(nei.eq.i.or.nei.gt.mesh%Nele.or.nei.le.0)then
                        head=(PMLinfo%Yb(i)-1)*pNp
                        PMLinfo%YvmapP(fhead1:ftail1)=&
                            mesh%vmapP(fhead:ftail)-i*pNp + head
                    elseif(PMLinfo%Yb(nei) .le. 0)then
                        head=(PMLinfo%Yb(i)-1)*pNp
                        PMLinfo%YvmapP(fhead1:ftail1)=&
                            mesh%vmapP(fhead:ftail)-i*pNp + head
                    else
                        head=(PMLinfo%Yb(nei)-1)*pNp
                        PMLinfo%YvmapP(fhead1:ftail1)=&
                            mesh%vmapP(fhead:ftail)-nei*pNp + head
                    endif
                enddo
            endif
        enddo
    endif
    if(PMLinfo%Nzb.gt.0)then
        allocate(PMLinfo%ZvmapP(PMLinfo%Nzb*Nfp*4))
        do i=1,mesh%Nele
            if(PMLinfo%Zb(i).gt.0)then
                do j=1,4
                    fhead=(i-1)*Nfp*4+(j-1)*Nfp+1
                    ftail=(i-1)*Nfp*4+j*Nfp
                    fhead1=(PMLinfo%Zb(i)-1)*Nfp*4+(j-1)*Nfp+1
                    ftail1=(PMLinfo%Zb(i)-1)*Nfp*4+j*Nfp
                    nei=(mesh%vmapP(fhead)-1)/pNp
                    if(nei.eq.i.or.nei.gt.mesh%Nele.or.nei.le.0)then
                        head=(PMLinfo%Zb(i)-1)*pNp
                        PMLinfo%ZvmapP(fhead1:ftail1)=&
                            mesh%vmapP(fhead:ftail)-i*pNp + head
                    elseif(PMLinfo%Zb(nei).le.0)then
                        head=(PMLinfo%Zb(i)-1)*pNp
                        PMLinfo%ZvmapP(fhead1:ftail1)=&
                            mesh%vmapP(fhead:ftail)-i*pNp + head
                    else
                        head=(PMLinfo%Zb(nei)-1)*pNp
                        PMLinfo%ZvmapP(fhead1:ftail1)=&
                            mesh%vmapP(fhead:ftail)-nei*pNp + head
                    endif
                enddo
            endif
        enddo
    endif

end subroutine init_PML

subroutine init_PML_variable(pNp,PMLinfo,PMLwf)
    integer :: pNp
    type(PML_geometry) :: PMLinfo
    type(PML_wavefield) :: PMLwf

    if(PMLinfo%Nxb.gt.0)then
        allocate( PMLwf%xv%x(pNp*PMLinfo%Nxb));PMLwf%xv%x=0d0
        allocate( PMLwf%xv%y(pNp*PMLinfo%Nxb));PMLwf%xv%y=0d0
        allocate( PMLwf%xv%z(pNp*PMLinfo%Nxb));PMLwf%xv%z=0d0
        allocate( PMLwf%xS%x(pNp*PMLinfo%Nxb));PMLwf%xS%x=0d0
        allocate( PMLwf%xS%y(pNp*PMLinfo%Nxb));PMLwf%xS%y=0d0
        allocate( PMLwf%xS%z(pNp*PMLinfo%Nxb));PMLwf%xS%z=0d0
    endif

    if(PMLinfo%Nyb.gt.0)then
        allocate( PMLwf%yv%x(pNp*PMLinfo%Nyb));PMLwf%yv%x=0d0
        allocate( PMLwf%yv%y(pNp*PMLinfo%Nyb));PMLwf%yv%y=0d0
        allocate( PMLwf%yv%z(pNp*PMLinfo%Nyb));PMLwf%yv%z=0d0
        allocate( PMLwf%yS%x(pNp*PMLinfo%Nyb));PMLwf%yS%x=0d0
        allocate( PMLwf%yS%y(pNp*PMLinfo%Nyb));PMLwf%yS%y=0d0
        allocate( PMLwf%yS%z(pNp*PMLinfo%Nyb));PMLwf%yS%z=0d0
    endif

    if(PMLinfo%Nzb.gt.0)then
        allocate( PMLwf%zv%x(pNp*PMLinfo%Nzb));PMLwf%zv%x=0d0
        allocate( PMLwf%zv%y(pNp*PMLinfo%Nzb));PMLwf%zv%y=0d0
        allocate( PMLwf%zv%z(pNp*PMLinfo%Nzb));PMLwf%zv%z=0d0
        allocate( PMLwf%zS%x(pNp*PMLinfo%Nzb));PMLwf%zS%x=0d0
        allocate( PMLwf%zS%y(pNp*PMLinfo%Nzb));PMLwf%zS%y=0d0
        allocate( PMLwf%zS%z(pNp*PMLinfo%Nzb));PMLwf%zS%z=0d0
    endif
end subroutine init_PML_variable

subroutine RHSElasticPML(v1,v2,v3,S1,S2,S3,DD,&
                 dv1,dv2,dv3,dS1,dS2,dS3,damp,alph,&
                 PMLv1, PMLv2, PMLv3, PMLS1, PMLS2, PMLS3,&
                rPMLv1,rPMLv2,rPMLv3,rPMLS1,rPMLS2,rPMLS3)

        real(kind=rkind) :: DD(pNp,pNp)
        real(kind=rkind) :: v1(pNp),v2(pNp),v3(pNp),&
                            S1(pNp),S2(pNp),S3(pNp)
        real(kind=rkind) :: dv1(pNp),dv2(pNp),dv3(pNp),&
                            dS1(pNp),dS2(pNp),dS3(pNp)
        real(kind=rkind) :: PMLv1(pNp),PMLv2(pNp),PMLv3(pNp),&
                            PMLS1(pNp),PMLS2(pNp),PMLS3(pNp)
        real(kind=rkind) :: rPMLv1(pNp),rPMLv2(pNp),rPMLv3(pNp),&
                            rPMLS1(pNp),rPMLS2(pNp),rPMLS3(pNp)
        real(kind=rkind) :: PMLrhs(pNp),damp(pNp),alph(pNp)

        rPMLS1 = (damp+PMLAlpha)*PMLS1 + damp*(matmul(DD,S1)+dS1) 

        rPMLS2 = (damp+PMLAlpha)*PMLS2 + damp*(matmul(DD,S2)+dS2)

        rPMLS3 = (damp+PMLAlpha)*PMLS3 + damp*(matmul(DD,S3)+dS3)

        rPMLv1 = (damp+PMLAlpha)*PMLv1 + damp*(matmul(DD,v1)+dv1)

        rPMLv2 = (damp+PMLAlpha)*PMLv2 + damp*(matmul(DD,v2)+dv2) 

        rPMLv3 = (damp+PMLAlpha)*PMLv3 + damp*(matmul(DD,v3)+dv3)

end subroutine RHSElasticPML

subroutine sponge(PMLinfo,V,E,Nele,PMLwf,amp)
    type(PML_geometry) :: PMLinfo
    type(PML_wavefield),intent(in) :: PMLwf
    type(vector_wavefield) :: V
    type(tensor_wavefield) :: E
    integer :: Nele,iele,head,ipmlx,ipmly,ipmlz
    real(kind=rkind) :: damp(pNp),amp
    real(kind=rkind),parameter :: pow=4d0
    
    do iele=1,Nele
        head=(iele-1)*pNp+1
        ipmlx=PMLinfo%Xb(iele)
        ipmly=PMLinfo%Yb(iele)
        ipmlz=PMLinfo%Zb(iele)
        if(ipmlx.gt.0)then
            damp=exp(-amp*(PMLinfo%alphx(:,ipmlx))**pow)
            call array_product(pNp,damp(1:),V%x(head:))
            call array_product(pNp,damp(1:),V%y(head:))
            call array_product(pNp,damp(1:),V%z(head:))
            call array_product(pNp,damp(1:),E%xx(head:))
            call array_product(pNp,damp(1:),E%yy(head:))
            call array_product(pNp,damp(1:),E%zz(head:))
            call array_product(pNp,damp(1:),E%yz(head:))
            call array_product(pNp,damp(1:),E%xz(head:))
            call array_product(pNp,damp(1:),E%xy(head:))
        endif
        if(ipmly.gt.0)then
            damp=exp(-amp*(PMLinfo%alphy(:,ipmly))**pow)
            call array_product(pNp,damp(1:),V%x(head:))
            call array_product(pNp,damp(1:),V%y(head:))
            call array_product(pNp,damp(1:),V%z(head:))
            call array_product(pNp,damp(1:),E%xx(head:))
            call array_product(pNp,damp(1:),E%yy(head:))
            call array_product(pNp,damp(1:),E%zz(head:))
            call array_product(pNp,damp(1:),E%yz(head:))
            call array_product(pNp,damp(1:),E%xz(head:))
            call array_product(pNp,damp(1:),E%xy(head:))
        endif
        if(ipmlz.gt.0)then
            damp=exp(-amp*(PMLinfo%alphz(:,ipmlz))**pow)
            call array_product(pNp,damp(1:),V%x(head:))
            call array_product(pNp,damp(1:),V%y(head:))
            call array_product(pNp,damp(1:),V%z(head:))
            call array_product(pNp,damp(1:),E%xx(head:))
            call array_product(pNp,damp(1:),E%yy(head:))
            call array_product(pNp,damp(1:),E%zz(head:))
            call array_product(pNp,damp(1:),E%yz(head:))
            call array_product(pNp,damp(1:),E%xz(head:))
            call array_product(pNp,damp(1:),E%xy(head:))
        endif
    enddo
end subroutine sponge

subroutine array_product(N,a,x)
        integer,intent(in) :: N
        real(kind=rkind),intent(in) :: a(N)
        real(kind=rkind),intent(inout) :: x(N)
        x=a*x
end subroutine array_product

subroutine PML_flux(pNp,Nfp,LIFT,n,k_media,fbctype,&
                        dS1,dS2,dS3,dv1,dv2,dv3,&
                        fS1,fS2,fS3,fv1,fv2,fv3)
    integer :: pNp,Nfp,k_media,fbctype(2)
    real(kind=rkind) :: LIFT(pNp,Nfp),n
    real(kind=rkind) :: dS1(Nfp),dS2(Nfp),dS3(Nfp)
    real(kind=rkind) :: dv1(Nfp),dv2(Nfp),dv3(Nfp)
    real(kind=rkind) :: fS1(pNp),fS2(pNp),fS3(pNp)
    real(kind=rkind) :: fv1(pNp),fv2(pNp),fv3(pNp)
    real(kind=rkind),parameter :: alpha=0.5d0

    call DGEMV('n',pNp,Nfp,n,LIFT,pNp,dS1,1,1d0,fS1,1)
    call DGEMV('n',pNp,Nfp,n,LIFT,pNp,dS2,1,1d0,fS2,1)
    call DGEMV('n',pNp,Nfp,n,LIFT,pNp,dS3,1,1d0,fS3,1)
    call DGEMV('n',pNp,Nfp,n,LIFT,pNp,dv1,1,1d0,fv1,1)
    call DGEMV('n',pNp,Nfp,n,LIFT,pNp,dv2,1,1d0,fv2,1)
    call DGEMV('n',pNp,Nfp,n,LIFT,pNp,dv3,1,1d0,fv3,1)

end subroutine PML_flux

subroutine PML_flux_acoustic(pNp,Nfp,LIFT,n,&
                        dP,dv,fP,fv)
    integer :: pNp,Nfp
    real(kind=rkind) :: LIFT(pNp,Nfp),n
    real(kind=rkind) :: dP(Nfp),dv(Nfp)
    real(kind=rkind) :: fP(pNp),fv(pNp)

    call DGEMV('n',pNp,Nfp,n,LIFT,pNp,dP,1,1d0,fP,1)
    call DGEMV('n',pNp,Nfp,n,LIFT,pNp,dv,1,1d0,fv,1)

end subroutine PML_flux_acoustic

subroutine PML_average(mesh,pNp,Nfp,Nb,Ib,IvmapP,&
        PMLwfv,rPMLwfv,PMLwfS,rPMLwfS)
    integer,intent(in) :: pNp,Nfp,Nb,Ib(1),IvmapP(1)
    type(tetmesh_geometry) :: mesh
    type(vector_wavefield) :: PMLwfv,rPMLwfv
    type(vector_wavefield) :: PMLwfS,rPMLwfS
    integer :: iele,nei,j,head1,head2,fhead,ftail
    real(kind=rkind) :: weigh(pNp),tmp(Nfp),nx,ny,nz
    integer :: mapM(Nfp),mapP(Nfp)

    if(Nb.le.0)return

    weigh=1d0
    do j=1,4
        mapM=mesh%vmapM((j-1)*Nfp+1:j*Nfp)
        weigh(mapM)=weigh(mapM)+1d0
    enddo

    rPMLwfv%x=PMLwfv%x
    rPMLwfv%y=PMLwfv%y
    rPMLwfv%z=PMLwfv%z
    rPMLwfS%x=PMLwfS%x
    rPMLwfS%y=PMLwfS%y
    rPMLwfS%z=PMLwfS%z
    do iele=1,mesh%Nhele
        if(Ib(iele).gt.0)then
            head1=(Ib(iele)-1)*pNp
            do j=1,4
                fhead=(iele-1)*Nfp*4+(j-1)*Nfp+1
                nei=(mesh%vmapP(fhead)-1)/pNp
                if(nei.eq.iele .or. &
                   nei.gt.mesh%Nhele .or. nei.le.0)cycle
                if(Ib(nei) .le. 0)cycle
                mapM=mesh%vmapM((j-1)*Nfp+1:j*Nfp)   + head1
                fhead=(Ib(iele)-1)*Nfp*4+(j-1)*Nfp+1
                ftail=(Ib(iele)-1)*Nfp*4+j*Nfp
                mapP=IvmapP(fhead:ftail)

                if(mesh%fbctype(1,j,iele).eq.0)then
                    if(mesh%fbctype(2,j,iele).eq.0)then
                rPMLwfv%x(mapM)=rPMLwfv%x(mapM)+PMLwfv%x(mapP)
                rPMLwfv%y(mapM)=rPMLwfv%y(mapM)+PMLwfv%y(mapP)
                rPMLwfv%z(mapM)=rPMLwfv%z(mapM)+PMLwfv%z(mapP)
                rPMLwfS%x(mapM)=rPMLwfS%x(mapM)+PMLwfS%x(mapP)
                rPMLwfS%y(mapM)=rPMLwfS%y(mapM)+PMLwfS%y(mapP)
                rPMLwfS%z(mapM)=rPMLwfS%z(mapM)+PMLwfS%z(mapP)
!                    elseif(mesh%fbctype(2,j,iele).eq.1)then ! solid-fluid
!                nx=mesh%nx(j,iele);ny=mesh%ny(j,iele);nz=mesh%nz(j,iele)
!                tmp=PMLwfv%x(mapP)*nx+PMLwfv%y(mapP)*ny+PMLwfv%z(mapP)*nz
!                rPMLwfv%x(mapM)=rPMLwfv%x(mapM)+tmp*nx
!                rPMLwfv%y(mapM)=rPMLwfv%y(mapM)+tmp*ny
!                rPMLwfv%z(mapM)=rPMLwfv%z(mapM)+tmp*nz
!                rPMLwfS%x(mapM)=rPMLwfS%x(mapM)+PMLwfS%x(mapP)
!                rPMLwfS%y(mapM)=rPMLwfS%y(mapM)+PMLwfS%y(mapP)
!                rPMLwfS%z(mapM)=rPMLwfS%z(mapM)+PMLwfS%z(mapP)
!                    else ! fluid-solid
!                nx=mesh%nx(j,iele);ny=mesh%ny(j,iele);nz=mesh%nz(j,iele)
!                tmp=PMLwfv%x(mapP)*nx+PMLwfv%y(mapP)*ny+PMLwfv%z(mapP)*nz
!                rPMLwfv%x(mapM)=rPMLwfv%x(mapM)+tmp*nx
!                rPMLwfv%y(mapM)=rPMLwfv%y(mapM)+tmp*ny
!                rPMLwfv%z(mapM)=rPMLwfv%z(mapM)+tmp*nz
!                tmp=PMLwfS%x(mapP)*nx+PMLwfS%y(mapP)*ny+PMLwfS%z(mapP)*nz
!                rPMLwfS%x(mapM)=rPMLwfS%x(mapM)+tmp*nx
!                rPMLwfS%y(mapM)=rPMLwfS%y(mapM)+tmp*ny
!                rPMLwfS%z(mapM)=rPMLwfS%z(mapM)+tmp*nz
                    endif
                endif
            enddo
        endif   
    enddo
    do iele=1,Nb
        fhead=(iele-1)*pNp+1;ftail=iele*pNp
        PMLwfv%x(fhead:ftail)=rPMLwfv%x(fhead:ftail)/weigh
        PMLwfv%y(fhead:ftail)=rPMLwfv%y(fhead:ftail)/weigh
        PMLwfv%z(fhead:ftail)=rPMLwfv%z(fhead:ftail)/weigh
        PMLwfS%x(fhead:ftail)=rPMLwfS%x(fhead:ftail)/weigh
        PMLwfS%y(fhead:ftail)=rPMLwfS%y(fhead:ftail)/weigh
        PMLwfS%z(fhead:ftail)=rPMLwfS%z(fhead:ftail)/weigh
    enddo

end subroutine PML_average

subroutine rotation_PML(PML_type,pNp,&
        PML_p,PMLH,PMLM,x,y,z,nx,ny,nz,flag,&
        fbctype,dtmp,PMLn,ipml)
    integer,intent(in) :: PML_type,pNp
    real(kind=rkind),intent(in) :: PML_p(3),PMLH,PMLM,&
        x(pNp),y(pNp),z(pNp),nx(4),ny(4),nz(4)
    logical,intent(inout) :: flag
    integer,intent(inout) :: fbctype(2,4)
    real(kind=rkind),intent(out) :: dtmp(pNp),PMLn(3)
    logical,intent(out) :: ipml
    real(kind=rkind) :: n1,n2,n3,rtmp
    integer :: j

    dtmp=0d0
    ipml=.false.
    if(PML_type .eq. 0)then ! plane
        dtmp=x*PML_p(1)+y*PML_p(2)+z*PML_p(3)
        if(any( dtmp .gt. PMLM-PMLH ))then
            ipml=.true.
            if(flag)then
                do j=1,4
                    if(fbctype(1,j) .ge. 1 .and. &
                        nx(j)*PML_p(1) + &
                        ny(j)*PML_p(2) + &
                        nz(j)*PML_p(3) .gt. 0.9d0)then
                            fbctype(1,j)=1
                            fbctype(2,j)=2
                    endif
                enddo
            else
                dtmp=(dtmp-PMLM+PMLH)/PMLH
                do j=1,pNp
                    if(dtmp(j).lt.0d0)dtmp(j)=0d0
                    if(dtmp(j).gt.1d0)dtmp(j)=1d0
                enddo
                PMLn(1)=-PML_p(1)
                PMLn(2)=-PML_p(2)
                PMLn(3)=-PML_p(3)
            endif
        else
            ipml=.false.
        endif
    elseif(PML_type .eq. 1)then ! sphere outside
        dtmp=(x-PML_p(1))**2+(y-PML_p(2))**2+(z-PML_p(3))**2
        if(any( dtmp .gt. (PMLM-PMLH)**2 ) )then
            ipml=.true.
            if(flag)then
                n1=sum(x)/pNp-PML_p(1)
                n2=sum(y)/pNp-PML_p(2)
                n3=sum(z)/pNp-PML_p(3)
                rtmp=sqrt(n1**2+n2**2+n3**2)
                n1=n1/rtmp;n2=n2/rtmp;n3=n3/rtmp
                do j=1,4
                    if(fbctype(1,j) .ge. 1 .and. &
                        n1*nx(j)+n2*ny(j)+n3*nz(j) .gt. 0.9d0 )then
                        fbctype(1,j)=1
                        fbctype(2,j)=2
                    endif
                enddo
            else
                do j=1,pNp
                    dtmp(j)=sqrt(dtmp(j))
                    dtmp(j)=(dtmp(j)-PMLM+PMLH)/PMLH
                    if(dtmp(j).lt.0d0)dtmp(j)=0d0
                    if(dtmp(j).gt.1d0)dtmp(j)=1d0
                enddo
                PMLn(1)=n1
                PMLn(2)=n2
                PMLn(3)=n3
            endif
        else
            ipml=.false.
        endif
    elseif(PML_type .eq. -1)then ! sphere inside
        dtmp=(x-PML_p(1))**2+(y-PML_p(2))**2+(z-PML_p(3))**2
        if(any( dtmp .lt. (PMLM+PMLH)**2 ) )then
            ipml=.true.
            if(flag)then
                n1=sum(x)/pNp-PML_p(1)
                n2=sum(y)/pNp-PML_p(2)
                n3=sum(z)/pNp-PML_p(3)
                rtmp=-sqrt(n1**2+n2**2+n3**2)
                n1=n1/rtmp;n2=n2/rtmp;n3=n3/rtmp
                do j=1,4
                    if(fbctype(1,j) .ge. 1 .and. &
                        n1*nx(j)+n2*ny(j)+n3*nz(j) .gt. 0.9d0 )then
                        fbctype(1,j)=1
                        fbctype(2,j)=2
                    endif
                enddo
            else
                do j=1,pNp
                    dtmp(j)=sqrt(dtmp(j))
                    dtmp(j)=(PMLM+PMLH-dtmp(j))/PMLH
                    if(dtmp(j).lt.0d0)dtmp(j)=0d0
                    if(dtmp(j).gt.1d0)dtmp(j)=1d0
                enddo
                PMLn(1)=-n1
                PMLn(2)=-n2
                PMLn(3)=-n3
            endif
        else
            ipml=.false.
        endif
    endif
    flag=ipml

end subroutine rotation_PML

end module PML_mod


