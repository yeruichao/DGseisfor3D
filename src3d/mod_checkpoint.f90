!*******************************************************************!
!*  This module contains checkpointing subroutines                 *!
!*                                                                 *!
!*  Author: Ruichao Ye    Email: Ruichao.ye@gmail.com              *!
!*******************************************************************!

!-------------------------------------------------------------------!
module checkpoint_mod
!-------------------------------------------------------------------!
    use datatype_mod, only: rkind,vector_array,tensor_array,&
                            tetmesh_geometry,tetmesh_material,&
                            tetmesh_Domain,PML_geometry,&
                            matrices,receivers,rupture,surface,&
                            send_recv_pack,wavefield,&
                            ini_tensor,del_tensor,&
                            ini_vector,reset_vector,del_vector
    use meshfile_mod, only: glob_Nele
    use geometry_mod, only: form_ref_matrices
    use para_mod,     only: Fnamelen,pNp,Nfp,verbose,logfid,&
                            startTime,job,def_pOrder,PMLM,friction,&
                            prestress,refer_T0
    use parallel_mod, only: init_SRtype
    use surface_mod,  only: rupture_malloc

    implicit none
!---------------------------------------------------------------------
    contains
!---------------------------------------------------------------------

subroutine subdomain_file(filename,mesh,subdomain,material,matrix,&
        Ndomain,k_inout)

! Save or load time-independent global informations:
!       glob_Nele,pNp,Nfp,Ndomain
! and subdomain informations:
!       mesh%(Nele,Nhele,Ninele,globID,nx,ny,nz,
!       Jac,invJac,sJac,detJ,Fscale,
!       vmapM,vmapP,fbctype,pOrder,
!       x,y,z,link)
!       subdomain%(inele,innod)

    character(len=*),intent(in)    :: filename,k_inout
    type(tetmesh_geometry)  :: mesh
    type(tetmesh_material)  :: material
    type(tetmesh_Domain)    :: subdomain
    type(matrices)          :: matrix
    integer,intent(in)      :: Ndomain
    character(len=Fnamelen) :: fname
    integer :: fid,itmp,pid
    integer :: bsize
    integer :: Nele,N_DD_Conn,maxN_host,maxN_ghost
    integer :: Ndof
    logical :: alive
    real(kind=rkind) :: r(pNp),s(pNp),t(pNp)
    integer,pointer :: Fmask(:,:)

    pid=subdomain%id-1

    if(k_inout .eq. 'save')then

        if(pid.eq.0)then
            call subdomain_ckpfname(fname,fid,filename,&
                'Para',-1,-1,'txt',alive)
            open(fid,file=trim(fname),status='replace')
            write(fid,*)glob_Nele
            write(fid,*)pNp
            write(fid,*)Nfp
            write(fid,*)Ndomain
            close(fid)
        endif

        if(verbose .and. logfid.gt.0)then
            write(logfid,*)&
                'Writing Subdomain information to file:'
            write(logfid,*)trim(fname)
            write(logfid,*)'Nele=',mesh%Nele
            write(logfid,*)'Nhele=',mesh%Nhele
            write(logfid,*)'Ninele=',mesh%Ninele
            write(logfid,*)'N_DD_Conn=',&
                subdomain%SRtype%host%N_DD_Conn,&
                subdomain%SRtype%gost%N_DD_Conn
            write(logfid,*)'DD_Conn_host=' ,&
                subdomain%SRtype%host%DD_Conn
            write(logfid,*)'DD_Conn_ghost=',&
                subdomain%SRtype%gost%DD_Conn
            write(logfid,*)'maxN_host=' ,&
                subdomain%SRtype%host%maxN_Conn
            write(logfid,*)'maxN_ghost=',&
                subdomain%SRtype%gost%maxN_Conn
        endif

        call subdomain_ckpfname(fname,fid,filename,&
            'Subdomain',pid,-1,'bin',alive)
        open(unit=fid,file=trim(fname),access='stream',&
            status='replace')
        call meshgeometry_file(fid,mesh,k_inout)
        call SRtype_file(fid,subdomain%SRtype,k_inout)
!        call surface_file(fid,mesh%surf,Nfp,k_inout)
        close(fid)

        if(mesh%PMLinfo%pml)then
            call subdomain_ckpfname(fname,fid,filename,&
                'PMLdamp',pid,-1,'bin',alive)
            open(unit=fid,file=trim(fname),access='stream',&
                status='replace')
            call PMLinfo_file(fid,mesh%Nele,pNp,mesh%PMLinfo,k_inout)
            call SRtype_file(fid,mesh%PMLinfo%SRtype,k_inout)
            close(fid)
        endif

        call vector_file(material%rho,mesh%Nele*pNp,filename,&
            'rho',-1,pid,'save')

        call subdomain_ckpfname(fname,fid,filename,&
            'Material',pid,-1,'bin',alive)
        open(unit=fid,file=trim(fname),access='stream',&
            status='replace')
        call material_file(fid,mesh%Nele,pNp,material,k_inout)
        close(fid)

        call subdomain_ckpfname(fname,fid,filename,&
            'Rupture',pid,-1,'bin',alive)
        open(unit=fid,file=trim(fname),access='stream',&
            status='replace')
        call rupture_file(fid,mesh%rupt,k_inout)
        close(fid)

    elseif(k_inout .eq. 'load')then

        call subdomain_ckpfname(fname,fid,filename,&
            'Para',-1,-1,'txt',alive)
        if(.not. alive)then
            print*,'file "',trim(fname),' " does not exist.'
            stop
        endif
        open(fid,file=trim(fname),status='old')
        read(fid,*)glob_Nele
        itmp=pNp
        read(fid,*)pNp
        if(pid.eq.0 .and. itmp.ne.pNp)then
            print*,'Polynomial order changed,', &
                ' running with pNp=',pNp
        endif
        read(fid,*)Nfp
        read(fid,*)itmp
        close(fid)
        if(itmp.ne.Ndomain)then
            print*,'Number of subdomains changed,', &
                   ' please run with ',itmp,' processors.'
            stop
        endif

        allocate(matrix%Fmask(Nfp,4));Fmask=>matrix%Fmask
        call form_ref_matrices(def_pOrder,pNp,Nfp,r,s,t,Fmask,matrix)

        if(verbose .and. logfid.gt.0)then
            write(logfid,*)&
                'Reading Subdomain information from file:'
            write(logfid,*)trim(fname)
            write(logfid,*)'Nele=',mesh%Nele
            write(logfid,*)'Nhele=',mesh%Nhele
            write(logfid,*)'Ninele=',mesh%Ninele
        endif

        call subdomain_ckpfname(fname,fid,filename,&
            'Subdomain',pid,-1,'bin',alive)
        if(.not. alive)then
            print*,'file "',trim(fname),' " does not exist.'
            stop
        endif
        open(unit=fid,file=trim(fname),access='stream',&
            status='old',position='rewind')
        call meshgeometry_file(fid,mesh,k_inout)
        call SRtype_file(fid,subdomain%SRtype,k_inout)
!        call surface_file(fid,mesh%surf,Nfp,k_inout)
        close(fid)

        call init_SRtype(subdomain%SRtype,pNp,mesh%Nele,9)

        if(verbose .and. logfid.gt.0)then
            write(logfid,*)&
                'Reading Subdomain information from file:'
            write(logfid,*)trim(fname)
            write(logfid,*)'Nele=',mesh%Nele
            write(logfid,*)'Nhele=',mesh%Nhele
            write(logfid,*)'Ninele=',mesh%Ninele
            write(logfid,*)'N_DD_Conn=',&
                subdomain%SRtype%host%N_DD_Conn,&
                subdomain%SRtype%gost%N_DD_Conn
            write(logfid,*)'DD_Conn_host=' ,&
                subdomain%SRtype%host%DD_Conn
            write(logfid,*)'DD_Conn_ghost=',&
                subdomain%SRtype%gost%DD_Conn
            write(logfid,*)'maxN_host=' ,&
                subdomain%SRtype%host%maxN_Conn
            write(logfid,*)'maxN_ghost=',&
                subdomain%SRtype%gost%maxN_Conn
        endif

        if(mesh%PMLinfo%pml)then
            call subdomain_ckpfname(fname,fid,filename,&
                'PMLdamp',pid,-1,'bin',alive)
            if(.not. alive)then
                print*,'file "',trim(fname),' " does not exist.'
                stop
            endif
            open(unit=fid,file=trim(fname),access='stream',&
                status='old',position='rewind')
            call PMLinfo_file(fid,mesh%Nele,pNp,mesh%PMLinfo,k_inout)
            call SRtype_file(fid,mesh%PMLinfo%SRtype,k_inout)
            close(fid)
            call init_SRtype(mesh%PMLinfo%SRtype,pNp,mesh%Nele,9)
        endif

        call subdomain_ckpfname(fname,fid,filename,&
            'Material',pid,-1,'bin',alive)
        if(.not. alive)then
            print*,'file "',trim(fname),' " does not exist.'
            stop
        endif
        open(unit=fid,file=trim(fname),access='stream',&
            status='old',position='rewind')
        call material_file(fid,mesh%Nele,pNp,material,k_inout)
        close(fid)

        call subdomain_ckpfname(fname,fid,filename,&
            'Rupture',pid,-1,'bin',alive)
        if(.not. alive)then
            print*,'file "',trim(fname),' " does not exist.'
            stop
        endif
        open(unit=fid,file=trim(fname),access='stream',&
            status='old')
        call rupture_file(fid,mesh%rupt,k_inout)
        close(fid)

    endif

end subroutine subdomain_file

subroutine time_subdomain_file(filename,mesh,&
        isnap,time,istep,irec,W,Ws,Wss,Vsss,pid,k_inout)

    character(len=*),intent(in) :: filename,k_inout
    type(tetmesh_geometry)      :: mesh
    type(wavefield)             :: W,Ws,Wss
    type(vector_array)          :: Vsss
    integer,intent(in)          :: pid
    integer          :: isnap,istep,irec
    real(kind=rkind) :: time
    character(len=Fnamelen) :: fname
    integer     :: fid,itmp
    integer     :: Nele,Ndof
    integer     :: iface,i,head,head1,head2
    logical     :: alive

    if(W%E%symmetric)then
        Ndof=mesh%Nele*pNp*9
    else
        Ndof=mesh%Nele*pNp*12
    endif

    if(k_inout .eq. 'save')then

        if(pid.eq.0)then
            call subdomain_ckpfname(fname,fid,filename,&
                'Para',-1,isnap,'bin',alive)
            open(unit=fid,file=trim(fname),access='stream',&
                status='replace')
            write(fid)time
            write(fid)istep
            write(fid)irec
            close(fid)
        endif

        call vector_file(W%array,Ndof,filename,&
            'Wave',isnap,pid,k_inout)

        call vector_file(W%U%array,mesh%Nele*pNp*3,filename,&
            'WaveU',isnap,pid,k_inout)

!        if(mesh%PMLinfo%pml)then
!            call subdomain_ckpfname(fname,fid,filename,&
!                'PMLt',pid,isnap,'bin',alive)
!            open(unit=fid,file=trim(fname),access='stream',&
!                status='replace')
!            write(fid)Ws%array,Wss%array,Vsss%array
!            close(fid)
!        endif

        if(mesh%rupt%Nface .gt.0)then
            call subdomain_ckpfname(fname,fid,filename,&
                'Rupt',pid,isnap,'bin',alive)
            open(unit=fid,file=trim(fname),access='stream',&
                status='replace')
            call time_rupture_file(fid,mesh%rupt,k_inout)
            close(fid)
        endif

    elseif(k_inout .eq. 'load')then
        call subdomain_ckpfname(fname,fid,filename,&
            'Para',-1,isnap,'bin',alive)
        if(.not. alive)then
            print*,'file "',trim(fname),' " does not exist.'
            stop
        endif
        open(unit=fid,file=trim(fname),access='stream',&
            status='old',position='rewind')
        read(fid)time
        read(fid)istep
        read(fid)irec
        close(fid)

        call vector_file(W%array,Ndof,filename,&
            'Wave',isnap,pid,k_inout)

        call vector_file(W%U%array,mesh%Nele*pNp*3,filename,&
            'WaveU',isnap,pid,k_inout)

!        if(mesh%PMLinfo%pml)then
!            call subdomain_ckpfname(fname,fid,filename,&
!                'PMLt',pid,isnap,'bin',alive)
!            if(.not. alive)then
!                print*,'file "',trim(fname),' " does not exist.'
!                stop
!            endif
!            open(unit=fid,file=trim(fname),access='stream',&
!                status='old')
!            read(fid)Ws%array,Wss%array,Vsss%array
!            close(fid)
!        endif

        if(mesh%rupt%Nface .gt.0)then
            call subdomain_ckpfname(fname,fid,filename,&
                'Rupt',pid,isnap,'bin',alive)
            if(.not. alive)then
                print*,'file "',trim(fname),' " does not exist.'
                stop
            endif
            open(unit=fid,file=trim(fname),access='stream',&
                status='old')
            call time_rupture_file(fid,mesh%rupt,k_inout)
            close(fid)
!            do iface=1,mesh%rupt%Nhface
!                head=(iface-1)*pNp
!                head1=(mesh%rupt%T2E(1,iface)-1)*pNp
!                head2=(mesh%rupt%T2E(2,iface)-1)*pNp
!                do i=1,pNp
!                    mesh%rupt%Em%xx(head1+i)=W%E%xx(head+i)
!                    mesh%rupt%Em%yx(head1+i)=W%E%yx(head+i)
!                    mesh%rupt%Em%zx(head1+i)=W%E%zx(head+i)
!                    mesh%rupt%Em%xy(head1+i)=W%E%xy(head+i)
!                    mesh%rupt%Em%yy(head1+i)=W%E%yy(head+i)
!                    mesh%rupt%Em%zy(head1+i)=W%E%zy(head+i)
!                    mesh%rupt%Em%xz(head1+i)=W%E%xz(head+i)
!                    mesh%rupt%Em%yz(head1+i)=W%E%yz(head+i)
!                    mesh%rupt%Em%zz(head1+i)=W%E%zz(head+i)
!                    mesh%rupt%Ep%xx(head2+i)=W%E%xx(head+i)
!                    mesh%rupt%Ep%yx(head2+i)=W%E%yx(head+i)
!                    mesh%rupt%Ep%zx(head2+i)=W%E%zx(head+i)
!                    mesh%rupt%Ep%xy(head2+i)=W%E%xy(head+i)
!                    mesh%rupt%Ep%yy(head2+i)=W%E%yy(head+i)
!                    mesh%rupt%Ep%zy(head2+i)=W%E%zy(head+i)
!                    mesh%rupt%Ep%xz(head2+i)=W%E%xz(head+i)
!                    mesh%rupt%Ep%yz(head2+i)=W%E%yz(head+i)
!                    mesh%rupt%Ep%zz(head2+i)=W%E%zz(head+i)
!                enddo
!            enddo
        endif

    endif

end subroutine time_subdomain_file

subroutine meshgeometry_file(fid,mesh,k_inout)
    integer,intent(in) :: fid
    type(tetmesh_geometry)  :: mesh
    character(len=*),intent(in) :: k_inout

    if(k_inout .eq. 'save')then
        write(fid)mesh%Nele,mesh%Nhele,mesh%Ninele
        write(fid)mesh%PMLinfo%pml
        write(fid)mesh%coord%array
        write(fid)mesh%globID
        write(fid)mesh%xmax,mesh%xmin
        write(fid)mesh%ymax,mesh%ymin
        write(fid)mesh%zmax,mesh%zmin
        write(fid)mesh%nx
        write(fid)mesh%ny
        write(fid)mesh%nz
        write(fid)mesh%Jac
        write(fid)mesh%invJac
        write(fid)mesh%detJ
        write(fid)mesh%Fscale
        write(fid)mesh%vmapM
        write(fid)mesh%vmapP
        write(fid)mesh%fbctype
    else
        read(fid)mesh%Nele,mesh%Nhele,mesh%Ninele
        read(fid)mesh%PMLinfo%pml

        call ini_vector      (pNp*mesh%Nele,mesh%coord)
        allocate(mesh%globID (    mesh%Nele  ))
        allocate(mesh%nx     (  4,mesh%Nele  ))
        allocate(mesh%ny     (  4,mesh%Nele  ))
        allocate(mesh%nz     (  4,mesh%Nele  ))
        allocate(mesh%Jac    (    mesh%Nele,9))
        allocate(mesh%invJac (    mesh%Nele,9))
        allocate(mesh%detJ   (    mesh%Nele  ))
        allocate(mesh%Fscale (  4,mesh%Nele  ))
        allocate(mesh%vmapM  (4*Nfp))
        allocate(mesh%vmapP  (4*Nfp*mesh%Nele))
        allocate(mesh%fbctype(  2,4,mesh%Nele))

        read(fid)mesh%coord%array
        read(fid)mesh%globID
        read(fid)mesh%xmax,mesh%xmin
        read(fid)mesh%ymax,mesh%ymin
        read(fid)mesh%zmax,mesh%zmin
        read(fid)mesh%nx
        read(fid)mesh%ny
        read(fid)mesh%nz
        read(fid)mesh%Jac
        read(fid)mesh%invJac
        read(fid)mesh%detJ
        read(fid)mesh%Fscale
        read(fid)mesh%vmapM
        read(fid)mesh%vmapP
        read(fid)mesh%fbctype
    endif

end subroutine meshgeometry_file

subroutine surface_file(fid,surf,Nfp,k_inout)
    integer,intent(in) :: fid
    type(surface) :: surf
    character(len=*),intent(in) :: k_inout
    integer :: Nfp,Nface

    if(k_inout .eq. 'save')then
        write(fid)surf%Nface
        write(fid)surf%coord%array
        write(fid)surf%nx
        write(fid)surf%ny
        write(fid)surf%nz
    else
        read(fid)Nface
        surf%Nface=Nface
        if(Nface.gt.0)then
            call ini_vector(Nfp*Nface,surf%coord)
            allocate(surf%nx(Nface))
            allocate(surf%ny(Nface))
            allocate(surf%nz(Nface))
            call ini_vector(Nfp*Nface,surf%Sn)
            read(fid)surf%coord%array
            read(fid)surf%nx
            read(fid)surf%ny
            read(fid)surf%nz
        endif
    endif
end subroutine surface_file

subroutine SRtype_file(fid,SRtype,k_inout)
    integer,intent(in) :: fid
    type(send_recv_pack) :: SRtype
    character(len=*),intent(in) :: k_inout
    integer :: i,j

    if(k_inout .eq. 'save')then
        write(fid)SRtype%host%N_DD_Conn
        if(SRtype%host%N_DD_Conn.gt.0)then
!            print*,SRtype%host%N_DD_Conn,SRtype%host%maxN_Conn
            write(fid)SRtype%host%DD_Conn
            write(fid)SRtype%host%N_Conn
            write(fid)SRtype%host%maxN_Conn
            do i=1,SRtype%host%N_DD_Conn
                do j=1,SRtype%host%maxN_Conn
                    write(fid)SRtype%host%Mask(j,i)
                enddo
            enddo
        endif
        write(fid)SRtype%gost%N_DD_Conn
        if(SRtype%gost%N_DD_Conn.gt.0)then
!            print*,SRtype%gost%N_DD_Conn,SRtype%gost%maxN_Conn
            write(fid)SRtype%gost%DD_Conn
            write(fid)SRtype%gost%N_Conn
            write(fid)SRtype%gost%maxN_Conn
            do i=1,SRtype%gost%N_DD_Conn
                do j=1,SRtype%gost%maxN_Conn
                    write(fid)SRtype%gost%Mask(j,i)
                enddo
            enddo
        endif
    else
        read(fid)SRtype%host%N_DD_Conn
        if(SRtype%host%N_DD_Conn.gt.0)then
            allocate(SRtype%host%DD_Conn &
                    (SRtype%host%N_DD_Conn ))
            allocate(SRtype%host%N_Conn &
                    (SRtype%host%N_DD_Conn ))
            read(fid)SRtype%host%DD_Conn
            read(fid)SRtype%host%N_Conn
            read(fid)SRtype%host%maxN_Conn
            allocate(SRtype%host%Mask&
                    (SRtype%host%maxN_Conn,&
                     SRtype%host%N_DD_Conn))
            read(fid)SRtype%host%Mask 
!            print*,SRtype%host%N_DD_Conn,SRtype%host%maxN_Conn
        endif
        read(fid)SRtype%gost%N_DD_Conn
        if(SRtype%gost%N_DD_Conn.gt.0)then
            allocate(SRtype%gost%DD_Conn&
                    (SRtype%gost%N_DD_Conn))
            allocate(SRtype%gost%N_Conn&
                    (SRtype%gost%N_DD_Conn))
            read(fid)SRtype%gost%DD_Conn
            read(fid)SRtype%gost%N_Conn
            read(fid)SRtype%gost%maxN_Conn
            allocate(SRtype%gost%Mask&
                    (SRtype%gost%maxN_Conn,&
                     SRtype%gost%N_DD_Conn))
            read(fid)SRtype%gost%Mask
!            print*,SRtype%gost%N_DD_Conn,SRtype%gost%maxN_Conn
        endif
    endif

end subroutine SRtype_file

subroutine PMLinfo_file(fid,Nele,pNp,PMLinfo,k_inout)
    integer,intent(in) :: fid,Nele,pNp
    type(PML_geometry) :: PMLinfo
    character(len=*),intent(in) :: k_inout
    integer :: Ndof

    Ndof=Nele*pNp
    if(k_inout .eq. 'save')then
        write(fid)PMLinfo%damp%array
    else
        call ini_vector(Ndof,PMLinfo%damp)
        read(fid)PMLinfo%damp%array
    endif
end subroutine PMLinfo_file

subroutine material_file(fid,Nele,pNp,material,k_inout)
    integer,intent(in) :: fid,Nele,pNp
    type(tetmesh_material)  :: material
    character(len=*),intent(in) :: k_inout
    integer :: Ndof

    Ndof=Nele*pNp
    if(k_inout .eq. 'save')then
        write(fid)material%job
        write(fid)material%rho
        if(job.eq.0)then
            write(fid)material%C(:,1)
        elseif(job.eq.1)then
            write(fid)material%C(:,1:2)
        elseif(job.eq.2)then
            write(fid)material%C(:,1:9)
        else
            write(fid)material%C(:,1:21)
        endif
        write(fid)material%k_media
        if(prestress.gt.0)write(fid)material%T0%array
    else
        read(fid)material%job
        allocate(material%rho    (Ndof))
        allocate(material%k_media(Nele))
        if(job.eq.0)then
            allocate(material%C(Ndof,1 ))
        elseif(job.eq.1)then
            allocate(material%C(Ndof,2 ))
        elseif(job.eq.2)then
            allocate(material%C(Ndof,9 ))
        else
            allocate(material%C(Ndof,21))
        endif
        read(fid)material%rho
        read(fid)material%C
        read(fid)material%k_media
        if(prestress.gt.0)then
            call ini_tensor(Ndof,material%T0,.true.)
!            read(fid)material%T0%array
            material%T0%xx=refer_T0(1)
            material%T0%yy=refer_T0(2)
            material%T0%zz=refer_T0(3)
            material%T0%yz=refer_T0(4)
            material%T0%xz=refer_T0(5)
            material%T0%xy=refer_T0(6)
        endif
    endif

end subroutine material_file

subroutine rupture_file(fid,rupt,k_inout)
    integer,intent(in) :: fid
    type(rupture) :: rupt
    character(len=*),intent(in) :: k_inout
    integer :: k

    if(k_inout .eq. 'save')then
        write(fid)rupt%globNface
        write(fid)rupt%Nface
        write(fid)rupt%Nhface
        if(rupt%Nface.gt.0)then
            write(fid)rupt%globID
            write(fid)rupt%coord%array
            write(fid)rupt%a
!            write(fid)rupt%sigma0
            write(fid)rupt%b
            write(fid)rupt%L
            write(fid)rupt%V0
            write(fid)rupt%f0
            write(fid)rupt%dVt
            write(fid)rupt%perm
            write(fid)rupt%T2E
            call SRtype_file(fid,rupt%SRtype,k_inout)
            write(fid)rupt%nmap
            write(fid)size(rupt%weigh)
            write(fid)rupt%weigh
        endif
    else
        read(fid)rupt%globNface
        read(fid)rupt%Nface
        read(fid)rupt%Nhface
        if(verbose .and. logfid.gt.0)then
            write(logfid,*)'rupt Nface=',rupt%Nface
            write(logfid,*)'rupt Nhface=',rupt%Nhface
        endif
        if(rupt%Nface.gt.0)then
            call rupture_malloc(rupt)
            read(fid)rupt%globID
            read(fid)rupt%coord%array
            read(fid)rupt%a
!            read(fid)rupt%sigma0
            read(fid)rupt%b
            read(fid)rupt%L
            read(fid)rupt%V0
            read(fid)rupt%f0
            read(fid)rupt%dVt
            read(fid)rupt%perm
            read(fid)rupt%T2E
            call SRtype_file(fid,rupt%SRtype,k_inout)
            allocate(rupt%nmap(Nfp*rupt%Nface))
            read(fid)rupt%nmap
            read(fid)k
            allocate(rupt%avg(k))
            allocate(rupt%weigh(k))
            read(fid)rupt%weigh

            call init_SRtype(rupt%SRtype,Nfp,rupt%Nface,1)

        endif
    endif

end subroutine rupture_file

subroutine time_rupture_file(fid,rupt,k_inout)
    integer,intent(in) :: fid
    type(rupture) :: rupt
    character(len=*),intent(in) :: k_inout

    if(k_inout .eq. 'save')then
        write(fid)rupt%sigma
        write(fid)rupt%tauf%array
        write(fid)rupt%Vt%array
        write(fid)rupt%dVt
        write(fid)rupt%f
        write(fid)rupt%psi
        write(fid)rupt%crack_t
        write(fid)rupt%Um%array
        write(fid)rupt%Up%array
        write(fid)rupt%Em%array
        write(fid)rupt%Ep%array
        write(fid)rupt%tau%array
    else
        read(fid)rupt%sigma
        read(fid)rupt%tauf%array
        read(fid)rupt%Vt%array
        read(fid)rupt%f
        read(fid)rupt%f
        read(fid)rupt%psi
        read(fid)rupt%crack_t
        read(fid)rupt%Um%array
        read(fid)rupt%Up%array
        read(fid)rupt%Em%array
        read(fid)rupt%Ep%array
    endif

end subroutine time_rupture_file

subroutine subdomain_ckpfname(fname,fid,snapfname,varname,&
        pid,ickp,ftype,alive)
    integer,intent(out):: fid
    logical,intent(out):: alive
    character(len=*),intent(out) :: fname
    character(len=*),intent(in)  :: snapfname,varname,ftype
    integer,intent(in) :: pid,ickp
    character(len=6)   :: cpid,cckp

    write(cckp,"(I6)")ickp
    cckp=adjustl(cckp)
    write(cpid,"(I6)")pid
    cpid=adjustl(cpid)

    fname=trim(snapfname)//'_'//trim(varname)
    if(ickp.ge.0)fname=trim(fname)//'_'//trim(cckp)
    if( pid.ge.0)fname=trim(fname)//'_'//trim(cpid)
    if(trim(ftype).eq.'txt')then
        fname=trim(fname)//'.txt'
    else
        fname=trim(fname)//'.dat'
    endif

    fid=3190+pid
    inquire(file=trim(fname),exist=alive)

end subroutine subdomain_ckpfname

subroutine vector_file(Vect,Ndof,filename,varname,ickp,pid,k_inout)

    character(len=*),intent(in) :: filename,varname,k_inout
    integer,intent(in)          :: Ndof,ickp,pid
    real(kind=rkind) :: Vect(Ndof)
    character(len=Fnamelen) :: fname
    integer :: fid
    logical :: alive

    call subdomain_ckpfname(fname,fid,filename,&
        varname,pid,ickp,'bin',alive)

    if(k_inout .eq. 'save')then

        if(verbose .and. logfid.gt.0)then
            write(logfid,*)'Writing into file:',trim(fname)
            write(logfid,*)'Filesize: ',Ndof*rkind
        endif
            
        open(unit=fid,file=trim(fname),access='stream',&
            status='replace')
        write(fid) Vect(1:Ndof)
        close(fid)

    elseif(k_inout .eq. 'load')then

        if(.not. alive)then
            print*,'file "',trim(fname),' " does not exist.'
            stop
        endif
        if(verbose .and. logfid.gt.0)then
            write(logfid,*)'Reading from file:',trim(fname)
        endif
            
        open(unit=fid,file=trim(fname),access='stream',&
        status='old',position='rewind')
        read(fid) Vect(1:Ndof)
        close(fid)

    endif

end subroutine vector_file

end module checkpoint_mod

