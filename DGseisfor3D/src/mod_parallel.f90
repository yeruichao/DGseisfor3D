!********************************************************************!
!*  This module contains parallel-related subroutines.              *!
!*                                                                  *!
!*  Author: Ruichao Ye    Email: Ruichao.ye@gmail.com               *!
!********************************************************************!

!---------------------------------------------------------------------
module parallel_mod
!---------------------------------------------------------------------
    use domain_mod,        only: domain_ele_count,t_domain
    use meshfile_mod,      only: nodes,tet,neigh,F2T,glob_Nface
    use datatype_mod,      only: rkind,&
                                 send_recv_info,send_recv_pack,&
                                 tetmesh_geometry,tetmesh_Domain,&
                                 vector_array,reset_vector
    use para_mod,          only: pNp,Nfp,job,verbose,logfid,Ndomain,&
                                 withrupt
    use Ssort_Module,      only: issort,Ibisearch

    implicit none

    include 'mpif.h'
    integer,parameter :: MDN=200
    character(len=6) :: ctmp
!---------------------------------------------------------------------
    contains
!---------------------------------------------------------------------

subroutine Domain_map3D(glob_Nele,pNp,N_domain,mesh,subdomain)

!input:
    integer :: glob_Nele
    ! global number of elements
    integer :: pNp
    ! number of nodes per element
    integer :: N_domain
    ! number of domains, same as number of processors
!output:
    type(tetmesh_geometry) :: mesh
    type(tetmesh_Domain)   :: subdomain
!pointer:
    integer local_Nele,Nele_this_domain
    integer N_DD_Conn
    integer,pointer :: DD_Conn(:),N_ghost(:),N_host(:)
    integer,pointer :: hostMask(:,:),ghostMask(:,:)
    integer,pointer :: local_tet(:,:),local_neigh(:,:)
    ! subdomain%{local_tet;local_neigh}
    integer,pointer :: local_globID(:)
    ! mesh%globID
!auxilary:
    logical,allocatable :: inele(:)
    integer,allocatable :: globID(:)
    integer,allocatable :: flag(:),Dcounter(:)
    integer :: Nghost,Nhost,maxN_ghost,maxN_host
    integer :: i,sk,skg
    integer :: j,l,f,fj,itmp,itmp2,key
    integer :: ierr
    integer :: dest,counter,tag,k,statu(MPI_STATUS_SIZE,MDN),ierror
    integer :: nid,ndomn

    allocate(flag(glob_Nele),globID(glob_Nele))
    allocate(Dcounter(N_domain))
    Nele_this_domain=0;Dcounter=0;Nghost=0;Nhost=1;flag=0
    do i=1,glob_Nele
        globID(i)=i
    enddo
    do i=1,glob_Nele
        if(t_domain(i).ne.subdomain%id)cycle
        Nele_this_domain=Nele_this_domain+1
        key=0

        do f=1,4
            nid=neigh(f,i)
            if(nid.le.0)cycle
            ndomn=t_domain(nid)
            if(ndomn.ne.subdomain%id)then
                if(flag(nid).eq.0)Dcounter(ndomn)=Dcounter(ndomn)+1
                flag(nid)=1
                if(key.eq.0)then
                    Nhost=Nhost+1
                    key=1
                endif
            endif
        enddo

    enddo

    N_DD_Conn=count(Dcounter .gt. 0)

    allocate(DD_Conn(N_DD_Conn))
    allocate(N_ghost(N_DD_Conn))
    allocate(N_host (N_DD_Conn)) 
    DD_Conn=0;N_ghost=0;N_host=0
    Nghost=maxval(Dcounter)
    itmp=1
    do i=1,N_domain
        if(Dcounter(i).gt.0)then
            if(itmp.gt.N_DD_Conn)&
                print*,'Error: Too many neighbour domains.'
            DD_Conn(itmp)=i
            itmp=itmp+1
        endif
    enddo
    if(itmp-1.ne.N_DD_Conn)&
                print*,'Error: neighbour domain count error.'

    do i=1,N_DD_Conn
        if(DD_Conn(i).le.0)exit
        N_ghost(i)=Dcounter(DD_Conn(i))
    enddo

    allocate(ghostMask(Nghost,N_DD_Conn))
    allocate( hostMask( Nhost,N_DD_Conn));

    if(Nele_this_domain.ne.domain_ele_count(subdomain%id))&
        print*,'Something wrong with Nele this domain.'
    local_Nele=Nele_this_domain+sum(Dcounter)

    if(verbose .and. logfid.gt.0)then
        write(logfid,*)'glob_Nele=',glob_Nele,'local_Nele=',local_Nele
        write(logfid,*)'Nele_this_domain=',Nele_this_domain
    endif

    allocate(subdomain%local_tet(4,local_Nele))
        local_tet=>subdomain%local_tet
    allocate(subdomain%local_neigh(4,local_Nele))
        local_neigh=>subdomain%local_neigh
    allocate(mesh%globID(local_Nele))
        local_globID=>mesh%globID
    allocate(inele(Nele_this_domain));inele=.true.

    sk=0;skg=Nele_this_domain;flag=0
    do i=1,glob_Nele
        if(t_domain(i).eq.subdomain%id)then
            sk=sk+1
            local_tet(:,sk)=tet(:,i)
            local_neigh(:,sk)=neigh(:,i)
            local_globID(sk)=i
        do f=1,4
            if(neigh(f,i).le.0)cycle
            itmp=neigh(f,i)
            if(t_domain(itmp).ne.subdomain%id &
            .and. flag(itmp).eq.0)then
                flag(itmp)=1
                skg=skg+1
                local_tet(:,skg)=tet(:,itmp)
                local_globID(skg)=globID(itmp)
                local_neigh(:,skg)=neigh(:,itmp)
            endif
        enddo

        endif
    enddo

    if(sk.ne.Nele_this_domain)print*,'Something wrong with host.'
    local_Nele=skg

    mesh%Nele=local_Nele;mesh%Nhele=Nele_this_domain
    
    ! Reorder ghost elements according to their global IDs
    call partial_reorder(local_globID(sk+1:skg),&
                 local_tet(:,sk+1:skg),&
                 local_neigh(:,sk+1:skg),skg-sk)

    ! Getting map: GlobalID->localID (Stored in globID(:)) 
    ! so as to configure local_neigh
    ! If connected to ele not belong to this domain, neigh=-1
    globID=-1 
    do i=1,local_Nele
        globID(local_globID(i))=i
    enddo

    ! Also forming hostMask and ghostMask simutaneously
    ! But I need map DomainID->DD_Conn. I store it in Dcounter
    Dcounter=0
    do i=1,N_DD_Conn
        Dcounter(DD_Conn(i))=i
    enddo

    N_host=0;N_ghost=0;hostMask=0;ghostMask=0
    do i=1,local_Nele
        do f=1,4
            if(local_neigh(f,i).le.0)cycle
            itmp=t_domain(local_neigh(f,i))
            local_neigh(f,i)=globID(local_neigh(f,i))
            if(i.le.Nele_this_domain .and. &
                itmp.ne.subdomain%id)then
    ! Is a host element
                itmp=Dcounter(itmp)
                N_host(itmp)=N_host(itmp)+1
                itmp2=N_host(itmp)
                hostMask(itmp2,itmp)=i
                inele(i)=.false.
                if(itmp2.eq.1)cycle
                if(hostMask(itmp2-1,itmp).eq.i)&
                    N_host(itmp)=N_host(itmp)-1
                    ! already in host list
            endif
        enddo
        if(i.gt.Nele_this_domain)then
    ! Is a ghost element
                itmp=t_domain(local_globID(i))
                itmp=Dcounter(itmp)
                N_ghost(itmp)=N_ghost(itmp)+1
                ghostMask(N_ghost(itmp),itmp)=i
        endif
    enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! permute the local elements in such ordering: 
!! interior eles, local hosts, local ghosts

    sk=0
    do i=1,Nele_this_domain
        if(inele(i))then
            sk=sk+1;globID(i)=sk
        endif
    enddo
    mesh%Ninele=sk
    do i=1,Nele_this_domain
        if(.not.(inele(i)))then
            sk=sk+1;globID(i)=sk
        endif
    enddo
!! the original ith element is currently globID(i)th element
    do j=1,N_DD_Conn
        do i=1,N_host(j)
            if(hostMask(i,j).le.0)cycle
            hostMask(i,j)=globID(hostMask(i,j))
        enddo
    enddo
    do i=1,local_Nele
        do j=1,4
            if(local_neigh(j,i).le.0 .or. &
               local_neigh(j,i).gt.Nele_this_domain )cycle
            local_neigh(j,i)=globID(local_neigh(j,i))
        enddo
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! save in Send Receive Type
    subdomain%SRtype%host%N_DD_Conn=N_DD_Conn
    subdomain%SRtype%gost%N_DD_Conn=N_DD_Conn
    allocate(subdomain%SRtype%host%DD_Conn(N_DD_Conn))
    allocate(subdomain%SRtype%gost%DD_Conn(N_DD_Conn))
    subdomain%SRtype%host%DD_Conn=DD_Conn
    subdomain%SRtype%gost%DD_Conn=DD_Conn
    allocate(subdomain%SRtype%host%N_Conn(N_DD_Conn))
    allocate(subdomain%SRtype%gost%N_Conn(N_DD_Conn))
    subdomain%SRtype%host%N_Conn=N_host 
    subdomain%SRtype%gost%N_Conn=N_ghost
    maxN_host =maxval(N_host )
    maxN_ghost=maxval(N_ghost)
    subdomain%SRtype%host%maxN_Conn=maxN_host 
    subdomain%SRtype%gost%maxN_Conn=maxN_ghost
    allocate(subdomain%SRtype%host%Mask(maxN_host ,N_DD_Conn))
    allocate(subdomain%SRtype%gost%Mask(maxN_ghost,N_DD_Conn))
    subdomain%SRtype%host%Mask=0
    subdomain%SRtype%gost%Mask=0
    do i=1,N_DD_Conn
        do j=1,maxN_host
            if(hostMask(j,i).le.0)exit
            subdomain%SRtype%host%Mask(j,i)=hostMask(j,i)
        enddo
        do j=1,maxN_ghost
            if(ghostMask(j,i).le.0)exit
            subdomain%SRtype%gost%Mask(j,i)=ghostMask(j,i)
        enddo
    enddo

    subdomain%SRtype1%host%N_DD_Conn=N_DD_Conn
    subdomain%SRtype1%gost%N_DD_Conn=N_DD_Conn
    allocate(subdomain%SRtype1%host%DD_Conn(N_DD_Conn))
    allocate(subdomain%SRtype1%gost%DD_Conn(N_DD_Conn))
    subdomain%SRtype1%host%DD_Conn=DD_Conn
    subdomain%SRtype1%gost%DD_Conn=DD_Conn
    allocate(subdomain%SRtype1%host%N_Conn(N_DD_Conn))
    allocate(subdomain%SRtype1%gost%N_Conn(N_DD_Conn))
    subdomain%SRtype1%host%N_Conn=N_host 
    subdomain%SRtype1%gost%N_Conn=N_ghost
    maxN_host =maxval(N_host )
    maxN_ghost=maxval(N_ghost)
    subdomain%SRtype1%host%maxN_Conn=maxN_host 
    subdomain%SRtype1%gost%maxN_Conn=maxN_ghost
    allocate(subdomain%SRtype1%host%Mask(maxN_host ,N_DD_Conn))
    allocate(subdomain%SRtype1%gost%Mask(maxN_ghost,N_DD_Conn))
    subdomain%SRtype1%host%Mask=0
    subdomain%SRtype1%gost%Mask=0
    do i=1,N_DD_Conn
        do j=1,maxN_host
            if(hostMask(j,i).le.0)exit
            subdomain%SRtype1%host%Mask(j,i)=hostMask(j,i)
        enddo
        do j=1,maxN_ghost
            if(ghostMask(j,i).le.0)exit
            subdomain%SRtype1%gost%Mask(j,i)=ghostMask(j,i)
        enddo
    enddo
    call init_SRtype(subdomain%SRtype1,1,mesh%Nele,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! save local-to-global map

    sk=0
    do i=1,Nele_this_domain
        if(inele(i))then
            sk=sk+1;globID(sk)=i
        endif
    enddo
    do i=1,Nele_this_domain
        if(.not.(inele(i)))then
            sk=sk+1;globID(sk)=i
        endif
    enddo
    do i=Nele_this_domain+1,local_Nele
        globID(i)=i
    enddo

!! the current ith element is originally globID(i)th element
!      local_tet(:,:)=local_tet(:,globID(1:local_Nele))
!    local_neigh(:,:)=local_neigh(:,globID(1:local_Nele))
!     local_globID(:)=local_globID(globID(1:local_Nele))
    do j=1,4
        do i=1,local_Nele
            flag(i)=local_tet(j,i)
        enddo
        do i=1,local_Nele
            local_tet(j,i)=flag(globID(i))
        enddo
        do i=1,local_Nele
            flag(i)=local_neigh(j,i)
        enddo
        do i=1,local_Nele
            local_neigh(j,i)=flag(globID(i))
        enddo
    enddo
    do i=1,local_Nele
        flag(i)=local_globID(i)
    enddo
    do i=1,local_Nele
        local_globID(i)=flag(globID(i))
    enddo

    deallocate(inele)

    globID=-1
    do i=1,local_Nele
        globID(local_globID(i))=i
    enddo
!!!! transform global T2E to local T2E
    if(withrupt)then
        do i=1,glob_Nface
            do j=1,2
                if(F2T(j,i).gt.0)then
                    F2T(j,i)=globID(F2T(j,i))
                else
                    F2T(j,i)=-1
                endif
                if(F2T(1,i).lt.0 .and. &
!                        F2T(2,i).gt.Nele_this_domain)then
                        F2T(2,i).gt.local_Nele)then
                    F2T(2,i)=-1
                elseif(F2T(2,i).lt.0 .and. &
!                        F2T(1,i).gt.Nele_this_domain)then
                        F2T(1,i).gt.local_Nele)then
                    F2T(1,i)=-1
                endif
            enddo
        enddo
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    deallocate(flag,globID,Dcounter)

    if(verbose .and. logfid.gt.0)then
        write(logfid,*)'Ninele=',mesh%Ninele
    endif

end subroutine Domain_map3D

subroutine init_SRtype(SRtype,Np,Nele,Ncomp)
    type(send_recv_pack) :: SRtype
    integer,intent(in) :: Np,Nele,Ncomp
    integer,allocatable  :: disp(:)
    integer,allocatable  :: block(:)
    integer :: i,j,k,l,m,n
    integer :: maxsize
    integer :: ierr

    maxsize=max(maxval(SRtype%host%N_Conn),maxval(SRtype%gost%N_Conn))
    allocate(SRtype%SStype(SRtype%host%N_DD_Conn))
    allocate(SRtype%RRtype(SRtype%gost%N_DD_Conn))
    allocate(disp (maxsize*Ncomp))
    allocate(block(maxsize*Ncomp))
    block=Np
    do j=1, SRtype%host%N_DD_Conn
        i = SRtype%host%N_Conn(j)
        if(i.le.0)cycle
        n=0
        do m=1,Ncomp
            do k=1,i
                l=SRtype%host%Mask(k,j)
                n=n+1
                disp(n)=(l-1)*Np+(m-1)*Nele*Np
            enddo
        enddo
        i=i*Ncomp
        call MPI_TYPE_INDEXED(i,block(1:i),disp(1:i),&
            MPI_DOUBLE_PRECISION,SRtype%SStype(j),ierr)
        call MPI_TYPE_COMMIT(SRtype%SStype(j),ierr)
    enddo
    do j=1, SRtype%gost%N_DD_Conn
        i = SRtype%gost%N_Conn(j)
        if(i.le.0)cycle
        n=0
        do m=1,Ncomp
            do k=1,i
                l=SRtype%gost%Mask(k,j)
                n=n+1
                disp(n)=(l-1)*Np+(m-1)*Nele*Np
            enddo
        enddo
        i=i*Ncomp
        call MPI_TYPE_INDEXED(i,block(1:i),disp(1:i),&
            MPI_DOUBLE_PRECISION,SRtype%RRtype(j),ierr)
        call MPI_TYPE_COMMIT(SRtype%RRtype(j),ierr)
    enddo
    deallocate(disp,block)
!    if(SRtype%host%N_DD_Conn .gt.0)&
    allocate(SRtype%send_req(SRtype%host%N_DD_Conn ))
!    if(SRtype%gost%N_DD_Conn .gt.0)&
    allocate(SRtype%recv_req(SRtype%gost%N_DD_Conn))
    i=max(SRtype%host%N_DD_Conn,SRtype%gost%N_DD_Conn)
    if(i.gt.0)allocate(SRtype%sta(MPI_STATUS_SIZE,i))
end subroutine init_SRtype

subroutine free_SRtype(SRtype)
    type(send_recv_pack) :: SRtype
    integer :: j,ierr
    if(SRtype%host%N_DD_Conn.gt.0)then
        do j=1,SRtype%host%N_DD_Conn
            call MPI_TYPE_FREE(SRtype%SStype(j),ierr)
        enddo
        deallocate(SRtype%host%N_Conn)
        deallocate(SRtype%host%DD_Conn)
        deallocate(SRtype%host%Mask)
        SRtype%host%N_DD_Conn=0
    endif
    if(SRtype%gost%N_DD_Conn.gt.0)then
        do j=1,SRtype%gost%N_DD_Conn
            call MPI_TYPE_FREE(SRtype%RRtype(j),ierr)
        enddo
        deallocate(SRtype%gost%DD_Conn)
        deallocate(SRtype%gost%N_Conn)
        deallocate(SRtype%gost%Mask)
        SRtype%gost%N_DD_Conn=0
    endif
end subroutine free_SRtype

subroutine sample_SRtype_info(flag,glb_host,loc_host)
    type(send_recv_info) :: glb_host,loc_host
    logical,intent(in) :: flag(:)
    integer,allocatable :: DD_Conn(:),N_conn(:),Mask(:,:)
    integer :: i,j,k,m,n
    m=glb_host%N_DD_Conn
    n=glb_host%maxN_Conn
    allocate(DD_Conn(m),N_Conn(m),Mask(n,m))
    Mask=0
    N_conn=0
    do i=1,m
        DD_conn(i)=glb_host%DD_conn(i)
        do j=1,glb_host%N_conn(i)
            k=glb_host%Mask(j,i)
            if(flag(k))then ! this element is sampled
                N_conn(i)=N_conn(i)+1
                Mask(N_conn(i),i)=k
            endif
        enddo
    enddo
    k=0
    do i=1,m
        if(N_conn(i).gt.0)k=k+1
    enddo
    loc_host%N_DD_Conn=k
    n=maxval(N_conn)
    loc_host%maxN_Conn=n
    if(k.gt.0)then
        allocate(loc_host%DD_Conn(k))
        allocate(loc_host%N_Conn(k))
        allocate(loc_host%Mask(n,k))
        loc_host%Mask=0
        k=0
        do i=1,m
            if(N_conn(i).gt.0)then
                k=k+1
                loc_host%DD_Conn(k)=DD_conn(i)
                loc_host%N_Conn(k)=N_conn(i)
                do j=1,N_conn(i)
                    loc_host%Mask(j,k)=Mask(j,i)
                enddo
            endif
        enddo
    endif
end subroutine sample_SRtype_info

subroutine copy_SRtype_info(host1,host2)
    type(send_recv_info) :: host1,host2
    integer :: m,n
    m=host1%N_DD_Conn
    n=host1%maxN_Conn
    host2%N_DD_Conn=m
    host2%maxN_Conn=n
    if(m.le.0)return
    allocate(host2%DD_Conn(m))
    allocate(host2%N_Conn(m))
    allocate(host2%Mask(n,m))
    host2%DD_Conn=host1%DD_Conn
    host2%N_Conn =host1%N_Conn 
    host2%Mask   =host1%Mask
end subroutine copy_SRtype_info

subroutine subdomain_reduce_variable(subdomain)
    type(tetmesh_Domain)   :: subdomain

    deallocate(subdomain%local_tet)
    deallocate(subdomain%local_neigh)
end subroutine subdomain_reduce_variable

subroutine partial_reorder(lglobID,ltet,lneigh,N)
    integer :: N,ltet(4,N),lneigh(4,N),i
    integer,allocatable :: itmp(:)
    integer :: lglobID(N)
    allocate(itmp(N))
    do i=1,N
        itmp(i)=i
    enddo
    call issort(lglobID,itmp)
    ltet(:,:)=ltet(:,itmp)
    lneigh(:,:)=lneigh(:,itmp)
    deallocate(itmp)
end subroutine partial_reorder

subroutine check_map(subdomain,SRtype,V,Vt,ndof,nhdof,ncomp)
    type(tetmesh_Domain)   :: subdomain
    type(send_recv_pack) :: SRtype
    integer :: ndof,nhdof,ncomp
    real(kind=rkind) :: V(1),Vt(1)
    integer:: i,j,k,l,m,ierr
    real(kind=rkind) :: rtmp,rerror
    integer :: N_DD_Conn

    do j=1,ncomp
        k=(j-1)*ndof
        do i=1,nhdof
            k=k+1
            Vt(k)=V(k)
        enddo
    enddo

    SRtype%send_req=MPI_REQUEST_NULL
    SRtype%recv_req=MPI_REQUEST_NULL
    call domain_exchange(subdomain,SRtype,Vt,&
            SRtype%send_req,SRtype%recv_req,100)
    N_DD_Conn=SRtype%host%N_DD_Conn
    if(N_DD_Conn.gt.0)&
    call MPI_Waitall(N_DD_Conn,SRtype%send_req,SRtype%sta,ierr)
    N_DD_Conn=SRtype%gost%N_DD_Conn
    if(N_DD_Conn.gt.0)&
    call MPI_Waitall(N_DD_Conn,SRtype%recv_req,SRtype%sta,ierr)

    rerror=0d0
    do j=1,ncomp
        k=(j-1)*ndof+nhdof
        do i=nhdof+1,ndof
            k=k+1
            rerror=rerror+((V(k)-Vt(k))**2)
        enddo
    enddo
    if(rerror.gt.1e-8)&
            print*,'Communication error:',rerror

end subroutine check_map

subroutine domain_exchange(subdomain,SRtype,x,send_req,recv_req,offset)

! inout
    type(tetmesh_Domain)   :: subdomain
    real(kind=rkind)       :: x(:)
! input
    type(send_recv_pack) :: SRtype
    integer :: send_req(:),recv_req(:)
    integer :: pid,poffset,offset
! auxilary
    integer :: i,j,k,l,ierr,dest,tag

!!!!!!!!!!! reveiving message !!!!!!!!!!!!

    poffset=subdomain%pid-subdomain%id
    pid=subdomain%pid

    do j=1,SRtype%gost%N_DD_Conn
        dest=SRtype%gost%DD_Conn(j)+poffset
        tag=dest+offset
        i=SRtype%gost%N_Conn(j)
        if(i.le.0)cycle
        call MPI_IRECV(x,1,SRtype%RRtype(j),&
            dest,tag,MPI_COMM_WORLD,&
            recv_req(j),ierr)
    enddo

!!!!!!!!!!!! sending message !!!!!!!!!!!

    do j=1,SRtype%host%N_DD_Conn
        dest=SRtype%host%DD_Conn(j)+poffset
        tag=pid+offset
        i=SRtype%host%N_Conn(j)
        if(i.le.0)cycle
        call MPI_ISEND(x,1,SRtype%SStype(j),&
            dest,tag,MPI_COMM_WORLD,&
            send_req(j),ierr)
    enddo

end subroutine domain_exchange

subroutine init_rupt_exchange_type(mesh,subdomain,rupt_SRtype)
    type(tetmesh_geometry) :: mesh
    type(tetmesh_Domain)   :: subdomain
    type(send_recv_pack) :: rupt_SRtype
    integer,allocatable :: ID(:),tID(:),Nghost(:),Nhost(:),&
               DD_Conn_ghost(:),DD_Conn_host(:)
    integer :: i,j,k,pid,Nface,Nhface,Nlface,ierr,dest,tag,NN
    integer :: Ng_DD_Conn,Nh_DD_Conn,N_ghost,N_host
    integer :: send_req(Ndomain),recv_req(Ndomain),&
               statu(MPI_STATUS_SIZE,Ndomain)
    integer,allocatable :: hostMask(:,:),ghostMask(:,:)
    Nface=mesh%rupt%globNface
    Nlface=mesh%rupt%Nface
    Nhface=mesh%rupt%Nhface
    allocate( ID(Nface))
    allocate(tID(Nface))
    pid=subdomain%id-1
    ID=-1;tID=-1
    if(Nhface.gt.0)then
        do i=1,Nhface
            ID(mesh%rupt%globID(i))=pid
        enddo
    endif
    call MPI_ALLREDUCE(ID,tID,Nface,MPI_INTEGER,MPI_MAX,&
            MPI_COMM_WORLD,ierr)
    if(minval(tID).lt.0)print*,'tid error'
    allocate(Nghost(Ndomain))
    Nghost=0;ID=0
    do i=Nhface+1,Nlface
        j=tID(mesh%rupt%globID(i))+1
        Nghost(j)=Nghost(j)+1
        ID(i-Nhface)=j
    enddo
    Ng_DD_Conn=0
    do i=1,Ndomain
        if(Nghost(i).gt.0)Ng_DD_Conn=Ng_DD_Conn+1
    enddo
    rupt_SRtype%gost%N_DD_Conn=Ng_DD_Conn
    if(Ng_DD_Conn.gt.0)then
        allocate(rupt_SRtype%gost%DD_Conn(Ng_DD_Conn))
        allocate(rupt_SRtype%gost%N_Conn(Ng_DD_Conn))
        allocate(DD_Conn_ghost(Ndomain))
        j=0;DD_Conn_ghost=0
        do i=1,Ndomain
            if(Nghost(i).gt.0)then
                j=j+1;DD_Conn_ghost(i)=j
                rupt_SRtype%gost%DD_Conn(j)=i
                rupt_SRtype%gost%N_Conn(j)=Nghost(i)
            endif
        enddo
        N_ghost=maxval(Nghost)
        rupt_SRtype%gost%maxN_Conn=N_ghost
        allocate(rupt_SRtype%gost%Mask(N_ghost,Ng_DD_Conn))
        Nghost=0
        do i=Nhface+1,Nlface
            j=tID(mesh%rupt%globID(i))+1
            Nghost(j)=Nghost(j)+1
            rupt_SRtype%gost%Mask(Nghost(j),DD_Conn_ghost(j))=i
        enddo
        allocate(ghostMask(N_ghost,Ng_DD_Conn));ghostMask=-1
    endif

    allocate(Nhost(Ndomain));Nhost=0
    send_req=MPI_REQUEST_NULL;recv_req=MPI_REQUEST_NULL
    do i=1,Ndomain
        if(i-1.eq.pid)cycle
        call MPI_IRECV(Nhost(i),1,MPI_INTEGER,&
            i-1,i-1,MPI_COMM_WORLD,recv_req(i),ierr)
        call MPI_ISEND(Nghost(i),1,MPI_INTEGER,&
            i-1,pid,MPI_COMM_WORLD,send_req(i),ierr)
    enddo
    call MPI_Waitall(Ndomain,recv_req,statu,ierr)
    call MPI_Waitall(Ndomain,send_req,statu,ierr)
    Nh_DD_Conn=0
    do i=1,Ndomain
        if(Nhost(i).gt.0)Nh_DD_Conn=Nh_DD_Conn+1
    enddo
    rupt_SRtype%host%N_DD_Conn=Nh_DD_Conn
    if(Nh_DD_Conn.gt.0)then
        allocate(rupt_SRtype%host%DD_Conn(Nh_DD_Conn))
        allocate(rupt_SRtype%host%N_Conn(Nh_DD_Conn))
        allocate(DD_Conn_host(Ndomain))
        j=0;DD_Conn_host=0
        do i=1,Ndomain
            if(Nhost(i).gt.0)then
                j=j+1;DD_Conn_host(i)=j
                rupt_SRtype%host%DD_Conn(j)=i
                rupt_SRtype%host%N_Conn(j)=Nhost(i)
            endif
        enddo
        N_host=maxval(Nhost)
        rupt_SRtype%host%maxN_Conn=N_host
        allocate(hostMask(N_host,Nh_DD_Conn));hostMask=-1
        allocate(rupt_SRtype%host%Mask(N_host,Nh_DD_Conn))
    endif
!    print*,N_host,N_ghost
    send_req=MPI_REQUEST_NULL;recv_req=MPI_REQUEST_NULL
!    print*,N_ghost,N_host
    do i=1,Nh_DD_Conn
        NN=rupt_SRtype%host%N_Conn(i)
        dest=rupt_SRtype%host%DD_Conn(i)-1
        tag=dest
        call MPI_IRECV(hostMask(1:NN,i),NN,MPI_INTEGER,dest,tag,&
                MPI_COMM_WORLD,recv_req(i),ierr)
    enddo
    do i=1,Ng_DD_Conn
        NN=rupt_SRtype%gost%N_Conn(i)
        do j=1,NN
            ghostMask(j,i)=mesh%rupt%globID(&
                rupt_SRtype%gost%Mask(j,i))
        enddo
        dest=rupt_SRtype%gost%DD_Conn(i)-1
        tag=pid
        call MPI_ISEND(ghostMask(1:NN,i),NN,MPI_INTEGER,dest,tag,&
                MPI_COMM_WORLD,send_req(i),ierr)
    enddo
    call MPI_Waitall(Ng_DD_Conn,send_req,statu,ierr)
    call MPI_Waitall(Nh_DD_Conn,recv_req,statu,ierr)
    do i=1,Nh_DD_Conn
        NN=rupt_SRtype%host%N_Conn(i)
        do j=1,NN
            tag=-1
            do k=1,mesh%rupt%Nhface
                if(mesh%rupt%globID(k).eq.hostMask(j,i))then
                        tag=k;exit
                endif
            enddo
            if(tag.le.0)print*,'node not found'
            rupt_SRtype%host%Mask(j,i)=tag
        enddo
    enddo
    call init_SRtype(rupt_SRtype,Nfp,mesh%rupt%Nface,1)
end subroutine init_rupt_exchange_type

function check_blowup(subdomain,x)
    logical :: check_blowup,local_blowup
    type(tetmesh_Domain) :: subdomain
    real(kind=rkind) :: x(:)
    real(kind=rkind),parameter :: tol=1d5
    integer :: ierr

    local_blowup=any(x.gt.tol) .or. any(x.lt.-tol) .or. any(x.ne.x)
    ! check inf, -inf and nan
    call MPI_ALLREDUCE(local_blowup,check_blowup,1,MPI_LOGICAL,&
        MPI_LOR, MPI_COMM_WORLD,ierr)
    return
end function check_blowup

subroutine My_barrier1(pid)
    integer,intent(in) :: pid
    integer :: sta(MPI_STATUS_SIZE),ierr
    integer :: i,j
    if(pid.gt.0 .and. pid.le.Ndomain-1)&
        call MPI_RECV(i,1,MPI_INTEGER,pid-1,pid-1,&
            MPI_COMM_WORLD,sta,ierr)
end subroutine My_barrier1

subroutine My_barrier2(pid)
    integer,intent(in) :: pid
    integer :: sta(MPI_STATUS_SIZE),ierr
    integer :: i,j
    if(pid.lt.Ndomain-1)&
        call MPI_SEND(i,1,MPI_INTEGER,pid+1,pid,&
            MPI_COMM_WORLD,ierr)
end subroutine My_barrier2

end module parallel_mod
