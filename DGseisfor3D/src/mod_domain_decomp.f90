!*******************************************************************!
!*  This module balancely decompose physical domain and rebuild    *!
!*  connection                                                     *!
!*                                                                 *!
!*  Author: Ruichao Ye    Email: Ruichao.ye@gmail.com              *!
!*******************************************************************!

!--------------------------------------------------------------------
module domain_mod
!--------------------------------------------------------------------
    use meshfile_mod, only: nodes,tet,neigh,glob_Nele
    use Ssort_Module, only: ssort
    use para_mod,     only: basename

    implicit none
    public :: domain_ele_count, t_domain,&
              domain_decomp3D, metis_decomp3D

    integer, allocatable :: t_domain(:)
    integer,allocatable :: domain_ele_count(:)
!--------------------------------------------------------------------
    contains
!--------------------------------------------------------------------

subroutine metis_decomp3D(fid,order)

    integer :: order,fid,error
    integer :: i

    allocate(t_domain(glob_Nele))

    do i=1,glob_Nele
            read(fid,*,iostat=error)t_domain(i)
    enddo

    t_domain=t_domain+1

    allocate(domain_ele_count(order))
    domain_ele_count=0

    do i=1,glob_Nele
            domain_ele_count(t_domain(i))=&
                    domain_ele_count(t_domain(i))+1
    enddo

end subroutine metis_decomp3D

subroutine domain_decomp3D(order,k_method)

    ! This subroutine do domain decomposition
    ! input: node, tet, neigh
    ! output: domain_ele_count, t_domain

    real,allocatable :: tx(:),ty(:),tz(:)
    integer,allocatable :: globID(:),idtmp(:)
    integer order
    character(len=*) :: k_method
    integer :: i
!    integer pid,ierr
!    include "mpif.h"

    allocate(tx(glob_Nele))
    allocate(ty(glob_Nele))
    allocate(tz(glob_Nele))
    allocate(globID(glob_Nele))
    allocate(idtmp(glob_Nele))
    allocate(domain_ele_count(2**order))
    allocate(t_domain(glob_Nele))

    do i=1,glob_Nele

            globID(i)=i

    tx(i)=(nodes(1,tet(1,i))+nodes(1,tet(2,i))&
            +nodes(1,tet(3,i))+nodes(1,tet(4,i)))/4;
    ty(i)=(nodes(2,tet(1,i))+nodes(2,tet(2,i))&
            +nodes(2,tet(3,i))+nodes(2,tet(4,i)))/4;
    tz(i)=(nodes(3,tet(1,i))+nodes(3,tet(2,i))&
            +nodes(3,tet(3,i))+nodes(3,tet(4,i)))/4;

            t_domain(i)=1

    enddo

    call domain_sort3D(tet,neigh,tx,ty,tz,globID,t_domain,&
            idtmp,glob_Nele,1,order,k_method)
    do i=1,glob_Nele
            idtmp(globID(i))=i
    enddo
    tet(:,:)=tet(:,idtmp)
    neigh(:,:)=neigh(:,idtmp)
    t_domain(:)=t_domain(idtmp)
    deallocate(tx,ty,tz,globID,idtmp)

!    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    end subroutine domain_decomp3D

    recursive subroutine domain_sort3D(t,n,tx,ty,tz,globID,&
            t_domain,idtmp,Nele,Nth,order,k_method)

    integer :: Nele
    integer :: Nth, order
    integer,intent(in out) :: t(4,Nele),n(4,Nele),&
            globID(Nele)
    integer,intent(in out) :: t_domain(Nele)
    integer :: idtmp(Nele)
    real,intent(in out):: tx(Nele),ty(Nele),tz(Nele)
    integer :: Nele1,Nele2,i
    character(len=*) :: k_method
    character :: kdir

    if(Nth .gt. order)then
        domain_ele_count(t_domain(1))=Nele
        return
    endif


    do i=1,Nele
        idtmp(i)=i
    enddo

    if(len(trim(k_method)).eq.3)then
        if(mod(Nth-1,3).eq.0)then     
            kdir=k_method(1:1)
        elseif(mod(Nth-1,3).eq.1)then 
            kdir=k_method(2:2)
        else                          
            kdir=k_method(3:3)
        endif
    elseif(len(k_method).eq.2)then
        if(mod(Nth-1,2).eq.0)then     
            kdir=k_method(1:1)
        else                          
            kdir=k_method(2:2)
        endif
    elseif(len(k_method).eq.1)then
        kdir=k_method(1:1)
    endif

    if(kdir.eq.'x' .or. kdir.eq.'X')then
        call ssort(tx,idtmp)
        ty=ty(idtmp)
        tz=tz(idtmp)
    elseif(kdir.eq.'y' .or. kdir.eq.'Y')then
        call ssort(ty,idtmp)
        tx=tx(idtmp)
        tz=tz(idtmp)
    else
        call ssort(tz,idtmp)
        tx=tx(idtmp)
        ty=ty(idtmp)
    endif

    Nele1=Nele/2
    Nele2=Nele-Nele1

    t(:,:)=t(:,idtmp)
    n(:,:)=n(:,idtmp)
    globID(:)=globID(idtmp)

    t_domain(Nele1+1:Nele)=t_domain(Nele1+1:Nele)+2**(Nth-1)

    call domain_sort3D(t(:,1:Nele1),n(:,1:Nele1),&
            tx(1:),ty(1:),tz(1:),&
            globID(1:),t_domain(1:),idtmp(1:),&
            Nele1,Nth+1,order,k_method)
    call domain_sort3D(t(:,Nele1+1:Nele),n(:,Nele1+1:Nele),&
            tx(Nele1+1:),ty(Nele1+1:),tz(Nele1+1:),&
            globID(Nele1+1:),t_domain(Nele1+1:),idtmp(Nele1+1:),&
            Nele2,Nth+1,order,k_method)

    return

    end subroutine domain_sort3D

end module domain_mod
