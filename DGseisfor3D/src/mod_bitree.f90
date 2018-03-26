module bitree_mod

implicit none

    type bitree
        integer :: Nnod=-1
        double precision,allocatable :: node(:,:)
        integer,allocatable :: lnod(:),rnod(:),pnod(:)
    end type bitree
    type(bitree) :: tree
    contains

    subroutine bitree_init(Ndim,Nvar)
        integer :: Ndim,Nvar
        allocate(tree%node(Ndim,Nvar))
        allocate(tree%lnod(Ndim));tree%lnod=-1
        allocate(tree%rnod(Ndim));tree%rnod=-1
        allocate(tree%pnod(Ndim));tree%pnod=-1
        tree%Nnod=0
    end subroutine bitree_init

    subroutine bitree_del(Nnod)
        integer :: Nnod
        deallocate(tree%node)
        deallocate(tree%lnod)
        deallocate(tree%rnod)
        deallocate(tree%pnod)
        Nnod=tree%Nnod
        tree%Nnod=-1
    end subroutine bitree_del

    subroutine bitree_add(var,Nvar,id,tol)
        integer :: pt,id,ivar,Nvar
        double precision :: var(Nvar),tol,lv
        if(tree%Nnod.eq.0)then
            do ivar=1,Nvar
                tree%node(1,ivar)=var(ivar)
            enddo
            tree%Nnod=1
            tree%pnod(1)=0
            id=1
            return
        endif
        pt=1;ivar=1;lv=tree%node(pt,ivar)
        do while(.true.)
            if(abs(var(ivar)-lv).le.tol)then
                if(ivar.eq.Nvar)then
                    id=pt;exit
                else
                    ivar=ivar+1
                    lv=tree%node(pt,ivar)
                    cycle
                endif
            elseif(var(ivar)+tol.lt.lv)then
                if(tree%lnod(pt).lt.0)then
                    tree%Nnod=tree%Nnod+1
                    id=tree%Nnod
                    tree%lnod(pt)=id
                    do ivar=1,Nvar
                        tree%node(id,ivar)=var(ivar)
                    enddo
                    tree%pnod(id)=pt
                    exit
                else
                    pt=tree%lnod(pt)
                    lv=tree%node(pt,ivar)
                endif
            else
                if(tree%rnod(pt).lt.0)then
                    tree%Nnod=tree%Nnod+1
                    id=tree%Nnod
                    tree%rnod(pt)=id
                    do ivar=1,Nvar
                        tree%node(id,ivar)=var(ivar)
                    enddo
                    tree%pnod(id)=pt
                    exit
                else
                    pt=tree%rnod(pt)
                    lv=tree%node(pt,ivar)
                endif
            endif
        enddo
    end subroutine bitree_add

    subroutine bitree_find(var,Nvar,id,tol)
        integer :: pt,id,ivar,Nvar
        double precision :: var(Nvar),tol,lv
        if(tree%Nnod.eq.0)then
            id=-1
            return
        endif
        pt=1;ivar=1;lv=tree%node(pt,ivar)
        do while(.true.)
            if(abs(var(ivar)-lv).le.tol)then
                if(ivar.eq.Nvar)then
                    id=pt;exit
                else
                    ivar=ivar+1
                    lv=tree%node(pt,ivar)
                    cycle
                endif
            elseif(var(ivar)+tol.lt.lv)then
                if(tree%lnod(pt).lt.0)then
                    id=-1
                    exit
                else
                    pt=tree%lnod(pt)
                    lv=tree%node(pt,ivar)
                endif
            else
                if(tree%rnod(pt).lt.0)then
                    id=-1
                    exit
                else
                    pt=tree%rnod(pt)
                    lv=tree%node(pt,ivar)
                endif
            endif
        enddo
    end subroutine bitree_find

    subroutine bitree_output(var,Ndim,Nvar)
        integer :: Ndim,Nvar,ivar,pt,i,j
        double precision :: var(Ndim,Nvar)
        logical,allocatable :: flag(:)
        if(tree%Nnod.le.0)then
            Ndim=tree%Nnod
            return
        endif
        allocate(flag(tree%Nnod))
        flag=.true.
        pt=1;i=0
        do while(.true.)
            if(tree%lnod(pt).gt.0)then
                j=pt
                pt=tree%lnod(pt)
                tree%lnod(j)=-1
                cycle
            endif
            if(flag(pt))then
                i=i+1
                do ivar=1,Nvar
                    var(i,ivar)=tree%node(pt,ivar)
                enddo
                flag(pt)=.false.
            endif
            if(tree%rnod(pt).gt.0)then
                j=pt
                pt=tree%rnod(pt)
                tree%rnod(j)=-1
                cycle
            endif
            if(tree%pnod(pt).gt.0)then
                pt=tree%pnod(pt)
                cycle
            else
                exit
            endif
        enddo
        call bitree_del(Ndim)
        deallocate(flag)
    end subroutine bitree_output

end module bitree_mod
