Module Ssort_Module

use datatype_mod, only: rkind

contains

subroutine Ssort(a,ID)

implicit none
real,intent(in out) :: a(:)
integer,intent(in out) :: ID(:)
integer :: N,i,j,k,gap,Na,Nb,itmp
real :: tmp

N=size(a)
k=1;gap=0
do while(gap.lt.N)
    gap=gap*3+1
!    print*,gap
enddo
do while(gap.gt.0)
    do i=gap,N
        tmp=a(i);itmp=ID(i)
        j=i-gap
!print*,'checkpoint',j
        do while(.true.)!(j .gt. 0 .and. a(j).gt.tmp)
            if(j .le. 0)exit
            if(a(j).le.tmp)exit
            a(j+gap)=a(j);ID(j+gap)=ID(j);j=j-gap
        enddo
        a(j+gap)=tmp;ID(j+gap)=itmp
    enddo
    gap=(gap-1)/3
enddo
end subroutine Ssort

subroutine ISsort(a,ID)

integer,intent(in out) :: a(:)
integer,intent(in out) :: ID(:)
integer :: N,k,gap,Na,Nb,itmp,tmp,i,j

N=size(a)
k=1;gap=0
do while(gap.lt.N)
    gap=gap*3+1
enddo
do while(gap.gt.0)
    do i=gap,N
        tmp=a(i);itmp=ID(i)
        j=i-gap
        do while(.true.)!(j .gt. 0 .and. a(j).gt.tmp)
            if(j .le. 0)exit
            if(a(j).le.tmp)exit
            a(j+gap)=a(j);ID(j+gap)=ID(j);j=j-gap
        enddo
        a(j+gap)=tmp;ID(j+gap)=itmp
    enddo
    gap=(gap-1)/3
enddo
end subroutine ISsort

subroutine Ibisearch(N,a,key,id)
    integer :: N,a(N),key,id
    integer :: i,j,m
    id=-1
    i=1;j=N
    do while (i .le. j) 
        m = i + (j - i) / 2
        if (a(m).lt.key)then
            i = m + 1
        elseif (a(m).gt.key)then
            j = m - 1
        else
            id=m;exit
        endif
    enddo
end subroutine Ibisearch

end Module Ssort_Module
