!*******************************************************************!
!*  This module declares and read in parameters variables          *!
!*                                                                 *!
!*  Author: Ruichao Ye    Email: Ruichao.ye@gmail.com              *!
!*******************************************************************!

!--------------------------------------------------------------------
module jacobi_mod
!--------------------------------------------------------------------

    use datatype_mod, only : rkind

!--------------------------------------------------------------------

    implicit none

    public :: JacobiGL,Vandermonde3D,Vandermonde2D

    real(kind=rkind),parameter :: TOL =1d-10

!--------------------------------------------------------------------
    contains
!--------------------------------------------------------------------

subroutine JacobiGL(gaussX,alpha,beta,N)
    integer :: N
    real(kind=rkind) :: alpha,beta,gaussX(N+1),alphat,betat
    real(kind=rkind),allocatable :: J(:,:)
    integer :: i, LDA, INFO, LWORK
    double precision ::  WORK(512)
    gaussX=0.0d0
    gaussX(1)=-1.0d0;gaussX(N+1)=1.0d0
    alphat=alpha+1d0; betat=beta+1d0
    if(N.eq.2)gaussX(2)=-(alphat-betat)/(alphat+betat+2d0)
    if(N.gt.2)then
        allocate(J(N-1,N-1))
        J=0d0
        J(1,1)=(alphat**2-betat**2)/(2d0+alphat+betat)/(alphat+betat)
        if(J(1,1).lt.1e-10)J(1,1)=0.0d0
        do i=1,N-2,1
            J(i+1,i+1)=-0.5d0*(alphat**2-betat**2) &
              /(2d0*dble(i+1)+alphat+betat)&
              /(2d0*dble(i)+alphat+betat)*2d0
            J(i,i+1)=2d0/(2d0*dble(i)+alphat+betat)*sqrt( &
              dble(i)*(dble(i)+alphat+betat)*&
              (dble(i)+alphat)*(dble(i)+betat) &
              /((2d0*dble(i)+alphat+betat)**2-1d0))
        enddo
        LDA=N-1
        LWORK=-1
! Computing eigenvalue of matrix J
        call DSYEV( 'N', 'U', N-1, J, LDA, gaussX(2:N), &
            WORK, LWORK, INFO )
        LWORK=min(512,int(WORK(1)))
        call DSYEV( 'N', 'U', N-1, J, LDA, gaussX(2:N), &
            WORK, LWORK, INFO )
    endif
end subroutine JacobiGL

function combination(n,alpha)
    integer n,alpha,beta,i
    real(kind=rkind) :: combination
    real(kind=rkind)::tmp
    beta=n-alpha
    tmp=1.0d0
    do i=1,beta
        tmp=tmp*dble(alpha+i)/dble(i)
    enddo
    combination=tmp
    return
end function combination

function JacobiP(x,alpha,beta,N)
    integer :: alpha,beta,N,i,hh
    real(kind=rkind) :: JacobiP,x,gamma0,gamma1,a1,a2,a3,p1,p2

    JacobiP=0.0d0
    gamma0=(2.0d0**(alpha+beta+1))/dble(alpha+beta+1) &
        /combination(alpha+beta,alpha)
    gamma1=dble(alpha+1)*dble(beta+1)/dble(alpha+beta+3)*gamma0
    if(N.eq.0)then
        JacobiP=1.0d0/sqrt(gamma0)
    elseif(N.eq.1)then
        JacobiP=(dble(alpha+beta+2)*x/2d0+dble(alpha-beta)/2d0)&
            /sqrt(gamma1)
    else
        p1=1.0d0/sqrt(gamma0)
        p2=(dble(alpha+beta+2)*x/2d0+dble(alpha-beta)/2d0)&
            /sqrt(gamma1)
        a1=2.0/dble(2+alpha+beta)*sqrt(dble(alpha+1)*dble(beta+1) &
            /dble(alpha+beta+3))
        do i=1,N-1,1
            hh=2*i+alpha+beta
            a2=2.0d0/dble(hh+2)*sqrt(dble((i+1)*(i+1+alpha+beta)&
              *(i+1+alpha)*(i+1+beta))/dble((hh+1)*(hh+3)))
            a3=-dble(alpha**2-beta**2)/dble(hh*(hh+2))
            JacobiP=1.0d0/a2*(-a1*p1+(x-a3)*p2)
            a1=a2;p1=p2;p2=JacobiP
        enddo
    endif
    return
end function JacobiP

subroutine Basis3D(p,r,s,t,xdim,i,j,k)
    integer i,j,k,xdim,n
    real(kind=rkind)::p(xdim),r(xdim),s(xdim),t(xdim),a,b,c,h1,h2,h3

    do n=1,xdim
        if(abs(s(n)+t(n)).le.TOL)then
            a=-1.0d0
        else
            a=2d0*(1.0d0+r(n))/(-s(n)-t(n))-1d0
        endif
        if(abs(t(n)-1).le.TOL)then
            b=-1d0
        else
            b=2d0*(1.0d0+s(n))/(1.0d0-t(n))-1d0
        endif
        c=t(n)
        h1=JacobiP(a,0,0,i)
        h2=JacobiP(b,2*i+1,0,j)
        h3=JacobiP(c,2*(i+j)+2,0,k)
        p(n)=2d0*sqrt(2.0d0)*h1&
            *((1.0d0-b)**i)*h2&
            *((1.0d0-c)**(i+j))*h3
    enddo
end subroutine Basis3D

subroutine Vandermonde3D(V3D,N,r,s,t,xdim)
    integer      :: N,xdim,i,j,k,sk
    real(kind=rkind) :: V3D(xdim,(N+1)*(N+2)*(N+3)/6)
    real(kind=rkind) :: r(xdim),s(xdim),t(xdim)
    sk=1;
    do i=0,N,1
        do j=0,N-i,1
            do k=0,N-i-j,1
                call Basis3D(V3D(:,sk),r,s,t,xdim,i,j,k)
                sk=sk+1
            enddo
        enddo
    enddo
end subroutine Vandermonde3D

function dJacobiP(x,alpha,beta,N)
    integer :: alpha,beta,N
    real(kind=rkind) :: dJacobiP,x
    if(N.eq.0)then
        dJacobiP=0.0d0
    else
        dJacobiP=sqrt(dble(N*(N+alpha+beta+1)))&
            *JacobiP(x,alpha+1,beta+1,N-1)
    endif
    return
end function dJacobiP

subroutine GradBasis3D(pr,ps,pt,r,s,t,xdim,i,j,k)
    integer i,j,k,xdim,n
    real(kind=rkind):: pr(xdim),ps(xdim),pt(xdim)
    real(kind=rkind)::  r(xdim), s(xdim), t(xdim)
    real(kind=rkind)::a,b,c,h1,h2,h3,dh1,dh2,dh3,tmp

    do n=1,xdim
        if(abs(s(n)+t(n)).le.TOL)then
            a=-1.0d0
        else
            a=2d0*(1.0d0+r(n))/(-s(n)-t(n))-1d0
        endif
        if(abs(t(n)-1).le.TOL)then
            b=-1d0
        else
            b=2d0*(1.0d0+s(n))/(1.0d0-t(n))-1d0
        endif
        c=t(n)

         h1= JacobiP(a,0,0,i)
        dh1=dJacobiP(a,0,0,i)
         h2= JacobiP(b,2*i+1,0,j)
        dh2=dJacobiP(b,2*i+1,0,j)
         h3= JacobiP(c,2*(i+j)+2,0,k)
        dh3=dJacobiP(c,2*(i+j)+2,0,k)

        pr(n)=dh1*h2*h3
        if(i.gt.1)pr(n)=pr(n)*((0.5d0*(1.0d0-b))**(i-1))
        if((i+j).gt.1)pr(n)=pr(n)*((0.5d0*(1.0d0-c))**(i+j-1))

        ps(n)=0.5d0*(1.0d0+a)*pr(n)
        tmp=dh2*((0.5d0*(1.0d0-b))**i)
        if(i.gt.0)tmp=tmp+(-0.5d0*i)*(h2*(0.5d0*(1.0d0-b))**(i-1))
        if((i+j).gt.1)tmp=tmp*((0.5d0*(1.0d0-c))**(i+j-1))
        tmp=tmp*h1*h3
        ps(n)=ps(n)+tmp

        pt(n)=0.5d0*(1.0d0+a)*pr(n)+0.5d0*(1.0d0+b)*tmp
        tmp=dh3*((0.5d0*(1.0d0-c))**(i+j))
        if((i+j).gt.0)tmp=tmp-0.5d0*(i+j)*(h3*((0.5d0*(1.0d0-c)) &
            **(i+j-1)))
        tmp=h1*h2*tmp*((0.5d0*(1.0d0-b))**i)
        pt(n)=pt(n)+tmp
    enddo

    tmp=2.0d0**(2*i+j+1.5d0)
    pr=pr*tmp
    ps=ps*tmp
    pt=pt*tmp

end subroutine GradBasis3D

subroutine GradVandermonde3D(V3Dr,V3Ds,V3Dt,N,r,s,t,xdim)
    integer      :: N,xdim,i,j,k,sk
    real(kind=rkind) :: V3Dr(xdim,(N+1)*(N+2)*(N+3)/6)
    real(kind=rkind) :: V3Ds(xdim,(N+1)*(N+2)*(N+3)/6)
    real(kind=rkind) :: V3Dt(xdim,(N+1)*(N+2)*(N+3)/6)
    real(kind=rkind) :: r(xdim),s(xdim),t(xdim)
    sk=1;
    do i=0,N,1
        do j=0,N-i,1
            do k=0,N-i-j,1
                call GradBasis3D(V3Dr(:,sk),V3Ds(:,sk),V3Dt(:,sk),&
                    r,s,t,xdim,i,j,k)
                sk=sk+1
            enddo
        enddo
    enddo
end subroutine GradVandermonde3D

subroutine Basis2D(p,r,s,xdim,i,j)
    integer i,j,xdim,n
    real(kind=rkind)::p(xdim),r(xdim),s(xdim),a,b,h1,h2

    do n=1,xdim,1
        if(abs(s(n)-1.0).le.TOL)then
            a=-1.0d0
        else
            a=2*(1.0d0+r(n))/(1.0d0-s(n))-1.0d0
        endif
        b=s(n)

        h1=JacobiP(a,0,0,i);h2=JacobiP(b,2*i+1,0,j)
        p(n)=sqrt(2.0d0)*h1*h2*(1.0d0-b)**i
    enddo
end subroutine Basis2D

subroutine Vandermonde2D(V2D,N,r,s,xdim)
    integer      :: N,xdim,i,j,sk
    real(kind=rkind) :: V2D(xdim,(N+1)*(N+2)/2)
    real(kind=rkind) :: r(xdim),s(xdim)
    sk=1
    do i=0,N,1
        do j=0,N-i,1
            call Basis2D(V2D(:,sk),r,s,xdim,i,j)
            sk=sk+1
        enddo
    enddo
end subroutine Vandermonde2D

subroutine GradBasis2D(pr,ps,r,s,xdim,i,j)
    integer i,j,xdim,n
    real(kind=rkind)::pr(xdim),ps(xdim),r(xdim),s(xdim)
    real(kind=rkind)::a,b,h1,h2,dh1,dh2,tmp

    do n=1,xdim
        if(abs(s(n)-1.0).le.TOL)then
                a=-1.0d0
        else
                a=2.0d0*(1.0d0+r(n))/(1.0d0-s(n))-1.0d0;
        endif
        b=s(n)

        h1=JacobiP(a,0,0,i);      dh1=dJacobiP(a,0,0,i)
        h2=JacobiP(b,2*i+1,0,j);  dh2=dJacobiP(b,2*i+1,0,j)

        pr(n)=dh1*h2
        if(i.gt.1)pr(n)=pr(n)*((0.5d0*(1.0d0-b))**(i-1))

        ps(n)=0.5d0*(1.0d0+a)*dh1*h2
        if(i.gt.1)ps(n)=ps(n)*((0.5d0*(1.0d0-b))**(i-1))

        tmp=dh2*((0.5d0*(1.0d0-b))**i)
        if(i.gt.0)tmp=tmp-0.5d0*i*h2*((0.5d0*(1.0d0-b))**(i-1))
        ps(n)=ps(n)+h1*tmp
    enddo

    tmp=2.0d0**(i+0.5d0)
    pr=pr*tmp
    ps=ps*tmp

end subroutine GradBasis2D

subroutine GradVandermonde2D(V2Dr,V2Ds,N,r,s,xdim)
    integer          :: N,xdim,i,j,sk
    real(kind=rkind) :: V2Dr(xdim,(N+1)*(N+2)/2)
    real(kind=rkind) :: V2Ds(xdim,(N+1)*(N+2)/2)
    real(kind=rkind) :: r(xdim),s(xdim)
    sk=1;
    do i=0,N,1
            do j=0,N-i,1
    call GradBasis2D(V2Dr(:,sk),V2Ds(:,sk),r,s,xdim,i,j)
    sk=sk+1
            enddo
    enddo
end subroutine GradVandermonde2D

subroutine Vandermonde1D(V1D,N,r,xdim)
    integer          :: N,xdim,i,j
    real(kind=rkind) :: V1D(xdim,N+1)
    real(kind=rkind) :: r(xdim)
    do i=0,N,1
        do j=1,xdim
            v1D(j,i+1)=JacobiP(r(j),0,0,i)
        enddo
    enddo
end subroutine Vandermonde1D

!-------------------------------------------------------------------

end module jacobi_mod

