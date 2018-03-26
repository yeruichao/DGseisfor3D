!*******************************************************************!
!*  This module contains friction laws for rupture problem         *!
!*                                                                 *!
!*  Author: Ruichao Ye    Email: Ruichao.ye@gmail.com              *!
!*******************************************************************!

!--------------------------------------------------------------------
module friction_mod
!--------------------------------------------------------------------

    use datatype_mod, only : rkind,matrices,sources,rupture

    implicit none
    real(kind=rkind),parameter :: eps=1d-12
    real(kind=rkind),parameter :: sig0=0.05d0,sig1=0.20d0
!--------------------------------------------------------------------
    contains
!--------------------------------------------------------------------

function RSfric_F(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_F,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0,g,x,sig
    g=exp((f0+b*log(s0*psi/L))/a)/(2d0*s0)
    x=sqrt(max(sigma-sig0,0d0)**2+eps)+sig0
    sig=sig1-sqrt(max(sig1-x,0d0)**2+eps)
    RSfric_F=a*sig*asinh(s*g)
end function RSfric_F

function RSfric_invF(tauf,sigma,s,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_invF,tauf,sigma,s
    real(kind=rkind) :: a,b,L,s0,f0,g,x,sig
    x=sqrt(max(sigma-sig0,0d0)**2+eps)+sig0
    sig=sig1-sqrt(max(sig1-x,0d0)**2+eps)
    g=a*log( 2*sinh( tauf/(a*sig) ) ) - f0 - a*log( s/s0 )
    RSfric_invF=L/s0*exp(g/b)
end function RSfric_invF

function RSfric_dFdP(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dFdP,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0,g,x,sig
    x=sqrt(max(sigma-sig0,0d0)**2+eps)+sig0
    sig=sig1-sqrt(max(sig1-x,0d0)**2+eps)
    g=exp((f0+b*log(s0*psi/L))/a)/(2d0*s0)
    RSfric_dFdP=a*sig/sqrt((s*g)**2+1d0)*s*g*b/(a*psi)
end function RSfric_dFdP

function RSfric_dFdN(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dFdN,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0,g,x,y,x1,x2,y1,y2,sig
    x1=max(sigma-sig0,0d0);x2=sqrt(x1**2+eps)
    x=sig0+x2
    y1=max(sig1-x,0d0);y2=sqrt(y1**2+eps)
    y=(x1*y1)/(x2*y2)
    g=exp((f0+b*log(s0*psi/L))/a)/(2d0*s0)
    RSfric_dFdN=a*asinh(s*g)*y
end function RSfric_dFdN

function RSfric_dFdS(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dFdS,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0,g,x,sig
    x=sqrt(max(sigma-sig0,0d0)**2+eps)+sig0
    sig=sig1-sqrt(max(sig1-x,0d0)**2+eps)
    g=exp((f0+b*log(s0*psi/L))/a)/(2d0*s0)
    RSfric_dFdS=a*sig*g/sqrt((g*s)**2+1d0)
end function RSfric_dFdS

function RSfric_G(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_G,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0
    RSfric_G=-1d0+s*psi/L
end function RSfric_G

function RSfric_dGdS(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dGdS,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0
    RSfric_dGdS=psi/L
end function RSfric_dGdS

function RSfric_dGdP(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dGdP,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0
    RSfric_dGdP=s/L
end function RSfric_dGdP

function RSfric_invG(sigma,s,psi,a,b,L,s0,f0,dt)
    real(kind=rkind) :: RSfric_invG,sigma,s,psi,dt
    real(kind=rkind) :: a,b,L,s0,f0
    RSfric_invG=(psi+dt)/(1d0+dt*(s/L))
end function RSfric_invG

end module friction_mod
