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
!--------------------------------------------------------------------
    contains
!--------------------------------------------------------------------

function RSfric_F(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_F,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0,g
    g=exp((f0+b*log(s0*psi/L))/a)/(2d0*s0)
    RSfric_F=a*sigma*asinh(s*g)
end function RSfric_F

function RSfric_dFdpsi(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dFdpsi,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0,g
    g=exp((f0+b*log(s0*psi/L))/a)/(2d0*s0)
    RSfric_dFdpsi=a*sigma/sqrt((s*g)**2+1d0)*s*g*b/(a*psi)
end function RSfric_dFdpsi

function RSfric_invF(tauf,sigma,s,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_invF,tauf,sigma,s
    real(kind=rkind) :: a,b,L,s0,f0,g
    g=a*log( 2*sinh( tauf/(a*sigma) ) ) - f0 - a*log( s/s0 )
    RSfric_invF=L/s0*exp(g/b)
end function RSfric_invF

function RSfric_invFs(tauf,sigma,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_invFs,tauf,sigma,psi
    real(kind=rkind) :: a,b,L,s0,f0,g
    g=exp((f0+b*log(s0*psi/L))/a)/(2d0*s0)
    RSfric_invFs=sinh(tauf/(sigma*a))/g
end function RSfric_invFs

function RSfric_dFdN(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dFdN,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0,g
    g=exp((f0+b*log(s0*psi/L))/a)/(2d0*s0)
    RSfric_dFdN=a*asinh(s*g)
end function RSfric_dFdN

function RSfric_dFds(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dFds,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0,g
    g=exp((f0+b*log(s0*psi/L))/a)/(2d0*s0)
    RSfric_dFds=a*sigma*g/sqrt((g*s)**2+1d0)
end function RSfric_dFds

function RSfric_dFdSN(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dFdSN,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0,g
    g=exp((f0+b*log(s0*psi/L))/a)/(2*s0)
    RSfric_dFdSN=a*g/sqrt((g*s)**2+1d0)
end function RSfric_dFdSN

function RSfric_dFdS2(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dFdS2,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0,g
    g=exp((f0+b*log(s0*psi/L))/a)/(2d0*s0)
    RSfric_dFdS2=-a*sigma*s*(g/sqrt((g*s)**2+1d0))**3
end function RSfric_dFdS2

function RSfric_G(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_G,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0
    RSfric_G=-1d0+s*psi/L
end function RSfric_G

function RSfric_dGds(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dGds,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0
    RSfric_dGds=psi/L
end function RSfric_dGds

function RSfric_dGdpsi(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dGdpsi,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0
    RSfric_dGdpsi=s/L
end function RSfric_dGdpsi

function RSfric_invG(sigma,s,psi,a,b,L,s0,f0,dt)
    real(kind=rkind) :: RSfric_invG,sigma,s,psi,dt
    real(kind=rkind) :: a,b,L,s0,f0
    RSfric_invG=(psi+dt)/(1d0+dt*(s/L))
end function RSfric_invG

end module friction_mod
