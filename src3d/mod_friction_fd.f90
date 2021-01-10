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
    real(kind=rkind) :: a,b,L,s0,f0
    RSfric_F=sigma*a*asinh(0.5d0*s/s0*exp(psi/a))
end function RSfric_F

function RSfric_invF(tauf,sigma,s,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_invF,tauf,sigma,s
    real(kind=rkind) :: a,b,L,s0,f0
    RSfric_invF=a*log( 2d0*s0/s*sinh(tauf/(sigma*a)) )
end function RSfric_invF

function RSfric_dFdN(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dFdN,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0
    RSfric_dFdN=a*asinh(0.5d0*s/s0*exp(psi/a))
end function RSfric_dFdN

function RSfric_dFds(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dFds,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0,g
    g=0.5d0/s0*exp(psi/a)
    RSfric_dFds=a*sigma*g/sqrt(1+(s*g)**2)
end function RSfric_dFds

function RSfric_dFdSN(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dFdSN,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0,g
    g=0.5d0/s0*exp(psi/a)
    RSfric_dFdSN=a*g/sqrt(1+(s*g)**2)
end function RSfric_dFdSN

function RSfric_dFdS2(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dFdS2,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0,g
    g=0.5d0/s0*exp(psi/a)
    RSfric_dFdS2=-a*sigma*s*((g/sqrt(1+(s*g)**2))**3)
end function RSfric_dFdS2

function RSfric_G(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_G,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0,f,fss
    f=a*asinh(0.5d0*s/s0*exp(psi/a))
    fss=f0-(b-a)*log(s/s0)
    RSfric_G=(s/L)*(f-fss)
end function RSfric_G

function RSfric_dGdS(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dGdS,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0,f,fss,df,dfss
    f=a*asinh(0.5d0*s/s0*exp(psi/a))
    fss=f0-(b-a)*log(s/s0)
    df=0.5d0*a/s0*exp(psi/a)/sqrt((0.5d0*s/s0*exp(psi/a))**2+1d0)
    dfss=-(b-a)/s
    RSfric_dGdS=(f-fss)/L+(s/L)*(df-dfss)
end function RSfric_dGdS

end module friction_mod
