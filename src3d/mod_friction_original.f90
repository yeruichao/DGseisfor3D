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
    RSfric_F=sigma*(a*log(s/s0)+psi)
end function RSfric_F

function RSfric_invF(tauf,sigma,s,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_invF,tauf,sigma,s
    real(kind=rkind) :: a,b,L,s0,f0
    RSfric_invF=tauf/sigma-a*log(s/s0)
end function RSfric_invF

function RSfric_invFs(tauf,sigma,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_invFs,tauf,sigma,psi
    real(kind=rkind) :: a,b,L,s0,f0
    RSfric_invFs=exp( (tauf/sigma-psi)/a )* s0
end function RSfric_invFs

function RSfric_dFdN(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dFdN,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0
    RSfric_dFdN=a*log(s/s0)+psi
end function RSfric_dFdN

function RSfric_dFds(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dFds,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0
    RSfric_dFds=a*sigma/s
end function RSfric_dFds

function RSfric_dFdSN(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dFdSN,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0
    RSfric_dFdSN=a/s
end function RSfric_dFdSN

function RSfric_dFdS2(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dFdS2,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0
    RSfric_dFdS2=-a*sigma/(s**2)
end function RSfric_dFdS2

function RSfric_G(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_G,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0
    RSfric_G=(s/L)*(psi-f0+b*log(s/s0))
end function RSfric_G

function RSfric_dGdS(sigma,s,psi,a,b,L,s0,f0)
    real(kind=rkind) :: RSfric_dGdS,sigma,s,psi
    real(kind=rkind) :: a,b,L,s0,f0
    RSfric_dGdS=(psi-f0+b*log(s/s0)+b)/L
end function RSfric_dGdS

end module friction_mod
