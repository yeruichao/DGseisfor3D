module PREM_mod 

    implicit none

! Layers:
!  1. Inner core            
!  2. Outer core            
!  3. D'' layer             
!  4. Lower Mantle          
!  5. Inner transition zone1
!  6. Inner transition zone2
!  7. Outer transition zone 
!  8. Low velocity zone lid 
!  9. Inner crust           
! 10. Outer crust           
! 11. Ocean                 
        
    double precision, parameter :: &
      h_PREM(12)=(/&
        0.0000d0, 1.2215d6, 3.4800d6, 3.6300d6, 5.7010d6, 5.7710d6,&
        5.9710d6, 6.1510d6, 6.3466d6, 6.3560d6, 6.3680d6, 6.3710d6 &
      /),&
      a_PREM(11)=(/&
        -2.1773d-10, -2.4123d-10, 0.0000d000,-3.0922d-11,&
         0.0000d000,  0.0000d000, 0.0000d000, 0.0000d000,&
         0.0000d000,  0.0000d000, 0.0000d000 &
      /),&
      b_PREM(11)=(/&
         1.9110d-8,  1.3976d-4, -5.0007d-4, -2.4441d-4,&
        -2.3286d-4, -1.2603d-3, -5.9706d-4,  1.0869d-4,&
         0.0000d00,  0.0000d00,  0.0000d00 &
      /),&
      c_PREM(11)=(/&
        1.3088d4, 1.2346d4, 7.3067d3, 6.7823d3, 5.3197d3, 1.1249d4,&
        7.1083d3, 2.6910d3, 2.9000d3, 2.6000d3, 1.0200d3 &
      /),&
      Grav_const=6.6743d-11,pi=3.1415926536

contains

subroutine prem_gravity(Ndim,x,y,z,x0,y0,z0,gx,gy,gz)
    integer :: Ndim,i
    double precision :: x(Ndim),y(Ndim),z(Ndim)
    double precision :: gx(Ndim),gy(Ndim),gz(Ndim)
    double precision :: x0,y0,z0,r,g,xx,yy,zz
    double precision :: M(11)
    do i=1,11
        M(i) = 4d0/5d0*pi*a_PREM(i)*(h_PREM(i+1)**5-h_PREM(i)**5) &
             +         pi*b_PREM(i)*(h_PREM(i+1)**4-h_PREM(i)**4) &
             + 4d0/3d0*pi*c_PREM(i)*(h_PREM(i+1)**3-h_PREM(i)**3)
    enddo
    do i=1,Ndim
        xx=(x0-x(i))*1d3
        yy=(y0-y(i))*1d3
        zz=(z0-z(i))*1d3
        r=sqrt(xx**2+yy**2+zz**2)
        g=Grav_const*prem_grav_eval(M,r)/(r**2)*1d-3
        gx(i)=g*xx/r;gy(i)=g*yy/r;gz(i)=g*zz/r
    enddo
end subroutine prem_gravity

function prem_grav_eval(M,r)
    double precision :: prem_grav_eval,M(11),r
    integer :: i
    prem_grav_eval=0d0
    do i=1,11
        prem_grav_eval=prem_grav_eval+&
            M(i)*h_region1(r,h_PREM(i+1))
    enddo
    do i=1,10
        prem_grav_eval=prem_grav_eval+&
            h_region2(r,h_PREM(i),h_PREM(i+1))*&
            ( 4d0/5d0*pi*a_PREM(i)*( r**5 - h_PREM(i)**5 ) &
            +         pi*b_PREM(i)*( r**4 - h_PREM(i)**4 ) &
            + 4d0/3d0*pi*c_PREM(i)*( r**3 - h_PREM(i)**3 ))
    enddo
end function prem_grav_eval

function h_region1(r,h)
    double precision :: h_region1,r,h
    h_region1=(sign(1d0,sign(1d0,r-h)-0.5d0)+1d0)/2d0
end function h_region1

function h_region2(r,h1,h2)
    double precision :: h_region2,r,h1,h2
    h_region2=(sign(1d0,sign(1d0,r-h1)-0.5d0)+1d0)*&
        (-sign(1d0,sign(1d0,r-h2)-0.5d0)+1d0)/4d0
end function h_region2

end module PREM_mod 
