MODULE utils
 USE nrtype
CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE calc_arg(f1,alpha1)
!------------------------------------------------------------------------------
 IMPLICIT NONE
 COMPLEX(DP), INTENT(IN) :: f1
 REAL(DP), INTENT(OUT) :: alpha1
 if(DIMAG(f1) .lt. 0.0d0)then
  alpha1 = 2.d0 * PI - dacos(REAL(f1,DP))
 else
  alpha1 = dacos(REAL(f1,DP))
 endif
END SUBROUTINE calc_arg
!------------------------------------------------------------------------------
SUBROUTINE bcor_to_vcor(theta_b,phi_b,theta_R, phi_R)
!------------------------------------------------------------------------------
! transforms the polar coordinates of the molecular axis from 
! the frame fixed to the projectile velocity vector (Oy axis) into the
! frame fixed on the impact parameter (Oz axis)   
!------------------------------------------------------------------------------
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: theta_b, phi_b
 REAL(DP), INTENT(OUT) :: theta_R, phi_R
 REAL(DP) :: costetar, sintetar, sinphir, cosphir
 costetar = dsin(theta_b) * dcos(phi_b)
 theta_r = dacos(costetar)
 sintetar = dsqrt(1 - costetar*costetar)
 if(sintetar .LT. 1.d-10)then
   phi_r = 0.0d0
 else
   sinphir = dcos(theta_b) / sintetar
   cosphir = dsin(theta_b) * dsin(phi_b) / sintetar
   if( sinphir .GT. 0.d0 )then
     phi_R = dacos(cosphir)
   else
     phi_R = 2.d0*PI - dacos(cosphir)
   endif
 endif
END SUBROUTINE
!------------------------------------------------------------------------------
SUBROUTINE factorial(n,fact)
 IMPLICIT NONE
 INTEGER(I4B), INTENT(IN) :: n
 REAL(DP), INTENT(OUT) :: fact
 INTEGER(I4B):: ii
 if(n .LT. 0)then
   write(6,*)"ERROR_ Factorial called with negative argument"
   fact = 0.0d0
 else
   fact = 1.0d0
   do ii = 1, n
      fact = fact * REAL(ii,DP)
   enddo
 endif
END SUBROUTINE factorial
!------------------------------------------------------------------------------
SUBROUTINE poshsym(x,n,posh)
 IMPLICIT NONE
 INTEGER(I4B), INTENT(IN) :: n
 REAL(DP), INTENT(IN) :: x
 REAL(DP), INTENT(OUT) :: posh
 INTEGER(I4B) :: ii
 posh = 1.d0
 do ii=1,n
   posh = posh * (x+ii-1)
 enddo
  
END SUBROUTINE poshsym
!------------------------------------------------------------------------------
!-----------calculates the polar coordinates from descart coordinates----------
!------------------------------------------------------------------------------
SUBROUTINE  descarttopolar(x0,y0,z0,theta0,phi0)
 USE nrtype
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: x0,y0,z0
 REAL(DP), INTENT(OUT) :: theta0, phi0
 REAL(DP) :: rr, sinphi, cosphi
 if(abs(x0) .LT. 1.d-5 .AND. abs(y0) .LT. 1.d-5 .AND. abs(z0) .LT. 1.d-5)then
   theta0 = 0.0d0
   phi0   = 0.0d0
 else
   rr = dsqrt(x0*x0 + y0*y0 + z0*z0)
   theta0 = acos(z0/rr)
   if(abs(x0) .LT. 1.d-5)then
     if(abs(y0) .LT. 1.d-5)then
       phi0 = 0.0d0
     else
       phi0 = SIGN(PI/2.d0,y0)
     endif     
   else
     if(abs(y0) .LT. 1.d-5)then
       phi0 = (-1.d0 + SIGN(1.d0,x0))*PI/2.d0
     else
       if(dsin(theta0) .LT. 1.d-5)then
         phi0 = 0.0d0
       else
         cosphi = x0 / rr / dsin(theta0)
         sinphi = y0 / rr / dsin(theta0)
         if(cosphi .GT. 0.d0)then
            phi0 = dasin(sinphi)
         else
            phi0 = PI - dasin(sinphi)
         endif
       endif
     endif
   endif
 endif
END SUBROUTINE descarttopolar
!------------------------------------------------------------------------------- 

END MODULE utils
