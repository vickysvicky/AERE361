! Basis function subroutine
SUBROUTINE Basis(x,z,nbasis)
IMPLICIT NONE
REAL*8::x,z
INTEGER:: nbasis
DIMENSION z(nbasis)

!!!!!!!!!!!!!!!!!           change basis functions here           !!!!!!!!!!!!!!!!!
!! Test Case 1
! z(1)=1
! z(2)=x
! z(3)=x**2.0

!! Test Case 2
! z(1)=x
! z(2)=exp(x)
! z(3)=cos(x)

!! EXAM 1
z(1)=1
z(2)=1/x
z(3)=1/(1+sqrt(x))

END