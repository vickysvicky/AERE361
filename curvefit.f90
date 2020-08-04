PROGRAM curvefit
IMPLICIT NONE

! Declare variables
INTEGER:: nbasis,ndata,irow,jcol,i,j,k,NDIM,N
REAL*8,ALLOCATABLE:: x(:), y(:), z(:), c(:), ZM(:,:), ZT(:,:), ZZT(:,:), Zy(:)
REAL*8:: COND, xinput
REAL*8,ALLOCATABLE:: A(:,:), B(:), WORK(:)
INTEGER*8,ALLOCATABLE:: IPVT(:)

! Get number of basis functions and data pairs
WRITE(*,*) "Enter number of basis function: "
READ(*,*) nbasis
WRITE(*,*) "Enter number of paired observations: "
READ(*,*) ndata
NDIM=nbasis
n=NDIM

! Allocate after getting ndata and nbasis
ALLOCATE(x(ndata), y(ndata), z(nbasis), c(ndata))
ALLOCATE(ZM(ndata, nbasis), ZT(nbasis, ndata), ZZT(nbasis, nbasis), Zy(nbasis))
ALLOCATE(A(n,n),B(n),WORK(n),IPVT(n))



! Open data file and read in x and y
!!!!!!!!!!!!!!!!!           change dat file name here              !!!!!!!!!!!!!!!!!
OPEN(unit=12, file="exam1.dat")


WRITE(*,*) "Reading data..."
DO irow = 1,ndata
		READ(12,*) x(irow),y(irow)
		WRITE(*,*) x(irow),y(irow)
END DO
WRITE(*,*) "...Done reading"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! Codes to create linear system of equations !!!!!!!!!!!!!!!!!!!!!!!
! Construct z matrix
write(*,*) "ZM:"
DO irow = 1,ndata
	xinput=x(irow)
     CALL Basis(xinput,z,nbasis)
     c ( irow ) = y ( irow )
     DO jcol = 1,nbasis
         ZM ( irow , jcol ) = z ( jcol )
     END DO
END DO
do i=1,ndata
	write(*,*) ZM(i,:)
end do

! Get z transpose
! write(*,*) "ZT:"
! DO irow = 1,nbasis
     ! DO jcol = 1,ndata
         ! ZT ( irow , jcol ) = ZM ( jcol , irow )
		 ! write(*,*) zt(irow,jcol)
     ! END DO
! END DO
ZT=transpose(ZM)
write(*,*) "ZT:"
do i=1,nbasis
	write(*,*) zt(i,:)
end do

! Solve for ZZT
write(*,*) "ZZT:"
DO i = 1 , nbasis
     DO j = 1 , nbasis
         ZZT ( i , j ) = 0.0
         DO k = 1 , ndata
             ZZT ( i , j ) = ZT ( i , k ) * ZM ( k , j ) + ZZT ( i , j )
         END DO
     END DO
END DO
do i=1,nbasis
	write(*,*) zzt(i,:)
end do

! Solve for Zy
write(*,*) "Zy:"
DO i = 1 , nbasis
     Zy ( i ) = 0.0
     DO k = 1 , ndata
         Zy ( i ) = ZT ( i , k ) * c ( k ) + Zy ( i )
     END DO
END DO
write(*,*) zy
!!!!!!!!!!!!!!!!!! Codes to create linear system of equations !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




! Call Decomp and Solve
CALL DECOMP (NDIM,N,COND,IPVT,WORK,ZZT)
write(*,*) "condition number is :",COND
write(*,*) "ipvt:"
write(*,*) IPVT
write(*,*) "Triangularized ZZT:"
do i=1,n
	write(*,*) zzt(i,:)
end do

CALL SOLVE (NDIM,N,Zy,IPVT,ZZT)
write(*,*) "Solution Vector:"
write(*,*) Zy

END PROGRAM