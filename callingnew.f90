! linear solver with allocatable variables

PROGRAM calling
IMPLICIT NONE
REAL*8:: COND
REAL*8,ALLOCATABLE:: A(:,:), B(:), WORK(:)
INTEGER*8:: NDIM, n, i, j, k
INTEGER*8,ALLOCATABLE:: IPVT(:)
OPEN(unit=12, file="data1.dat")
READ(12,*)n
ALLOCATE(A(n,n),B(n),WORK(n),IPVT(n))
NDIM=n


	DO i = 1,n
		DO j = 1,n
			READ(12,*) A(i,j)
		END DO
	END DO
	DO k=1,n
		READ(12,*) B(k)
	END DO
	
CALL DECOMP(NDIM,N,COND,IPVT,WORK,A)
	CALL SOLVE (NDIM,N,B,IPVT,A)
	WRITE(*,*) "condition number = ",COND
	DO I=1,N
		WRITE(*,*) B(I)
	END DO


END PROGRAM