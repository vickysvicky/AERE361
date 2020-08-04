! *****************************************************
!
! TEST PROGRAM FOR DECOMP AND SOLVE ROUTINES
! CREATED 1/30/03
!
! *****************************************************

	PROGRAM LINEAR_SOLVER

	REAL(KIND=8) :: A(20,20),WORK(20),B(20),COND
	INTEGER :: NDIM,N,I,IPVT(20)
!
! DIMENSION OF A,WORK,B,IPVT ARE ALL 20
! HENCE NDIM = 20
!
	NDIM=20
!
! N IS THE NUMBER OF EQUATIONS
!
	N=3
!
! ENTERING ELEMENTS OF MATRRIX A
!
	A(1,1) = 1.0
	A(1,2) = 4.0
	A(1,3) = 9.0
	A(2,1) = 4.0
	A(2,2) = 9.0
	A(2,3) = 16.0
	A(3,1) = 9.0
	A(3,2) = 16.0
	A(3,3) = 25.0
! 
!ENTERING ELEMENTS OF VECTOR B
!

	B(1) = 14.0
	B(2) = 29.0
	B(3) = 50.0
	
	CALL DECOMP(NDIM,N,COND,IPVT,WORK,A)
	CALL SOLVE (NDIM,N,B,IPVT,A)
!
! WRITING CONDITION NUMBER
!

	WRITE(*,*) "condition number = ",COND
!
! WRITING SOLUTION STORED IN B
!	
	DO 5 I=1,N
		WRITE(*,*) B(I)
5	END DO


	END PROGRAM LINEAR_SOLVER