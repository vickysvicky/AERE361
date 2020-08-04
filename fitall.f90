! fit support and coefficient matrices together

SUBROUTINE fitAll(support, coef, length, nsupport, xend, A)
IMPLICIT NONE
INTEGER::i, j, nsupport, supnum
REAL*8::support(nsupport,3,5), coef(4,2,5), xend(2), length, A(nsupport+4,nsupport+4), x2(nsupport)



!WRITE(*,*) "elements of A for each calculation"
! A matrix row 1 and 2
DO i=1,2
	DO j=1,4
		IF (coef(j,2,2)>=0) THEN
			A(i,j) = coef(j,1,2) * (xend(i))**coef(j,2,2)
		ELSE
			A(i,j) = 0
		END IF
	END DO
	DO j=5,(nsupport+4)
		supnum = j - 4
		IF ((support(supnum,3,2)>=0).AND.(xend(i)>=support(supnum,2,2))) THEN
			A(i,j) = support(supnum,1,2) * (xend(i)-support(supnum,2,2))**support(supnum,3,2)
		ELSE
			A(i,j) = 0
		END IF
	END DO
END DO

! A matrix row 3 and 4
DO i=1, 2
	DO j=1,4
		IF (coef(j,2,3)>=0) THEN
			A(i+2,j) = coef(j,1,3) * (xend(i))**coef(j,2,3)
		ELSE
			A(i+2,j) = 0
		END IF
	END DO
	DO j=5,(nsupport+4)
		supnum = j - 4
		IF ((support(supnum,3,3)>=0).AND.(xend(i)>=support(supnum,2,3))) THEN
			A(i+2,j) = support(supnum,1,3) * (xend(i)-support(supnum,2,3))**support(supnum,3,3)
		ELSE
			A(i+2,j) = 0
		END IF
	END DO
END DO

! Positions of support
DO j=1,nsupport
	x2(j) = support(j,2,5)
END DO


! A matrix row 5 onwards
DO i=1, nsupport
	DO j=1,4
		IF (coef(j,2,5)>=0) THEN
			A(i+4,j) = coef(j,1,5) * (x2(i))**coef(j,2,5)
		ELSE
			A(i+4,j) = 0
		END IF
	END DO
	DO j=5,(nsupport+4)
		supnum = j - 4
		IF ((support(supnum,3,5)>=0).AND.(x2(i)>=support(supnum,2,5))) THEN
			A(i+4,j) = support(supnum,1,5) * (x2(i)-support(supnum,2,5))**support(supnum,3,5)
		ELSE
			A(i+4,j) = 0
		END IF
	END DO
END DO

END