PROGRAM bending
IMPLICIT NONE

! Declare variables
INTEGER::nload, nsupport, i, k, kc, j, NDIM, N
REAL*8:: coef(4,2,5), power, length, xend(2), COND, x, V, M, y, EI
REAL*8,ALLOCATABLE:: load(:,:,:), support(:,:,:), A(:,:), RHS(:), WORK(:)
INTEGER, ALLOCATABLE::IPVT(:)

! Hard coded Coefficient matrix ----------------------------------------------------------------------
coef(1,1,5)=(1.0/6.0)
coef(1,2,5)=3
coef(2,1,5)=(1.0/2.0)
coef(2,2,5)=2
coef(3,1,5)=1
coef(3,2,5)=1
coef(4,1,5)=1
coef(4,2,5)=0.0

! Ask user for load and support information in dat file ----------------------------------------------
! Takes in load data first, then support
! Data in position of magnitude, position, power
! Ask user for number of loads and support
WRITE(*,*) "Please save load data as beambend.dat"
WRITE(*,*) "Units of data have to be consistant (eg in, lb, m, N)"
WRITE(*,*) "How many loads are there?"
READ(*,*) nload
WRITE(*,*) "How many supports are there?"
READ(*,*) nsupport
ALLOCATE(load(nload,3,5), support(nsupport,3,5), A((nsupport+4),(nsupport+4)), RHS(nsupport+4))
NDIM=(nsupport+4)
N=NDIM
ALLOCATE(WORK(N),IPVT(N))

! Ask user for beam length
WRITE(*,*) "What is the length of beam?"
READ(*,*) length
length=length+1.0d-8
xend(1)=-1.0d-8
xend(2)=length

! Ask user for EI
WRITE(*,*) "What is EI of beam?"
READ(*,*) EI

! ! Ask user for point of interest
! WRITE(*,*) "Point of interest:"
! READ(*,*) x
! For this test, x=length/3
x=length/3

!!!!! Read in dat file -------------------------------------------------------------------------------
OPEN(unit=12, file="beambend.dat")
! Layer 1: load data
! load layer 1
DO i = 1,nload
	READ(12,*) load(i,1,1), load(i,2,1), load(i,3,1)
END DO
! support layer 1
DO i = 1,nsupport
	READ(12,*) support(i,1,1), support(i,2,1), support(i,3,1)
END DO

!!!!! Get layer 2-5 -----------------------------------------------------------------------------------
! Layer 2: shear data
! Layer 3: moment data
! Layer 4: dy/dx
! Layer 5: y
DO k=2,5
	! load layer
	DO i=1,nload
		load(i,3,k) = load(i,3,k-1) + 1
		power=load(i,3,k)
		IF (power/=0) THEN
			load(i,1,k) = load(i,1,k-1) / abs(power)
		ELSE
			load(i,1,k) = load(i,1,k-1)
		END IF
		load(i,2,k) = load(i,2,k-1)
	END DO
	! support layer
	DO i=1,nsupport
		support(i,3,k) = support(i,3,k-1) + 1
		power=support(i,3,k)
		IF (power/=0) THEN
			support(i,1,k) = support(i,1,k-1) / abs(power)
		ELSE
			support(i,1,k) = support(i,1,k-1)
		END IF
		support(i,2,k) = support(i,2,k-1)
	END DO
	! coefficient layer
	DO i=1,4
		kc=6-k
		power=coef(i,2,kc+1)
		coef(i,1,kc) = coef(i,1,kc+1) * abs(power)
		IF (power/=0) THEN
			coef(i,2,kc) = coef(i,2,kc+1) - 1
		ELSE
			coef(i,2,kc) = coef(i,2,kc+1)
		END IF
	END DO
END DO

!!!!! Show calculated layers --------------------------------------------------------------------------
DO k=1,5
	WRITE(*,*) "Load layer:",k
	DO i=1,nload
		WRITE(*,*) load(i,:,k)
	END DO
END DO
DO k=1,5
	WRITE(*,*) "Support Layer:",k
	DO i=1,nsupport
		WRITE(*,*) support(i,:,k)
	END DO
END DO
DO k=1,5
	WRITE(*,*) "Coefficient Layer:",k
	DO i=1,4
		WRITE(*,*) coef(i,:,k)
	END DO
END DO

! Show beam ends
WRITE(*,*) "Length"
WRITE(*,*) xend


!!!!! Fitting them all together ----------------------------------------------------------------------
! call function
CALL fitAll(support, coef, length, nsupport, xend, A)

! Show the final matrix
WRITE(*,*) "A matrix"
DO i=1,(nsupport+4)
	WRITE(*,*) A(i,:) 
END DO

!!!!! RHS matrix -------------------------------------------------------------------------------------
! RHS row 1 and 2
DO i=1,2
	RHS(i)=0.0d0
	DO j=1,nload
		IF ((load(j,3,2) >=0).AND.(xend(i)>=load(j,2,2))) THEN
			RHS(i) = RHS(i) - load(j,1,2)*(xend(i)-load(j,2,2))**load(j,3,2)
		ELSE 
			RHS(i) = RHS(i)-0
		END IF
	END DO
END DO

! RHS row 3 and 4
DO i=1,2
	RHS(i+2)=0.0d0
	DO j=1,nload
		IF ((load(j,3,3) >=0).AND.(xend(i)>=load(j,2,3))) THEN
			RHS(i+2) = RHS(i+2) - load(j,1,3)*(xend(i)-load(j,2,3))**load(j,3,3)
		ELSE
			RHS(i+2)= RHS(i+2)-0
		END IF
	END DO
END DO

! RHS row 5 onwards
DO i=1,nsupport
	RHS(i+4)=0.0d0
	DO j=1,nload
		IF ((load(j,3,5) >=0).AND.(support(i,2,5)>=load(j,2,5))) THEN
			RHS(i+4) = RHS(i+4) - load(j,1,5)*(support(i,2,5)-load(j,2,5))**load(j,3,5)
		ELSE
			RHS(i+4)= RHS(i+4)-0
		END IF
	END DO
END DO

! Show the final matrix
WRITE(*,*) "RHS"
DO i=1,(nsupport+4)
	WRITE(*,*) RHS(i) 
END DO

!!!!! Decomp and solve ----------------------------------------------------------------------------
CALL DECOMP (NDIM,N,COND,IPVT,WORK,A)
write(*,*) "condition number is :",COND
write(*,*) "ipvt:", IPVT
write(*,*) "Triangularized A:"
do i=1,n
	write(*,*) A(i,:)
end do

CALL SOLVE (NDIM,N,RHS,IPVT,A)
write(*,*) "Solution Vector:"
write(*,*) RHS

!!!! Find shear, moment, and deflection at given point--------------------------------------------
WRITE(*,*) "x:",x
V=0.0d0
M=0.0d0
y=0.0d0

! Calculate shear
DO j=1,nload
	IF ((load(j,3,2) >=0).AND.(x>=load(j,2,2))) THEN
		V = V + load(j,1,2)*(x-load(j,2,2))**load(j,3,2)
	ELSE
		V= V
	END IF
END DO
DO j=1,nsupport
	IF ((support(j,3,2)>=0).AND.(x>=support(j,2,2))) THEN
		V = V + RHS(j+4)*support(j,1,2) * (x-support(j,2,2))**support(j,3,2)
	ELSE
		V = V
	END IF
END DO
DO j=1,4
	IF (coef(j,2,2)>=0) THEN
		V = V + RHS(j)*coef(j,1,2) * x**coef(j,2,2)
	ELSE
		V = V
	END IF
END DO

write(*,*) "V:", V

! Calculate moment
DO j=1,nload
	IF ((load(j,3,3) >=0).AND.(x>=load(j,2,3))) THEN
		M = M + load(j,1,3)*(x-load(j,2,3))**load(j,3,3)
	ELSE
		M = M
	END IF
END DO
DO j=1,nsupport
	IF ((support(j,3,3)>=0).AND.(x>=support(j,2,3))) THEN
		M = M + RHS(j+4)*support(j,1,3) * (x-support(j,2,3))**support(j,3,3)
	ELSE
		M = M
	END IF
END DO
DO j=1,4
	IF (coef(j,2,3)>=0) THEN
		M = M + RHS(j)*coef(j,1,3) * x**coef(j,2,3)
	ELSE
		M = M
	END IF
END DO

write(*,*) "M:", M

! Calculate deflection
DO j=1,nload
	IF ((load(j,3,5) >=0).AND.(x>=load(j,2,5))) THEN
		y = y + load(j,1,5)*(x-load(j,2,5))**load(j,3,5)
	ELSE
		y = y
	END IF
END DO
DO j=1,nsupport
	IF ((support(j,3,5)>=0).AND.(x>=support(j,2,5))) THEN
		y = y + RHS(j+4)*support(j,1,5) * (x-support(j,2,5))**support(j,3,5)
	ELSE
		y = y
	END IF
END DO
DO j=1,4
	IF (coef(j,2,5)>=0) THEN
		y = y + RHS(j)*coef(j,1,5) * x**coef(j,2,5)
	ELSE
		y = y
	END IF
END DO
y=y/EI

write(*,*) "y:", y

END PROGRAM