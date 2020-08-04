! could replace hard coded data in calling.f90
! so it could read data files and generalize data abd make data allocatable
! it's a part of a program so there is no program and end

OPEN(unit=12, file="data.dat")
READ(12,*)n
	DO i=1,n
		DO j=1,n
			READ(12,*) A(i,j)
		END DO
	END DO
	DO k=1,n
		READ(12,*) B(n)
	END DO