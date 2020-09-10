MODULE STRING_IO

USE STRING_COMMONS

CONTAINS

!---------------------------------------------------------------------------------------------------
!This routine reads in the input minima
!---------------------------------------------------------------------------------------------------
SUBROUTINE INITIALISE_MINIMA()
IMPLICIT NONE

!1)Initialise the working arrays
ALLOCATE(MIN_1(N_DOF), MIN_2(N_DOF))

!2)Read in the minima
OPEN(1,FILE='data/minimum_1')
READ(1,'(F20.10)') MIN_1
CLOSE(1)
OPEN(1,FILE='data/minimum_2')
READ(1,'(F20.10)') MIN_2
CLOSE(1)

PRINT *, 'INITIALISE_MINIMA> Initialisation complete'

END SUBROUTINE INITIALISE_MINIMA
!---------------------------------------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------
!This routine outputs the string data
!---------------------------------------------------------------------------------------------------
SUBROUTINE OUTPUT_STRING(STEP_NUMBER)
IMPLICIT NONE
INTEGER :: STEP_NUMBER, J1, J2
CHARACTER(len=30) :: FILENAME
CHARACTER(len=80) :: STRNAME, ENNAME, POSNAME

write(FILENAME,'(I20)') STEP_NUMBER

WRITE(STRNAME,'(A,A,A)') 'outputs/string', trim(adjustl(FILENAME)), '.out'
WRITE(ENNAME,'(A,A,A)') 'outputs/energy', trim(adjustl(FILENAME)), '.out'
WRITE(POSNAME,'(A,A,A)') 'outputs/position', trim(adjustl(FILENAME)), '.out'



OPEN(1,FILE=STRNAME)
DO J1=1,N_IMAGES
	DO J2=1,N_DOF
		WRITE(1,*) IMAGES(J1,J2)
	ENDDO
ENDDO
CLOSE(1)

OPEN(1,FILE=ENNAME)
DO J1=1,N_IMAGES
	WRITE(1,*) E_IMAGES(J1)
ENDDO
CLOSE(1)

OPEN(1,FILE=POSNAME)
DO J1=1,N_IMAGES
	WRITE(1,*) POS_IMAGES(J1)
ENDDO
CLOSE(1)

OPEN(1,FILE='outputs/maxe.out', ACCESS='append')
	WRITE(1,*) STEP_NUMBER, MAXVAL(E_IMAGES)
CLOSE(1)

OPEN(1,FILE='outputs/angles.out')
	WRITE(1,*) STRING_ANGLES
CLOSE(1)

END SUBROUTINE OUTPUT_STRING
!---------------------------------------------------------------------------------------------------






END MODULE STRING_IO
