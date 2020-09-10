MODULE STRING_ROUTINES

USE STRING_COMMONS
USE GLOBAL_VARIABLES
USE BFGS_LOGIC
USE STRING_IO
USE POTENTIAL

IMPLICIT NONE

INTEGER :: CLIM_NO, STEP
DOUBLE PRECISION, ALLOCATABLE :: TAU_ZERO(:), DOWNHILL_IMAGE(:), TS_P(:), TS_M(:)
DOUBLE PRECISION :: E_TSP, E_TSM, E_TS

CONTAINS


!---------------------------------------------------------------------------------------------------
!Intialise the l-bfgs algorithm and potential
!---------------------------------------------------------------------------------------------------
SUBROUTINE INITIALISE_LBFGS()
IMPLICIT NONE

!1) Read in data from the data.f90 file
call set_defaults()
!2) Allocate the lbfgs arrays
call allocate_arrays()
!3) Initialise the potential
call init()


PRINT *, 'INITIALISE_LBFGS> Initialisation of L-BFGS complete'

END SUBROUTINE INITIALISE_LBFGS
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
!Define the initial string interpolation
!---------------------------------------------------------------------------------------------------
SUBROUTINE INITIALISE_STRING()
IMPLICIT NONE
INTEGER :: J1, J2, K1, K2, K3, CUR

 CLIMBSTAT=.FALSE.
 EIGENSTAT=.FALSE.

!1) Initialise the array of images
ALLOCATE(IMAGES(N_IMAGES,N_DOF), E_IMAGES(N_IMAGES), POS_IMAGES(N_IMAGES), STRING_ANGLES(N_IMAGES))
ALLOCATE(TAU_ZERO(N_DOF), DOWNHILL_IMAGE(N_DOF), TS_P(N_DOF), TS_M(N_DOF))

IMAGES(1,:)=MIN_1
IMAGES(N_IMAGES,:)=MIN_2

!2) Form the initial string
! Simple linear interpolation
DO J1=2,N_IMAGES-1
	IMAGES(J1,:) = MIN_1 + (J1-1.0)/(N_IMAGES-1.0)*(MIN_2-MIN_1)
ENDDO

IMAGES(1,:)=MIN_1
IMAGES(N_IMAGES,:)=MIN_2
PRINT *, 'INITIALISE_STRING> Initialisation of', N_IMAGES, ' image string complete'

END SUBROUTINE INITIALISE_STRING
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Evolve the string using the L-BFGS algorithm
!---------------------------------------------------------------------------------------------------
SUBROUTINE EVOLVE_STRING()
IMPLICIT NONE
INTEGER :: K
DOUBLE PRECISION :: STR_DIST_CHANGE, OLD_IMAGES(N_IMAGES,N_DOF)

DO STEP=1,MAXSTEP
	!1) Compute the energy of the string
	CALL STRING_STEP(100, 1.0D10, 0.25D0,.FALSE.)
	print *, 'step number = ', STEP

	!2) If called for, output the string, energy, and position
	IF (MOD(STEP-1,OUTPUTSTEP)==0) THEN
		PRINT *, 'EVOLVE_STRING> output at step', STEP-1
		CALL OUTPUT_STRING(0)
	ENDIF

	OLD_IMAGES=IMAGES
	!3) Perform a partial relaxation
	IF (STEP < 20) THEN
		CALL STRING_STEP(100, 1.0D-6, 0.25D0,.FALSE.)
	ELSEIF (STEP < 40) THEN
		CALL STRING_STEP(100, 1.0D-6, 0.25D-1,.FALSE.)
	ELSEIF (STEP < 60) THEN
		CALL STRING_STEP(100, 1.0D-6, 0.25D-2,.FALSE.)
	ELSEIF (STEP < 200) THEN
		CALL STRING_STEP(100, 1.0D-6, 0.25D-3,.FALSE.)
	ENDIF

	!!3) If called for, output the string, energy, and position
	!IF (STEP==10) THEN
	!	CALL OUTPUT_STRING(0)
	!	STOP
	!ENDIF
	
	!4) Reparameterise the string
	CALL REPARAMETERISE(OLD_IMAGES,STR_DIST_CHANGE)
	CALL STRING_STEP(100, 1.0D10, 0.25D0,.FALSE.)
	!CALL OUTPUT_STRING(0)
	!PRINT *, STR_DIST_CHANGE

	!5) Compute the angle of the tangent relative to the gardient
	CALL STRING_STEP(100, 1.0D10, 0.25D0,.TRUE.)
	print *, 'EVOLVE_STRING> Step ', STEP, ' average angle = ', SUM(STRING_ANGLES(2:N_IMAGES-1))/(N_IMAGES-2)

ENDDO

print *, 'EVOLVE_STRING> Complete'


END SUBROUTINE EVOLVE_STRING
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Evolve each image downhill and return the energy and position of each image
!---------------------------------------------------------------------------------------------------
SUBROUTINE STRING_STEP(im_iter,im_conv,im_stepsize,ANGLE_STAT)
IMPLICIT NONE
INTEGER :: J1, im_iter
DOUBLE PRECISION :: im_conv, im_stepsize, LOC_TAN(N_DOF), PI, GVEC(N_IMAGES,N_DOF)
LOGICAL :: ANGLE_STAT
PI=4*ATAN(1.0)

G=0
GVEC=0.0
max_iterations  = im_iter
convergence_rms = im_conv
max_step_size   = im_stepsize

DO J1=1,N_IMAGES
	!1) Relax the images accoring to the convergence inputs
	X=IMAGES(J1,:)
	E=0 
	CALL QUENCH()
	E_IMAGES(J1)=E
	IMAGES(J1,:)=X
	GVEC(J1,:)=G
ENDDO
DO J1=1,N_IMAGES
	!2) If required, compute the angle between the local tangent and the gradient
	IF (ANGLE_STAT) THEN
		!1) Compute tangent vector based on the downhill direction. If at a local max or min, average tangent on both sides
		IF (J1==1) THEN
			LOC_TAN=IMAGES(2,:)-IMAGES(1,:)
		ELSEIF (J1==N_IMAGES) THEN
			LOC_TAN=IMAGES(N_IMAGES-1,:)-IMAGES(N_IMAGES,:)
		ELSEIF ( (E_IMAGES(J1+1) > E_IMAGES(J1)) .AND. (E_IMAGES(J1) > E_IMAGES(J1-1)) ) THEN
			LOC_TAN=IMAGES(J1,:)-IMAGES(J1-1,:)
		ELSEIF ( (E_IMAGES(J1+1) < E_IMAGES(J1)) .AND. (E_IMAGES(J1) < E_IMAGES(J1-1)) ) THEN
			LOC_TAN=IMAGES(J1,:)-IMAGES(J1+1,:)
		ELSEIF ( (E_IMAGES(J1+1) > E_IMAGES(J1)) .AND. (E_IMAGES(J1) < E_IMAGES(J1-1)) ) THEN
			LOC_TAN=0.5*(IMAGES(J1+1,:)+IMAGES(J1-1,:)-2*IMAGES(J1,:))
		ELSEIF ( (E_IMAGES(J1+1) < E_IMAGES(J1)) .AND. (E_IMAGES(J1) > E_IMAGES(J1-1)) ) THEN
			LOC_TAN=-0.5*(IMAGES(J1+1,:)+IMAGES(J1-1,:)-2*IMAGES(J1,:))
		ENDIF

		STRING_ANGLES(J1) = DOT_PRODUCT(LOC_TAN,-GVEC(J1,:))
		STRING_ANGLES(J1) = -STRING_ANGLES(J1)/(NORM2(LOC_TAN)*NORM2(GVEC(J1,:)))
		STRING_ANGLES(J1) = ACOS(STRING_ANGLES(J1))*180/PI
		!PRINT *, J1, STRING_ANGLES(J1)
	ENDIF
ENDDO

POS_IMAGES=0.0
DO J1=2,N_IMAGES
	POS_IMAGES(J1) = POS_IMAGES(J1-1) + NORM2( IMAGES(J1,:)-IMAGES(J1-1,:) )
ENDDO

END SUBROUTINE STRING_STEP
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Redistribute the images along the string
!---------------------------------------------------------------------------------------------------
SUBROUTINE REPARAMETERISE(OLD_IMAGES,STR_DIST_CHANGE)
IMPLICIT NONE
INTEGER :: J1, J2
DOUBLE PRECISION :: AV_LEN, POS_J(N_IMAGES), NEW_IMAGES(N_IMAGES,N_DOF), STR_DIST_CHANGE, OLD_IMAGES(N_IMAGES,N_DOF)
DOUBLE PRECISION :: MIN_ENERGY, MAX_ENERGY, ENERGY_LENGTH, W_EXP

!For simple linear interpolation scheme, set exponent WEIGHT_EXP=0


!1) Get the total path length and average separation
PATH_LENGTH=POS_IMAGES(N_IMAGES)
AV_LEN=PATH_LENGTH/(N_IMAGES-1)


!1.1) Get the total energy length
MIN_ENERGY=MINVAL(E_IMAGES)
MAX_ENERGY=MAXVAL(E_IMAGES)
ENERGY_LENGTH=0.0

!1.2) Ensure a bunched-up initialisation becomes spread out, before applying image density weighting
IF (STEP<2) THEN
 	W_EXP=0.0
ELSE
	W_EXP=WEIGHT_EXP
ENDIF

!1.3) Apply image density weighting
DO J1=2,N_IMAGES
	ENERGY_LENGTH=ENERGY_LENGTH + ( (E_IMAGES(J1)-0.99*MIN_ENERGY)/MAX_ENERGY )**W_EXP
ENDDO

POS_J(1)=0.0
DO J1=2,N_IMAGES
	POS_J(J1)=POS_J(J1-1) + PATH_LENGTH/ENERGY_LENGTH*( (E_IMAGES(J1)-0.99*MIN_ENERGY)/MAX_ENERGY )**W_EXP
ENDDO

!print *, E_IMAGES
!print *, POS_J


!2) Redistribute the images equally along the piecewise-linearly interpolated path
NEW_IMAGES=0.0
NEW_IMAGES(1,:)=IMAGES(1,:)
NEW_IMAGES(N_IMAGES,:)=IMAGES(N_IMAGES,:)
DO J1=2,N_IMAGES-1


	DO J2=2,N_IMAGES
		IF ( (POS_IMAGES(J2) >= POS_J(J1)) .AND. (POS_IMAGES(J2-1) <= POS_J(J1)) ) THEN
			NEW_IMAGES(J1,:)=IMAGES(J2-1,:) + (POS_J(J1) - POS_IMAGES(J2-1))/(POS_IMAGES(J2)-POS_IMAGES(J2-1))*(IMAGES(J2,:)-IMAGES(J2-1,:))
		ENDIF
	ENDDO  

ENDDO

!3) Compute the total distance the string has moved
STR_DIST_CHANGE=0.0
DO J1=1,N_IMAGES
	STR_DIST_CHANGE=STR_DIST_CHANGE + NORM2( NEW_IMAGES(J1,:)-OLD_IMAGES(J1,:) )
	!print *, J1, NORM2( NEW_IMAGES(J1,:)-OLD_IMAGES(J1,:) )
ENDDO

IMAGES=NEW_IMAGES

END SUBROUTINE REPARAMETERISE
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Climbing image method
!---------------------------------------------------------------------------------------------------
SUBROUTINE CLIMBING_IMAGE()
IMPLICIT NONE
INTEGER :: J1

 ALLOCATE(CLIM_IMAGE(N_DOF))
 CLIM_IMAGE=0.0 

 !1) Get the highest energy image, and local tangent vector
 CLIM_NO=MAXLOC(E_IMAGES,1)
 CLIM_IMAGE = IMAGES(CLIM_NO,:)
 !TAU_ZERO = 0.5* ( IMAGES(CLIM_NO+1,:) + IMAGES(CLIM_NO-1,:) -2*IMAGES(CLIM_NO,:))
 TAU_ZERO=IMAGES(CLIM_NO,:) - IMAGES(CLIM_NO-1,:)
 TAU_ZERO = TAU_ZERO/NORM2(TAU_ZERO)

 !2) Call the l-bfgs to perform the climbing image
 CLIMBSTAT=.TRUE.
 G=0
 max_iterations  = 1D5
 convergence_rms = 1D-6
 max_step_size   = 0.25D0

 !3) Relax the images accoring to the convergence inputs
 X=CLIM_IMAGE
 CALL QUENCH()
 CLIM_IMAGE=X
 E_TS=E


 OPEN(1,FILE='outputs/climbing_image.out')
 DO J1=1,N_DOF
 	WRITE(1,*) CLIM_IMAGE(J1)
 ENDDO
 CLOSE(1)

 OPEN(1,FILE='outputs/climb_energy.out')
 WRITE(1,*) E_TS, NORM2(G)
 CLOSE(1)

print *, 'CLIMBING_IMAGE> Complete'

 CLIMBSTAT=.FALSE.

END SUBROUTINE CLIMBING_IMAGE
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Gradient modification for climbing image
!---------------------------------------------------------------------------------------------------
SUBROUTINE CLIMBING_GRADIENT()
IMPLICIT NONE

 G = G - T_SPEED*DOT_PRODUCT(G,TAU_ZERO)*TAU_ZERO
 PRINT *, 'CLIMBING_GRADIENT>',  NORM2(G)


END SUBROUTINE CLIMBING_GRADIENT
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Find the eigenvector and negative eigenvalue at the saddle point
!---------------------------------------------------------------------------------------------------
SUBROUTINE FIND_EIGEN()
IMPLICIT NONE
DOUBLE PRECISION :: EIGEN_P(N_DOF), EIGEN_M(N_DOF), CONVP, CONVM

EIGENSTAT=.TRUE.
ALLOCATE(EIGENVEC(N_DOF))

!1) Positive direction
DOWNHILL_IMAGE=IMAGES(CLIM_NO+1,:)
X=DOWNHILL_IMAGE

 G=0
 max_iterations  = 1D4
 convergence_rms = 1D-6
 max_step_size   = 0.25D0
 CALL QUENCH()
 TS_P=X
 E_TSP=E
 CONVP=NORM2(G)

!2) Negative direction
DOWNHILL_IMAGE=2*CLIM_IMAGE-TS_P
X=DOWNHILL_IMAGE
 G=0
 max_iterations  = 1D4
 convergence_rms = 1D-6
 max_step_size   = 0.25D0
 CALL QUENCH()
 TS_M=X
 E_TSM=E
 CONVM=NORM2(G)

!3) Compute the eigenvector
EIGEN_P=(TS_P-CLIM_IMAGE)/NORM2(TS_P-CLIM_IMAGE)
EIGEN_M=(CLIM_IMAGE-TS_M)/NORM2(CLIM_IMAGE-TS_M)

EIGENVEC=0.5*(EIGEN_P + EIGEN_M)
EIGENVAL=(E_TSP+E_TSM-2*E_TS)/(EIGEN_LENGTH**2)

OPEN(1,FILE='outputs/eigenconv.out')
WRITE(1,*) 'PLUS_SIDE', E_TSP, CONVP
WRITE(1,*) 'MINUS_SIDE', E_TSM, CONVM
CLOSE(1)

OPEN(1,FILE='outputs/eigenvec.out')
WRITE(1,*) EIGENVEC
CLOSE(1)

OPEN(1,FILE='outputs/eigenval.out')
WRITE(1,*) EIGENVAL
CLOSE(1)

print *, 'FIND_EIGEN> Complete'

EIGENSTAT=.FALSE.

END SUBROUTINE FIND_EIGEN
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Gradient modification for eigen finding
!---------------------------------------------------------------------------------------------------
SUBROUTINE EIGEN_GRADIENT()
IMPLICIT NONE
DOUBLE PRECISION :: SEPARATION
 
 !Apply a soft constraint to evolve an image in the downhill direction (but not too far)
 SEPARATION = NORM2(X-CLIM_IMAGE)
 E = E + EIGEN_SPRING*(SEPARATION-EIGEN_LENGTH)**2
 G = G + 2.0*EIGEN_SPRING*(1.0-EIGEN_LENGTH/SEPARATION)*(X-CLIM_IMAGE)


END SUBROUTINE EIGEN_GRADIENT
!---------------------------------------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------
! Peform a downhill descent after TS convergence to get MEP
!---------------------------------------------------------------------------------------------------
SUBROUTINE DOWNHILL_MEP()
IMPLICIT NONE
INTEGER :: J1, J2
DOUBLE PRECISION :: EMIN_P, EMIN_M, DOWN_IMAGE(N_DOF)
DOUBLE PRECISION :: DOWN_KICK

ALLOCATE(MEP(2*DOWN_SAVES+1,N_DOF), MEP_DATA(2*DOWN_SAVES+1,2))

DOWN_KICK=1.0D-3

!1) Step off the TS in both directions to get the connected minima

!1.1) Do a downhill run, stepping off in the positive direction, to get the energy change
PRINT *, 'DOWNHILL_MEP>', 'Connecting minima, plus-side'
DOWN_IMAGE=CLIM_IMAGE+DOWN_KICK*EIGENVEC
CALL DOWNHILL_RUN(DOWN_IMAGE,EMIN_P,0,.FALSE.)

!1.2) Do a downhill run, stepping off in the negative direction, to get the energy change
PRINT *, 'DOWNHILL_MEP>', 'Connecting minima, minus-side' 
DOWN_IMAGE=CLIM_IMAGE-DOWN_KICK*EIGENVEC
CALL DOWNHILL_RUN(DOWN_IMAGE,EMIN_M,0,.FALSE.)


!2) Repeat the downhill descent, saving DOWNSAVES images separated approximately by the same energy

!2.1) Save the TS and TS data
MEP(DOWN_SAVES+1,:)=CLIM_IMAGE
MEP_DATA(DOWN_SAVES+1,1)=0.0
MEP_DATA(DOWN_SAVES+1,2)=E_TS

!2.2) Do a downhill run in the positive direction
PRINT *, 'DOWNHILL_MEP>', 'Forming plus-side MEP' 
DOWN_IMAGE=CLIM_IMAGE+DOWN_KICK*EIGENVEC
CALL DOWNHILL_RUN(DOWN_IMAGE,EMIN_P,1,.TRUE.)

!2.3) Do a downhill run in the negative direction
PRINT *, 'DOWNHILL_MEP>', 'Forming minus-side MEP' 
DOWN_IMAGE=CLIM_IMAGE-DOWN_KICK*EIGENVEC
CALL DOWNHILL_RUN(DOWN_IMAGE,EMIN_M,-1,.TRUE.)

OPEN(10,FILE='outputs/mep.data')
DO J1=1,2*DOWN_SAVES+1
	WRITE(10,'(2f20.10)') MEP_DATA(J1,:)
ENDDO
CLOSE(10)

OPEN(10,FILE='outputs/mep.out')
DO J1=1,2*DOWN_SAVES+1
	DO J2=1,N_DOF
		WRITE(10,*) MEP(J1,J2)
	ENDDO
ENDDO
CLOSE(10)

END SUBROUTINE DOWNHILL_MEP
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Gradient descent for MEP calculation
!---------------------------------------------------------------------------------------------------
SUBROUTINE DOWNHILL_RUN(DOWN_INIT_IMAGE,DOWN_MIN_E,SAVE_DIR,DOWN_SAVESTAT)
IMPLICIT NONE
INTEGER :: J1, J_SAVE, SAVE_DIR
DOUBLE PRECISION :: DOWN_INIT_IMAGE(N_DOF), DOWN_MIN_E, E_OLD, SAVE_E
LOGICAL :: DOWN_SAVESTAT

INTEGER :: DOWN_FREQ, DOWN_STEPS
DOUBLE PRECISION :: dt, DOWN_RMS

DOWN_FREQ=10
dt=1.0D-3
DOWN_STEPS=1.0D4
DOWN_RMS=1.0D-6


max_iterations  = 1
convergence_rms = 1.0D60
max_step_size   = 0.25D-1

X=DOWN_INIT_IMAGE
!OPEN(10,FILE='downhill.out')
J_SAVE=1
E_OLD=E_TS
DO J1=1,DOWN_STEPS

	E=0
	CALL QUENCH()


	IF (DOWN_SAVESTAT) THEN
		SAVE_E=E_TS-(E_TS-DOWN_MIN_E)/(DOWN_SAVES)*(J_SAVE)
		IF ( (E_OLD > SAVE_E) .AND. (E <= SAVE_E) ) THEN
			MEP_DATA(DOWN_SAVES+1 + J_SAVE*SAVE_DIR,1)=SAVE_DIR*NORM2(CLIM_IMAGE-X)
			MEP_DATA(DOWN_SAVES+1 + J_SAVE*SAVE_DIR,2)=E
			MEP(DOWN_SAVES+1 + J_SAVE*SAVE_DIR,:)=X
			J_SAVE=J_SAVE+1
		ENDIF
	ENDIF


	X=X-G*DT

	!WRITE(10,*) J1, E, NORM2(V_DOWN)

	IF (NORM2(G)<=DOWN_RMS) THEN
		DOWN_MIN_E=E
		PRINT *, 'DOWNHILL_MEP>', 'Converged to E = ', DOWN_MIN_E, 'in', J1, 'steps' 
		EXIT
	ENDIF
	
ENDDO

IF (J1==DOWN_STEPS) THEN
	PRINT *, 'DOWNHILL_MEP>', 'Did not converge in', DOWN_STEPS, 'steps'
ENDIF

!CLOSE(10)
END SUBROUTINE DOWNHILL_RUN




END MODULE STRING_ROUTINES
