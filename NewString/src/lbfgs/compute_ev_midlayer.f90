SUBROUTINE COMPUTE_EV_MIDLAYER()

USE STRING_COMMONS
USE GLOBAL_VARIABLES
USE POTENTIAL
USE STRING_ROUTINES
IMPLICIT NONE

!---------------------------------------------------------------------------------------------------
!This subroutine is sandwiched between bfgs_logic and potential. It's purpose is to enable 
!modifications to the gradient and outputs of potential, without requiring the user to make
!any changes to the potential.f90 file from its GMIN version.
!---------------------------------------------------------------------------------------------------

call calc_energy_gradient()
!PRINT *, E, NORM2(G)
IF (CLIMBSTAT) THEN
	CALL CLIMBING_GRADIENT
ENDIF

IF (EIGENSTAT) THEN
	CALL EIGEN_GRADIENT
ENDIF

END SUBROUTINE COMPUTE_EV_MIDLAYER
