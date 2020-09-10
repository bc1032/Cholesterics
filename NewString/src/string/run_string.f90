PROGRAM RUN_STRING

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!Control the string program
!------------------------------------------------------------------------
!------------------------------------------------------------------------

USE STRING_COMMONS
USE STRING_IO
USE STRING_ROUTINES

IMPLICIT NONE

!1) Read in the data
CALL STRING_DATA()

!2) Read in the end points
CALL INITIALISE_MINIMA()

!3) Initialise the l-bfgs algorithm and potential data
CALL INITIALISE_LBFGS()

!4) Initialise the string
CALL INITIALISE_STRING()

!5) Run the simplified string algorithm
CALL EVOLVE_STRING()

!6) Perform the climbing image method
!CALL CLIMBING_IMAGE()

!7) Find the negatve eigenvector and eigenvalue at the TS
!CALL FIND_EIGEN()

!8) Finally, compute the accurate MEP, setting off downhill from the saddle point
!CALL DOWNHILL_MEP()

END PROGRAM RUN_STRING
