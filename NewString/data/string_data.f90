SUBROUTINE STRING_DATA
USE STRING_COMMONS

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!This file specifies the user input data for the string algorithm
!------------------------------------------------------------------------
!------------------------------------------------------------------------

!String method parameters
N_IMAGES=20 !The total number of images in the string (including the ends)
N_DOF=5*30*30!The number of elements in the image coordinate vectors
MAXSTEP=40 !The maximum number of string evolution steps to try
OUTPUTSTEP=1 !The step interval to outout string data during convergence
WEIGHT_EXP=-0.0 !Weighting term for image concentration near TS (should be = 0 for linear, or less than zero)

!Climbing image parameters
T_SPEED=2 !The climb speed parameter nu in 'Ren, Vanden-Eijnden, (2013), A climbing string method for saddle point search, J. Chem. Phys.'

!Negative eigenvector finding method
EIGEN_LENGTH=1 !The distance of the image away from the saddle point to estimate the EV
EIGEN_SPRING=1.0D4 !The spring constant of the soft length constraint

!MEP saving
DOWN_SAVES=100 !The number of images to save at equal-energy intervals each side of the TS (2*DOWNSAVES+1 images saved in total)

END SUBROUTINE STRING_DATA
