MODULE potential
  use global_variables
  implicit none
  DOUBLE PRECISION :: TOYX, TOYY, CM1

  ! Functions to override
  logical, parameter :: run_init = .true.
  logical, parameter :: perturb_coordinates_override = .true.

  ! Declare model globals here.
  ! This replaces the commons file.


CONTAINS
  !Compute the energy, energy derivatives and gradient (if called for)
  subroutine calc_energy_gradient()
    implicit none
    ! Wrapper to phase field model
    call IMPLEMENT_POTENTIAL(X, G, E, .true.)
  end subroutine

  subroutine init()
    ! Wrapped for code that needs to be called
    write(*,*) "InitWetting"
    call INIT_SYSTEM()
    write(*,*) "InitWetting end"
  end subroutine

  subroutine perturb_coordinates()
    ! Wrapper
    IMPLICIT NONE
  end subroutine perturb_coordinates

!Compute the energy, energy derivatives and gradient (if called for)
Subroutine IMPLEMENT_POTENTIAL(COORDS2,V,E,GTEST)
  IMPLICIT NONE
  LOGICAL GTEST
  DOUBLE PRECISION COORDS2(N), E
  DOUBLE PRECISION :: V(N)

  E = 0.0
  V(:) = 0.0

  TOYX=COORDS2(1)
  TOYY=COORDS2(2)

  CALL COMPUTE_ENERGY(E)

  IF (GTEST) THEN
   	CALL COMPUTE_GRAD(V)
  ENDIF

  V = V * CM1
  E = E * CM1

END SUBROUTINE IMPLEMENT_POTENTIAL


SUBROUTINE COMPUTE_ENERGY(E)
  IMPLICIT NONE
  DOUBLE PRECISION :: E
  E = 1.0D-7*(TOYX**4+TOYY**4) - 2.0D-4*(TOYX**3+TOYY**3) + 0.1442*(TOYX**2+TOYY**2) + 10924.6 &
     & -47.18*TOYX - 46.58*TOYY + 5.0D-3*TOYX*TOYY
END SUBROUTINE COMPUTE_ENERGY


SUBROUTINE COMPUTE_GRAD(V)
  IMPLICIT NONE
  DOUBLE PRECISION :: V(N)
  V(1) = TOYY/200.0 - 29.0*TOYX/2500.0 + (TOYX-500.0)**3/(2.5D6) + 141.0/50.0
  V(2) = TOYX/200.0 - 29.0*TOYY/2500.0 + (TOYY-500.0)**3/(2.5D6) + 171.0/50.0
END SUBROUTINE COMPUTE_GRAD


SUBROUTINE INIT_SYSTEM()
  IMPLICIT NONE
  INTEGER J1, J2, J3, Cur, S1
  DOUBLE PRECISION PI
  character(len=10)       :: datechar,timechar,zonechar
  integer                 :: values(8),itime1
  CALL DATE_AND_TIME(datechar,timechar,zonechar,values)
  itime1= values(7)*39 + values(8)
  CALL SDPRND(itime1)
END SUBROUTINE INIT_SYSTEM


END MODULE potential
