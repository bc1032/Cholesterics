MODULE potential
  use global_variables
  implicit none
  DOUBLE PRECISION :: TOYX, TOYY, CM1, gridsize

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
  INTEGER :: I, J, K, Cur, lx, ly, lz
  DOUBLE PRECISION COORDS2(N), E, BULK, SPLAY, SURFACE, q_0, WS, TWIST, A, B, C, ks, kt, gridsize, s
  DOUBLE PRECISION :: V(N)
  DOUBLE PRECISION, ALLOCATABLE :: Q1(:,:), Q2(:,:), Q3(:,:), weight(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Q4(:,:), Q5(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Qt1(:,:), Qt2(:,:), Qt3(:,:), Qt4(:,:), Qt5(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: GRADXQ1(:,:), GRADXQ2(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: GRADXQ3(:,:), GRADXQ4(:,:), GRADXQ5(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: GRADZQ1(:,:), GRADZQ2(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: GRADZQ3(:,:), GRADZQ4(:,:), GRADZQ5(:,:)

  !Allocate array memory for all dimensions; length of lx, lz
  ALLOCATE(weight(lz,lx))
  ALLOCATE(Q1(lz,lx))
  ALLOCATE(Q2(lz,lx))
  ALLOCATE(Q3(lz,lx))
  ALLOCATE(Q4(lz,lx))
  ALLOCATE(Q5(lz,lx))
  ALLOCATE(Qt1(lz+1,lx+1))
  ALLOCATE(Qt2(lz+1,lx+1))
  ALLOCATE(Qt3(lz+1,lx+1))
  ALLOCATE(Qt4(lz+1,lx+1))
  ALLOCATE(Qt5(lz+1,lx+1))
  ALLOCATE(GRADXQ1(lz,lx))
  ALLOCATE(GRADXQ2(lz,lx))
  ALLOCATE(GRADXQ3(lz,lx))
  ALLOCATE(GRADXQ4(lz,lx))
  ALLOCATE(GRADXQ5(lz,lx))
  ALLOCATE(GRADZQ1(lz,lx))
  ALLOCATE(GRADZQ2(lz,lx))
  ALLOCATE(GRADZQ3(lz,lx))
  ALLOCATE(GRADZQ4(lz,lx))
  ALLOCATE(GRADZQ5(lz,lx))
  !Initialise array values for completeness
  weight(:,:) = 0.0D0
  GRADZQ1(:,:) = 0.0D0
  GRADZQ2(:,:) = 0.0D0
  GRADZQ3(:,:) = 0.0D0
  GRADZQ4(:,:) = 0.0D0
  GRADZQ5(:,:) = 0.0D0
  GRADXQ1(:,:) = 0.0D0
  GRADXQ2(:,:) = 0.0D0
  GRADXQ3(:,:) = 0.0D0
  GRADXQ4(:,:) = 0.0D0
  GRADXQ5(:,:) = 0.0D0
  Q1(:,:) = 0.0D0
  Q2(:,:) = 0.0D0
  Q3(:,:) = 0.0D0
  Q4(:,:) = 0.0D0
  Q5(:,:) = 0.0D0
  Qt1(:,:) = 0.0D0
  Qt2(:,:) = 0.0D0
  Qt3(:,:) = 0.0D0
  Qt4(:,:) = 0.0D0
  Qt5(:,:) = 0.0D0

  ws = 0.01D0
  bulk = 0.0D0
  splay = 0.0D0
  twist = 0.0D0
  surface = 0.0D0
  E = 0.0
  V(:) = 0.0
  lx = 30
  ly = 1
  lz = 30
  gridsize = 1.0D0 / (lz)
  A = 0.3D0
  B = 0.2D0
  C = 1.0D0
  ks = 1.0D0
  kt = 1.0D0
  ! q_0 = (((WINDNUMin+WINDNUMfin)*PI))/(2.0D0*(lz))
  ! a_0 = (((WINDNUMin)*PI))/((lz))
  ! bulkan = -((b**4) + 36*a*c*(b**2) + 216*(a**2)*(c**2) + (b**3)*SQRT((b**2) + 24*a*c) + 24*a*b*c*SQRT((b**2) + 24*a*c))&
  !      / (864 * (c**3) )
  ! twistan = (ly)*(kt*s**2/3.0D0)*(3*a_0**2-6*a_0*q_0+4*q_0**2)
  ! surfacean = ((b + SQRT(b**2 + 24*a*c))**2 * (lx)*(ly)*(SIN(pi*WINDNUMin))**2 )/ (8*c**2)

  s = (B + DSQRT(B**2 + 24D0*A*C))/(4.0D0*C)

!In the following loop we use central difference formulae (explicit method) for finite differences.
  do k = 1,lz
    do i = 1,lx
      Q1(k,i)=COORDS2(1 + ((i-1)*5) + ((k-1)*lx*5))
      Q2(k,i)=COORDS2(2 + ((i-1)*5) + ((k-1)*lx*5))
      Q3(k,i) = COORDS2(3 + ((i-1)*5) + ((k-1)*lx*5))
      Q4(k,i) = COORDS2(4 + ((i-1)*5) + ((k-1)*lx*5))
      Q5(k,i) = COORDS2(5 + ((i-1)*5) + ((k-1)*lx*5))

      Qt1(k,i) = (2.0D0*s/3.0D0)
      Qt2(k,i) = 0.0D0
      Qt3(k,i) = 0.0D0
      Qt4(k,i) = -s/3.0D0
      Qt5(k,i) = 0.0D0

      IF ( (K .EQ. lz) .OR. (K .EQ. 1) ) THEN
        !WEIGHT(K,I) = 0.5 * PFGRIDSIZE**3
        WEIGHT(K,I) = 0.5D0*gridsize**2
      ELSE
        !WEIGHT(K,I) = 1.0 * PFGRIDSIZE**3
        WEIGHT(K,I) = 1.0D0*gridsize**2
      end if

      bulk = bulk + weight(k,i)*(-a*(Q1(k,i)*Q4(k,i) + Q1(k,i)**2 + Q2(k,i)**2 + Q3(k,i)**2 + Q4(k,i)**2 + Q5(k,i)**2) &
            + b*(Q1(k,i)*(-Q2(k,i)**2 + Q4(k,i)**2 + Q5(k,i)**2) + (Q1(k,i)**2)*Q4(k,i) - 2*Q2(k,i)*Q3(k,i)*Q5(k,i) &
            - (Q2(k,i)**2)*Q4(k,i) + (Q3(k,i)**2)*Q4(k,i)) &
            + c*(Q1(k,i)*Q4(k,i) + Q1(k,i)**2 + Q2(k,i)**2 + Q3(k,i)**2 + Q4(k,i)**2 + Q5(k,i)**2)**2 )

      IF (k .EQ. lz) THEN
        GRADZQ1(k,i)   =  (WS/gridsize) * (3*Q1(k,i)-4*Q1(k-1,i)+Q1(k-2,i)) / (2.0D0)
        GRADZQ2(k,i)   =  (WS/gridsize) * (3*Q2(k,i)-4*Q2(k-1,i)+Q2(k-2,i)) / (2.0D0)
        GRADZQ3(k,i)   =  (WS/gridsize) * (3*Q3(k,i)-4*Q3(k-1,i)+Q3(k-2,i)) / (2.0D0)
        GRADZQ4(k,i)   =  (WS/gridsize) * (3*Q4(k,i)-4*Q4(k-1,i)+Q4(k-2,i)) / (2.0D0)
        GRADZQ5(k,i)   =  (WS/gridsize) * (3*Q5(k,i)-4*Q5(k-1,i)+Q5(k-2,i)) / (2.0D0)

        surface = surface + (weight(k,i)/gridsize)*0.5D0*Ws*( (Q1(k,i)-Qt1(k,i))**2 + 2*(Q2(k,i)-Qt2(k,i))**2 &
                + 2*(Q3(k,i)-Qt3(k,i))**2+(Qt1(k,i) + Qt4(k,i) - Q1(k,i) - Q4(k,i))**2 &
                + (Q4(k,i) - Qt4(k,i))**2 +2*(Q5(k,i) - Qt5(k,i))**2 )
      ELSEIF (k .EQ. 1) THEN
        GRADZQ1(k,i)   =  (WS/gridsize) * (-3*Q1(k,i)+4*Q1(k+1,i)-Q1(k+2,i)) / (2.0D0)
        GRADZQ2(k,i)   =  (WS/gridsize) * (-3*Q2(k,i)+4*Q2(k+1,i)-Q2(k+2,i)) / (2.0D0)
        GRADZQ3(k,i)   =  (WS/gridsize) * (-3*Q3(k,i)+4*Q3(k+1,i)-Q3(k+2,i)) / (2.0D0)
        GRADZQ4(k,i)   =  (WS/gridsize) * (-3*Q4(k,i)+4*Q4(k+1,i)-Q4(k+2,i)) / (2.0D0)
        GRADZQ5(k,i)   =  (WS/gridsize) * (-3*Q5(k,i)+4*Q5(k+1,i)-Q5(k+2,i)) / (2.0D0)
        surface = surface + (weight(k,i)/gridsize)*0.5D0*Ws*( (Q1(k,i)-Qt1(k,i))**2 + 2*(Q2(k,i)-Qt2(k,i))**2 &
                + 2*(Q3(k,i)-Qt3(k,i))**2+(Qt1(k,i) + Qt4(k,i) - Q1(k,i) - Q4(k,i))**2 &
                + (Q4(k,i) - Qt4(k,i))**2 +2*(Q5(k,i) - Qt5(k,i))**2 )
      ELSE
        !Calculate dz components.
        GRADZQ1(k,i) = (Q1(k+1,i) - Q1(k-1,i)) / (2.0D0)
        GRADZQ2(k,i) = (Q2(k+1,i) - Q2(k-1,i)) / (2.0D0)
        GRADZQ3(k,i) = (Q3(k+1,i) - Q3(k-1,i)) / (2.0D0)
        GRADZQ4(k,i) = (Q4(k+1,i) - Q4(k-1,i)) / (2.0D0)
        GRADZQ5(k,i) = (Q5(k+1,i) - Q5(k-1,i)) / (2.0D0)
      ENDIF
      IF (i == lx) THEN
        !Here we enforce periodicity in the x-direction (vertical).
        GRADXQ1(k,i) = (Q1(k,1) - Q1(k,i-1)) / (2.0D0 )
        GRADXQ2(k,i) = (Q2(k,1) - Q2(k,i-1)) / (2.0D0 )
        GRADXQ3(k,i) = (Q3(k,1) - Q3(k,i-1)) / (2.0D0 )
        GRADXQ4(k,i) = (Q4(k,1) - Q4(k,i-1)) / (2.0D0 )
        GRADXQ5(k,i) = (Q5(k,1) - Q5(k,i-1)) / (2.0D0 )
      ELSEIF (i == 1) THEN
        GRADXQ1(k,i) = (Q1(k,i+1) - Q1(k,lx+1)) / (2.0D0 )
        GRADXQ2(k,i) = (Q2(k,i+1) - Q2(k,lx+1)) / (2.0D0 )
        GRADXQ3(k,i) = (Q3(k,i+1) - Q3(k,lx+1)) / (2.0D0 )
        GRADXQ4(k,i) = (Q4(k,i+1) - Q4(k,lx+1)) / (2.0D0 )
        GRADXQ5(k,i) = (Q5(k,i+1) - Q5(k,lx+1)) / (2.0D0 )
      ELSE
        !Calculate dx components for divergence squared (splay) term
        GRADXQ1(k,i) = GRADXQ1(k,i) + (Q1(k,i+1) - Q1(k,i-1)) / (2.0D0)
        GRADXQ2(k,i) = GRADXQ2(k,i) + (Q2(k,i+1) - Q2(k,i-1)) / (2.0D0)
        GRADXQ3(k,i) = GRADXQ3(k,i) + (Q3(k,i+1) - Q3(k,i-1)) / (2.0D0)
        GRADXQ4(k,i) = GRADXQ4(k,i) + (Q4(k,i+1) - Q4(k,i-1)) / (2.0D0)
        GRADXQ5(k,i) = GRADXQ5(k,i) + (Q5(k,i+1) - Q5(k,i-1)) / (2.0D0)

      ENDIF

      twist = twist + (weight(k,i)/gridsize)*0.5D0*kt*( (-2*q_0*Q1(k,i)+gradzQ2(k,i))**2 &
            + 2*(2*q_0*Q3(k,i) - gradzQ5(k,i))*(2*q_0*Q3(k,i)+gradxQ2(k,i)) &
            + 2*(2*q_0*Q2(k,i)-gradzQ4(k,i))*(2*q_0*Q2(k,i)+gradzQ1(k,i) - gradxQ3(k,i)) &
            + 2*(2*q_0*Q5(k,i)+gradxQ4(k,i))*(2*q_0*Q5(k,i)+gradzQ3(k,i)+gradxQ1(k,i)-gradxQ4(k,i)) &
            + (2*q_0*Q4(k,i)+gradzQ2(k,i)-GRADXQ5(k,i))**2 + (-2*q_0*(Q1(k,i)+Q4(k,i)) + GRADXQ5(k,i))**2)

      splay = splay + 0.5D0*(weight(k,i)/gridsize)*ks*((gradzQ3(k,i)+gradxQ1(k,i))**2 &
            + ( gradzQ5(k,i) + gradxQ2(k,i) )**2 + (gradzQ1(k,i) + gradzQ4(k,i) - gradxQ3(k,i))**2 )
    enddo
  enddo

  CALL COMPUTE_ENERGY(E, bulk, twist, splay, surface)

  IF (GTEST) THEN
   	CALL COMPUTE_GRAD(V)
  ENDIF

  V = V * CM1
  E = E * CM1

END SUBROUTINE IMPLEMENT_POTENTIAL


SUBROUTINE COMPUTE_ENERGY(E, bulk, twist, splay, surface)
  IMPLICIT NONE
  DOUBLE PRECISION :: E, bulk, twist, splay, surface
  E = BULK + TWIST + SPLAY + SURFACE
  !write(*,*) E
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
