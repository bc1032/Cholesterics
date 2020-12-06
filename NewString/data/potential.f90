MODULE potential
  use global_variables
  implicit none
  DOUBLE PRECISION :: TOYX, TOYY, CM1, gridsize
  INTEGER :: lx, ly, lz

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
  INTEGER :: I, J, K, Cur, lx, ly, lz, w_0,w_1, r1,l1
  DOUBLE PRECISION COORDS2(N), E, BULK, SPLAY, SURFACE, q_0, WS, TWIST, a, b, c, ks, kt, gridsize, s
  DOUBLE PRECISION :: V(N), pi
  DOUBLE PRECISION, ALLOCATABLE :: Q1(:,:), Q2(:,:), Q3(:,:), weight(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Q4(:,:), Q5(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Qt1(:,:), Qt2(:,:), Qt3(:,:), Qt4(:,:), Qt5(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: GRADXQ1(:,:), GRADXQ2(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: GRADXQ3(:,:), GRADXQ4(:,:), GRADXQ5(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: GRADZQ1(:,:), GRADZQ2(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: GRADZQ3(:,:), GRADZQ4(:,:), GRADZQ5(:,:)

  !Allocate array memory for all dimensions; length of lx, lz
  lx = 100
  ly = 1
  lz = 100
  ALLOCATE(weight(lz,lx))
  ALLOCATE(Q1(lz,lx))
  ALLOCATE(Q2(lz,lx))
  ALLOCATE(Q3(lz,lx))
  ALLOCATE(Q4(lz,lx))
  ALLOCATE(Q5(lz,lx))
  ALLOCATE(Qt1(lz,lx))
  ALLOCATE(Qt2(lz,lx))
  ALLOCATE(Qt3(lz,lx))
  ALLOCATE(Qt4(lz,lx))
  ALLOCATE(Qt5(lz,lx))
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
  PI = 4.0D0*ATAN(1.0D0)
  w_0 = 5
  w_1 = 2
  gridsize = 1.0D0 / (lz)
  a = 0.3D0
  b = 0.2D0
  c = 1.0D0
  ks = 1.0D0
  kt = 1.0D0
  q_0 = (((w_0+w_1)*PI))/(2.0D0*(lz))
  ! a_0 = (((WINDNUMin)*PI))/((lz))
  ! bulkan = -((b**4) + 36*a*c*(b**2) + 216*(a**2)*(c**2) + (b**3)*SQRT((b**2) + 24*a*c) + 24*a*b*c*SQRT((b**2) + 24*a*c))&
  !      / (864 * (c**3) )
  ! twistan = (ly)*(kt*s**2/3.0D0)*(3*a_0**2-6*a_0*q_0+4*q_0**2)
  ! surfacean = ((b + SQRT(b**2 + 24*a*c))**2 * (lx)*(ly)*(SIN(pi*WINDNUMin))**2 )/ (8*c**2)

  s = (B + DSQRT(B**2 + 24D0*A*C))/(4.0D0*C)

!In the following loop we use central difference formulae (explicit method) for finite differences.
  l1 = 0
  r1 = 0

  do k = 1,lz
    do i = 1,lx
      if ( i == 1 ) THEN
        l1 = lx!i-1
        r1 = i+1
      else if ( i == 2 ) THEN
        l1 = 1
        r1 = i+1
                ! write(*,*) "i=2"
      else if ( i == lx-1 ) THEN
        l1 = i-1
        r1 = lx!i+1
        ! write(*,*) "i=lx-1"
      else if ( i == lx ) THEN
        l1 = i-1
        r1 = 1!i+1
        ! write(*,*) "i=lx"
      else
        l1 = i-1
        r1 = i+1
      end if
      Q1(k,i) = COORDS2(1 + ((i-1)*5) + ((k-1)*lx*5))
      Q2(k,i) = COORDS2(2 + ((i-1)*5) + ((k-1)*lx*5))
      Q3(k,i) = COORDS2(3 + ((i-1)*5) + ((k-1)*lx*5))
      Q4(k,i) = COORDS2(4 + ((i-1)*5) + ((k-1)*lx*5))
      Q5(k,i) = COORDS2(5 + ((i-1)*5) + ((k-1)*lx*5))

      Qt1(k,i) = (2.0D0*s/3.0D0)
      Qt2(k,i) = 0.0D0
      Qt3(k,i) = 0.0D0
      Qt4(k,i) = -s/3.0D0
      Qt5(k,i) = 0.0D0

      IF ( (k .EQ. lz) .OR. (k .EQ. 1) ) THEN
        !WEIGHT(K,I) = 0.5 * PFGRIDSIZE**3
        weight(k,i) = 1.0D0*gridsize**2!change from 0.5D0 to 1.0D0
      ELSE
        !WEIGHT(K,I) = 1.0 * PFGRIDSIZE**3
        weight(k,i) = 1.0D0*gridsize**2
      end if
      !write(*,*) Q1(k,i)
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
      ! IF (i == lx) THEN
      !   !Here we enforce periodicity in the x-direction (vertical).
      !   GRADXQ1(k,i) = (Q1(k,1) - Q1(k,l1)) / (2.0D0 )
      !   GRADXQ2(k,i) = (Q2(k,1) - Q2(k,l1)) / (2.0D0 )
      !   GRADXQ3(k,i) = (Q3(k,1) - Q3(k,l1)) / (2.0D0 )
      !   GRADXQ4(k,i) = (Q4(k,1) - Q4(k,l1)) / (2.0D0 )
      !   GRADXQ5(k,i) = (Q5(k,1) - Q5(k,l1)) / (2.0D0 )
      ! ELSEIF (i == 1) THEN
      !   GRADXQ1(k,i) = (Q1(k,r1) - Q1(k,lx)) / (2.0D0 )
      !   GRADXQ2(k,i) = (Q2(k,r1) - Q2(k,lx)) / (2.0D0 )
      !   GRADXQ3(k,i) = (Q3(k,r1) - Q3(k,lx)) / (2.0D0 )
      !   GRADXQ4(k,i) = (Q4(k,r1) - Q4(k,lx)) / (2.0D0 )
      !   GRADXQ5(k,i) = (Q5(k,r1) - Q5(k,lx)) / (2.0D0 )
      ! ELSE
        !Calculate dx components for divergence squared (splay) term
        GRADXQ1(k,i) = GRADXQ1(k,i) + (Q1(k,r1) - Q1(k,l1)) / (2.0D0)
        GRADXQ2(k,i) = GRADXQ2(k,i) + (Q2(k,r1) - Q2(k,l1)) / (2.0D0)
        GRADXQ3(k,i) = GRADXQ3(k,i) + (Q3(k,r1) - Q3(k,l1)) / (2.0D0)
        GRADXQ4(k,i) = GRADXQ4(k,i) + (Q4(k,r1) - Q4(k,l1)) / (2.0D0)
        GRADXQ5(k,i) = GRADXQ5(k,i) + (Q5(k,r1) - Q5(k,l1)) / (2.0D0)

      ! ENDIF

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
   	CALL COMPUTE_GRAD(V,a,b,c,ws,Q1,Q2,Q3,Q4,Q5,Qt1,Qt2,Qt3,Qt4,Qt5,weight,lx,lz,ks,kt,q_0)
  ENDIF

  V = V * CM1
  E = E * CM1

END SUBROUTINE IMPLEMENT_POTENTIAL


SUBROUTINE COMPUTE_ENERGY(E, bulk, twist, splay, surface)
  IMPLICIT NONE
  DOUBLE PRECISION :: E, bulk, twist, splay, surface
  E = bulk + twist + splay + surface
END SUBROUTINE COMPUTE_ENERGY


SUBROUTINE COMPUTE_GRAD(V,a,b,c,ws,Q1,Q2,Q3,Q4,Q5,Qt1,Qt2,Qt3,Qt4,Qt5,weight,lx,lz,ks,kt,q_0)
  IMPLICIT NONE
  DOUBLE PRECISION :: V(N), a, b, c, ws, gridsize,ks,kt,q_0, dz, dx
  INTEGER :: i, k, lx, lz, coxl2, coxl1, coxr1, coxr2
  DOUBLE PRECISION, ALLOCATABLE :: Q1(:,:), Q2(:,:), Q3(:,:), weight(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Q4(:,:), Q5(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Qt1(:,:), Qt2(:,:), Qt3(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Qt4(:,:), Qt5(:,:)
  gridsize = 1.0D0/lz
  dx = 1.0D0/lx
  dz = 1.0D0/lz
  V(:) = 0.0
  coxl2 = 0
  coxl1 = 0
  coxr1 = 0
  coxr2 = 0

  do k = 1,lz
    do i = 1,lx
      if ( i == 1 ) THEN
        coxl2 = lx-1!i-2
        coxl1 = lx!i-1
        coxr1 = i+1
        coxr2 = i+2
        ! write(*,*) "i=1"
      else if ( i == 2 ) THEN
        coxl2 = lx!i - 2
        coxl1 = 1
        coxr1 = i+1
        coxr2 = i+2
                ! write(*,*) "i=2"
      else if ( i == lx-1 ) THEN
        coxl1 = i-1
        coxl2 = i-2
        coxr1 = lx!i+1
        coxr2 = 1!i+2
        ! write(*,*) "i=lx-1"
      else if ( i == lx ) THEN
        coxl1 = i-1
        coxl2 = i-2
        coxr1 = 1!i+1
        coxr2 = 2!i+2
        ! write(*,*) "i=lx"
      else
        coxl2 = i-2
        coxl1 = i-1
        coxr1 = i+1
        coxr2 = i+2
        ! write(*,*) "bulk"
      endif

      !Surface energy
      if ( k == 1 .or. k == lz ) THEN
        V(1 + ((i-1)*5) + ((k-1)*lx*5)) = V(1 + ((i-1)*5) + ((k-1)*lx*5)) &
        + (weight(k,i)/gridsize)*(ws/2.0D0)*(2*(Q1(k,i)-Qt1(k,i))-2*(Qt1(k,i)+Qt4(k,i)-Q1(k,i)-Q4(k,i)))

        V(2 + ((i-1)*5) + ((k-1)*lx*5)) = V(2 + ((i-1)*5) + ((k-1)*lx*5)) &
        + (weight(k,i)/gridsize)*2.0D0*ws*(Q2(k,i)-Qt2(k,i))

        V(3 + ((i-1)*5) + ((k-1)*lx*5)) = V(3 + ((i-1)*5) + ((k-1)*lx*5)) &
        + (weight(k,i)/gridsize)*2.0D0*ws*(Q3(k,i)-Qt3(k,i))

        V(4 + ((i-1)*5) + ((k-1)*lx*5)) = V(4 + ((i-1)*5) + ((k-1)*lx*5))&
        + (weight(k,i)/gridsize)*(ws/2.0D0)*(2*(Q4(k,i)-Qt4(k,i))-2*(Qt1(k,i)+Qt4(k,i)-Q1(k,i)-Q4(k,i)))

        V(5 + ((i-1)*5) + ((k-1)*lx*5)) = V(5 + ((i-1)*5) + ((k-1)*lx*5)) &
        + (weight(k,i)/gridsize)*2.0D0*ws*(Q5(k,i)-Qt5(k,i))

      endif

      if ( (k == 1) ) THEN
        !Elastic ENERGY
        V(1 + ((i-1)*5) + ((k-1)*lx*5)) = V(1 + ((i-1)*5) + ((k-1)*lx*5)) &
        - (1.0D0/(4.0D0*dx*dz))*( (ks + kt)*(dz**2)*Q1(k,coxl2) &
        - 2.0D0*Q1(k,i)*(5.0D0*(ks+kt)*(dx**2) + (dz**2)*(ks+kt+16.0D0*kt*(q_0**2)*(dx**2)*dz**2)) &
        + 12.0D0*(ks+kt)*(dx**2)*Q1(k+1,i) - 2.0D0*(ks+kt)*(dx**2)*Q1(k+2,i) &
        + (ks+kt)*(dz**2)*Q1(k,coxr2) + 20.0D0*kt*q_0*(dx**2)*dz*Q2(k+1,i) &
        - 4.0D0*kt*q_0*(dx**2)*dz*Q2(k+2,i) + 6.0D0*dx*dz*Q3(k,coxl1)*(ks+kt) - 3.0D0*dx*dz*Q3(k+1,coxl1)*(ks+kt) &
        + dx*dz*Q3(k+2,coxl1)*(ks+kt) - 6.0D0*dx*dz*Q3(k,coxr1)*(ks+kt) + 3.0D0*dx*dz*Q3(k+1,coxr1)*(ks+kt) &
        - dx*dz*Q3(k+2,coxr1)*(ks+kt) + kt*(dz**2)*Q4(k,coxl2) &
        - 10.0D0*ks*(dx**2)*Q4(k,i) - 2.0D0*kt*(dz**2)*Q4(k,i) - 16.0D0*kt*(q_0**2)*(dx**2)*(dz**2)*Q4(k,i) &
        + 12.0D0*ks*(dx**2)*Q4(k+1,i)- 2.0*ks*(dx**2)*Q4(k+2,i) + kt*(dz**2)*Q4(k,coxr2)&
        + 8.0D0*kt*q_0*dx*(dz**2)*(Q5(k,coxr1)-Q5(k,coxl1)) )

        V(2 + ((i-1)*5) + ((k-1)*lx*5)) = V(2 + ((i-1)*5) + ((k-1)*lx*5)) &
        + (1.0D0/4.0D0)*(20.0D0*kt*q_0*dx*Q1(2,i) - 4.0D0*kt*q_0*dx*Q1(3,i) - (1.0/dx)*(ks+kt)*dz*Q2(1,coxl2) &
        + 20.0D0*kt*dx*(1.0D0/dz)*Q2(1,i) + 2.0D0*(ks+kt)*dz*(1.0D0/dx)*Q2(1,i) + 32.0D0*kt*(q_0**2)*dx*dz*Q2(1,i) &
        - 24.0D0*kt*dx*(1.0d0/dz)*Q2(2,i) + 4.0d0*kt*dx*(1.0d0/dz)*Q2(3,i) - (ks+kt)*dz*Q2(1,coxr2)*(1.0d0/dx) &
        + 8.0d0*kt*q_0*dz*(Q3(1,coxl1)-Q3(1,coxr1))  - 20.0d0*kt*q_0*dx*Q4(2,i) +4.0d0*kt*q_0*dx*Q4(3,i) &
        - 3.0d0*(ks+kt)*Q5(1, coxl1) + 4.0*ks*Q5(2,coxl1) - kt*Q5(2,coxl1) - ks*Q5(3,coxl1) &
        + 3.0d0*(ks+kt)*Q5(1,coxr1) - 4.0d0*ks*Q5(2,coxl1) + kt*Q5(2,coxr1) + ks*Q5(3,coxr1) )

        V(3 + ((i-1)*5) + ((k-1)*lx*5)) = V(3 + ((i-1)*5) + ((k-1)*lx*5)) &
        - (1.0D0/(4.0D0*dx*dz))*(-6.0D0*(ks+kt)*dx*dz*Q1(1,coxl1) &
        + 3.0D0*(ks+kt)*dx*dz*Q1(2,coxl1) - (ks+kt)*dx*dz*Q1(3,coxl1) &
        + 6.0D0*(ks+kt)*dx*dz*Q1(1,coxr1) - 3.0*(ks+kt)*dx*dz*Q1(2,coxr1) &
        + (ks+kt)*dx*dz*Q1(3,coxr1) + 8.0D0*kt*q_0*dx*(dz**2)*(Q2(1,coxl1)-Q2(1,coxr1)) &
        + (ks+kt)*(dz**2)*Q3(1,coxl2) - 10.0D0*(ks+kt)*(dx**2)*Q3(1,i) - 2.0D0*(ks+kt)*(dz**2)*Q3(1,i) &
        - 32.0D0*kt*(q_0**2)*(dx**2)*(dz**2)*Q3(1,i) + 12.0D0*(ks+kt)*(dx**2)*Q3(2,i) &
        - 2.0D0*(ks+kt)*(dx**2)*Q3(3,i) + (ks+kt)*(dz**2)*Q3(1,coxr2) - 3.0D0*(ks+kt)*dx*dz*Q4(1,coxl1) &
        + 4.0d0*ks*dx*dz*Q4(2,coxl1) - kt*dx*dz*Q4(2,coxl1) - ks*dx*dz*Q4(3,coxl1) + 3.0d0*(ks+kt)*dx*dz*Q4(1,coxr1) &
        - 4.0d0*ks*dx*dz*Q4(2,coxr1) + kt*dx*dz*Q4(2,coxr1) + ks*dx*dz*Q4(3,coxr1) &
        + 20.0d0*kt*q_0*(dx**2)*dz*Q5(2,i) - 4.0d0*kt*q_0*(dx**2)*dz*Q5(3,i) )

        V(4 + ((i-1)*5) + ((k-1)*lx*5)) = V(4 + ((i-1)*5) + ((k-1)*lx*5))&
        + (1.0D0/(4.0d0*dx*dz))*(-kt*(dz**2)*Q1(1,coxl2) + 2.0d0*Q1(1,i)*(5*ks*(dx**2) &
        + kt*(dz**2)*(1.0 + 8.0d0*(q_0**2)*(dx**2))) - 12.0d0*ks*(dx**2)*Q1(2,i) + 2.0d0*ks*(dx**2)*Q1(3,i) &
        - kt*(dz**2)*Q1(1,coxr2) + 20.0d0*kt*q_0*(dx**2)*dz*Q2(2,i) - 4.0d0*kt*q_0*(dx**2)*dz*Q2(3,i) &
        - 3.0d0*(ks+kt)*dx*dz*Q3(1,coxl1) - ks*dx*dz*Q3(2,coxl1) + 4.0d0*kt*dx*dz*Q3(2,coxl1) - kt*dx*dz*Q3(3,coxl1)&
        + 3.0d0*(ks+kt)*dx*dz*Q3(1,coxr1) + ks*dx*dz*Q3(2,coxr1) - 4.0d0*kt*dx*dz*Q3(2,coxr1) + kt*dx*dz*Q3(3,coxr1)&
        - 2.0d0*kt*(dz**2)*Q4(1,coxl2) + 10.0d0*(ks+kt)*(dx**2)*Q4(1,i) + 4.0d0*kt*(dz**2)*Q1(1,i)&
        + 32.0d0*kt*(q_0**2)*(dx**2)*(dz**2)*Q4(1,i) - 12.0d0*(ks+kt)*(dx**2)*Q4(2,i) + 2.0d0*(ks+kt)*(dx**2)*Q4(3,i)&
        + 2.0d0*kt*(dz**2)*Q4(1,coxr2) + 16.0d0*kt*q_0*dx*(dz**2)*Q5(1,coxl1) - 32.0d0*kt*q_0*dx*(dz**2)*Q5(1,coxr1) )

        V(5 + ((i-1)*5) + ((k-1)*lx*5)) = V(5 + ((i-1)*5) + ((k-1)*lx*5)) &
        + (1.0d0/4.0d0)*( -8.0d0*kt*q_0*dz*Q1(1,coxl1) + 8.0d0*kt*q_0*dz*Q1(1,coxr1) &
        + 3.0d0*kt*Q2(1,coxl1) - 4.0d0*kt*Q2(2,coxl1) + kt*Q2(3,coxl1) - 3.0d0*kt*Q2(1,coxr1) &
        + ks*(3.0d0*Q2(1,coxl1) + Q2(2,coxl1) - 3.0d0*Q2(1,coxr1) - Q2(2,coxr1)) + 4.0d0*kt*Q2(2,coxr1) &
        - kt*Q2(3,coxr1) + 20.0d0*kt*q_0*dx*Q3(2,i) - 4.0d0*kt*q_0*dx*Q3(3,i) &
        + 16.0d0*kt*q_0*dz*(Q4(1,coxr1) - Q4(1,coxl1)) - 2.0d0*(kt/dx)*dz*Q5(1,coxl2) + 10.0d0*(ks+kt)*(dx/dz)*Q5(1,i) &
        + 4.0d0*kt*(dz/dx)*Q5(1,i) + 32.0d0*kt*(q_0**2)*dx*dz*Q5(1,i) - 12.0d0*(ks+kt)*(dx/dz)*Q5(2,i) &
        + 2.0d0*(ks+kt)*(dx/dz)*Q5(3,i) - 2.0d0*kt*(dz/dx)*Q5(1,coxr2) )

      else if ( k  == 2 ) THEN
        !Elastic ENERGY
        V(1 + ((i-1)*5) + ((k-1)*lx*5)) = V(1 + ((i-1)*5) + ((k-1)*lx*5)) &
        + (-1.0D0/(4.0d0*dx*dz) )*( (ks+kt)*(dz**2)**Q1(2,coxl2) + 12.0d0*(ks+kt)*(dx**2)*Q1(1,i)&
        - 17.0d0*(ks+kt)*(dz**2)*Q1(2,i)- 2.0d0*(ks+kt)*Q1(2,i) - 32.0d0*kt*(q_0**2)*(dx**2)*(dz**2)*Q1(2,i)&
        + 4.0d0*(ks+kt)*(dx**2)*Q1(3,i) + (ks+kt)*(dx**2)*Q1(4,i) + (ks+kt)*(dz**2)*Q1(2,coxr2)&
        - 20.0d0*kt*q_0*(dx**2)*dz*Q2(1,i) + 8.0d0*kt*q_0*(dx**2)*dz*Q2(3,i) - 3.0d0*(ks+kt)*dx*dz*Q3(1,coxl1)&
        + 3.0d0*(ks+kt)*dx*dz*Q3(1,coxr1) + kt*(dz**2)*Q4(2,coxl2) + 12.0d0*ks*(dx**2)*Q4(1,i) &
        - 17.0d0*ks*(dx**2)*Q4(2,i) - 2.0d0*kt*(dz**2)*Q4(2,i) - 16.0d0*kt*(q_0**2)*(dx**2)*(dz**2)*Q4(2,i)&
        + 4.0d0*ks*(dx**2)*Q4(3,i) + ks*(dx**2)*Q4(4,i) + kt*(dz**2)*Q4(2,coxr2) - 8.0d0*kt*q_0*dx*(dz**2)*Q5(2,coxl1)&
        + 8.0d0*kt*q_0*dx*(dz**2)*Q5(2,coxr1) )

        V(2 + ((i-1)*5) + ((k-1)*lx*5)) = V(2 + ((i-1)*5) + ((k-1)*lx*5)) &
        + (1.0d0/4.0d0)*( -20.0d0*kt*q_0*dx*Q1(1,i) + 8.0d0*kt*q_0*dx*Q1(3,i) - (ks+kt)*(dz/dx)*Q2(2,coxl2)&
        - 24.0d0*kt*(dx/dz)*Q2(1,i) + 34.0d0*kt*(dx/dz)*Q2(2,i) + 2.0d0*(ks+kt)*(dz/dx)*Q2(2,i) + 32.0d0*kt*(q_0**2)*dx*dz*Q2(2,i)&
        - 8.0d0*kt*(dx/dz)*Q2(3,i) - 2.0d0*kt*(dx/dz)*Q2(4,i) - (ks+kt)*(dz/dx)*Q2(2,coxr2)&
        + 8.0d0*kt*q_0*dz*(Q3(2,coxl2)-Q3(2,coxr2))&
        + 20.0d0*kt*q_0*dx*Q4(1,i) - 8.0*kt*q_0*dx*Q4(3,i) - ks*Q5(1,coxl1) + 4.0d0*kt*Q5(1,coxl1)&
        + (ks-kt)*Q5(3,coxl1) + ks*Q5(1,coxr1) - 4.0d0*kt*Q5(1,coxr1) + (kt-ks)*Q5(3,coxr1) )

        V(3 + ((i-1)*5) + ((k-1)*lx*5)) = V(3 + ((i-1)*5) + ((k-1)*lx*5)) &
        + (-1.0d0/(4.0d0*dx*dz))*( 3.0d0*(ks+kt)*dx*dz*Q1(1,coxl1) - 3.0d0*(ks+kt)*dx*dz*Q1(1,coxr1)&
        + 8.0d0*kt*q_0*dx*(dz**2)*(Q2(2,coxl1) - Q2(2,coxr1)) + (ks+kt)*(dz**2)*Q3(2,coxl2)&
        + 12.0d0*(kt+ks)*(dx**2)*Q3(1,i) - 17.0d0*(ks+kt)*(dx**2)*Q3(2,i) - 2.0d0*(ks+kt)*(dz**2)*Q3(2,i)&
        - 32.0d0*kt*(q_0**2)*(dx**2)*(dz**2)*Q3(2,i) + 4.0d0*(ks+kt)*(dx**2)*Q3(3,i) + (ks+kt)*(dx**2)*Q3(4,i)&
        + (ks+kt)*(dz**2)*Q3(2,coxr2) - ks*dx*dz*Q4(1,coxl1) + 4.0d0*kt*dx*dz*Q4(1,coxl1) + (ks-kt)*dx*dz*Q4(1,coxr1)&
        + 4.0d0*kt*dx*dz*Q4(1,coxl1) + (ks-kt)*dx*dz*Q4(3,coxl1) + ks*dx*dz*Q4(1,coxr1) - 4.0d0*kt*dx*dz*Q4(1,coxr1)&
        + (kt-ks)*dx*dz*Q4(3,coxr1) - 20.0d0*kt*q_0*(dx**2)*dz*Q5(1,i) + 8.0d0*kt*q_0*(dx**2)*dz*Q5(3,i) )

        V(4 + ((i-1)*5) + ((k-1)*lx*5)) = V(4 + ((i-1)*5) + ((k-1)*lx*5)) &
        + (-1.0d0/(4.0d0*dx*dz))*( kt*(dz**2)*Q1(2,coxl2) + ks*(dx**2)*(12.0d0*Q1(1,i) - 17.0d0*Q1(2,i))&
        - 2.0*kt*(1.0d0 + 8.0d0*((q_0*dx)**2))*(dz**2)*Q1(2,i) + 4.0d0*ks*(dx**2)*Q1(3,i) + ks*(dx**2)*Q1(4,i)&
        + kt*(dz**2)*Q1(2,coxr2) + 20.0d0*kt*q_0*(dx**2)*dz*Q2(1,i) - 8.0d0*kt*q_0*(dx**2)*dz*Q2(3,i)&
        - 4.0d0*ks*dx*dz*Q3(1,coxl1) + (kt+ks)*dx*dz*Q3(1,coxl1) - kt*dx*dz*Q3(3,coxl1) + 4.0d0*ks*dx*dz*Q3(1,coxr1)&
        - (kt-ks)*dx*dz*Q3(1,coxr1) + kt*dx*dz*Q3(3,coxr1) + 2.0d0*kt*(dz**2)*Q4(2,coxl2)&
        + 12.0d0*(ks+kt)*(dx**2)*Q4(1,i) - 17.0d0*(ks+kt)*(dx**2)*Q4(2,i) - 4.0d0*kt*(dz**2)*Q4(2,i)&
        - 32.0d0*kt*((q_0*dx*dz)**2)*Q4(2,i) + 4.0d0*(kt+ks)*(dx**2)*Q4(3,i) + (ks+kt)*(dx**2)*Q4(4,i)&
        + 2.0d0*kt*(dz**2)*Q4(2,coxr2) + 16.0d0*kt*q_0*dx*(dz**2)*(Q5(2,coxr1)-Q5(2,coxl1)) )

        V(5 + ((i-1)*5) + ((k-1)*lx*5)) = V(5 + ((i-1)*5) + ((k-1)*lx*5))&
        + (1.0d0/4.0d0)*( 8.0d0*kt*q_0*dz*(Q1(2,coxr1) - Q1(2,coxl1))&
        + kt*(Q2(1,coxl1) - Q2(3,coxl1) - Q2(1,coxr1)) + ks*(-4.0d0*Q2(1,coxl1) + Q2(3,coxl1) + 4.0d0*Q2(1,coxr1)&
        - Q2(3,coxr1)) + kt*Q2(3,coxr1) - 20.0d0*kt*q_0*dx*Q3(1,i) + 8.0d0*kt*q_0*dx*Q3(3,i)&
        + 16.0d0*kt*q_0*dz*(Q4(2,coxr1)-Q4(2,coxl1)) - 2.0d0*kt*(dz/dx)*Q5(2,coxl2) - 12.0d0*(kt+ks)*(dx/dz)*Q5(1,i)&
        + 17.0d0*ks*(dx/dz)*Q5(2,i) + 4.0d0*kt*(dz/dx)*Q5(2,i) + 32.0d0*kt*(q_0**2)*dx*dz*Q5(2,i)&
        - 4.0d0*(ks+kt)*(dx/dz)*Q5(3,i) - (ks+kt)*(dx/dz)*Q5(4,i) - 2.0d0*kt*(dz/dx)*Q5(2,coxr2) )

      else if ( k == 3 ) THEN
        V(1 + ((i-1)*5) + ((k-1)*lx*5)) = V(1 + ((i-1)*5) + ((k-1)*lx*5)) &
        - (1.0d0/(4.0d0*dx*dz))*( (ks+kt)*(dz**2)*Q1(3,coxl2) - 2.0d0*(ks+kt)*(dx**2)*Q1(1,i)&
        + 4.0d0*(ks+kt)*(dx**2)*Q1(2,i) - 3.0d0*(ks+kt)*(dx**2)*Q1(3,i) - 2.0d0*(ks+kt)*(dz**2)*Q1(3,i)&
        - 32.0d0*kt*((q_0*dx*dz)**2)*Q1(3,i) + (ks+kt)*(dx**2)*Q1(5,i) + (ks+kt)*(dz**2)*Q1(3,coxr2)&
        + 4.0d0*kt*q_0*(dx**2)*dz*Q2(1,i) - 8.0d0*kt*q_0*(dx**2)*dz*Q2(2,i) + 8.0d0*kt*q_0*(dx**2)*dz*Q2(4,i)&
        + (ks+kt)*dx*dz*(Q3(1,coxl1) - Q3(1,coxr1)) + kt*(dz**2)*Q4(3,coxl2) - 2.0d0*ks*(dx**2)*Q4(1,i)&
        + 4.0d0*ks*(dx**2)*Q4(2,i) - 3.0d0*ks*(dx**2)*Q3(3,i) - 2.0d0*kt*(dz**2)*Q4(3,i)&
        - 16.0d0*kt*((q_0*dx*dz)**2)*Q4(3,i) + ks*(dx**2)*Q4(5,i) + kt*(dz**2)*Q4(3,coxr2)&
        + 8.0d0*kt*q_0*dx*(dz**2)*(Q5(3,coxr1)-Q5(3,coxl1)) )

        V(2 + ((i-1)*5) + ((k-1)*lx*5)) = V(2 + ((i-1)*5) + ((k-1)*lx*5)) &
        + (1.0d0/4.0d0)*( 4.0d0*kt*q_0*dx*Q1(1,i) + 8.0d0*kt*q_0*dx*(Q1(4,i)-Q1(2,i))&
        - (ks+kt)*(dz/dx)*Q2(3,coxl2) + 4.0d0*kt*(dx/dz)*Q2(1,i) - 8.0d0*kt*(dx/dz)*Q2(2,i)&
        + 6.0d0*kt*(dx/dz)*Q2(3,i) + 2.0d0*(ks+kt)*(dz/dx)*Q2(3,i) + 32.0d0*kt*(q_0**2)*dx*dz*Q2(3,i)&
        - 2.0d0*kt*(dx/dz)*Q2(5,i) - (ks+kt)*(dz/dx)*Q2(3,coxr2) + 8.0d0*kt*q_0*dz*(Q3(3,coxl1) - Q3(3,coxr1)) &
        - 4.0d0*kt*q_0*dx*Q4(1,i) + 8.0d0*kt*q_0*dx*(Q4(2,i) - Q4(4,i)) - kt*Q5(1,coxl1)&
        + (kt-ks)*Q5(2,coxl1) + kt*Q5(2,coxl1) + (ks-kt)*Q5(4,coxl1) + kt*Q5(2,coxl1) + (ks-kt)*Q5(4,coxl1)&
        + kt*Q5(1,coxr1) + (ks-kt)*Q5(2,coxr1) + (kt-ks)*Q5(4,coxr1) )

        V(3 + ((i-1)*5) + ((k-1)*lx*5)) = V(3 + ((i-1)*5) + ((k-1)*lx*5)) &
        + (-1.0d0/(4.0d0*dx*dz))*( -(ks+kt)*dx*dz*Q1(1,coxl1) + (ks+kt)*dx*dz*Q1(1,coxr1)&
        + 8.0d0*kt*q_0*dx*(dz**2)*(Q2(3,coxl1) - Q2(3,coxr1)) + (ks+kt)*(dz**2)*Q3(3,coxl2)&
        - 2.0d0*(ks+kt)*(dx**2)*Q3(1,i) + 4.0d0*(ks+kt)*(dx**2)*Q3(2,i) - 3.0d0*(ks+kt)*(dx**2)*Q3(3,i)&
        - 2.0d0*(ks+kt)*(dz**2)*Q3(3,i) - 32.0d0*kt*((q_0*dx*dz)**2)*Q3(3,i) + (ks+kt)*(dx**2)*Q3(5,i)&
        + (ks+kt)*(dz**2)*Q3(3,coxr2) - kt*dx*dz*Q4(1,coxl1) + (kt-ks)*dx*dz*Q4(2,coxl1)&
        + (ks-kt)*dx*dz*Q4(4,coxl1) + kt*dx*dz*Q4(1,coxr1) + (ks-kt)*dx*dz*Q4(2,coxr1) + (kt-ks)*dx*dz*Q4(4,coxr1)&
        + 4.0d0*kt*q_0*(dx**2)*dz*Q5(1,i) + 8.0d0*kt*q_0*(dx**2)*dz*(Q5(4,i) - Q5(2,i)) )

        V(4 + ((i-1)*5) + ((k-1)*lx*5)) = V(4 + ((i-1)*5) + ((k-1)*lx*5))&
        - (1.0d0/(4.0d0*dx*dz))*( kt*(dz**2)*Q1(3,coxl2) + ks*(dx**2)*(-2.0d0*Q1(1,i) + 4.0d0*Q1(2,i) -3.0d0*Q1(3,i))&
        + kt*(dz**2)*Q1(3,coxr2) - 4.0d0*kt*q_0*(dx**2)*dz*Q2(1,i) + 8.0d0*kt*q_0*(dx**2)*dz*(Q2(2,i) - Q2(4,i))&
        + ks*dx*dz*Q3(1,coxl1) + (kt-ks)*dx*dz*Q3(2,coxl1) + (ks-kt)*dx*dz*Q3(4,coxl1) - ks*dx*dz*Q3(1,coxr1)&
        + (ks-kt)*dx*dz*Q3(2,coxr1) + (kt-ks)*dx*dz*Q3(4,coxl1) - 2.0d0*kt*(dz**2)*Q4(3,coxl2)&
        - 2.0d0*(ks+kt)*(dx**2)*Q4(1,i) - 4.0d0*(ks+kt)*(dx**2)*Q4(2,i) - 3.0d0*(ks+kt)*(dx**2)*Q4(3,i)&
        - 4.0d0*kt*(dz**2)*Q4(3,i) - 32.0d0*kt*((q_0*dx*dz)**2)*Q4(3,i) + (ks+kt)*(dx**2)*Q4(5,i)&
        + 2.0d0*kt*(dz**2)*Q4(3,coxr2) + 16.0d0*kt*q_0*dx*(dz**2)*(Q5(3,coxr1) - Q5(3,coxl1)) )

        V(5 + ((i-1)*5) + ((k-1)*lx*5)) = V(5 + ((i-1)*5) + ((k-1)*lx*5)) &
        -(1.0d0/(4.d0*dx*dz))*( 8.0d0*kt*q_0*dx*(dz**2)*(Q1(3,coxl1) - Q1(3,coxr1)) - ks*dx*dz*Q2(1,coxl1)&
        + (ks-kt)*dx*dz*Q2(2,coxl1) + (kt-ks)*dx*dz*Q2(4,coxl1) + ks*dx*dz*Q2(1,coxl1)&
        + (kt-ks)*dx*dz*Q2(2,coxr1) + (ks-kt)*dx*dz*Q2(4,coxr1) - 4.0d0*kt*q_0*(dx**2)*dz*Q3(1,i)&
        + 8.0d0*kt*q_0*(dx**2)*dz*(Q3(2,i) - Q3(4,i)) + 16.0d0*kt*q_0*dx*(dz**2)*(Q4(3,coxl1) - Q4(3,coxr1))&
        + 2.0d0*(kt*(dz**2)*Q5(3,coxl2) - (ks+kt)*(dx**2)*Q5(1,i)) + 4.0d0*(ks+kt)*(dx**2)*Q5(2,i)&
        - 3.0d0*(ks+kt)*(dx**2)*Q5(3,i) - 4.0d0*kt*(dz**2)*Q5(3,i) - 32.0d0*kt*((q_0*dx*dz)**2)*Q5(3,i)&
        + (ks+kt)*(dx**2)*Q5(5,i) + 2.0d0*kt*(dz**2)*Q5(3,coxr2) )

      else if ( k == lz ) THEN
        !!Surface Energy & Bulk Energy with backward difference
        !!!!SURFACE ENERGY!!!
        V(1 + ((i-1)*5) + ((k-1)*lx*5)) = V(1 + ((i-1)*5) + ((k-1)*lx*5)) &
        - (1.0d0/(4.0d0*dx*dz))*( (ks+kt)*(dz**2)*Q1(lz,coxl2) - 2.0d0*(ks+kt)*(dx**2)*Q1(lz-2,i)&
        + 12.0d0*(kt+ks)*(dx**2)*Q1(lx-1,i)  - 10.0d0*(kt+ks)*(dx**2)*Q1(lz,i) - 2.0d0*(kt+ks)*(dz**2)*Q1(lz,i)&
        - 32.0d0*kt*((q_0*dx*dz)**2)*Q1(lz,i) + (ks+kt)*(dz**2)*Q1(lz,coxr2)&
        + kt*(q_0)*(dx**2)*dz*(4.0d0*Q2(lz-2,i) - 20.0*Q2(lz-1,i)) - (ks+kt)*dx*dz*Q3(lz-2,coxl1)&
        + 3.0d0*(ks+kt)*dx*dz*Q3(lx-1,coxl1) - 6.0d0*(kt+ks)*dx*dz*Q3(lz,coxl1) + (kt+ks)*dx*dz*Q3(lz-2,coxr1)&
        - 3.0d0*(ks+kt)*dx*dz*Q3(lz-1,coxr1) + 6.0d0*(kt+ks)*dx*dz*Q3(lz,coxr1) + kt*(dz**2)*Q4(lz,coxl2)&
        - 2.0d0*ks*(dx**2)*Q4(lz-2,i) + 12.0d0*ks*(dx**2)*Q4(lz-1,i) - 10.0d0*ks*(dx**2)*Q4(lz,i)&
        - 2.0d0*kt*(dz**2)*Q4(lz,i) - 16.0d0*kt*((q_0*dx*dz)**2)*Q4(lz,i) + kt*(dz**2)*Q4(lz,coxr2)&
        + 8.0d0*kt*q_0*dx*(dz**2)*(Q5(lz,coxr1) - Q5(lz,coxl1)) )

        V(2 + ((i-1)*5) + ((k-1)*lx*5)) = V(2 + ((i-1)*5) + ((k-1)*lx*5)) &
        + (1.0d0/4.0d0)*( kt*q_0*dx*(4.0d0*Q1(lz-2,i) - 20.0d0*Q1(lz-1,i)) - (ks+kt)*(dz/dx)*Q2(lz,coxl2)&
        + 4.0d0*kt*(dx/dz)*Q2(lz-2,i) - 24.0d0*kt*(dx/dz)*Q2(lz-1,i) + 20.0d0*kt*(dx/dz)*Q2(lz,i)&
        + 2.0d0*(kt+ks)*(dz/dx)*Q2(lz,i) + 32.0d0*kt*(q_0**2)*dx*dz*Q2(lz,i) - (ks+kt)*(dz/dx)*Q2(lz,coxr2)&
        + 8.0d0*kt*q_0*dz*(Q3(lz,coxl1) - Q3(lz,coxr1)) - 4.0d0*kt*q_0*dx*Q4(lz-2,i) + 20.0d0*kt*q_0*dx*Q4(lz-1,i)&
        + ks*Q5(lz-2,coxl1) + 4.0d0*ks*Q5(lz-1,coxl1) + kt*Q5(lz-1,coxl1) + 3.0d0*(ks+kt)*Q5(lz,coxl1)&
        - ks*Q5(lz-2,coxr1) + 4.0d0*ks*Q5(lz-1,coxr1) - kt*Q5(lz-1,coxr1) - 3.0d0*(ks+kt)*Q5(lz,coxr1) )

        V(3 + ((i-1)*5) + ((k-1)*lx*5)) = V(3 + ((i-1)*5) + ((k-1)*lx*5)) &
        - (1.0d0/(4.0d0*dx*dz))*( (ks+kt)*dx*dz*Q1(lz-2,coxl1) - 3.0d0*(ks+kt)*dx*dz*Q1(lz-1,coxl1)&
        + 6.0d0*(ks+kt)*dx*dz*Q1(lz,coxl1) - (ks+kt)*dx*dz*Q1(lz-2,coxr1) + 3.0d0*(ks+kt)*dx*dz*Q1(lz-1,coxr1)&
        - 6.0d0*(ks+kt)*dx*dz*Q1(lz,coxr1) + 8.0d0*kt*q_0*dx*(dz**2)*(Q2(lz,coxl1)-Q2(lz,coxr1))&
        + (ks+kt)*(dz**2)*Q3(lz,coxl2) - 2.0d0*(ks+kt)*(dx**2)*Q3(lz-2,i) + 12.0d0*(ks+kt)*(dx**2)*Q3(lz-1,i)&
        - 10.0d0*(ks+kt)*(dx**2)*Q3(lz,i) - 2.0d0*(ks+kt)*(dz**2)*Q3(lz,i) - 32.0d0*kt*((q_0*dx*dz)**2)*Q3(lz,i)&
        + (ks+kt)*(dz**2)*Q3(lz,coxr2) + ks*dx*dz*Q4(lz-2,coxl1) - 4.0d0*ks*dx*dz*Q4(lz-1,coxl1)&
        + kt*dx*dz*Q4(lz-1,coxl1) + 3.0d0*(ks+kt)*dx*dz*Q4(lz,coxl1) - ks*dx*dz*Q4(lz-2,coxr1)&
        + 4.0d0*ks*dx*dz*Q4(lz-1,coxr1) - kt*dx*dz*Q4(lz-1,coxr1) - 3.0d0*(ks+kt)*dx*dz*Q4(lz,coxr1)&
        + 4.0d0*kt*q_0*(dx**2)*dz*Q5(lz-2,i) - 20.0d0*kt*q_0*(dx**2)*dz*Q5(lz-1,i) )

        V(4 + ((i-1)*5) + ((k-1)*lx*5)) = V(4 + ((i-1)*5) + ((k-1)*lx*5))&
        + (1.0d0/(4.0d0*dx*dz))*( -kt*(dz**2)*Q1(lz,coxl2) + 2.0d0*kt*(1.0d0 + 8.0d0*((q_0*dx)**2))*(dz**2)*Q1(lz,i)&
        + 2.0d0*ks*(dx**2)*(Q1(lz-2,i) - 6.0d0*Q1(lz-1,i) + 5.0d0*Q1(lz,i)) - kt*(dz**2)*Q1(lz,coxr2)&
        + 4.0d0*kt*q_0*(dx**2)*dz*Q2(lz-2,i) - 20.0d0*kt*q_0*(dx**2)*dz*Q2(lz-1,i) + kt*dx*dz*Q3(lz-2,coxl1)&
        + ks*dx*dz*Q3(lz-1,coxl1) - 4.0d0*kt*dx*dz*Q3(lz-1,coxl1) + 3.0d0*(kt+ks)*dx*dz*Q3(lz,coxl1)&
        - kt*dx*dz*Q3(lz-2,coxr1) + (4.0d0*kt-ks)*dx*dz*Q3(lz-1,coxr1) - 3.0d0*(ks+kt)*dx*dz*Q3(lz,coxr1)&
        - 2.0d0*kt*(dz**2)*Q4(lz,coxl2) + 2.0d0*(ks+kt)*(dx**2)*Q4(lz-2,i) - 12.0d0*(ks+kt)*(dx**2)*Q4(lz-1,i)&
        + 10.0d0*(ks+kt)*(dx**2)*Q4(lz,i) + 4.0d0*kt*(dz**2)*Q4(lz,i) + 32.0d0*kt*((q_0*dx*dz)**2)*Q4(lz,i)&
        - 2.0d0*kt*(dz**2)*Q4(lz,coxr2) + 16.0d0*kt*q_0*dx*(dz**2)*(Q5(lz,coxl1) - Q5(lz,coxr1)) )

        V(5 + ((i-1)*5) + ((k-1)*lx*5)) = V(5 + ((i-1)*5) + ((k-1)*lx*5)) &
        + (1.0d0/4.0d0)*( 8.0d0*kt*q_0*dz*(Q1(lz,coxr1) - Q1(lz,coxl1)) - kt*Q2(lz-2,coxl1) + 4.0d0*kt*Q2(lz-1,coxl1)&
        - 3.0d0*kt*Q2(lz,coxl1) + kt*Q2(lz-2,coxr1) - 4.0d0*kt*Q2(lz-1,coxr1) + 3.0d0*kt*Q2(lz,coxr1)&
        + ks*(-Q2(lz-1,coxl1) - 3.0d0*Q2(lz,coxl1) + Q2(lz-1,coxr1) + 3.0d0*Q2(lz,coxr1)) + 4.0d0*kt*q_0*dx*Q3(lz-2,i)&
        - 20.0d0*kt*q_0*dx*Q3(lz-1,i) + 16.0d0*kt*q_0*dz*(Q4(lz,coxr1) - Q4(lz,coxl1))&
        - 2.0d0*kt*(dz/dx)*Q5(lz,coxl2) + 2.0d0*(ks+kt)*(dx/dz)*Q5(lz-2,i) - 12.0d0*(ks+kt)*(dx/dz)*Q5(lz-1,i)&
        + 10.0d0*(ks+kt)*(dx/dz)*Q5(lz-1,i) + 4.0d0*kt*(dz/dx)*Q5(lz,i) + 32.0d0*kt*(q_0**2)*dx*dz*Q5(lz,i)&
        - 2.0d0*kt*(dz/dx)*Q5(lz,coxr2) )

      else if ( k == lz - 1 ) then
        V(1 + ((i-1)*5) + ((k-1)*lx*5)) = V(1 + ((i-1)*5) + ((k-1)*lx*5)) &
        - (1.0d0/(4.0d0*dx*dz))*( (ks+kt)*(dz**2)*Q1(lz-1,coxl2) + (ks+kt)*(dx**2)*Q1(lz-3,i)&
        + 4.0d0*(ks+kt)*(dx**2)*Q1(lz-2,i) - 17.0d0*(ks+kt)*Q1(lz-1,i) - 2.0d0*(ks+kt)*(dz**2)*Q1(lz-1,i)&
        - 32.0d0*kt*((q_0*dx*dz)**2)*Q1(lz-1,i) + 12.0d0*(ks+kt)*(dx**2)*Q1(lz,i) + (ks+kt)*(dz**2)*Q1(lz-1,coxr2)&
        - 8.0d0*kt*q_0*(dx**2)*dz*Q2(lz-2,i) + 20.0d0*kt*q_0*(dx**2)*dz*Q2(lz,i) + 3.0d0*(kt+ks)*dx*dz*Q3(lz,coxl1)&
        - 3.0d0*(kt+ks)*dx*dz*Q3(lz,coxr1) + kt*(dz**2)*Q4(lz-1,coxl2) + ks*(dx**2)*Q4(lz-3,i)&
        + 4.0d0*ks*(dx**2)*Q4(lz-2,i) - 17.0d0*ks*(dx**2)*Q4(lz-1,i) - 2.0d0*kt*(dz**2)*Q4(lz-1,i)&
        - 16.0d0*kt*((q_0*dx*dz)**2)*Q4(lz-1,i) + 12.0d0*ks*(dx**2)*Q4(lz,i) + kt*(dz**2)*Q4(lz-1,coxr2)&
        + 8.0d0*kt*q_0*dx*(dz**2)*(Q5(lz-1,coxr1) - Q5(lz-1,coxl1)) )

        V(2 + ((i-1)*5) + ((k-1)*lx*5)) = V(2 + ((i-1)*5) + ((k-1)*lx*5)) &
        + (1.0d0/4.0d0)*( - 8.0d0*kt*q_0*dx*Q1(lz-2,i) + 20.0d0*kt*q_0*dx*Q1(lz,i)&
        - (ks+kt)*(dz/dx)*Q2(lz-1,coxl2) + kt*(dx/dz)*( - 2.0d0*Q2(lz-3,i) - 8.0d0*Q2(lz-2,i) + 34.0d0*Q2(lz-1,i))&
        + 2.0d0*(ks+kt)*(dz/dx)*Q2(lz-1,i) + 32.0d0*kt*(q_0**2)*dx*dz*Q2(lz-1,i) - 24.0d0*kt*(dx/dz)*Q2(lz,i)&
        - (ks+kt)*(dz/dx)*Q2(lz-2,coxr2) + 8.0d0*kt*q_0*dz*(Q3(lz-1,coxl1) - Q3(lz-1,coxr1))&
        + 8.0d0*kt*dx*Q4(lz-2,i) - 20.0d0*kt*q_0*dx*Q4(lz,i) + (kt-ks)*Q5(lz-2,coxl1) + (ks - 4.0d0*kt)*Q5(lz,coxl1)&
        + (ks-kt)*Q5(lz-2,coxr1) + (4.0d0*kt - ks)*Q5(lz,coxr1) )

        V(3 + ((i-1)*5) + ((k-1)*lx*5)) = V(3 + ((i-1)*5) + ((k-1)*lx*5)) &
        - (1.0d0/(4.0d0*dx*dz))*( - 3.0d0*(ks+kt)*dx*dz*Q1(lz,coxl1) + 3.0d0*(ks+kt)*dx*dz*Q1(lz,coxr1)&
        + 8.0d0*kt*q_0*dx*(dz**2)*(Q2(lz-1,coxl1) - Q2(lz-1,coxr1)) + (ks+kt)*(dz**2)*Q3(lz-1,coxl2)&
        + (ks+kt)*(dx**2)*Q3(lz-3,i) + 4.0d0*(kt+ks)*(dx**2)*Q3(lz-2,i) - 17.0d0*(ks+kt)*(dx**2)*Q3(lz-1,i)&
        - 2.0d0*(ks+kt)*(dz**2)*Q3(lz-1,i) - 32.0d0*kt*((q_0*dx*dz)**2)*Q3(lz-1,i) + 12.0d0*(ks+kt)*(dx**2)*Q3(lz,i)&
        + (ks+kt)*(dz**2)*Q3(lz-1,coxr2) + (kt-ks)*dx*dz*Q4(lz-2,coxl1) + (ks - 4.0d0*kt)*dx*dz*Q4(lz,coxl1)&
        + (ks-kt)*dx*dz*Q4(lz-2,coxr1) + (4.0d0*kt - ks)*dx*dz*Q4(lz,coxr1) - 8.0d0*kt*q_0*(dx**2)*dz*Q5(lz-2,i)&
        + 20.0d0*kt*q_0*(dx**2)*dz*Q5(lz,i) )

        V(4 + ((i-1)*5) + ((k-1)*lx*5)) = V(4 + ((i-1)*5) + ((k-1)*lx*5)) &
        -(1.0D0/(4.0d0*dx*dz))*( kt*(dz**2)*Q1(lz-1,coxl2) + ks*(dx**2)*(Q1(lz-3,i) + 4.0d0*Q1(lz-2,i) - 17.0d0*Q1(lz-1,i))&
        - 2.0d0*kt*(1.0d0 + 8.0d0*((q_0*dx)**2))*(dz**2)*Q1(lz-1,coxr2)&
        + kt*q_0*(dx**2)*dz*(8.0d0*Q2(lz-2,i) - 20.0d0*Q2(lz,i))&
        + (kt-ks)*dx*dz*Q3(lz-2,coxl1) + (4.0d0*ks-kt)*dx*dz*Q3(lz,coxl1)&
        + 12.0d0*ks*(dx**2)*Q1(lz,i) + kt*(dz**2)*Q1(lz-1,coxr2)&
        + (ks-kt)*dx*dz*Q3(lz-2,coxr2) + (kt-4.0d0*ks)*dx*dz*Q3(lz,coxr1) + 2.0d0*kt*(dz**2)*Q4(lz-1,coxl2)&
        + (ks+kt)*(dx**2)*Q4(lz-3,i) + 4.0d0*(ks+kt)*(dx**2)*Q4(lz-2,i) - 17.0d0*(ks+kt)*(dx**2)*Q4(lz-1,i)&
        - kt*Q4(lz-1,i)*( 4.0d0*(dz**2) + 32.0d0*kt*((q_0*dx*dz)**2) ) + 12.0d0*(ks+kt)*(dx**2)*Q4(lz,i)&
        + 2.0d0*kt*(dz**2)*Q4(lz-1,coxr2) + 16.0d0*kt*q_0*dx*(dz**2)*(Q5(lz-1,coxr1) - Q5(lz-1,coxl1)) )

        V(5 + ((i-1)*5) + ((k-1)*lx*5)) = V(5 + ((i-1)*5) + ((k-1)*lx*5)) &
        + (1.0d0/4.0d0)*( 8*kt*q_0*dz*(Q1(lz-1,coxr1) - Q1(lz-1,coxl1)) + kt*Q2(lz-2,coxl1)&
        + kt*(Q2(lz-2,coxl1) - Q2(lz-2,coxr1)) + ks*(-Q2(lz-2,coxl1) + Q2(lz-2,coxr1)&
        + 4.0d0*Q2(lz,coxl1) - 4.0d0*Q2(lz,coxr1)) + kt*Q2(lz,coxr1) - 8.0d0*kt*q_0*dx*Q3(lz-2,i)&
        + 20.0d0*kt*q_0*dx*Q3(lz,i) + 16.0d0*kt*q_0*dz*(Q4(lz-1,coxr1)-Q4(lz-1,coxl1)) - 2.0d0*kt*(dz/dx)*Q5(lz-1,coxl2)&
        - (kt+ks)*(dx/dz)*Q5(lz-3,i) - 4.0d0*(kt+ks)*(dx/dz)*Q5(lz-2,i) + 17.0d0*(ks+kt)*(dx/dz)*Q5(lz-1,i)&
        + 4.0d0*kt*(dz/dx)*Q5(lz-1,i) + 32.0d0*kt*(q_0**2)*dx*dz*Q5(lz-1,i) - 12.0d0*(kt+ks)*(dx/dz)*Q5(lz,i)&
        - 2.0d0*kt*(dz/dx)*Q5(lz-1,coxr2) )

      else if (k == lz-2) then
        V(1 + ((i-1)*5) + ((k-1)*lx*5)) = V(1 + ((i-1)*5) + ((k-1)*lx*5))&
        -(1.0d0/(4.0d0*dx*dz))*( (ks+kt)*(dz**2)*Q1(lz-2,coxl2) + (ks+kt)*(dx**2)*Q1(lz-4,i) - 3.0d0*(ks+kt)*(dx**2)*Q1(lz-2,i)&
        - 2.0d0*(ks+kt)*(dz**2)*Q1(lz-2,i) - 32.0d0*kt*((q_0*dx*dz)**2)*Q1(lz-2,i) + 4.0d0*(ks+kt)*(dx**2)*Q1(lz-1,i)&
        - 2.0d0*(ks+kt)*(dx**2)*Q1(lz,i) + (ks+kt)*(dz**2)*Q1(lz-2,coxr2) + 8.0d0*kt*q_0*(dx**2)*dz*(Q2(lz-1,i)-Q2(lz-3,i))&
        - 4.0d0*kt*q_0*(dx**2)*dz*Q2(lz,i) - (kt+ks)*dx*dz*Q3(lz,coxl1) + (kt+ks)*dx*dz*Q3(lz,coxr1) + kt*(dz**2)*Q4(lz-2,coxl2)&
        + ks*(dx**2)*Q4(lz-4,i) - 3.0d0*ks*(dx**2)*Q4(lz-2,i) - 2.0d0*kt*(dz**2)*Q4(lz-2,i) - 16.0d0*kt*((q_0*dx*dz)**2)*Q4(lz-2,i)&
        + 4.0d0*ks*(dx**2)*Q4(lz-1,i) - 2.0d0*ks*(dx**2)*Q4(lz,i) + kt*(dz**2)*Q4(lz-2,coxr2) &
        + 8.0d0*kt*q_0*dx*(dz**2)*(Q5(lz-2,coxr1) - Q5(lz-2,coxl1)) )

        V(2 + ((i-1)*5) + ((k-1)*lx*5)) = V(2 + ((i-1)*5) + ((k-1)*lx*5)) &
        + (1.0d0/4.0d0)*( 8.0d0*kt*q_0*dx*(Q1(lz-1,i)-Q1(lz-3,i)) - 4.0d0*kt*q_0*dx*Q1(lz,i)&
        - (ks+kt)*(dz/dx)*Q2(lz-2,coxl2) - 2.0d0*kt*(dx/dz)*Q2(lz-4,i) + 6.0d0*kt*(dx/dz)*Q2(lz-2,i)&
        + 2.0d0*(ks+kt)*(dz/dx)*Q2(lz-2,i) + 32.0d0*kt*(q_0**2)*dx*dz*Q2(lz-2,i) - kt*(dx/dz)*(4.0d0*Q2(lz,i) - 8.0d0*Q2(lz-1,i))&
        - (ks+kt)*(dz/dx)*Q2(lz-2,coxr2) + 8.0d0*kt*q_0*dz*(Q3(lz-2,coxr1) - Q3(lz-2,coxr1))&
        + 8.0d0*q_0*kt*dx*(Q4(lz-3,i) - Q4(lz-1,i)) + 4.0d0*kt*q_0*dx*Q4(lz,i)&
        + (kt-ks)*Q5(lz-3,coxl1) + (ks-kt)*Q5(lz-1,coxl1) + kt*Q5(lz,coxl1) + (ks-kt)*Q5(lz-2,coxr1)&
        + (kt-ks)*Q5(lz-1,coxr1) - kt*Q5(lz,coxr1)  )

        V(3 + ((i-1)*5) + ((k-1)*lx*5)) = V(3 + ((i-1)*5) + ((k-1)*lx*5)) &
        - (1.0d0/(4.0d0*dx*dz))*( (ks+kt)*dx*dz*(Q1(lz,coxl1) - Q1(lz,coxr1))&
        + 8.0d0*kt*q_0*dx*(dz**2)*(Q2(lz-2,coxl1) - Q2(lz-2,coxr1)) + (ks+kt)*(dz**2)*Q3(lz-2,coxl2)&
        + (ks+kt)*(dx**2)*Q3(lz-4,i) - 3.0d0*(kt+ks)*(dx**2)*Q3(lz-2,i) - 2.0d0*(kt+ks)*(dz**2)*Q3(lz-2,i)&
        - 32.0d0*kt*((q_0*dx*dz)**2)*Q3(lz-2,i) + 4.0d0*(ks+kt)*(dx**2)*Q3(lz-1,i) - 2.0d0*(ks+kt)*(dx**2)*Q3(lz,i)&
        + (ks+kt)*(dz**2)*Q3(lz-2,coxr2) + (kt-ks)*dx*dz*Q4(lz-3,coxl1) + (ks-kt)*dx*dz*Q4(lz-1,coxl1)&
        + kt*dx*dz*Q4(lz,coxl1) + (ks-kt)*dx*dz*Q4(lz-3,coxr1) + (kt-ks)*dx*dz*Q4(lz-1,coxr1)&
        - kt*dx*dz*Q4(lz,coxr1) - 8.0d0*kt*q_0*(dx**2)*dz*(Q5(lz-1,i)-Q5(lz-3,i)) - 4.0d0*kt*q_0*(dx**2)*dz*Q5(lz,i) )

        V(4 + ((i-1)*5) + ((k-1)*lx*5)) = V(4 + ((i-1)*5) + ((k-1)*lx*5)) &
        - (1.0d0/(4.0d0*dx*dz))*( kt*(dz**2)*Q1(lz-2,coxl2) + ks*(dx**2)*(Q1(lz-4,i) - 3.0d0*Q1(lz-2,i))&
        - 2.0d0*kt*(1.0d0 + 8.0d0*((q_0*dx)**2))*(dz**2)*Q1(lz-2,i) + 4.0d0*ks*(dx**2)*Q1(lz-1,i) - 2.0d0*ks*(dx**2)*Q1(lz,i)&
        + kt*(dz**2)*Q1(lz-2,coxr2) + 8.0d0*kt*q_0*(dx**2)*dz*(Q2(lz-3,i) - Q2(lz-1,i)) + 4.0d0*kt*q_0*(dx**2)*Q2(lz,i)&
        + (kt-ks)*dx*dz*Q3(lz-3,coxl1) + (ks-kt)*dx*dz*Q3(lz-1,coxl1) - ks*dx*dz*Q3(lz,coxl1)&
        + (ks-kt)*dx*dz*Q3(lz-3,coxr1) + (kt-ks)*dx*dz*Q3(lz-1,coxr1) + ks*dx*dz*Q3(lz,coxr1)&
        + 2.0d0*kt*(dz**2)*Q4(lz-2,coxl2) + (ks+kt)*(dx**2)*Q4(lz-4,i) - 3.0d0*(ks+kt)*(dx**2)*Q4(lz-2,i)&
        - 4.0d0*kt*(dz**2)*Q4(lz-2,i) - 32.0d0*kt*((q_0*dx*dz)**2)*Q4(lz-2,i) + 4.0d0*(ks+kt)*(dx**2)*Q4(lz-1,i)&
        - 2.0d0*(ks+kt)*(dx**2)*Q4(lz,i) + 2.0d0*kt*(dz**2)*Q4(lz-2,coxr2)&
        + 16.0d0*kt*q_0*dx*(dz**2)*(Q5(lz-2,coxr1)-Q5(lz-2,coxl1)) )

        V(5 + ((i-1)*5) + ((k-1)*lx*5)) = V(5 + ((i-1)*5) + ((k-1)*lx*5)) &
        - (1.0d0/(4.0d0*dx*dz))*( 8.0d0*kt*q_0*dx*(dz**2)*(Q1(lz-2, coxl1) - Q1(lz-2,coxr1))&
        + (ks-kt)*dx*dz*Q2(lz-3,coxl1) + (kt-ks)*dx*dz*Q2(lz-1,coxl1) + ks*dx*dz*Q2(lz,coxl1)&
        + (kt-ks)*dx*dz*Q2(lz-3,coxr1) + (ks-kt)*dx*dz*Q2(lz-1,coxr1) - ks*dx*dz*Q2(lz,coxr1)&
        + 8.0d0*kt*q_0*(dx**2)*dz*(Q3(lz-3,i) - Q3(lz-1,i)) + 4.0d0*kt*q_0*(dx**2)*dz*Q3(lz,i)&
        + 16.0d0*kt*q_0*dx*(dz**2)*(Q4(lz-2,coxl1)-Q4(lz-2,coxr1)) + 2.0d0*kt*(dz**2)*Q5(lz-2,coxl2)&
        + (ks+kt)*(dx**2)*Q5(lz-4,i) - 3.0d0*(ks+kt)*(dx**2)*Q5(lz-2,i) - 4.0d0*kt*(dz**2)*Q5(lz-2,i)&
        - 32.0d0*kt*((q_0*dx*dz)**2)*Q5(lz-2,i) + 4.0d0*(ks+kt)*(dx**2)*Q5(lz-1,i)&
        - 2.0d0*(ks+kt)*(dx**2)*Q5(lz,i) + 2.0d0*kt*(dz**2)*Q5(lz-2,coxr2) )

      else
        V(1 + ((i-1)*5) + ((k-1)*lx*5)) = V(1 + ((i-1)*5) + ((k-1)*lx*5))&
        - (1.0d0/(4.0d0*dx*dz))*( (ks+kt)*(dz**2)*Q1(k,coxl2) + (ks+kt)*(dx**2)*Q1(k-2,i) - 2.0d0*(ks+kt)*(dx**2)*Q1(k,i)&
        - 2.0d0*(ks+kt)*(dz**2)*Q1(k,i) - 32.0d0*kt*((q_0*dx*dz)**2)*Q1(k,i) + (ks+kt)*(dx**2)*Q1(k+2,i)&
        + (ks+kt)*(dz**2)*Q1(k,coxr2) + 8.0d0*kt*q_0*(dx**2)*dz*(Q2(k+1,i)-Q2(k-1,i))&
        + kt*(dz**2)*Q4(k,coxl2) + ks*(dx**2)*Q4(k-2,i) - 2.0d0*(dx**2)*ks*Q4(k,i) - 2.0d0*(dz**2)*Q4(k,i)&
        - 16.0d0*kt*((q_0*dx*dz)**2)*Q4(k,i) + ks*(dx**2)*Q4(k+2,i) + kt*(dz**2)*Q4(k,coxr2)&
        + 8.0d0*kt*q_0*dx*(dz**2)*(Q5(k,coxr1) - Q5(k,coxl1)) )

        V(2 + ((i-1)*5) + ((k-1)*lx*5)) = V(2 + ((i-1)*5) + ((k-1)*lx*5)) &
        - (1.0d0/(4.0d0*dx*dz))*( 8.0d0*kt*q_0*(dx**2)*dz*(Q1(k-1,i) - Q1(k+1,i)) + (ks+kt)*(dz**2)*Q2(k,coxl2)&
        + 2.0d0*kt*(dx**2)*(Q2(k-2,i) - 2.0d0*Q3(k,i)) - 2.0d0*(ks+kt)*(dz**2)*Q2(k,i) - 32.0d0*kt*((q_0*dx*dz)**2)*Q2(k,i)&
        + 2.0d0*kt*(dx**2)*Q2(k+2,i) + (ks+kt)*(dz**2)*Q2(k,coxr2) + 8.0d0*kt*q_0*dx*(dz**2)*(Q3(k,coxr1)-Q3(k,coxl1))&
        + 8.0d0*kt*(dx**2)*dz*(Q4(k+1,i) - Q4(k-1,i)) + (ks-kt)*dx*dz*(Q5(k-1,coxl1)) + (ks-ks)*dx*dz*Q5(k+1,coxl1)&
        + (ks-ks)*dx*dz*Q5(k-1,coxr1) + (ks-kt)*dx*dz*Q5(k+1,coxr1) )

        V(3 + ((i-1)*5) + ((k-1)*lx*5)) = V(3 + ((i-1)*5) + ((k-1)*lx*5)) &
        - (1.0d0/(4.0d0*dx*dz))*( 8.0d0*kt*q_0*dx*(dz**2)*(Q2(k,coxl1) - Q2(k,coxr1)) + (ks+kt)*(dz**2)*Q3(k,coxl2)&
        - 2.0d0*(ks+kt)*(dx**2)*Q3(k,i) - 2.0d0*(ks+kt)*(dz**2)*Q3(k,i) - 32.0d0*kt*((q_0*dx*dz)**2)*Q3(k,i)&
        + (ks+kt)*(dx**2)*Q3(k+2,i) + (ks+kt)*(dz**2)*Q3(k,coxr2) + (ks-kt)*dx*dz*Q4(k-1,coxl1)&
        + (ks-kt)*dx*dz*Q4(k+1,coxl1) + (ks-kt)*dx*dz*Q4(k-1,coxr1) + (kt-ks)*dx*dz*Q4(k+1,coxr1)&
        + 8.0d0*kt*q_0*(dx**2)*dz*(Q5(k+1,i) - Q5(k-1,i)) )

        V(4 + ((i-1)*5) + ((k-1)*lx*5)) = V(4 + ((i-1)*5) + ((k-1)*lx*5)) &
        - (1.0d0/(4.0d0*dx*dz))*( kt*(dz**2)*Q1(k,coxl2) + ks*(dx**2)*(Q1(k-2,i) - 2.0d0*Q1(k,i))&
        - 2.0d0*kt*(1.0d0+8.0d0*((q_0*dx)**2))*(dz**2)*Q1(k,i) + ks*(dx**2)*Q1(k+2,i) + kt*(dz**2)*Q1(k,coxr2)&
        + 8.0d0*kt*q_0*(dx**2)*dz*(Q2(k-1,i) - Q2(k+1,i)) + (kt-ks)*dx*dz*Q2(k-1,coxl1) + (ks-kt)*dx*dz*Q3(k+1,coxl1)&
        + (ks-kt)*dx*dz*Q3(k-1,coxr1) + (kt-ks)*dx*dz*Q3(k+1,coxr1) + 2.0d0*kt*(dz**2)*Q4(k,coxl2)&
        + (ks+kt)*(dx**2)*Q4(k-2,i) - 2.0d0*(ks+kt)*(dx**2)*Q4(k,i) - 4.0d0*kt*(dz**2)*Q4(k,i)&
        - 32.0d0*kt*((q_0*dx*dz)**2)*Q4(k,i) + (ks+kt)*(dx**2)*Q4(k+2,i) + 2.0d0*kt*(dz**2)*Q4(k,coxr2)&
        + 16.0d0*kt*q_0*dx*(dz**2)*(Q5(k,coxr1) - Q5(k,coxl1))  )

        V(5 + ((i-1)*5) + ((k-1)*lx*5)) = V(5 + ((i-1)*5) + ((k-1)*lx*5)) &
        + (1.0d0/4.0d0)*( 8.0d0*kt*q_0*dz*(Q1(k,coxr1) - Q1(k,coxl1)) + kt*(Q2(k-1,coxl1) - Q2(k+1,coxl1) - Q2(k-1,coxr1))&
        + ks*(-Q2(k-1,coxl1) + Q2(k+1,coxl1) + Q2(k-1,coxr1) - Q2(k+1,coxr1)) + kt*Q2(k+1,coxr1)&
        + 8.0d0*kt*q_0*dx*(Q3(k+1,i) - Q3(k-1,i)) + 16.0d0*kt*q_0*dz*(Q4(k,coxr1)-Q4(k,coxl1))&
        - 2.0d0*kt*(dz/dx)*(Q5(k,coxl2)) - (ks+kt)*(dx/dz)*Q5(k-2,i) + 2.0d0*(kt+ks)*(dx/dz)*Q5(k,i)&
        + 4.0d0*kt*(dz/dx)*Q5(k,i) + 32.0d0*kt*(q_0**2)*dx*dz*Q5(k,i) - (ks+kt)*(dx/dz)*Q5(k+2,i)&
        - 2.0d0*kt*(dz/dx)*Q5(k,coxr2)  )

      end if

      !BULK
      V(1 + ((i-1)*5) + ((k-1)*lx*5)) = V(1 + ((i-1)*5) + ((k-1)*lx*5)) &
      + weight(k,i)*(-a*(2*Q1(k,i) + Q4(k,i)) &
      + b*(-(Q2(k,i)**2) + 2*Q1(k,i)*Q4(k,i) + Q4(k,i)**2 + Q5(k,i)**2) &
      + 2.0D0*c*(2*Q1(k,i)+Q4(k,i))*(Q1(k,i)**2 + Q2(k,i)**2 + Q3(k,i)**2 + Q4(k,i)**2 + Q5(k,i)**2 + Q1(k,i)*Q4(k,i)))

      V(2 + ((i-1)*5) + ((k-1)*lx*5)) = V(2 + ((i-1)*5) + ((k-1)*lx*5)) &
      + weight(k,i)*(-2*a*Q2(k,i) &
      + 2*b*(-Q1(k,i)*Q2(k,i)-Q2(k,i)*Q4(k,i)-Q3(k,i)*Q5(k,i)) & !changed from -*- tp -*+
      + 4.0D0*c*Q2(k,i)*(Q1(k,i)**2 + Q2(k,i)**2 + Q3(k,i)**2 + Q4(k,i)**2 + Q5(k,i)**2 + Q1(k,i)*Q4(k,i)))

      V(3 + ((i-1)*5) + ((k-1)*lx*5)) = V(3 + ((i-1)*5) + ((k-1)*lx*5)) &
      + weight(k,i)*(-2*a*Q3(k,i) &
      + 2*b*(Q3(k,i)*Q4(k,i)-Q2(k,i)*Q5(k,i)) &
      + 4*c*Q3(k,i)*(Q1(k,i)**2 + Q2(k,i)**2 + Q3(k,i)**2 + Q4(k,i)**2 + Q5(k,i)**2 + Q1(k,i)*Q4(k,i)))

      V(4 + ((i-1)*5) + ((k-1)*lx*5)) = V(4 + ((i-1)*5) + ((k-1)*lx*5)) &
      + weight(k,i)*(-a*(Q1(k,i) + 2*Q4(k,i)) &
      + b*(-(Q2(k,i)**2)+2*Q1(k,i)*Q4(k,i) + Q1(k,i)**2 + Q3(k,i)**2) &
      + 2*c*(Q1(k,i) + 2*Q4(k,i))*(Q1(k,i)**2 + Q2(k,i)**2 + Q3(k,i)**2 + Q4(k,i)**2 + Q5(k,i)**2 + Q1(k,i)*Q4(k,i)))

      V(5 + ((i-1)*5) + ((k-1)*lx*5)) = V(5 + ((i-1)*5) + ((k-1)*lx*5)) &
      + weight(k,i)*(-2*a*Q5(k,i) &
      + 2*b*(Q1(k,i)*Q5(k,i)-Q2(k,i)*Q3(k,i)) &
      + 4*c*Q5(k,i)*(Q1(k,i)**2 + Q2(k,i)**2 + Q3(k,i)**2 + Q4(k,i)**2 + Q5(k,i)**2 + Q1(k,i)*Q4(k,i)))

    end do
  end do

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
