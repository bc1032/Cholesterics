!TOY_2D
subroutine set_defaults()
  use potential
  !#################!
  ! User Parameters !
  !#################!

  CM1 = .1

  lz = 100
  lx = 100
  ly = 1

  N = 5*(lx)*(lz)
  M = 6

  !###################!
  ! L-BFGS Parameters !
  !###################!

  max_iterations  = 1000000
  convergence_rms = 1d-6
  max_step_size   = 0.25d0
  H0init          = 1d1
  relative_energy_check = .true.
  dE_max          = 1d-2
  ! Print convergence information to GMIN_OUT
  debug = .true.

  !##########################!
  ! Basin Hopping Parameters !
  !##########################!

  ! Number of basin hopping steps to take
  max_runs          = 1
  max_similarity    = 1d-4
  temperature       = 1d10
  binary_io         = .false.
  ! Save only the lowest i minima. 0 Saves all
  minima_save_limit = 0
end subroutine set_defaults
