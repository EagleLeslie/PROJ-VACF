PROGRAM main

  USE element_table_mod, only: get_element, atomic_element

  IMPLICIT NONE

  INTEGER :: i,j,k,l,n
  INTEGER :: ntimesteps, nions, totatoms
  INTEGER, ALLOCATABLE, DIMENSION(:) :: natom
  REAL*8 :: start, finish
  REAL*8, ALLOCATABLE, DIMENSION(:) :: volume, wmass
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: acell, xyz, vacf, deriv
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: poly
  CHARACTER(len=3), ALLOCATABLE, DIMENSION(:) :: atype
  COMPLEX*8, ALLOCATABLE, DIMENSION(:,:) :: proj_vacf
  COMPLEX*8, ALLOCATABLE, DIMENSION(:,:,:) :: wvel
  LOGICAL :: sim_type
  TYPE(atomic_element) :: element

  ! Time how long program runs
  CALL CPU_TIME(start)

  ! Find parameters: number of timesteps, number of ions in system, and total atoms
  CALL find_params(ntimesteps, nions, totatoms)

  ! Allocate arrays according to parameters
  ALLOCATE(natom(nions))
  ALLOCATE(atype(nions))
  ALLOCATE(wmass(nions))
  ALLOCATE(volume(ntimesteps))
  ALLOCATE(acell(ntimesteps,3,3))
  ALLOCATE(xyz(ntimesteps,totatoms,3))
  ALLOCATE(vacf(nions,ntimesteps,3))
  ALLOCATE(deriv(ntimesteps,totatoms,3))
  ALLOCATE(wvel(ntimesteps,totatoms,3))
  ALLOCATE(proj_vacf(ntimesteps,3))
  ALLOCATE(poly(ntimesteps,totatoms,3,3))

  ! Read in XDATCAR
  CALL read_xdat(ntimesteps,nions,totatoms,acell,volume,natom,xyz,atype, sim_type)

  ! True == NPT or NPH simulation
  ! False == NVT simulation
  PRINT*, sim_type

  ! Get atomic masses of the system
  DO i=1,nions
    element = get_element(atype(i))
    wmass(i) = element%mass
    WRITE(*,*) "Atomic Mass: ", atype(i), " = ", wmass(i)
  END DO 

  ! Calculate velocity autocorrelation function
  CALL calc_vacf(ntimesteps, nions, natom, xyz, atype, acell, poly, vacf, deriv)

  ! Calculate projected velocity autocorrleation function
  CALL calc_qvel(ntimesteps, nions, natom, xyz, acell, deriv, wmass, wvel, proj_vacf)

  ! Deallocate arrays
  DEALLOCATE(natom, atype, wmass, volume, acell, xyz, vacf, deriv, wvel, proj_vacf, poly)

  CALL CPU_TIME(finish)
  PRINT*, "Time = ", finish-start, "seconds."

END PROGRAM main