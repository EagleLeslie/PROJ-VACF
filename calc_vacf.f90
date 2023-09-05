SUBROUTINE calc_vacf(ntimesteps, nions, natom, xyz, atype, acell, poly, vacf, deriv)

    IMPLICIT NONE

    INTEGER :: i,j,k,l,m,tau
    INTEGER, ALLOCATABLE :: totatoms
    INTEGER, INTENT(IN) :: ntimesteps, nions
    INTEGER, DIMENSION(:) :: ntype(nions+1)
    INTEGER, DIMENSION(:), INTENT(IN) :: natom(nions)
    INTEGER, PARAMETER :: mint=3
    REAL*8 :: vx, vy, vz, vel, x, a, b, c, poop1, poop2, poop3, poop4
    REAL*8, DIMENSION(:) :: cofx(3), cofy(3), cofz(3), time(ntimesteps)
    REAL*8, DIMENSION(:,:,:), INTENT(IN) :: xyz(ntimesteps,SUM(natom),3), acell(ntimesteps,3,3)
    REAL*8, DIMENSION(:,:,:), INTENT(OUT) :: vacf(nions,ntimesteps,3), deriv(ntimesteps,SUM(natom),3)
    REAL*8, DIMENSION(:,:,:,:), INTENT(OUT) :: poly(ntimesteps,SUM(natom),3,3)
    CHARACTER(len=3), DIMENSION(:) :: aname(nions+1)
    CHARACTER(len=3), DIMENSION(:), INTENT(IN) :: atype(nions)

    WRITE(*,*) 'BEGINNING VACF'
    WRITE(*,*)

    ! Initialize data
    ntype(:) = 0
    cofx(:) = 0
    cofy(:) = 0
    cofz(:) = 0

    ! Make time array into a real number
    DO i = 1,ntimesteps
        time(i) = REAL(i)
    END DO

    DO i = 1,nions
        ntype(i+1) = ntype(i+1) + ntype(i) + natom(i)
        aname(i+1) = atype(i)
        WRITE(*,*) ntype(i+1), aname(i+1)
    END DO

    ! Assumes orthogonal cell
    a = acell(1,1,1)
    b = acell(1,2,2)
    c = acell(1,3,3)
    
    totatoms = SUM(natom)
    print*, 'TOTAL ATOMS: ',totatoms
    print*, 'TOTAL TIME: ', ntimesteps

    DO i = 1,ntimesteps-1
        DO j = 1,totatoms
        ! Polynomial interpolation and calculate polynomial interpolation coefficients
        ! polint(xa,ya,n,x,y,dy)
        ! REAL :: dy,x,y,xa(n),ya(n)
        ! INTEGER :: n, NMAX
        ! Given arrays xa and ya, each of length n, and given a value of x, polint
        ! returns a value y, and an error estiamte dy
        !CALL polint(time(i), xyz(i,j,1), mint, x, y, dy)
        CALL polcof(time(i), xyz(i,j,1), mint, cofx)

        !CALL polint(time(i), xyz(i,j,2), mint, x, y, dy)
        CALL polcof(time(i), xyz(i,j,2), mint, cofy(:))

        !CALL polint(time(i), xyz(i,j,3), mint, x, y, dy)
        CALL polcof(time(i), xyz(i,j,3), mint, cofz(:))

        poly(i,j,1,:) = cofx
        poly(i,j,2,:) = cofy
        poly(i,j,3,:) = cofz

        ! ! Calculate velocities; derivative of polynomial interpolation output
        deriv(i,j,1) = poly(i,j,1,2) + (2.*time(i)*poly(i,j,1,3))
        deriv(i,j,2) = poly(i,j,2,2) + (2.*time(i)*poly(i,j,2,3))
        deriv(i,j,3) = poly(i,j,3,2) + (2.*time(i)*poly(i,j,3,3))

        END DO
    END DO

    OPEN(UNIT=40,FILE='vacfout')
    OPEN(UNIT=50,FILE='unnorm-vacf')

    m = 1
    DO l = 1,nions
        WRITE(*,*) "CALCULATING: ", atype(l), atype(l)
        DO tau = 1, ntimesteps-1
        ! i is the position's "origin"
        DO i = 1,ntimesteps-tau
            vel = 0.0
            ! calculate v(t) * v(t + tau) for each individual atom
            DO j = m,ntype(l+1)
            vx = deriv(i+tau-1,j,1) * deriv(i,j,1)
            vy = deriv(i+tau-1,j,2) * deriv(i,j,2)
            vz = deriv(i+tau-1,j,3) * deriv(i,j,3)
            ! sum displacement for ALL atom type 1
            vel = vel + (((a**2.)*vx) + ((b**2.)*vy) + ((c**2.)*vz))
            END DO
            ! average over all atom type l
            vacf(l,tau,1) = vacf(l,tau,1) + (vel/ntype(l+1))
        END DO
        ! VAC average over all time intervals
        vacf(l,tau,2) = vacf(l,tau,1)/(ntimesteps-tau) ! un-normalized vacf
        vacf(l,tau,3) = vacf(l,tau,2)/vacf(l,1,2) ! normalized vacf
        ! TIMESTEP PRINT OUT; FROM FRED !
        IF (mod(i,100) .eq. 0) write(6,'(a25,i7,a10)',ADVANCE='NO') 'COMPUTING VAC ',tau,'         '//CHAR(13)
        END DO
        m = 1 + ntype(l+1)
        WRITE(*,*)
    END DO

    ! Write out VACF to vacfout
    DO tau = 1, ntimesteps-1
        WRITE(40,*) tau, vacf(1,tau,3), vacf(2,tau,3), vacf(3,tau,3) ! normalized vacf
        WRITE(50,*) tau, vacf(1,tau,2), vacf(2,tau,2), vacf(2,tau,3) ! un-normalized vacf
    END DO

    WRITE(*,*)
    WRITE(*,*) "FINISHED COMPUTING VACF"
    WRITE(*,*)

    CLOSE(40)
    CLOSE(50)

END SUBROUTINE calc_vacf