SUBROUTINE calc_qvel(ntimesteps, nions, natom, xyz, acell, deriv, wmass, wvel, proj_vacf)

! Subroutine that calculates the q projected velocity autocorrelation function
! q is the Brillion zone wave vector of the primitive cell sampled in the MD simulation
! q = (2 pi/ (a,b,c)) (h,k,l)

! deriv = velocities calculated from calc_vacf subroutine
! output: 
!   wvel = weighted velocity
!   proj_vacf = projected velocity autocorrelation 

IMPLICIT NONE

INTEGER :: i,j,n,h,k,l,tau
INTEGER, INTENT(IN) :: ntimesteps, nions
INTEGER, DIMENSION(:), INTENT(IN) ::  natom(nions)
REAL*8 :: qrx, qry, qrz, a, b, c
REAL*8, DIMENSION(:), INTENT(IN) :: wmass(nions)
REAL*8, DIMENSION(:,:,:), INTENT(IN) :: deriv(ntimesteps,SUM(natom),3), xyz(ntimesteps,SUM(natom),3), acell(ntimesteps,3,3)
REAL*8, PARAMETER :: pi = 4.0*ATAN(1.0)
COMPLEX*8, DIMENSION(:,:,:), INTENT(OUT) :: wvel(ntimesteps,SUM(natom),3)
COMPLEX*8 :: vel, vx, vy, vz
COMPLEX*8, DIMENSION(:,:), INTENT(OUT) :: proj_vacf(ntimesteps,3)
COMPLEX*8, PARAMETER :: imag = (0,1) ! sqrt(-1)


WRITE(*,*) "BEGINNING MODE-PROJECTED VELOCITY AUTOCORRELATION"

! Initialize (h,k,l) indices here
h = 0
k = 1
l = 0

! Initialize a, b, and c vectors
! a in the q vector will vary depending on the direction
! i.e. q = 2pi/a (1,0,0) or q = 2pi/c (0,0,1) 
 a = acell(1,1,1)/3 ! primitive vector for 3x3x3 supercell
 b = acell(1,2,2)/3 ! primitive vector for 3x3x3 supercell
 c = acell(1,3,3)/3 ! primitive vector for 3x3x3 supercell


! Calculate weighted velocities and exp(iq dot R)
DO n = 1,nions
    ! WRITE(*,*) wmass(n)
    DO i = 1, ntimesteps-1
        DO j = 1, natom(n)
        ! Calculate q dot R
        ! Note: fixing acell here is limiting this calculation to an orthorhombic or cubic system
        qrx = (2*pi/b) * (h*xyz(i,j,1)*a) 
        qry = (2*pi/b) * (k*xyz(i,j,2)*b)
        qrz = (2*pi/b) * (l*xyz(i,j,3)*c)

        wvel(i,j,1) = deriv(i,j,1) * SQRT(wmass(n)) * EXP(imag * qrx)
        wvel(i,j,2) = deriv(i,j,2) * SQRT(wmass(n)) * EXP(imag * qry)
        wvel(i,j,3) = deriv(i,j,3) * SQRT(wmass(n)) * EXP(imag * qrz)
        ! WRITE(*,*) (wvel(i,j,:))

        END DO
    END DO
END DO

! Calculated weighted velocity autocorrelation function
DO tau = 1,ntimesteps-1
    vel = 0
    DO i = 1,ntimesteps-tau
        DO j = 1,SUM(natom)
            vx = wvel(i+tau-1,j,1) * wvel(i,j,1)
            vy = wvel(i+tau-1,j,2) * wvel(i,j,2)
            vz = wvel(i+tau-1,j,3) * wvel(i,j,3)

            vel = vel + (((a**2.)*vx) + ((b**2.)*vy) + ((c**2.)*vz))
        END DO
        proj_vacf(tau,1) = proj_vacf(tau,1) + vel/SUM(natom)
    END DO
    proj_vacf(tau,2) = proj_vacf(tau,1)/(ntimesteps-tau)
    proj_vacf(tau,3) = proj_vacf(tau,2)/proj_vacf(1,2)
    IF (mod(i,100) .eq. 0) write(6,'(a25,i7,a10)',ADVANCE='NO') 'COMPUTING MODE-PROJ VAC ',tau,'         '//CHAR(13)
    WRITE(21,*) tau, REAL(proj_vacf(tau,2)), AIMAG(proj_vacf(tau,2))
    WRITE(22,*) tau, REAL(proj_vacf(tau,3)), AIMAG(proj_vacf(tau,3))
END DO

WRITE(*,*)
WRITE(*,*) "FINISHED COMPUTING MODE-PROJ VAC"

! ! Calculate the fourier transform of real part
! CALL four1(REAL(proj_vacf(:,3)),ntimesteps-1,1)

END SUBROUTINE calc_qvel