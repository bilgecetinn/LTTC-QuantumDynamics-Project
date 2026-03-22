!-----------------!
 module parameters
!-----------------!

 real(8), parameter :: length=5.12d0   ! length of the box (in Bohr)
 real(8), parameter :: mass = 1822.88839d0 != 1 atomic mass unit (=1 g/mol)
 real(8), parameter :: pi=3.141592653589793d0
 real(8), parameter :: au2kcalmol=627.509d0
 real(8), parameter :: fs2au=41.341373336561d0
 real(8)            :: angfreq, barrier
 integer            :: nmax
 character(10)      :: potentialtype
 
 end module parameters
!---------------------!

!-----------------!
 program propagate
!-----------------!
 use parameters
 implicit none

 integer                        :: npoints,ntime,snapshot,i,n,v
 real(8)                        :: alpha,dt,t,dx,x0,x
 real(8), allocatable           :: pot(:),kin(:),psisquare(:)
 complex(8), allocatable        :: psi(:),psi0(:),exppot(:),expkin(:)
 !real(8),allocatable            :: c(:)
 real(8)                        :: c(0:10)

!allocate(c(npoints))
 c = 0.0d0

 open(unit=10,file='wavepacket')
   read(10,*) npoints               !Number of lattice points
   read(10,*) x0                    !Initial position
   read(10,*) alpha                 !Governs the initial width of the wave packet
   read(10,*) dt                    !Propagation time step
   read(10,*) ntime                 !Number of propagation steps
   read(10,*) snapshot              !snapshot frequency 
 close(10)

 open(unit=11,file='potential')
   read(11,*) potentialtype         !harmonic or double well potential
   read(11,*) angfreq               !Angular frequency for harmonic potential
   read(11,*) barrier               !Height of barrier in double well potential (in kcal/mol)
 close(11)

 open(12, file='coefficients')
 do n = 0, 4
    read(12, *) c(n)
 end do
 close(12)

 dt=dt*fs2au                        !convert femtoseconds to atomic units
 angfreq=angfreq/fs2au              !convert femtoseconds to atomic units

 allocate(psi(npoints),psi0(npoints))
 allocate(pot(npoints),exppot(npoints))
 allocate(kin(npoints),expkin(npoints))
 allocate(psisquare(npoints))

 dx=length/dble(npoints)
 nmax = 4     !# of hermite polynomials

 call initpsi(npoints,dx,alpha,x0,psi0,angfreq,nmax,c) 	!Obtain initial wavepacket psi0
 call fourier(0,npoints,psi0)                          	!Initialize the FFT
 call operators(npoints,dx,dt,pot,kin,exppot,expkin)   	!Calculate the kinetic and potential operators

 psi=psi0                                            	!Set the wavepacket psi at t=0 equal psi0
 do i=0,ntime                                        	!Start propagation
    t=i*dt
    if (i>0) then
       psi=psi*exppot                                	!Multiply psi with exp(-i*dt*potential)
       call fourier(1,npoints,psi)                  	!Forward FFT to momentum space
       psi=psi*expkin                               	!Multiply psi with the exp(-i*dt*kinetic operator)
       call fourier(-1,npoints,psi)                  	!Backward FFT to position space
   
    endif
    if (mod(i,snapshot)==0) then                     !Take a snapshot if the remainder of i/snapshot equals 0
       call initgraph(i/snapshot,t)                  !Initialize graph
       psisquare=(abs(psi))**2
       call graphpot(dx,npoints)                     !Plot the potential
       call graphpsi(dx,npoints,psisquare)           !Plot |psi|^2
    endif
 end do                                              !End propagation

 deallocate(psi,psi0)
 deallocate(pot,exppot)
 deallocate(kin,expkin)
 deallocate(psisquare)
 !deallocate(c)

 end program propagate
!---------------------!

!---------------------
subroutine initpsi(npoints, dx, alpha, x0, psi0, angfreq, nmax, c)
!---------------------!
 !-----------------------------------------------------------------------------
 ! This subroutine initializes the wavefunction for a quantum harmonic oscillator
 ! using a superposition of eigenstates. The wavefunction is computed on a 
 ! discrete grid and normalized.
 !
 ! INPUTS:
 !   npoints   - Number of grid points
 !   dx        - Grid spacing
 !   alpha     - Parameter controlling width of wavefunction (not used here)
 !   x0        - Initial position (not used here)
 !   angfreq   - Angular frequency of the harmonic potential
 !   nmax      - Maximum quantum number considered
 !   c(0:10)   - Coefficients for constructing the wave packet
 !
 ! OUTPUTS:
 !   psi0(npoints) - Initialized wavefunction (complex)
 !
 ! EXTERNAL FUNCTIONS:
 !   factorial - Computes factorial of a given integer
 !   hermite_poly - Computes Hermite polynomials needed for wavefunction
 !-----------------------------------------------------------------------------

 ! Use parameters and precision settings
 implicit none

 ! Physical constants
 real(8), parameter                    :: mass = 1822.88839d0        ! Atomic mass unit in a.u.
 real(8), parameter                    :: pi = 3.141592653589793d0  

 ! Input parameters
 integer, intent(in)                   :: npoints, nmax
 real(8), intent(in)                   :: alpha, x0, dx
 complex(8), intent(out)               :: psi0(npoints)
 real(8), intent(in)                   :: angfreq                ! Angular frequency of harmonic oscillator
 real(8), external                     :: factorial             ! Declare factorial function
 real(8)                               :: c(0:10)               ! Coefficients for superposition of states

 ! Local variables
 integer :: i, j, n
 real(8) :: x, norm, prefactor, hbar, y, fact
 real(8) :: H(0:nmax), psi_real(npoints)  ! Arrays for Hermite polynomials and real part of psi

 ! Set reduced Planck's constant to 1 (atomic units)
 hbar = 1.0d0

 ! Initialize the wavefunction array to zero
 psi0 = (0.0d0, 0.0d0)

 ! Loop over grid points to compute the wavefunction
 do i = -npoints/2 + 1, npoints/2
     x = dble(i) * dx  ! Convert grid index to real space coordinate

     ! Ensure positive index mapping for Fortran arrays
     if (i > 0) then
         j = i
     else
         j = i + npoints
     endif

     ! Parameter to be passed to hermite polynomial 
     y = x * dsqrt(mass * angfreq / hbar)

     ! Compute Hermite polynomials up to H_nmax(y)
     call hermite_poly(y, H, nmax)

     ! Initialize the real part of the wavefunction
     psi_real(j) = 0.0d0

     ! Construct wave packet as a superposition of eigenstates (up to n=4)
     do n = 0, 4
         ! Retrieve n!:
         fact = factorial(n)

         ! Prefactor for the Hermite polynomial contribution (assumes hbar=1)
         prefactor = 1.0d0 / (dsqrt(dble(2**n) * fact)) * (mass * angfreq / pi)**0.25d0
         
         ! Compute wavefunction component using Hermite polynomial H(n)
         psi_real(j) = psi_real(j) + c(n) * prefactor * H(n) * dexp(-0.5d0 * mass * angfreq * x**2 / hbar)
     end do

     ! Convert real wavefunction to complex form (imaginary part is zero)
     psi0(j) = dcmplx(psi_real(j), 0.0d0)
 end do

 ! Normalize the wavefunction (task4)
  norm = 0.0d0
  do i = 1, npoints
      norm = norm + abs(psi0(i))**2 * dx
  end do
  norm = 1.0d0 / dsqrt(norm)
  do i = 1, npoints
      psi0(i) = psi0(i) * norm
  end do


end subroutine initpsi
!----------------------!

!------------------------------------------------------!
subroutine operators(npoints,dx,dt,pot,kin,exppot,expkin)
!------------------------------------------------------!
 ! Computes potential and kinetic energy operators for a quantum system.
 ! Also evaluates their exponential forms for time evolution.

 use parameters
 implicit none

 ! Input/Output variables
 integer :: i, j, npoints
 real(8) :: x, p, dt, dx, dp, pot(npoints), kin(npoints)
 complex(8) :: exppot(npoints), expkin(npoints)

 ! Define momentum step size using periodic boundary conditions
 dp = 2.d0 * pi / length

 ! Loop over grid points to define operators
 do i = -npoints/2 + 1, npoints/2
     x = dble(i) * dx  ! Position coordinate
     p = dble(i - 1) * dp  ! Discretized momentum

     ! Ensure positive index mapping for Fortran arrays
     if (i > 0) then
         j = i
     else
         j = i + npoints
     endif

     ! Define potential energy depending on chosen type
     if (potentialtype == 'harmonic') then
         pot(j) = 0.5d0 * mass * angfreq**2 * x**2  ! Harmonic oscillator potential
     elseif (potentialtype == 'doublewell') then
         pot(j) = barrier * (16.d0 * x**4 - 8.d0 * x**2 + 1.d0) / au2kcalmol  ! Double-well potential
     endif

     ! Define kinetic energy operator (p²/2m)
     kin(j) = 0.5d0 * p**2 / mass

     ! Compute exponentials for imaginary-time evolution
     exppot(j) = exp(-dt * (0,1) * pot(j))  ! exp(-iV dt)
     expkin(j) = exp(-dt * (0,1) * kin(j))  ! exp(-iT dt)
    end do

end subroutine operators
!------------------------------------------------------!

!-----------------------------------!
subroutine fourier(dir, npoints, psi)
!-----------------------------------!
 ! Performs Fast Fourier Transform (FFT) using external FFTPACK routines.
 ! Supports forward and inverse Fourier transforms.

 implicit none

 ! Input/Output variables
 integer :: i, npoints, dir
 real(8) :: nr
 complex(8) :: psi(npoints)
 real(8), allocatable, save :: wsave(:)  ! FFT workspace array

 if (dir == 1) then
     ! Forward FFT: Transform to momentum space
     call dcfftf(npoints, psi, wsave)
     nr = 1.d0 / dble(npoints)  ! Normalize after FFT
     do i = 1, npoints
         psi(i) = psi(i) * nr
     end do
 elseif (dir == -1) then
     ! Inverse FFT: Transform back to real space
     call dcfftb(npoints, psi, wsave)
 elseif (dir == 0) then
     ! Initialize FFT workspace
     if (allocated(wsave)) deallocate(wsave)
     allocate(wsave(4 * npoints + 20))
     call dcffti(npoints, wsave)
 endif

end subroutine fourier
!-----------------------------------!


!----------------------!
subroutine hermite_poly(y, H, nmax)
!----------------------!
 ! Computes Hermite polynomials H_n(y) up to order nmax using recursion.
 
 implicit none
 integer, intent(in)  :: nmax      ! Maximum order of Hermite polynomials
 real(8), intent(in)  :: y         ! Input variable y (scaled coordinate)
 real(8), intent(out) :: H(0:nmax) ! Array storing computed polynomials
 integer              :: n         ! Loop variable

 ! Base cases
 H(0) = 1.0d0
 if (nmax >= 1) then
     H(1) = 2.0d0 * y
 endif

 ! Loop from n=2 to nmax - recursion relation
 do n = 2, nmax
     H(n) = 2.0d0 * y * H(n-1) - 2.0d0 * (n-1) * H(n-2)
 end do

end subroutine hermite_poly
!----------------------!

!----------------------!
function factorial(n) result(fact) 
!----------------------!
 ! Computes factorial of an integer n using iterative multiplication
 
 implicit none
 integer, intent(in) :: n
 integer             :: i
 real(8)             :: fact   ! Using real(8) to maintain consistency

 ! Initialize factorial value
 fact = 1.0d0  
 if (n > 0) then
     do i = 1, n
         fact = fact * i
     end do
 endif

end function factorial
!----------------------!

