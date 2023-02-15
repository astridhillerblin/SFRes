!======================================================================!
! MODULE THAT DECLARES PHYSICAL CONSTANTS
!======================================================================!
	module phys_consts
	implicit none
	real(8), parameter :: mn = 0.93827
	real(8), parameter :: mpi = 0.135
	real(8), parameter :: meta = 0.547853
	real(8), parameter :: alpha = 1./137.0366
	end module phys_consts
!======================================================================!
! MODULE THAT DECLARES MATHEMATICAL CONSTANTS
!======================================================================!
	module math_consts
	implicit none
	real(8), parameter :: pi = 4.d0*datan(1.d0)
	real(8), parameter :: zero = 0.d0
	real(8), parameter :: ex = exp(1.d0)
	end module math_consts 