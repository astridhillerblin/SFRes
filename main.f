!======================================================================!
! Run with: gfortran SF.f
!======================================================================!
!======================================================================!
! Includes all external fortran files
!======================================================================!
        include 'Dependencies/useful_constants.f'
        include 'Dependencies/res_arrays.f'
	include 'Dependencies/model.f'
        include 'Dependencies/intx.f'
	include 'Dependencies/resonance_files.f'
	include 'Dependencies/coupling_mc.f'
	include 'Dependencies/computer.f'
!======================================================================!
! MAIN
!======================================================================!
	program resonances
	use phys_consts
	use res_decl
	implicit none
!======================================================================!
! Start program
!======================================================================!
	real :: start, finish
	real(8) :: qsq,ebeam
	integer :: single,nmc
	character*1 :: trunc,bands
	common/singlenr/single
	common/nsamp/nmc
!======================================================================!
! Read input from terminal
!======================================================================!
	write(*,*) 'Please enter the Q^2 value in GeV^2.'
	read(*,*) qsq
	write(*,*) 'Please enter the beam energy in GeV.'
	read(*,*) ebeam
	write(*,*) 'Calculate truncated moments? y/n'
	write(*,*) '(Slightly long computation time.)'
	read(*,*) trunc
	write(*,*) 'Calculate uncertainty bands? y/n'
	write(*,*) '(Even longer computation time.)'
	read(*,*) bands
!======================================================================!
! Call initializing subroutines
!======================================================================!
	call res_params
	call cpu_time(start)
	print '("Time = ",f10.1," seconds.")',start
	write(*,*) 'Getting started'
!======================================================================!
! Initializes bootstrap
!======================================================================!
	nmc = 10000
	call readcouplings
	call init_random_seed2 ()
	write(*,*)'Writing single resonance and (in)coherent sum files'
	call write_singres(qsq,ebeam)
	call cpu_time(finish)
	print '("Time = ",f10.1," seconds.")',finish-start
	if (trunc=='y') then
	write(*,*)'Now writing truncated moment files'
	call write_trunc_mom
	call cpu_time(finish)
	print '("Time = ",f10.1," seconds.")',finish-start
	endif
	if (bands=='y') then
	write(*,*) 'Generating MC for uncertainty propagation'
	call write_ressamp(qsq,ebeam)
	call cpu_time(finish)
	print '("Time = ",f10.1," seconds.")',finish-start
	write(*,*) 'Creating uncertainty band files'
	call mc_meanstd
	call cpu_time(finish)
	print '("Time = ",f10.1," seconds.")',finish-start
	endif
	write(*,*) 'Finished all tasks'
!======================================================================!
! End of program
!======================================================================!
	return
	end program resonances
