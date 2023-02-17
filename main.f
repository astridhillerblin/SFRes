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
	include 'Dependencies/resonance_plots.f'
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
! Starts program
!======================================================================!
	real :: start, finish
	real(8) :: qsq
	integer :: single,nmc
	common/singlenr/single
	common/nsamp/nmc
	call cpu_time(start)
	print '("Time = ",f10.1," seconds.")',start
	write(*,*) 'Getting started'
	nmc = 10000
	qsq = 1.3d0
	call res_params
	call readcouplings
	call init_random_seed2 ()
	write(*,*) 'Now writing single resonance files'
	call write_singres(qsq)
	call cpu_time(finish)
	print '("Time = ",f10.1," seconds.")',finish-start
!	write(*,*) 'Now writing truncated moments'
!	call write_trunc_mom
!	call cpu_time(finish)
!	print '("Time = ",f10.1," seconds.")',finish-start
!	write(*,*) 'Generating MC for uncertainty propagation'
!	call write_ressamp(qsq)
!	call cpu_time(finish)
!	print '("Time = ",f10.1," seconds.")',finish-start
!	write(*,*) 'Creating observables with uncertainty propagation'
!	call mc_meanstd
!	call cpu_time(finish)
!	print '("Time = ",f10.1," seconds.")',finish-start
	write(*,*) 'Finished all tasks'
!======================================================================!
! End of program
!======================================================================!
	return
	end program resonances
