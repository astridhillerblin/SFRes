!======================================================================!
! Write plot files
!======================================================================!
	subroutine write_ressamp(qsq)
	use phys_consts
	implicit none
	character*60, dimension(2) :: res_file
	real(8) qsq,wsq
	real(8), dimension (19) :: flag
	integer iloop,jloop,floop
	integer iwrt
	integer isup, samp
	integer nmc
	real(8) f1tot,f2tot,fltot,g1tot,g2tot
	real(8) h12tot,h32tot
	real(8) a1tot,a2tot
	real(8) interf
	common/sf/f1tot,f2tot,fltot,g1tot,g2tot
	common/hel/h12tot,h32tot
	common/asym/a1tot,a2tot
	common/sing/flag
	common/isupov/isup
	common/nsamp/nmc
	common/swinterf/interf
	res_file(1)='Output/res_samp.dat'
	res_file(2)='Output/resinterf_samp.dat'
	isup = 0
	do floop=1,19
	flag(floop)=1.d0
	enddo
	iwrt = 101
	do iloop=1,2
	interf=(iloop-1)*1.d0
	open (iwrt, file = res_file(iloop))
	do samp = 1,nmc
	do jloop = 0,112
	wsq = (mn+mpi+0.01d0*jloop)**2
	call res_calc(wsq,qsq)
	write (iwrt,"(F10.2,9(F10.4))")
     >		dsqrt(wsq),f1tot,f2tot,fltot,g1tot,g2tot
     >		,h12tot,h32tot,a1tot,a2tot
	enddo
	enddo
	close (iwrt)
	enddo
	end subroutine write_ressamp
!======================================================================!
! Write single resonances
!======================================================================!
	subroutine write_singres(qsq)
	use phys_consts
	implicit none
	character*60, dimension (21) :: sing_file
	real(8), dimension (19) :: flag
	character*1 chint1
	character*2 chint2
	real(8) wsq,qsq
	integer iloop,jloop,floop,rloop
	integer iwrt
	integer isup
	real(8) f1tot,f2tot,fltot,g1tot,g2tot
	real(8) h12tot,h32tot
	real(8) a1tot,a2tot
	real(8) interf
	common/sf/f1tot,f2tot,fltot,g1tot,g2tot
	common/hel/h12tot,h32tot
	common/asym/a1tot,a2tot
	common/sing/flag
	common/isupov/isup
	common/swinterf/interf
	isup = 1
	iwrt = 101
	do rloop = 1,9
		write(chint1,'(i1)') rloop
		sing_file(rloop) = 
     >		'Output/sing_res_'//chint1//
     >		'.dat'
	enddo
	do rloop = 10,19
		write(chint2,'(i2)') rloop
		sing_file(rloop) = 
     >		'Output/sing_res_'//chint2//
     >		'.dat'
	enddo
	sing_file(20) = 'Output/sum_res.dat'
	sing_file(21) = 'Output/sum_res_interf.dat'
!Now the single resonance contributions are written for each resonance.
	do rloop=1,19
	open (iwrt, file = sing_file(rloop))
	do floop = 1,19
	flag(floop) = 0.d0	
	enddo
	flag(rloop) = 1.d0
	do jloop = 0,112
	wsq = (mn+mpi+0.01d0*jloop)**2
	call res_calc(wsq,qsq)
	write (iwrt,"(F10.2,9(F10.4))")
     >		dsqrt(wsq),f1tot,f2tot,fltot,g1tot,g2tot
     >		,h12tot,h32tot,a1tot,a2tot
	enddo
	close (iwrt)
	enddo
!Now the total resonance contributions are written for each resonance.
	do floop = 1,19
	flag(floop) = 1.d0	
	enddo
!First without, then with interferences.
	do iloop=1,2
	open (iwrt, file = sing_file(19+iloop))
	interf = (iloop-1)*1.d0
	do jloop = 0,112
	wsq = (mn+mpi+0.01d0*jloop)**2
	call res_calc(wsq,qsq)
	write (iwrt,"(F10.2,9(F10.4))")
     >		dsqrt(wsq),f1tot,f2tot,fltot,g1tot,g2tot
     >		,h12tot,h32tot,a1tot,a2tot
	enddo
	close (iwrt)
	enddo
	end subroutine write_singres
!======================================================================!
! Write truncated moments
!======================================================================!
	subroutine write_trunc_mom
	use phys_consts
	implicit none
	real(8),dimension(20000)::xrd
	real(8),dimension(20000)::integrand1,integrand2,integrand3,
     >		integrand4,integrand5
	real(8), dimension (19) :: flag
	real(8) wsq,qsq,wsqrmax,numax,x0,nu
	real(8) wsqrmax2,numax2,x02,wsqrmax3,numax3,x03,
     >		wsqrthr,nuthr,xthr
	integer floop,qloop
	integer iwrt1,iwrt2,iwrt3,iwrt4
	integer isup
	integer ninteg,n1p,i1
	integer, parameter :: ndat_f2clinterHC = 10794
	real(8) f1tot,f2tot,fltot,g1tot,g2tot
	real(8) f1trunc,f2trunc,fltrunc,g1trunc,g2trunc
	real(8) f1trunc2,f2trunc2,fltrunc2,g1trunc2,g2trunc2
	real(8) f1trunc3,f2trunc3,fltrunc3,g1trunc3,g2trunc3
	real(8) f1truncall,f2truncall,fltruncall,g1truncall,g2truncall
	real(8) interf
	common/swinterf/interf
	common/sf/f1tot,f2tot,fltot,g1tot,g2tot
	common/sing/flag
	common/isupov/isup
	interf=1.d0
	isup = 1
	ninteg = 10
	iwrt1 = 101
	iwrt2 = 102
	iwrt3 = 103
	iwrt4 = 104
	do floop = 1,19
	flag(floop) = 1.d0	
	enddo
	open (iwrt1, file = './Output/Trunc_Reg1.dat')
	open (iwrt2, file = './Output/Trunc_Reg2.dat')
	open (iwrt3, file = './Output/Trunc_Reg3.dat')
	open (iwrt4, file = './Output/Trunc_RegTot.dat')
	do qloop = 10,35
	qsq = qloop/10.d0
	wsqrthr=1.25d0
	wsqrmax=1.9d0
	wsqrmax2=2.5d0
	wsqrmax3=3.1d0
	nuthr=(wsqrthr-mn**2+qsq)/(2.d0*mn)
	numax=(wsqrmax-mn**2+qsq)/(2.d0*mn)
	numax2=(wsqrmax2-mn**2+qsq)/(2.d0*mn)
	numax3=(wsqrmax3-mn**2+qsq)/(2.d0*mn)
	xthr = qsq/(2.d0*mn*nuthr)
	x0 = qsq/(2.d0*mn*numax)
	x02 = qsq/(2.d0*mn*numax2)
	x03 = qsq/(2.d0*mn*numax3)
	call dsg20r(x0,xthr,ninteg,xrd,n1p)
	do i1=1,n1p
                nu=qsq/(2.d0*mn*xrd(i1))
                wsq=mn**2-qsq+2.d0*mn*nu
		call res_calc(wsq,qsq)
		integrand1(i1)=f1tot
		integrand2(i1)=f2tot
		integrand3(i1)=fltot
		integrand4(i1)=g1tot
		integrand5(i1)=g2tot
	enddo
	call DRG20r(x0,xthr,ninteg,integrand1,f1trunc)
	call DRG20r(x0,xthr,ninteg,integrand2,f2trunc)
	call DRG20r(x0,xthr,ninteg,integrand3,fltrunc)
	call DRG20r(x0,xthr,ninteg,integrand4,g1trunc)
	call DRG20r(x0,xthr,ninteg,integrand5,g2trunc)
	call dsg20r(x02,x0,ninteg,xrd,n1p)
	do i1=1,n1p
                nu=qsq/(2.d0*mn*xrd(i1))
                wsq=mn**2-qsq+2.d0*mn*nu
		call res_calc(wsq,qsq)
		integrand1(i1)=f1tot
		integrand2(i1)=f2tot
		integrand3(i1)=fltot
		integrand4(i1)=g1tot
		integrand5(i1)=g2tot
	enddo
	call DRG20r(x02,x0,ninteg,integrand1,f1trunc2)
	call DRG20r(x02,x0,ninteg,integrand2,f2trunc2)
	call DRG20r(x02,x0,ninteg,integrand3,fltrunc2)
	call DRG20r(x02,x0,ninteg,integrand4,g1trunc2)
	call DRG20r(x02,x0,ninteg,integrand5,g2trunc2)
	call dsg20r(x03,x02,ninteg,xrd,n1p)
	do i1=1,n1p
                nu=qsq/(2.d0*mn*xrd(i1))
                wsq=mn**2-qsq+2.d0*mn*nu
		call res_calc(wsq,qsq)
		integrand1(i1)=f1tot
		integrand2(i1)=f2tot
		integrand3(i1)=fltot
		integrand4(i1)=g1tot
		integrand5(i1)=g2tot
	enddo
	call DRG20r(x03,x02,ninteg,integrand1,f1trunc3)
	call DRG20r(x03,x02,ninteg,integrand2,f2trunc3)
	call DRG20r(x03,x02,ninteg,integrand3,fltrunc3)
	call DRG20r(x03,x02,ninteg,integrand4,g1trunc3)
	call DRG20r(x03,x02,ninteg,integrand5,g2trunc3)
	call dsg20r(x03,xthr,ninteg,xrd,n1p)
	do i1=1,n1p
                nu=qsq/(2.d0*mn*xrd(i1))
                wsq=mn**2-qsq+2.d0*mn*nu
		call res_calc(wsq,qsq)
		integrand1(i1)=f1tot
		integrand2(i1)=f2tot
		integrand3(i1)=fltot
		integrand4(i1)=g1tot
		integrand5(i1)=g2tot
	enddo
	call DRG20r(x03,xthr,ninteg,integrand1,f1truncall)
	call DRG20r(x03,xthr,ninteg,integrand2,f2truncall)
	call DRG20r(x03,xthr,ninteg,integrand3,fltruncall)
	call DRG20r(x03,xthr,ninteg,integrand4,g1truncall)
	call DRG20r(x03,xthr,ninteg,integrand5,g2truncall)
	write (iwrt1,"(F10.2,F10.4,F10.4,F10.4,F10.4,F10.4)") qsq,
     >		f1trunc,f2trunc,fltrunc,g1trunc,g2trunc
	write (iwrt2,"(F10.2,F10.4,F10.4,F10.4,F10.4,F10.4)") qsq,
     >		f1trunc2,f2trunc2,fltrunc2,g1trunc2,g2trunc2
	write (iwrt3,"(F10.2,F10.4,F10.4,F10.4,F10.4,F10.4)") qsq,
     >		f1trunc3,f2trunc3,fltrunc3,g1trunc3,g2trunc3
	write (iwrt4,"(F10.2,F10.4,F10.4,F10.4,F10.4,F10.4)") qsq,
     >		f1truncall,f2truncall,fltruncall
     >		,g1truncall,g2truncall
	enddo
	close (iwrt1)
	close (iwrt2)
	close (iwrt3)
	close (iwrt4)
	end subroutine write_trunc_mom
