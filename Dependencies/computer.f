!======================================================================!
! SUBROUTINE: PRODUCE MEAN AND STD FROM MC SAMPLED OBSERVABLES
!======================================================================!
!======================================================================!
! SUBROUTINE: READ MC with interferences for polarized functions
!======================================================================!
	subroutine mc_meanstd
	implicit none
	integer iloop, jloop, readwhich, readinit
	integer nmc, sig1, sig2, sig3
	real (8), dimension (nmc) :: f1intmc,f2intmc,flintmc,
     >		g1intmc,g2intmc,h12intmc,h32intmc,a1intmc,a2intmc
	real (8), dimension (nmc) :: dismcf1,dismcf2,dismcfl,
     >		dismcg1,dismcg2,dismch12,dismch32,dismca1,dismca2
	real (8) bob, wmc
	real (8) averagef1int,averagef2int,averageflint,
     >		averageg1int,averageg2int,averageh12int,averageh32int,
     >		averagea1int,averagea2int
	common/nsamp/nmc
!======================================================================!
! 0.6827,0.9545,0.9973 are 1,2,3 sigma
!======================================================================!
!This gives the integer number of samples within 1,2,3 sigma
	sig1 = 6827*nmc/10000
!	sig2 = 9545*nmc/10000
!	sig3 = 9973*nmc/10000
	averagef1int = 0.d0
	averagef2int = 0.d0
	averageflint = 0.d0
	averageg1int = 0.d0
	averageg2int = 0.d0
	averageh12int = 0.d0
	averageh32int = 0.d0
	averagea1int = 0.d0
	averagea2int = 0.d0
! Computations are to be done with and without interferences:
	open (101,file = 'Output/res_meanstd.dat')
	open (102,file = 'Output/resinterf_meanstd.dat')
	do iloop=1,2
	write (100+iloop,"(19(A10))") "W [GeV]","<F1>","DF1",
     >	"<F2>","DF2","<FL>","DFL","<g1>","Dg1","<g2>","Dg2",
     >	"<H1/2>","DH1/2","<H3/2>","DH3/2","<A1>","DA1","<A2>","DA2"
	do readinit=1,113
	readwhich=readinit
	open (201,file='Output/res_samp.dat')
	open (202,file='Output/resinterf_samp.dat')
! Reads the sampling files that have been generated.
! For each value of W, sums over all samples to get the mean:
	do jloop =1,nmc*113
	if (modulo(jloop+113-readinit,113).eq.0) then
	read (200+iloop,*) wmc,
     >	f1intmc(readwhich-readinit+1),f2intmc(readwhich-readinit+1),
     >	flintmc(readwhich-readinit+1),
     >	g1intmc(readwhich-readinit+1),g2intmc(readwhich-readinit+1),
     >	h12intmc(readwhich-readinit+1),h32intmc(readwhich-readinit+1),
     >	a1intmc(readwhich-readinit+1),a2intmc(readwhich-readinit+1)
	averagef1int=averagef1int+f1intmc(readwhich-readinit+1)
	averagef2int=averagef2int+f2intmc(readwhich-readinit+1)
	averageflint=averageflint+flintmc(readwhich-readinit+1)
	averageg1int=averageg1int+g1intmc(readwhich-readinit+1)
	averageg2int=averageg2int+g2intmc(readwhich-readinit+1)
	averageh12int=averageh12int+h12intmc(readwhich-readinit+1)
	averageh32int=averageh32int+h32intmc(readwhich-readinit+1)
	averagea1int=averagea1int+a1intmc(readwhich-readinit+1)
	averagea2int=averagea2int+a2intmc(readwhich-readinit+1)
	readwhich=readwhich+1
	else
	read(200+iloop,*) bob,bob,bob,bob,bob,bob,bob,bob,bob,bob
	endif
	enddo
	close (200+iloop)
! After the sum has been done, gets the mean:
	averagef1int=averagef1int/nmc
	averagef2int=averagef2int/nmc
	averageflint=averageflint/nmc
	averageg1int=averageg1int/nmc
	averageg2int=averageg2int/nmc
	averageh12int=averageh12int/nmc
	averageh32int=averageh32int/nmc
	averagea1int=averagea1int/nmc
	averagea2int=averagea2int/nmc
! Gets the distances between sampling values and mean:
	do jloop = 1,nmc
		dismcf1(jloop)=abs(averagef1int-f1intmc(jloop))
		dismcf2(jloop)=abs(averagef2int-f2intmc(jloop))
		dismcfl(jloop)=abs(averageflint-flintmc(jloop))
		dismcg1(jloop)=abs(averageg1int-g1intmc(jloop))
		dismcg2(jloop)=abs(averageg2int-g2intmc(jloop))
		dismch12(jloop)=abs(averageg1int-g1intmc(jloop))
		dismch32(jloop)=abs(averageg2int-g2intmc(jloop))
		dismca1(jloop)=abs(averagea1int-a1intmc(jloop))
		dismca2(jloop)=abs(averagea2int-a2intmc(jloop))
	enddo
! Sorts the sampled observables by their distances to the mean:
	call hpsort(nmc,dismcf1,f1intmc)
	call hpsort(nmc,dismcf2,f2intmc)
	call hpsort(nmc,dismcfl,flintmc)
	call hpsort(nmc,dismcg1,g1intmc)
	call hpsort(nmc,dismcg2,g2intmc)
	call hpsort(nmc,dismch12,h12intmc)
	call hpsort(nmc,dismch32,h32intmc)
	call hpsort(nmc,dismca1,a1intmc)
	call hpsort(nmc,dismca2,a2intmc)
! Stores average and 1-sigma deviation
! (by getting the element at the sig1-th index):
	write(100+iloop,"(F10.2,18(F10.4))") wmc,
     >		averagef1int,abs(averagef1int-f1intmc(sig1)),
     >		averagef2int,abs(averagef2int-f2intmc(sig1)),
     >		averageflint,abs(averageflint-flintmc(sig1)),
     >		averageg1int,abs(averageg1int-g1intmc(sig1)),
     >		averageg2int,abs(averageg2int-g2intmc(sig1)),
     >		averageh12int,abs(averageh12int-h12intmc(sig1)),
     >		averageh32int,abs(averageh32int-h32intmc(sig1)),
     >		averagea1int,abs(averagea1int-a1intmc(sig1)),
     >		averagea2int,abs(averagea2int-a2intmc(sig1))
	enddo
	enddo
	close (101)
	close (102)
	end subroutine mc_meanstd
!======================================================================!
! SUBROUTINE: SORTING ROUTINE
!======================================================================!
	SUBROUTINE HPSORT(N,RA,LI)
	real(8) RA(N), LI(N)
	L=N/2+1
	IR=N
10 	continue
	if(L > 1)then
		L=L-1
		RRA=RA(L)
		RLI=LI(L)
	else
		RRA=RA(IR)
		RA(IR)=RA(1)
		RLI=LI(IR)
		LI(IR)=LI(1)
		IR=IR-1
		if(IR.eq.1)then
			RA(1)=RRA
			LI(1)=RLI
			return
		end if
	end if
	I=L
	J=L+L
20 	if(J.le.IR)then
	if(J < IR)then
	if(RA(J) < RA(J+1))  J=J+1
	end if
	if(RRA < RA(J))then
	RA(I)=RA(J)
	LI(I)=LI(J)
	I=J; J=J+J
	else
	J=IR+1
	end if
	goto 20
	end if
	RA(I)=RRA
	LI(I)=RLI
	goto 10
	END
