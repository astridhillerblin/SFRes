!======================================================================!
! SUBROUTINE FOR RESONANT CONTRIBUTIONS TO ALL OBSERVABLES
!======================================================================!
	subroutine res_calc(wsq,qsq,ebeam)
	use res_decl
	use phys_consts
	use math_consts
	implicit none
	real(8) interf
	real(8) wsq,qsq,ebeam
	real(8) f1tot,f2tot,fltot,g1tot,g2tot
	real(8) sigt,sigl,dxsdq2dw
	real(8) h12tot,h32tot
	real(8) a1tot,a2tot
	real(8), dimension (19) :: coupsqt,coupsql,coupt1,coupt2,coupl
	real(8), dimension (19) :: flag
        complex(16) gpm0r,gpm0rdel
        complex(16) gpsumsqr,gpwgtsumsqr,gmsumsqr,g0sumsqr
	integer ind
	integer isup
	real(8) xbj,nu
	real(8) eps,gamv
	common/sf/f1tot,f2tot,fltot,g1tot,g2tot
	common/hel/h12tot,h32tot
	common/asym/a1tot,a2tot
	common/xs/sigt,sigl,dxsdq2dw
	common/sing/flag
	common/isupov/isup
	common/swinterf/interf
! Bjorken x and nu definitions
	xbj = qsq/(wsq-mn**2+qsq)
	nu = (wsq-mn**2+qsq)/(2.d0*mn)
! Initializes |G+|^2,|G0|^2,|G-|^2,(G+)(G0*)
        gpsumsqr=complex(0.d0,0.d0)
        gpwgtsumsqr=complex(0.d0,0.d0)
        gmsumsqr=complex(0.d0,0.d0)
        g0sumsqr=complex(0.d0,0.d0)
! Either call central electrocoupling values, or random MC for uncertainty bands
	if (isup.eq.0) then
		call listecoupread(wsq,qsq,coupsqt,coupsql,coupt1,coupt2,coupl)
	else if (isup.eq.1) then
		call listecoup(wsq,qsq,coupsqt,coupsql,coupt1,coupt2,coupl)
	endif
! Calculates |G+|^2,|G0|^2,|G-|^2,(G+)(G0*) without interferences, by summing over
! all resonances (except for Delta(1232)3/2+ handled differently right below)
	do ind=1,19
	if (ind.ne.12) then
                gpsumsqr=gpsumsqr+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*
     >		conjg(gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1))*flag(ind)
                gmsumsqr=gmsumsqr+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*
     >		conjg(gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1))*flag(ind)
                g0sumsqr=g0sumsqr+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*
     >		conjg(gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0))*flag(ind)
                gpwgtsumsqr=gpwgtsumsqr+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)*
! Note that for the mixed (G+)(G0*) term the above factor is needed,
! see also Eqs.10 in our work 2212.11952
     >		conjg(gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0))
	endif
	enddo
! The Delta contributions are handled here
! The reason for separating them is due to the different treatment of 
! the hadronic decay width 
        gpsumsqr=gpsumsqr+gpm0rdel(wsq,qsq,coupt1(12),1)*
     >		conjg(gpm0rdel(wsq,qsq,coupt1(12),1))*flag(12)
        gmsumsqr=gmsumsqr+gpm0rdel(wsq,qsq,coupt2(12),-1)*
     >		conjg(gpm0rdel(wsq,qsq,coupt2(12),-1))*flag(12)
        g0sumsqr=g0sumsqr+gpm0rdel(wsq,qsq,coupl(12),0)*
     >		conjg(gpm0rdel(wsq,qsq,coupl(12),0))*flag(12)
        gpwgtsumsqr=gpwgtsumsqr+gpm0rdel(wsq,qsq,coupt1(12),1)
     >			*flag(12)*(-1)**(jres(12)-0.5d0)*parres(12)*
     >		conjg(gpm0rdel(wsq,qsq,coupl(12),0))
!======================================================================!
! Now the interference terms are included.
! This is done by taking into account the interferences between 
! resonances with the same quantum numbers.
! Ultimately, this affects the interference between:
! 1) N(1720)3/2+ and N'(1720)3/2+ (indices 8 and 19)
! 1) N(1440)1/2+ and N(1710)1/2+ (indices 1 and 7)
! 1) N(1535)1/2- and N(1650)1/2+ (indices 3 and 4)
! Note also that the interference between N(1720) and N'(1720)
! is taken into account differently than the ohers, by reducing their
! interference by a factor 1.72.
!======================================================================!
! As explained just above, interference terms in |G+|^2:
        gpsumsqr=gpsumsqr+interf*(
     >		(gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupt1(8),jres(8),parres(8),1)*
     >		conjg(gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupt1(19),jres(19),parres(19),1))
     >		+conjg(gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupt1(8),jres(8),parres(8),1))*
     >		gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupt1(19),jres(19),parres(19),1))
     >		*flag(19)*flag(8)/1.72d0+
     >		(gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupt1(1),jres(1),parres(1),1)*
     >		conjg(gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupt1(7),jres(7),parres(7),1))
     >		+conjg(gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupt1(1),jres(1),parres(1),1))*
     >		gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupt1(7),jres(7),parres(7),1))
     >		*flag(7)*flag(1)+
     >		(gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupt1(3),jres(3),parres(3),1)*
     >		conjg(gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupt1(4),jres(4),parres(4),1))
     >		+conjg(gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupt1(3),jres(3),parres(3),1))*
     >		gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupt1(4),jres(4),parres(4),1))
     >		*flag(3)*flag(4)
     >		)
! As explained above, interference terms in |G-|^2:
        gmsumsqr=gmsumsqr+interf*(
     >		(gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupt2(8),jres(8),parres(8),-1)*
     >		conjg(gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupt2(19),jres(19),parres(19),-1))
     >		+conjg(gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupt2(8),jres(8),parres(8),-1))*
     >		gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupt2(19),jres(19),parres(19),-1))
     >		*flag(19)*flag(8)/1.72d0+
     >		(gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupt2(1),jres(1),parres(1),-1)*
     >		conjg(gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupt2(7),jres(7),parres(7),-1))
     >		+conjg(gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupt2(1),jres(1),parres(1),-1))*
     >		gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupt2(7),jres(7),parres(7),-1))
     >		*flag(7)*flag(1)+
     >		(gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupt2(3),jres(3),parres(3),-1)*
     >		conjg(gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupt2(4),jres(4),parres(4),-1))
     >		+conjg(gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupt2(3),jres(3),parres(3),-1))*
     >		gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupt2(4),jres(4),parres(4),-1))
     >		*flag(3)*flag(4)
     >		)
! As explained above, interference terms in |G0|^2:
        g0sumsqr=g0sumsqr+interf*(
     >		(gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupl(8),jres(8),parres(8),0)*
     >		conjg(gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupl(19),jres(19),parres(19),0))
     >		+conjg(gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupl(8),jres(8),parres(8),0))*
     >		gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupl(19),jres(19),parres(19),0))
     >		*flag(19)*flag(8)/1.72d0+
     >		(gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupl(1),jres(1),parres(1),0)*
     >		conjg(gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupl(7),jres(7),parres(7),0))
     >		+conjg(gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupl(1),jres(1),parres(1),0))*
     >		gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupl(7),jres(7),parres(7),0))
     >		*flag(7)*flag(1)+
     >		(gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupl(3),jres(3),parres(3),0)*
     >		conjg(gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupl(4),jres(4),parres(4),0))
     >		+conjg(gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupl(3),jres(3),parres(3),0))*
     >		gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupl(4),jres(4),parres(4),0))
     >		*flag(3)*flag(4)
     >		)
! As explained above, interference terms in (G+)(G0*):
        gpwgtsumsqr=gpwgtsumsqr+interf*(
     >		(gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupt1(8),jres(8),parres(8),1)*
     >		(-1)**(jres(8)-0.5d0)*parres(8)*
     >		conjg(gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupl(19),jres(19),parres(19),0))
     >		+conjg(gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupl(8),jres(8),parres(8),0))*
     >		gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupt1(19),jres(19),parres(19),1)*
     >		(-1)**(jres(19)-0.5d0)*parres(19))
     >		*flag(19)*flag(8)/1.72d0+
     >		(gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupt1(1),jres(1),parres(1),1)*
     >		(-1)**(jres(1)-0.5d0)*parres(1)*
     >		conjg(gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupl(7),jres(7),parres(7),0))
     >		+conjg(gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupl(1),jres(1),parres(1),0))*
     >		gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupt1(7),jres(7),parres(7),1)*
     >		(-1)**(jres(7)-0.5d0)*parres(7))
     >		*flag(7)*flag(1)+
     >		(gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupt1(3),jres(3),parres(3),1)*
     >		(-1)**(jres(3)-0.5d0)*parres(3)*
     >		conjg(gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupl(4),jres(4),parres(4),0))
     >		+conjg(gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupl(3),jres(3),parres(3),0))*
     >		gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupt1(4),jres(4),parres(4),1)*
     >		(-1)**(jres(4)-0.5d0)*parres(4))
     >		*flag(3)*flag(4)
     >		)
! Finally, the observables can be calculated from the above amplitudes
        g1tot = real(mn**2/(1.d0+qsq/nu**2)*(gpsumsqr-gmsumsqr+
     >			gpwgtsumsqr*dsqrt(2.d0*qsq)/nu))
        g2tot = real(-mn**2/(1.d0+qsq/nu**2)*(gpsumsqr-gmsumsqr-
     >			gpwgtsumsqr*dsqrt(2.d0/qsq)*nu))
        f1tot = real(mn**2*(gpsumsqr+gmsumsqr))
        f2tot = real(mn*nu/(1.d0+nu**2/qsq)*(gpsumsqr+
     >	gmsumsqr+2.d0*g0sumsqr))
        fltot = (1.d0+4.d0*mn**2*xbj**2/qsq)*f2tot
     >		-2.d0*xbj*f1tot
	h12tot=f1tot+g1tot-qsq/nu**2*g2tot
	h32tot=f1tot-g1tot+qsq/nu**2*g2tot
	a1tot=(g1tot-4.d0*mn**2*xbj**2/qsq*g2tot)/f1tot
	a2tot=sqrt(qsq)/nu*(g1tot+g2tot)/f1tot
	sigt=8.d0*pi**2*alpha/(wsq-mn**2)*f1tot*389.379
	sigl=8.d0*pi**2*alpha/(wsq-mn**2)*
     >		((1.d0+qsq/nu**2)/(2.d0*xbj)*f2tot-f1tot)*389.379
	dxsdq2dw=1000.d0*
     >		(sigt+eps(wsq,qsq,ebeam)*sigl)*gamv(wsq,qsq,ebeam)
	endsubroutine res_calc
!======================================================================!
! G+, G0, G- helicity amplitudes
!======================================================================!
	complex(16) function gpm0r(wsq,qsq,gres,mres,xbc,lres,
     >		bfpi,bf2pi,bfeta,ecoup,jres,parres,plusminus)
	use phys_consts
	use math_consts
	implicit none
        integer plusminus
	real(8) wsq,qsq,gres,mres,xbc,bfpi,bf2pi,bfeta,ecoup,jres
	integer lres,parres
	real(8) ghad
        real(8) p1232,qgam
        p1232 = parres*(-1.d0)**(jres-0.5d0)
	gpm0r = dsqrt((mres**2-mn**2)/mres)*ecoup/(2.d0*mn)
        gpm0r = gpm0r/dsqrt(pi)*mres
     >		*ghad(wsq,gres,mres,xbc,lres,bfpi,bf2pi,bfeta)/
     >		(mres**2-wsq-complex(0.d0,1.d0)*mres*
     >		ghad(wsq,gres,mres,xbc,lres,bfpi,bf2pi,bfeta))
	gpm0r=gpm0r*qgam(mres**2,qsq)/qgam(wsq,qsq)*
     >		dsqrt(2.d0*dsqrt(wsq)*mn*(wsq-mn**2)/
     >		(2.d0*dsqrt(wsq)*pi*alpha*(mres**2-mn**2))/
     >		ghad(wsq,gres,mres,xbc,lres,bfpi,bf2pi,bfeta))
        if (plusminus.eq.-1) then
          gpm0r=gpm0r*p1232
        else if (plusminus.eq.0) then
          gpm0r=gpm0r*p1232/dsqrt(2.d0)
        endif
	end function gpm0r
!======================================================================!
! G+, G0, G- helicity amplitudes for Delta(1232)3/2+
!======================================================================!
	complex(16) function gpm0rdel(wsq,qsq,ecoup,plusminus)
	use res_decl
	use phys_consts
	use math_consts
	implicit none
        integer plusminus
	real(8) wsq,qsq,ecoup
	real(8) ghaddel
	real(8) mompi
	real(8) b1j,b1y,b1jr,b1yr
        real(8) p1232,qgam
        p1232 = parres(12)*(-1.d0)**(jres(12)-0.5d0)
	b1jr = 0.488866d0
	b1yr = -0.648204d0
	b1j = dbesj1(5.07614d0*mompi(wsq))
	b1y = dbesy1(5.07614d0*mompi(wsq))
	ghaddel=gres(12)*mres(12)/dsqrt(wsq)*(b1jr**2+b1yr**2)
     >		/(b1j**2+b1y**2)
	gpm0rdel = dsqrt((mres(12)**2-mn**2)/mres(12))
     >	*ecoup/(2.d0*mn)
        gpm0rdel = gpm0rdel/dsqrt(pi)*mres(12)
     >		*ghaddel/(mres(12)**2-wsq-
     >		complex(0.d0,1.d0)*mres(12)*ghaddel)
	gpm0rdel=gpm0rdel*qgam(mres(12)**2,qsq)/qgam(wsq,qsq)*
     >		dsqrt(2.d0*dsqrt(wsq)*mn*(wsq-mn**2)/
     >		(2.d0*dsqrt(wsq)*pi*alpha*(mres(12)**2-mn**2))/
     >		ghaddel)
        if (plusminus.eq.-1) then
          gpm0rdel=gpm0rdel*p1232
        else if (plusminus.eq.0) then
          gpm0rdel=gpm0rdel*p1232/dsqrt(2.d0)
        endif
	end function gpm0rdel
!======================================================================!
! FUNCTION OF THE ELECTROMAGNETIC WIDTH
!======================================================================!
	real(8) function qgam(wsq,qsq)
	use phys_consts
	implicit none
	real(8) wsq,qsq
	real(8) egam
	egam = (wsq-qsq-mn**2)/(2.d0*dsqrt(wsq))
	qgam = dsqrt(qsq+egam**2)
	end function qgam
!======================================================================!
! SUBROUTINE FOR ELECTROCOUPLINGS
! All the functions, except where explicitly noted, have been fitted to
! electrocoupling data by E. Isupov.
!======================================================================!
	subroutine listecoup(wsq,qsq,listt,listl,listt1,listt2,listl1)
	use phys_consts
	implicit none
	integer ind
	real(8) wsq,qsq,kcm
	real(8) qgam
	real(8), dimension (19) :: listt,listl,listt1,listt2,listl1
	if (qsq.eq.0.d0) then
		do ind=1,19
			listl(ind) = 0.d0
		enddo
	else
	listl(1) = (31.19227d0+3.53338d0*qsq)/
     >		(1.d0-0.278265d0*qsq**2+0.3677575d0*qsq**2*dsqrt(qsq))
	listl(2) = (-67.32d0)/
     >	(1.d0+1.73d0*qsq-2.8d0*qsq**2+2.91d0*qsq**2*dsqrt(qsq))
	listl(3) = (-9.758811d0-4.231412d0*qsq)/
     >		(1.d0-0.7341952d0*qsq**2+0.5087887d0*qsq**2*dsqrt(qsq))
	listl(4) = (-2.67d0)/
     >		(1.d0-2.82d0*(qsq+0.1d0)+
     >		2.d0*(qsq+0.1d0)*dsqrt(qsq+0.1d0))
	listl(5) = (-213.d0-60.25d0*qsq)/
     >		(1.d0-16.376d0*qsq**2+68.87d0*qsq*dsqrt(qsq))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This electrocoupling has been refitted by Astrid in October 2021,
! since we want to use only the CLAS values from references [7,8]
! on the electrocouplings CLAS website.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	listl(6) = (-50.32d0)/
     >		(1.d0+2.55d0*qsq**4)! Re-fit October 2021
	listl(7) = (-7.25d0)/
     >		(1.d0+0.0733d0*qsq**(1.5d0))
!=====================================================================!
! The following is because we do not trust the parameterization of
! P13(1720) below Q2 < 0.95 GeV^2, so we take the interpolation between
! the data point at 0.65 GeV^2 (-32.0) and the function at 0.95 GeV^2.
! In addition, note the usual 0.1 shift in Q^2.
!=====================================================================!
	if (qsq+0.1d0 .lt. 0.95d0) then
	listl(8) = -23.1975d0-(-23.1975d0+32.d0)*(0.95d0-(qsq+0.1d0))
     >		/(0.95d0-0.65d0)
	else
	listl(8) = (-3.09d0)/
     >		(1.d0-3.7d0*(qsq+0.1d0)+
     >		2.86d0*(qsq+0.1d0)*dsqrt(qsq+0.1d0))
	endif
	listl(9) = 0.d0
	listl(10) = 0.d0
	listl(11) = 0.d0
	listl(12) = (21.64d0-12.41d0*qsq+1.909d0*qsq**2)/
     >		(1.d0-0.4583d0*qsq+0.1422d0*qsq**2
     >		-0.0525d0*qsq**3+0.00931d0*qsq**4)
	listl(13) = -(25.0d0)/
     >		((qsq+0.1d0)*dsqrt(qsq+0.1d0))
	listl(14) = (24.5d0)/
     >		(1.d0+0.15d0*(qsq+0.1d0)*dsqrt(qsq+0.1d0))
	listl(15) = 0.d0
	listl(16) = 0.d0
	listl(17) = 0.d0
	listl(18) = 0.d0
	if (qsq+0.1d0 .le. 2.d0) then
	listl(19) = (93.056d0)/
     >		(1.d0+1.61d0*(qsq+0.1d0)*dsqrt(qsq+0.1d0))
	else if ((qsq+0.1d0 .gt. 2.d0).and.(qsq+0.1d0 .le. 4.5d0)) then
	listl(19) = (93.056d0)/
     >		(1.d0+1.61d0*(qsq+0.1d0)*dsqrt(qsq+0.1d0))
     >		*(1.d0-0.24d0*(qsq+0.1d0-2.d0))
	else if (qsq+0.1d0 .gt. 4.5d0) then
	listl(19) = (93.056d0)/
     >		(1.d0+1.61d0*(qsq+0.1d0)*dsqrt(qsq+0.1d0))*
     >		(0.4d0-0.24d0/5.d0*(qsq+0.1d0-4.5d0))
	endif
!======================================================================!
! The second endif is the question whether Q^2 = 0.
!======================================================================!
	endif
	listt1(1) = (-68.7866d0+21.3966d0*qsq+79.8415d0*dsqrt(qsq))/
     >		(1.d0-0.7178d0*qsq**2+0.5663d0*qsq**2*dsqrt(qsq))
	listt1(2) = 0.9d0*(-23.357d0-151.199533d0*qsq)/
     >	(1.d0+2.01489898d0*qsq**2-0.2654327d0*qsq**2*dsqrt(qsq))
	listt1(3) = (92.5029d0+1.45023d0*qsq)/
     >		(1.d0+0.1095d0*qsq**2-0.000322d0*qsq**2*dsqrt(qsq))
	listt1(4) = (47.4d0-19.6d0*(qsq+0.1d0))/
     >		(1.d0-1.46d0*(qsq+0.1d0)*dsqrt(qsq+0.1d0)
     >		+1.17d0*(qsq+0.1d0)**3)
	listt1(5) = (14.85d0-2.49d0*qsq)/
     >		(1.d0 - 0.597d0*qsq*dsqrt(qsq)+0.368d0*qsq**2)
	listt1(6) = -20.d0*(1.d0+3.04d0*qsq+0.0305d0*qsq**2)/
     >		((1.d0+0.034d0*(qsq-2.d0)**2)*(1.d0+0.09d0*qsq**2)
     >		*(1.d0+0.765d0*qsq*dsqrt(qsq)))
	listt1(7) = (27.36d0)/
     >		(1.d0+0.77d0*qsq**(1.5d0))
	listt1(8) = (90.d0)/
     >		(1.d0+3.16d0*(qsq+0.1d0)*dsqrt(qsq+0.1d0))
	listt1(9) = 0.d0
	listt1(10) = 0.d0
	listt1(11) = 0.d0
	listt1(12) = -(178.45d0)/
     >		((1.d0+qsq)*(1.d0+0.3457d0*qsq**2
     >		-0.087d0*qsq**3+0.00806d0*qsq**4))
	listt1(13) = (47.2d0)/
     >		(1.d0+3.71d0*(qsq+0.1d0)*dsqrt(qsq+0.1d0))
	listt1(14) = (118.28d0)/
     >		(1.d0+0.72d0*(qsq+0.1d0)**2+3.26d0*(qsq+0.1d0)**2.5d0)
	listt1(15) = 0.d0
	listt1(16) = 0.d0
	listt1(17) = 0.d0
	listt1(18) = 0.d0
	if (qsq+0.1d0 .le. 2.d0) then
	listt1(19) = (36.5d0+1044.05d0*(qsq+0.1d0))/
     >		(1.d0+19.5d0*(qsq+0.1d0)*dsqrt(qsq+0.1d0))
	else if ((qsq+0.1d0 .gt. 2.d0).and.(qsq+0.1d0 .le. 4.5d0)) then
	listt1(19) = (36.5d0+1044.05d0*(qsq+0.1d0))/
     >		(1.d0+19.5d0*(qsq+0.1d0)*dsqrt(qsq+0.1d0))*
     >		(1.d0-0.24d0*(qsq+0.1d0-2.d0))
	else if (qsq+0.1d0 .gt. 4.5d0) then
	listt1(19) = (36.5d0+1044.05d0*(qsq+0.1d0))/
     >		(1.d0+19.5d0*(qsq+0.1d0)*dsqrt(qsq+0.1d0))*
     >		(0.4d0-0.24d0/5.d0*(qsq+0.1d0-4.5d0))
	endif
	listt2(1) = 0.d0
	listt2(2) = 0.9d0*(162.458285d0)/
     >	(1.d0+3.322979d0*qsq-2.0339966d0*qsq**2+
     >	1.622563d0*qsq**2*dsqrt(qsq))
	listt2(3) = 0.d0
	listt2(4) = 0.d0
	listt2(5) = (14.95d0-17.8d0*dsqrt(qsq)+4.75d0*qsq)/
     >		(1.d0-0.78d0*qsq*dsqrt(qsq)+0.405d0*qsq**2)
	listt2(6) = (134.2d0)/
     >		(1.d0+4.83d0*qsq-2.68d0*qsq**2+1.1d0*qsq**2*dsqrt(qsq))
	listt2(7) = 0.d0
	listt2(8) = (-35.87d0-6.85d0*(qsq+0.1d0))/
     >		(1.d0+0.118d0*(qsq+0.1d0)**4)
	listt2(9) = 0.d0
	listt2(10) = 0.d0
	listt2(11) = 0.d0
	listt2(12) = -(339.06d0)/
     >		((1.d0+qsq)*(1.d0+0.3481d0*qsq**2-0.0854d0*qsq**3 +
     >		0.00758d0*qsq**4))
	listt2(13) = 0.d0
	listt2(14) = (101.d0)/
     >		(1.d0+16.6d0*(qsq+0.1d0)**2-5.0d0*(qsq+0.1d0)**2
     >		*dsqrt(qsq+0.1d0))
	listt2(15) = 0.d0
	listt2(16) = 0.d0
	listt2(17) = 0.d0
	listt2(18) = 0.d0
!======================================================================!
! The following is because we do not trust the parameterization of
! the new state N'(1720) below Q2 < 1.0 GeV^2, so we take a constant
! value below that
!======================================================================!
	if ((qsq+0.1d0) .le. 1.d0) then
	listt2(19) = (-37.9d0)/
     >		(1.d0+0.455d0)
	else if (((qsq+0.1d0).gt.1.d0).and.((qsq+0.1d0).le.2.d0)) then
	listt2(19) = (-37.9d0)/
     >		(1.d0+0.455d0*(qsq+0.1d0)**2*dsqrt(qsq+0.1d0))
	else if (((qsq+0.1d0).gt.2.d0).and.((qsq+0.1d0).le.4.5d0)) then
	listt2(19) = (-37.9d0)/
     >		(1.d0+0.455d0*(qsq+0.1d0)**2*dsqrt(qsq+0.1d0))*
     >		(1.d0-0.24d0*(qsq+0.1d0-2.d0))
	else if ((qsq+0.1d0).gt.(4.5d0)) then
	listt2(19) = (-37.9d0)/
     >		(1.d0+0.455d0*(qsq+0.1d0)**2*dsqrt(qsq+0.1d0))*
     >		(0.4d0-0.24d0/5.d0*(qsq+0.1d0-4.5d0))
	endif
! Overall order of magnitude of the units.
	do ind=1,19
		listl1(ind) = dsqrt(2.d0)*10.d0**(-3)*listl(ind)
		listl(ind) = 2.d0*10.d0**(-6)*listl(ind)**2
		listt(ind) =10.d0**(-6)*(listt1(ind)**2+listt2(ind)**2)
		listt1(ind) = 10.d0**(-3)*listt1(ind)
		listt2(ind) = 10.d0**(-3)*listt2(ind)
	enddo
	kcm = (wsq-mn**2)/(2.d0*dsqrt(wsq))
! Resonances 8, 13, 14, 19 need to be weighted with other kinematic
! factors than the rest, to stay consistent with the way they
! were extracted from data.
	listl(8) = listl(8)/2.d0*qgam(wsq,qsq)/kcm
	listl(13) = listl(13)/2.d0*qgam(wsq,qsq)/kcm
	listl(14) = listl(14)/2.d0*qgam(wsq,qsq)/kcm
	listl(19) = listl(19)/2.d0*qgam(wsq,qsq)/kcm
	listt(8) = listt(8)*qgam(wsq,qsq)/kcm
	listt(13) = listt(13)*qgam(wsq,qsq)/kcm
	listt(14) = listt(14)*qgam(wsq,qsq)/kcm
	listt(19) = listt(19)*qgam(wsq,qsq)/kcm
	listt1(8) = listt1(8)*dsqrt(qgam(wsq,qsq)/kcm)
	listt1(13) = listt1(13)*dsqrt(qgam(wsq,qsq)/kcm)
	listt1(14) = listt1(14)*dsqrt(qgam(wsq,qsq)/kcm)
	listt1(19) = listt1(19)*dsqrt(qgam(wsq,qsq)/kcm)
	listt2(8) = listt2(8)*dsqrt(qgam(wsq,qsq)/kcm)
	listt2(13) = listt2(13)*dsqrt(qgam(wsq,qsq)/kcm)
	listt2(14) = listt2(14)*dsqrt(qgam(wsq,qsq)/kcm)
	listt2(19) = listt2(19)*dsqrt(qgam(wsq,qsq)/kcm)
	listl1(8) = listl1(8)*dsqrt(qgam(wsq,qsq)/kcm/2.d0)
	listl1(13) = listl1(13)*dsqrt(qgam(wsq,qsq)/kcm/2.d0)
	listl1(14) = listl1(14)*dsqrt(qgam(wsq,qsq)/kcm/2.d0)
	listl1(19) = listl1(19)*dsqrt(qgam(wsq,qsq)/kcm/2.d0)
	endsubroutine listecoup
!======================================================================!
! SUBROUTINE FOR USING RANDOMLY GENERATED ELECTROCOUPLINGS
!======================================================================!
	subroutine listecoupread(wsq,qsq,listt,listl,listt1,listt2,
     >		listl1)
	use phys_consts
	implicit none
	integer ind
	real(8) wsq,qsq,kcm
	real(8) qgam
	real(8), dimension (19) :: listt,listl,listt1,listt2,listl1
	do ind=1,19
	call mcamplitudes(ind,qsq,listt1(ind),listt2(ind),listl(ind))
	enddo
	do ind=1,19
		listl1(ind) = dsqrt(2.d0)*10.d0**(-3)*listl(ind)
		listl(ind) = 2.d0*10.d0**(-6)*listl(ind)**2
		listt(ind) =10.d0**(-6)*(listt1(ind)**2+listt2(ind)**2)
		listt1(ind) = 10.d0**(-3)*listt1(ind)
		listt2(ind) = 10.d0**(-3)*listt2(ind)
	enddo
	kcm = (wsq-mn**2)/(2.d0*dsqrt(wsq))
	listl(8) = listl(8)/2.d0*qgam(wsq,qsq)/kcm
	listl(13) = listl(13)/2.d0*qgam(wsq,qsq)/kcm
	listl(14) = listl(14)/2.d0*qgam(wsq,qsq)/kcm
	listl(19) = listl(19)/2.d0*qgam(wsq,qsq)/kcm
	listt(8) = listt(8)*qgam(wsq,qsq)/kcm
	listt(13) = listt(13)*qgam(wsq,qsq)/kcm
	listt(14) = listt(14)*qgam(wsq,qsq)/kcm
	listt(19) = listt(19)*qgam(wsq,qsq)/kcm
	listt1(8) = listt1(8)*dsqrt(qgam(wsq,qsq)/kcm)
	listt1(13) = listt1(13)*dsqrt(qgam(wsq,qsq)/kcm)
	listt1(14) = listt1(14)*dsqrt(qgam(wsq,qsq)/kcm)
	listt1(19) = listt1(19)*dsqrt(qgam(wsq,qsq)/kcm)
	listt2(8) = listt2(8)*dsqrt(qgam(wsq,qsq)/kcm)
	listt2(13) = listt2(13)*dsqrt(qgam(wsq,qsq)/kcm)
	listt2(14) = listt2(14)*dsqrt(qgam(wsq,qsq)/kcm)
	listt2(19) = listt2(19)*dsqrt(qgam(wsq,qsq)/kcm)
	listl1(8) = listl1(8)*dsqrt(qgam(wsq,qsq)/kcm/2.d0)
	listl1(13) = listl1(13)*dsqrt(qgam(wsq,qsq)/kcm/2.d0)
	listl1(14) = listl1(14)*dsqrt(qgam(wsq,qsq)/kcm/2.d0)
	listl1(19) = listl1(19)*dsqrt(qgam(wsq,qsq)/kcm/2.d0)
	endsubroutine listecoupread
!======================================================================!
! FUNCTION OF THE HADRONIC WIDTH INTO 2PI
!======================================================================!
	real(8) function g2pi(wsq,gres,mres,xbc,lres)
	use phys_consts
	implicit none
	real(8) wsq
	real(8) gres
	real(8) mres
	real(8) xbc
	integer lres
	real(8) mom2pi
	g2pi = gres*(mom2pi(wsq)/mom2pi(mres**2))**(2*lres+4)*
     >		((mom2pi(mres**2)**2+xbc**2)/
     >		(mom2pi(wsq)**2+xbc**2))**(lres+2)
	end function g2pi
!======================================================================!
! FUNCTION OF THE HADRONIC WIDTH INTO PI
!======================================================================!
	real(8) function gpi(wsq,gres,mres,xbc,lres)
	use phys_consts
	implicit none
	real(8) wsq
	real(8) gres
	real(8) mres
	real(8) xbc
	integer lres
	real(8) mompi
	gpi = gres*(mompi(wsq)/mompi(mres**2))**(2*lres+1)*
     >		((mompi(mres**2)**2+xbc**2)/
     >		(mompi(wsq)**2+xbc**2))**lres
	end function gpi
!======================================================================!
! FUNCTION OF THE HADRONIC WIDTH INTO ETA
!======================================================================!
	real(8) function geta(wsq,gres,mres,xbc,lres)
	use phys_consts
	implicit none
	real(8) wsq
	real(8) gres
	real(8) mres
	real(8) xbc
	integer lres
	real(8) mometa
	geta = gres*(mometa(wsq)/mometa(mres**2))**(2*lres+1)*
     >		((mometa(mres**2)**2+xbc**2)/
     >		(mometa(wsq)**2+xbc**2))**lres
	end function geta
!======================================================================!
! FUNCTION THAT CALCULATES THE DECAY PION MOMENTUM
!======================================================================!
	real(8) function mompi(wsq)
	use phys_consts
	implicit none
	real(8) wsq
	if (wsq < (mn+mpi)**2) then
		mompi = 0.
	else
		mompi = 1./(2.*dsqrt(wsq))*
     >			dsqrt((wsq-(mn+mpi)**2)*(wsq-(mn-mpi)**2))
	endif
	end function mompi
!======================================================================!
! FUNCTION THAT CALCULATES THE DOUBLE-PION SYSTEM MOMENTUM
!======================================================================!
	real(8) function mom2pi(wsq)
	use phys_consts
	implicit none
	real(8) wsq
	if (wsq < (mn+2.*mpi)**2) then
		mom2pi = 0.
	else
		mom2pi = 1./(2.*dsqrt(wsq))*
     >			dsqrt((wsq-(mn+2.*mpi)**2)*(wsq-(mn-2.*mpi)**2))
	endif
	end function mom2pi
!======================================================================!
! FUNCTION THAT CALCULATES THE DECAY ETA MOMENTUM
!======================================================================!
	real(8) function mometa(wsq)
	use phys_consts
	implicit none
	real(8) wsq
	if (wsq < (mn+meta)**2) then
		mometa = 0.
	else
		mometa = 1./(2.*dsqrt(wsq))*
     >			dsqrt((wsq-(mn+meta)**2)*(wsq-(mn-meta)**2))
	endif
	end function mometa
!======================================================================!
! FUNCTION OF THE FULL HADRONIC WIDTH
!======================================================================!
	real(8) function ghad(wsq,gres,mres,xbc,lres,bfpi,bf2pi,bfeta)
	use phys_consts
	implicit none
	real(8) wsq
	real(8) gres
	real(8) mres
	real(8) xbc
	integer lres
	real(8) bfpi, bf2pi, bfeta
	real(8) geta, gpi, g2pi
! The following cases check whether meson production thresholds have
! already been crossed:
	if (wsq < (mn + mpi)**2 .or. mres < mn + mpi) then
		ghad = 0.d0
	else if (wsq < (mn + 2.d0*mpi)**2 .or. mres < mn+2.d0*mpi) then
		ghad = gpi(wsq,gres,mres,xbc,lres)*bfpi
	else if (wsq < (mn + meta)**2 .or. mres < mn + meta) then
		ghad = gpi(wsq,gres,mres,xbc,lres)*bfpi+
     >			g2pi(wsq,gres,mres,xbc,lres)*bf2pi
	else
		ghad = gpi(wsq,gres,mres,xbc,lres)*bfpi+
     >			g2pi(wsq,gres,mres,xbc,lres)*bf2pi+
     >			geta(wsq,gres,mres,xbc,lres)*bfeta
	endif
	end function ghad
!======================================================================!
! FUNCTION THAT CALCULATES THE VIRTUAL PHOTON TRANSVERSE POLARIZATION
!======================================================================!
	real(8) function eps(wsq,qsq,ebeam)
	use phys_consts
	implicit none
	real(8) wsq,qsq,ebeam,nu,sin2th
	nu = (wsq-mn**2+qsq)/(2.*mn)
	sin2th = qsq/(4.d0*ebeam*(ebeam-nu))
	eps = 1./(1.+2.*(1.+nu**2/qsq)*sin2th/(1.d0-sin2th))
	end function eps
!======================================================================!
! FUNCTION THAT CALCULATES THE VIRTUAL PHOTON FLUX GAMMA
!======================================================================!
	real(8) function gamv(wsq,qsq,ebeam)
	use phys_consts
	use math_consts
	implicit none
	real(8) wsq,qsq,ebeam,nu,eps
	nu = (wsq-mn**2+qsq)/(2.*mn)
	gamv = alpha/(4.d0*pi)*dsqrt(wsq)/(ebeam*mn)**2*(wsq-mn**2)/
     >		(qsq*(1-eps(wsq,qsq,ebeam)))
	end function gamv
