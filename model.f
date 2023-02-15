!======================================================================!
! SUBROUTINE FOR RESONANT CONTRIBUTIONS TO OBSERVABLES
!======================================================================!
	subroutine res_calc(wsq,qsq)
	use res_decl
	use phys_consts
	use math_consts
	implicit none
	real(8) interf
	real(8) wsq,qsq
	real(8) f1tot,f2tot,fltot,g1tot,g2tot
	real(8) h12tot,h32tot
	real(8) a1tot,a2tot
	real(8), dimension (19) :: coupsqt,coupsql,coupt1,coupt2,coupl
	real(8), dimension (19) :: flag
        complex(16) gpm0r,gpm0rdel
        complex(16) gpsumsqr,gpwgtsumsqr,gmsumsqr,g0sumsqr
        complex(16) gpsum12pn,gpwgtsum12pn,gmsum12pn,g0sum12pn
        complex(16) gpsum32pn,gpwgtsum32pn,gmsum32pn,g0sum32pn
        complex(16) gpsum52pn,gpwgtsum52pn,gmsum52pn,g0sum52pn
        complex(16) gpsum12pd,gpwgtsum12pd,gmsum12pd,g0sum12pd
        complex(16) gpsum32pd,gpwgtsum32pd,gmsum32pd,g0sum32pd
        complex(16) gpsum12mn,gpwgtsum12mn,gmsum12mn,g0sum12mn
        complex(16) gpsum32mn,gpwgtsum32mn,gmsum32mn,g0sum32mn
        complex(16) gpsum52mn,gpwgtsum52mn,gmsum52mn,g0sum52mn
        complex(16) gpsum12md,gpwgtsum12md,gmsum12md,g0sum12md
        complex(16) gpsum32md,gpwgtsum32md,gmsum32md,g0sum32md
	integer ind
	integer isup
	real(8) xbj,nu
	common/sf/f1tot,f2tot,fltot,g1tot,g2tot
	common/hel/h12tot,h32tot
	common/asym/a1tot,a2tot
	common/sing/flag
	common/isupov/isup
	common/swinterf/interf
	xbj = qsq/(wsq-mn**2+qsq)
	nu = (wsq-mn**2+qsq)/(2.d0*mn)
        gpsumsqr=complex(0.d0,0.d0)
        gpwgtsumsqr=complex(0.d0,0.d0)
        gmsumsqr=complex(0.d0,0.d0)
        g0sumsqr=complex(0.d0,0.d0)
        gpsum12pn=complex(0.d0,0.d0)
        gpwgtsum12pn=complex(0.d0,0.d0)
        gmsum12pn=complex(0.d0,0.d0)
        g0sum12pn=complex(0.d0,0.d0)
        gpsum52pn=complex(0.d0,0.d0)
        gpwgtsum52pn=complex(0.d0,0.d0)
        gmsum52pn=complex(0.d0,0.d0)
        g0sum52pn=complex(0.d0,0.d0)
        gpsum12pd=complex(0.d0,0.d0)
        gpwgtsum12pd=complex(0.d0,0.d0)
        gmsum12pd=complex(0.d0,0.d0)
        g0sum12pd=complex(0.d0,0.d0)
        gpsum32pd=complex(0.d0,0.d0)
        gpwgtsum32pd=complex(0.d0,0.d0)
        gmsum32pd=complex(0.d0,0.d0)
        g0sum32pd=complex(0.d0,0.d0)
        gpsum12mn=complex(0.d0,0.d0)
        gpwgtsum12mn=complex(0.d0,0.d0)
        gmsum12mn=complex(0.d0,0.d0)
        g0sum12mn=complex(0.d0,0.d0)
        gpsum32mn=complex(0.d0,0.d0)
        gpwgtsum32mn=complex(0.d0,0.d0)
        gmsum32mn=complex(0.d0,0.d0)
        g0sum32mn=complex(0.d0,0.d0)
        gpsum52mn=complex(0.d0,0.d0)
        gpwgtsum52mn=complex(0.d0,0.d0)
        gmsum52mn=complex(0.d0,0.d0)
        g0sum52mn=complex(0.d0,0.d0)
        gpsum12md=complex(0.d0,0.d0)
        gpwgtsum12md=complex(0.d0,0.d0)
        gmsum12md=complex(0.d0,0.d0)
        g0sum12md=complex(0.d0,0.d0)
        gpsum32md=complex(0.d0,0.d0)
        gpwgtsum32md=complex(0.d0,0.d0)
        gmsum32md=complex(0.d0,0.d0)
        g0sum32md=complex(0.d0,0.d0)
	if (isup.eq.0) then
		call listecoupread(wsq,qsq,coupsqt,coupsql,coupt1,coupt2,coupl)
	else if (isup.eq.1) then
		call listecoup(wsq,qsq,coupsqt,coupsql,coupt1,coupt2,coupl)
	endif
	do ind=1,11
		if (spiso(ind).eq."12n".and.parres(ind).eq.1) then
                gpsum12pn=gpsum12pn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum12pn=gmsum12pn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum12pn=g0sum12pn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum12pn=gpwgtsum12pn+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		elseif (spiso(ind).eq."32n".and.parres(ind).eq.1) then
                gpsum32pn=gpsum32pn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum32pn=gmsum32pn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum32pn=g0sum32pn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum32pn=gpwgtsum32pn+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		elseif (spiso(ind).eq."52n".and.parres(ind).eq.1) then
                gpsum52pn=gpsum52pn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum52pn=gmsum52pn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum52pn=g0sum52pn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum52pn=gpwgtsum52pn+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		elseif (spiso(ind).eq."12d".and.parres(ind).eq.1) then
                gpsum12pd=gpsum12pd+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum12pd=gmsum12pd+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum12pd=g0sum12pd+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum12pd=gpwgtsum12pd+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		elseif (spiso(ind).eq."32d".and.parres(ind).eq.1) then
                gpsum32pd=gpsum32pd+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum32pd=gmsum32pd+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum32pd=g0sum32pd+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum32pd=gpwgtsum32pd+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		elseif (spiso(ind).eq."12n".and.parres(ind).eq.-1) then
                gpsum12mn=gpsum12mn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum12mn=gmsum12mn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum12mn=g0sum12mn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum12mn=gpwgtsum12mn+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		elseif (spiso(ind).eq."32n".and.parres(ind).eq.-1) then
                gpsum32mn=gpsum32mn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum32mn=gmsum32mn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum32mn=g0sum32mn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum32mn=gpwgtsum32mn+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		elseif (spiso(ind).eq."52n".and.parres(ind).eq.-1) then
                gpsum52mn=gpsum52mn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum52mn=gmsum52mn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum52mn=g0sum52mn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum52mn=gpwgtsum52mn+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		elseif (spiso(ind).eq."12d".and.parres(ind).eq.-1) then
                gpsum12md=gpsum12md+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum12md=gmsum12md+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum12md=g0sum12md+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum12md=gpwgtsum12md+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		elseif (spiso(ind).eq."32d".and.parres(ind).eq.-1) then
                gpsum32md=gpsum32md+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum32md=gmsum32md+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum32md=g0sum32md+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum32md=gpwgtsum32md+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		endif
	enddo
        gpsum32pd=gpsum32pd+gpm0rdel(wsq,qsq,coupt1(12),1)*flag(12)
        gmsum32pd=gmsum32pd+gpm0rdel(wsq,qsq,coupt2(12),-1)*flag(12)
        g0sum32pd=g0sum32pd+gpm0rdel(wsq,qsq,coupl(12),0)*flag(12)
        gpwgtsum32pd=gpwgtsum32pd+gpm0rdel(wsq,qsq,coupt1(12),1)
     >			*flag(12)*(-1)**(jres(ind)-0.5d0)*parres(ind)
	do ind=13,19
		if (spiso(ind).eq."12n".and.parres(ind).eq.1) then
                gpsum12pn=gpsum12pn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum12pn=gmsum12pn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum12pn=g0sum12pn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum12pn=gpwgtsum12pn+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		elseif (spiso(ind).eq."32n".and.parres(ind).eq.1) then
                gpsum32pn=gpsum32pn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum32pn=gmsum32pn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum32pn=g0sum32pn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum32pn=gpwgtsum32pn+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		elseif (spiso(ind).eq."52n".and.parres(ind).eq.1) then
                gpsum52pn=gpsum52pn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum52pn=gmsum52pn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum52pn=g0sum52pn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum52pn=gpwgtsum52pn+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		elseif (spiso(ind).eq."12d".and.parres(ind).eq.1) then
                gpsum12pd=gpsum12pd+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum12pd=gmsum12pd+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum12pd=g0sum12pd+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum12pd=gpwgtsum12pd+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		elseif (spiso(ind).eq."32d".and.parres(ind).eq.1) then
                gpsum32pd=gpsum32pd+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum32pd=gmsum32pd+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum32pd=g0sum32pd+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum32pd=gpwgtsum32pd+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		elseif (spiso(ind).eq."12n".and.parres(ind).eq.-1) then
                gpsum12mn=gpsum12mn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum12mn=gmsum12mn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum12mn=g0sum12mn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum12mn=gpwgtsum12mn+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		elseif (spiso(ind).eq."32n".and.parres(ind).eq.-1) then
                gpsum32mn=gpsum32mn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum32mn=gmsum32mn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum32mn=g0sum32mn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum32mn=gpwgtsum32mn+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		elseif (spiso(ind).eq."52n".and.parres(ind).eq.-1) then
                gpsum52mn=gpsum52mn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum52mn=gmsum52mn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum52mn=g0sum52mn+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum52mn=gpwgtsum52mn+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		elseif (spiso(ind).eq."12d".and.parres(ind).eq.-1) then
                gpsum12md=gpsum12md+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum12md=gmsum12md+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum12md=g0sum12md+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum12md=gpwgtsum12md+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		elseif (spiso(ind).eq."32d".and.parres(ind).eq.-1) then
                gpsum32md=gpsum32md+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
                gmsum32md=gmsum32md+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt2(ind),jres(ind),parres(ind),-1)*flag(ind)
                g0sum32md=g0sum32md+gpm0r(wsq,qsq,gres(ind),mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupl(ind),jres(ind),parres(ind),0)*flag(ind)
                gpwgtsum32md=gpwgtsum32md+gpm0r(wsq,qsq,gres(ind),
     >			mres(ind),
     >		xbc(ind),lres(ind),bfpi(ind),bf2pi(ind),bfeta(ind),
     >		coupt1(ind),jres(ind),parres(ind),1)*flag(ind)
     >			*(-1)**(jres(ind)-0.5d0)*parres(ind)
		endif
	enddo
!======================================================================!
! Note that for all the sums squared the interference between N(1720)
! and N'(1720) is taken into account differently, by reducing their
! interference by a factor 1.72. This is done by removing the 32pn
! sum above and handling it differently below!
!======================================================================!
	gpsumsqr=!gpsum12pn*conjg(gpsum12pn)+
!     >		gpsum32pn*conjg(gpsum32pn)+
     >		gpsum52pn*conjg(gpsum52pn)+gpsum12pd*conjg(gpsum12pd)+
     >		gpsum32pd*conjg(gpsum32pd)+
!     >		gpsum12mn*conjg(gpsum12mn)+
     >		gpsum32mn*conjg(gpsum32mn)+
     >		gpsum52mn*conjg(gpsum52mn)+gpsum12md*conjg(gpsum12md)+
     >		gpsum32md*conjg(gpsum32md)
        gpsumsqr=gpsumsqr+gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupt1(8),jres(8),parres(8),1)*
     >		conjg(gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupt1(8),jres(8),parres(8),1))*flag(8)
     >		+gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupt1(19),jres(19),parres(19),1)*
     >		conjg(gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupt1(19),jres(19),parres(19),1))*flag(19)
     >		+1.d0/1.72d0*gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupt1(8),jres(8),parres(8),1)*
     >		conjg(gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupt1(19),jres(19),parres(19),1))*flag(19)*flag(8)
     >		*interf
     >		+1.d0/1.72d0*conjg(gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupt1(8),jres(8),parres(8),1))*
     >		gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupt1(19),jres(19),parres(19),1)*flag(19)*flag(8)
     >		*interf
        gpsumsqr=gpsumsqr+gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupt1(1),jres(1),parres(1),1)*
     >		conjg(gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupt1(1),jres(1),parres(1),1))*flag(1)
     >		+gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupt1(7),jres(7),parres(7),1)*
     >		conjg(gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupt1(7),jres(7),parres(7),1))*flag(7)
     >		+gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupt1(1),jres(1),parres(1),1)*
     >		conjg(gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupt1(7),jres(7),parres(7),1))*flag(7)*flag(1)
     >		*interf
     >		+conjg(gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupt1(1),jres(1),parres(1),1))*
     >		gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupt1(7),jres(7),parres(7),1)*flag(7)*flag(1)
     >		*interf
        gpsumsqr=gpsumsqr+gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupt1(3),jres(3),parres(3),1)*
     >		conjg(gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupt1(3),jres(3),parres(3),1))*flag(3)
     >		+gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupt1(4),jres(4),parres(4),1)*
     >		conjg(gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupt1(4),jres(4),parres(4),1))*flag(4)
     >		+gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupt1(3),jres(3),parres(3),1)*
     >		conjg(gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupt1(4),jres(4),parres(4),1))*flag(4)*flag(3)
     >		*interf
     >		+conjg(gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupt1(3),jres(3),parres(3),1))*
     >		gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupt1(4),jres(4),parres(4),1)*flag(4)*flag(3)
     >		*interf
	gmsumsqr=!gmsum12pn*conjg(gmsum12pn)+
!     >		gmsum32pn*conjg(gmsum32pn)+
     >		gmsum52pn*conjg(gmsum52pn)+gmsum12pd*conjg(gmsum12pd)+
     >		gmsum32pd*conjg(gmsum32pd)+
!     >		gmsum12mn*conjg(gmsum12mn)+
     >		gmsum32mn*conjg(gmsum32mn)+
     >		gmsum52mn*conjg(gmsum52mn)+gmsum12md*conjg(gmsum12md)+
     >		gmsum32md*conjg(gmsum32md)
        gmsumsqr=gmsumsqr+gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupt2(8),jres(8),parres(8),-1)*
     >		conjg(gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupt2(8),jres(8),parres(8),-1))*flag(8)
     >		+gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupt2(19),jres(19),parres(19),-1)*
     >		conjg(gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupt2(19),jres(19),parres(19),-1))*flag(19)
     >		+1.d0/1.72d0*gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupt2(8),jres(8),parres(8),-1)*
     >		conjg(gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupt2(19),jres(19),parres(19),-1))*flag(19)*flag(8)
     >		*interf
     >		+1.d0/1.72d0*conjg(gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupt2(8),jres(8),parres(8),-1))*
     >		gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupt2(19),jres(19),parres(19),-1)*flag(19)*flag(8)
     >		*interf
        gmsumsqr=gmsumsqr+gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupt2(1),jres(1),parres(1),-1)*
     >		conjg(gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupt2(1),jres(1),parres(1),-1))*flag(1)
     >		+gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupt2(7),jres(7),parres(7),-1)*
     >		conjg(gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupt2(7),jres(7),parres(7),-1))*flag(7)
     >		+gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupt2(1),jres(1),parres(1),-1)*
     >		conjg(gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupt2(7),jres(7),parres(7),-1))*flag(7)*flag(1)
     >		*interf
     >		+conjg(gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupt2(1),jres(1),parres(1),-1))*
     >		gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupt2(7),jres(7),parres(7),-1)*flag(7)*flag(1)
     >		*interf
        gmsumsqr=gmsumsqr+gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupt2(3),jres(3),parres(3),-1)*
     >		conjg(gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupt2(3),jres(3),parres(3),-1))*flag(3)
     >		+gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupt2(4),jres(4),parres(4),-1)*
     >		conjg(gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupt2(4),jres(4),parres(4),-1))*flag(4)
     >		+gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupt2(3),jres(3),parres(3),-1)*
     >		conjg(gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupt2(4),jres(4),parres(4),-1))*flag(4)*flag(3)
     >		*interf
     >		+conjg(gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupt2(3),jres(3),parres(3),-1))*
     >		gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupt2(4),jres(4),parres(4),-1)*flag(4)*flag(3)
     >		*interf
	g0sumsqr=!g0sum12pn*conjg(g0sum12pn)+
!     >		g0sum32pn*conjg(g0sum32pn)+
     >		g0sum52pn*conjg(g0sum52pn)+g0sum12pd*conjg(g0sum12pd)+
     >		g0sum32pd*conjg(g0sum32pd)+
!     >		g0sum12mn*conjg(g0sum12mn)+
     >		g0sum32mn*conjg(g0sum32mn)+
     >		g0sum52mn*conjg(g0sum52mn)+g0sum12md*conjg(g0sum12md)+
     >		g0sum32md*conjg(g0sum32md)
        g0sumsqr=g0sumsqr+gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupl(8),jres(8),parres(8),0)*
     >		conjg(gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupl(8),jres(8),parres(8),0))*flag(8)
     >		+gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupl(19),jres(19),parres(19),0)*
     >		conjg(gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupl(19),jres(19),parres(19),0))*flag(19)
     >		+1.d0/1.72d0*gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupl(8),jres(8),parres(8),0)*
     >		conjg(gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupl(19),jres(19),parres(19),0))*flag(19)*flag(8)
     >		*interf
     >		+1.d0/1.72d0*conjg(gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupl(8),jres(8),parres(8),0))*
     >		gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupl(19),jres(19),parres(19),0)*flag(19)*flag(8)
     >		*interf
        g0sumsqr=g0sumsqr+gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupl(1),jres(1),parres(1),0)*
     >		conjg(gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupl(1),jres(1),parres(1),0))*flag(1)
     >		+gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupl(7),jres(7),parres(7),0)*
     >		conjg(gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupl(7),jres(7),parres(7),0))*flag(7)
     >		+gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupl(1),jres(1),parres(1),0)*
     >		conjg(gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupl(7),jres(7),parres(7),0))*flag(7)*flag(1)
     >		*interf
     >		+conjg(gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupl(1),jres(1),parres(1),0))*
     >		gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupl(7),jres(7),parres(7),0)*flag(7)*flag(1)
     >		*interf
        g0sumsqr=g0sumsqr+gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupl(3),jres(3),parres(3),0)*
     >		conjg(gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupl(3),jres(3),parres(3),0))*flag(3)
     >		+gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupl(4),jres(4),parres(4),0)*
     >		conjg(gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupl(4),jres(4),parres(4),0))*flag(4)
     >		+gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupl(3),jres(3),parres(3),0)*
     >		conjg(gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupl(4),jres(4),parres(4),0))*flag(4)*flag(3)
     >		*interf
     >		+conjg(gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupl(3),jres(3),parres(3),0))*
     >		gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupl(4),jres(4),parres(4),0)*flag(4)*flag(3)
     >		*interf
	gpwgtsumsqr=!gpwgtsum12pn*conjg(g0sum12pn)+
!     >		gpwgtsum32pn*conjg(g0sum32pn)+
     >		gpwgtsum52pn*conjg(g0sum52pn)+
     >		gpwgtsum12pd*conjg(g0sum12pd)+
     >		gpwgtsum32pd*conjg(g0sum32pd)+
!     >		gpwgtsum12mn*conjg(g0sum12mn)+
     >		gpwgtsum32mn*conjg(g0sum32mn)+
     >		gpwgtsum52mn*conjg(g0sum52mn)+
     >		gpwgtsum12md*conjg(g0sum12md)+
     >		gpwgtsum32md*conjg(g0sum32md)
        gpwgtsumsqr=gpwgtsumsqr+gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupt1(8),jres(8),parres(8),1)*
     >		(-1)**(jres(8)-0.5d0)*parres(8)*
     >		conjg(gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupl(8),jres(8),parres(8),0))*flag(8)
     >		+gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupt1(19),jres(19),parres(19),1)*
     >		(-1)**(jres(19)-0.5d0)*parres(19)*
     >		conjg(gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupl(19),jres(19),parres(19),0))*flag(19)
     >		+1.d0/1.72d0*gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupt1(8),jres(8),parres(8),1)*
     >		(-1)**(jres(8)-0.5d0)*parres(8)*
     >		conjg(gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupl(19),jres(19),parres(19),0))*flag(19)*flag(8)
     >		*interf
     >		+1.d0/1.72d0*conjg(gpm0r(wsq,qsq,gres(8),mres(8),
     >		xbc(8),lres(8),bfpi(8),bf2pi(8),bfeta(8),
     >		coupl(8),jres(8),parres(8),0))*
     >		gpm0r(wsq,qsq,gres(19),mres(19),
     >		xbc(19),lres(19),bfpi(19),bf2pi(19),bfeta(19),
     >		coupt1(19),jres(19),parres(19),1)*
     >		(-1)**(jres(19)-0.5d0)*parres(19)*flag(19)*flag(8)
     >		*interf
        gpwgtsumsqr=gpwgtsumsqr+gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupt1(1),jres(1),parres(1),1)*
     >		(-1)**(jres(1)-0.5d0)*parres(1)*
     >		conjg(gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupl(1),jres(1),parres(1),0))*flag(1)
     >		+gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupt1(7),jres(7),parres(7),1)*
     >		(-1)**(jres(7)-0.5d0)*parres(7)*
     >		conjg(gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupl(7),jres(7),parres(7),0))*flag(7)
     >		+gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupt1(1),jres(1),parres(1),1)*
     >		(-1)**(jres(1)-0.5d0)*parres(1)*
     >		conjg(gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupl(7),jres(7),parres(7),0))*flag(7)*flag(1)
     >		*interf
     >		+conjg(gpm0r(wsq,qsq,gres(1),mres(1),
     >		xbc(1),lres(1),bfpi(1),bf2pi(1),bfeta(1),
     >		coupl(1),jres(1),parres(1),0))*
     >		gpm0r(wsq,qsq,gres(7),mres(7),
     >		xbc(7),lres(7),bfpi(7),bf2pi(7),bfeta(7),
     >		coupt1(7),jres(7),parres(7),1)*
     >		(-1)**(jres(7)-0.5d0)*parres(7)*flag(7)*flag(1)
     >		*interf
        gpwgtsumsqr=gpwgtsumsqr+gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupt1(3),jres(3),parres(3),1)*
     >		(-1)**(jres(3)-0.5d0)*parres(3)*
     >		conjg(gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupl(3),jres(3),parres(3),0))*flag(3)
     >		+gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupt1(4),jres(4),parres(4),1)*
     >		(-1)**(jres(4)-0.5d0)*parres(4)*
     >		conjg(gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupl(4),jres(4),parres(4),0))*flag(4)
     >		+gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupt1(3),jres(3),parres(3),1)*
     >		(-1)**(jres(3)-0.5d0)*parres(3)*
     >		conjg(gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupl(4),jres(4),parres(4),0))*flag(4)*flag(3)
     >		*interf
     >		+conjg(gpm0r(wsq,qsq,gres(3),mres(3),
     >		xbc(3),lres(3),bfpi(3),bf2pi(3),bfeta(3),
     >		coupl(3),jres(3),parres(3),0))*
     >		gpm0r(wsq,qsq,gres(4),mres(4),
     >		xbc(4),lres(4),bfpi(4),bf2pi(4),bfeta(4),
     >		coupt1(4),jres(4),parres(4),1)*
     >		(-1)**(jres(4)-0.5d0)*parres(4)*flag(4)*flag(3)
     >		*interf
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
          gpm0rdel=gpm0rdel*p1232!!!! Corrected (-1.d0)**p1232 to p1232 in May 2022!!
        else if (plusminus.eq.0) then
          gpm0rdel=gpm0rdel*p1232/dsqrt(2.d0)!!!! Corrected (-1.d0)**p1232 to p1232 in May 2022!!
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
	endsubroutine listecoup
!======================================================================!
! SUBROUTINE FOR RANDOM GENERATED ELECTROCOUPLINGS
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
