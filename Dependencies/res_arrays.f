!======================================================================!
! MODULE THAT DECLARES THE RESONANCE ARRAYS: 11 N*'S AND 8 DELTAS
!======================================================================!
	module res_decl
	implicit none
	real(8), dimension (19) :: mres
	real(8), dimension (19) :: gres
	integer, dimension (19) :: lres
	real(8), dimension (19) :: jres
	integer, dimension (19) :: parres
	real(8), dimension (19) :: bfpi
	real(8), dimension (19) :: bf2pi
	real(8), dimension (19) :: bfeta
	real(8), dimension (19) :: xbc
	character (3), dimension (19) :: spiso
	end module res_decl
!======================================================================!
! SUBROUTINE THAT EXPLICITLY DECLARES ARRAY ELEMENTS
!======================================================================!
	subroutine res_params
	use res_decl
	mres(1) = 1.430
	mres(2) = 1.515
	mres(3) = 1.535
	mres(4) = 1.655
	mres(5) = 1.675
	mres(6) = 1.685
	mres(7) = 1.710
	mres(8) = 1.748
	mres(9) = 2.190
	mres(10) = 2.250
	mres(11) = 2.280
	mres(12) = 1.232
	mres(13) = 1.630
	mres(14) = 1.700
	mres(15) = 1.880
	mres(16) = 1.890
	mres(17) = 1.930
	mres(18) = 2.420
	mres(19) = 1.725
	spiso(1) = "12n"
	spiso(2) = "32n"
	spiso(3) = "12n"
	spiso(4) = "12n"
	spiso(5) = "52n"
	spiso(6) = "52n"
	spiso(7) = "12n"
	spiso(8) = "32n"
	spiso(9) = "000"
	spiso(10) = "000"
	spiso(11) = "000"
	spiso(12) = "32d"
	spiso(13) = "12d"
	spiso(14) = "32d"
	spiso(15) = "000"
	spiso(16) = "000"
	spiso(17) = "000"
	spiso(18) = "000"
	spiso(19) = "32n"
	gres(1) = 0.350
	gres(2) = 0.115
	gres(3) = 0.150
	gres(4) = 0.140
	gres(5) = 0.150
	gres(6) = 0.130
	gres(7) = 0.100
	gres(8) = 0.114
	gres(9) = 0.500
	gres(10) = 0.400
	gres(11) = 0.500
	gres(12) = 0.117
	gres(13) = 0.140
	gres(14) = 0.293
	gres(15) = 0.330
	gres(16) = 0.280
	gres(17) = 0.285
	gres(18) = 0.400
	gres(19) = 0.120
	lres(1) = 1
	lres(2) = 2
	lres(3) = 0
	lres(4) = 0
	lres(5) = 2
	lres(6) = 3
	lres(7) = 1
	lres(8) = 1
	lres(9) = 4
	lres(10) = 5
	lres(11) = 4
	lres(12) = 1
	lres(13) = 0
	lres(14) = 2
	lres(15) = 3
	lres(16) = 1
	lres(17) = 3
	lres(18) = 5
	lres(19) = 1
	parres(1) = 1
	parres(2) = -1
	parres(3) = -1
	parres(4) = -1
	parres(5) = -1
	parres(6) = 1
	parres(7) = 1
	parres(8) = 1
	parres(9) = 0
	parres(10) = 0
	parres(11) = 0
	parres(12) = 1
	parres(13) = 1
	parres(14) = -1
	parres(15) = 0
	parres(16) = 0
	parres(17) = 0
	parres(18) = 0
	parres(19) = 1
	jres(1) = 0.5
	jres(2) = 1.5
	jres(3) = 0.5
	jres(4) = 0.5
	jres(5) = 2.5
	jres(6) = 2.5
	jres(7) = 0.5
	jres(8) = 1.5
	jres(9) = 3.5
	jres(10) = 4.5
	jres(11) = 4.5
	jres(12) = 1.5
	jres(13) = 0.5
	jres(14) = 1.5
	jres(15) = 2.5
	jres(16) = 0.5
	jres(17) = 3.5
	jres(18) = 5.5
	jres(19) = 1.5
	bfpi(1) = 0.65
	bfpi(2) = 0.6
	bfpi(3) = 0.45
	bfpi(4) = 0.6
	bfpi(5) = 0.4
	bfpi(6) = 0.68
	bfpi(7) = 0.13
	bfpi(8) = 0.14
	bfpi(9) = 0.15
	bfpi(10) = 0.2
	bfpi(11) = 0.1
	bfpi(12) = 1.
	bfpi(13) = 0.25
	bfpi(14) = 0.10
	bfpi(15) = 0.12
	bfpi(16) = 0.23
	bfpi(17) = 0.4
	bfpi(18) = 0.1
	bfpi(19) = 0.38
	bf2pi(1) = 0.35
	bf2pi(2) = 0.4
	bf2pi(3) = 0.13
	bf2pi(4) = 0.22
	bf2pi(5) = 0.6
	bf2pi(6) = 0.32
	bf2pi(7) = 0.57
	bf2pi(8) = 0.82
	bf2pi(9) = 0.85
	bf2pi(10) = 0.8
	bf2pi(11) = 0.9
	bf2pi(12) = 0.
	bf2pi(13) = 0.75
	bf2pi(14) = 0.9
	bf2pi(15) = 0.88
	bf2pi(16) = 0.77
	bf2pi(17) = 0.6
	bf2pi(18) = 0.9
	bf2pi(19) = 0.62
	bfeta(1) = 0.
	bfeta(2) = 0.
	bfeta(3) = 0.42
	bfeta(4) = 0.18
	bfeta(5) = 0.
	bfeta(6) = 0.
	bfeta(7) = 0.3
	bfeta(8) = 0.04
	bfeta(9) = 0.
	bfeta(10) = 0.
	bfeta(11) = 0.
	bfeta(12) = 0.
	bfeta(13) = 0.
	bfeta(14) = 0.
	bfeta(15) = 0.
	bfeta(16) = 0.
	bfeta(17) = 0.
	bfeta(18) = 0.
	bfeta(19) = 0.
	xbc(1) = 0.3
	xbc(2) = 0.1
	xbc(3) = 0.5
	xbc(4) = 0.5
	xbc(5) = 0.5
	xbc(6) = 0.2
	xbc(7) = 0.5
	xbc(8) = 0.5
	xbc(9) = 0.5
	xbc(10) = 0.5
	xbc(11) = 0.5
!======================================================================!
! For the xbc values we use the results from Aznauryan's model; the
! Delta(1232) resonance is an exception treated differently, so no xbc
! value is needed.
!======================================================================!
	xbc(13) = 0.5
	xbc(14) = 0.22
	xbc(15) = 0.5
	xbc(16) = 0.5
	xbc(17) = 0.5
	xbc(18) = 0.5
	xbc(19) = 0.5
	end subroutine res_params
