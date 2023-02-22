!     ----------------------------------------
      module numbers
        integer, parameter :: ione = 1, itwo = 2, ithree = 3
        real (8), parameter :: zero = 0.d0
        real (8), parameter :: one = 1.d0, two = 2.d0, three = 3.d0,
     +    four = 4.d0, five = 5.d0, six = 6.d0
      end module numbers
!     ----------------------------------------
!   resonance info table
!    1.N(1440)  2.N(1520)   3.N(1535)  4.N(1650)   5.N(1675)
!    6.N(1680)  7.N(1710)   8.N(1720)  9.N(2190)  10.N(2220)
!   11.N(2250) 12.D(1232)  13.D(1620) 14.D(1700)  15.D(1905)
!   16.D(1910) 17.D(1950)  18.D(2420) 19.N'(1720)
      module couplingsandinterpolation
        character (10), parameter :: namefolder = './resbars/' ! folder name storing the amplitudes
        integer, parameter :: sizearray = 701 ! size of the files storing the amplitudes
        integer, parameter :: nres = 19 ! number of resonances
        integer, parameter :: ncol = 6
        real (8), dimension (sizearray) :: Q2
        real (8), dimension (sizearray,ncol) :: Res1, Res1coef
        real (8), dimension (sizearray,ncol) :: Res2, Res2coef
        real (8), dimension (sizearray,ncol) :: Res3, Res3coef
        real (8), dimension (sizearray,ncol) :: Res4, Res4coef
        real (8), dimension (sizearray,ncol) :: Res5, Res5coef
        real (8), dimension (sizearray,ncol) :: Res6, Res6coef
        real (8), dimension (sizearray,ncol) :: Res7, Res7coef
        real (8), dimension (sizearray,ncol) :: Res8, Res8coef
        real (8), dimension (sizearray,ncol) :: Res9, Res9coef
        real (8), dimension (sizearray,ncol) :: Res10, Res10coef
        real (8), dimension (sizearray,ncol) :: Res11, Res11coef
        real (8), dimension (sizearray,ncol) :: Res12, Res12coef
        real (8), dimension (sizearray,ncol) :: Res13, Res13coef
        real (8), dimension (sizearray,ncol) :: Res14, Res14coef
        real (8), dimension (sizearray,ncol) :: Res15, Res15coef
        real (8), dimension (sizearray,ncol) :: Res16, Res16coef
        real (8), dimension (sizearray,ncol) :: Res17, Res17coef
        real (8), dimension (sizearray,ncol) :: Res18, Res18coef
        real (8), dimension (sizearray,ncol) :: Res19, Res19coef
        integer :: npoints
        save npoints
        save Q2
        save Res1, Res1coef, Res2, Res2coef, Res3, Res3coef
        save Res4, Res4coef, Res5, Res5coef, Res6, Res6coef
        save Res7, Res7coef, Res8, Res8coef, Res9, Res9coef
        save Res10, Res10coef, Res11, Res11coef, Res12, Res12coef
        save Res13, Res13coef, Res14, Res14coef, Res15, Res15coef
        save Res16, Res16coef, Res17, Res17coef, Res18, Res18coef
        save Res19, Res19coef
      end module couplingsandinterpolation!     ----------------------------------------
!     input
!     ires: resonance id numbers
!     s:    q^2 value
!     output
!     a12: A_{1/2} amplitude resampled
!     a32: A_{3/2} amplitude resampled
!     s12: S_{1/2} amplitude resampled
      subroutine mcamplitudes ( ires, s, a12, a32, s12 )
        use couplingsandinterpolation
        implicit none
        ! input variables
        integer  :: ires ! selects the resonance
        real (8) :: s   ! q2 value
        ! internal variables
        real (8), dimension (6) :: x
        real (8) :: mean1, mean2, mean3, error1, error2, error3
        ! output
        real (8) :: a12, a32, s12 ! resampled amplitudes
        ! external functions
        real (8), external :: random_gaussian_number
!       code
        call random_number (x)
        select case (ires)
          case (1)
            call splint (Q2,Res1(1,1),Res1coef(1,1),npoints,s,mean1)
            call splint (Q2,Res1(1,2),Res1coef(1,2),npoints,s,error1)
            call splint (Q2,Res1(1,3),Res1coef(1,3),npoints,s,mean2)
            call splint (Q2,Res1(1,4),Res1coef(1,4),npoints,s,error2)
            call splint (Q2,Res1(1,5),Res1coef(1,5),npoints,s,mean3)
            call splint (Q2,Res1(1,6),Res1coef(1,6),npoints,s,error3)
          case (2)
            call splint (Q2,Res2(1,1),Res2coef(1,1),npoints,s,mean1)
            call splint (Q2,Res2(1,2),Res2coef(1,2),npoints,s,error1)
            call splint (Q2,Res2(1,3),Res2coef(1,3),npoints,s,mean2)
            call splint (Q2,Res2(1,4),Res2coef(1,4),npoints,s,error2)
            call splint (Q2,Res2(1,5),Res2coef(1,5),npoints,s,mean3)
            call splint (Q2,Res2(1,6),Res2coef(1,6),npoints,s,error3)
          case (3)
            call splint (Q2,Res3(1,1),Res3coef(1,1),npoints,s,mean1)
            call splint (Q2,Res3(1,2),Res3coef(1,2),npoints,s,error1)
            call splint (Q2,Res3(1,3),Res3coef(1,3),npoints,s,mean2)
            call splint (Q2,Res3(1,4),Res3coef(1,4),npoints,s,error2)
            call splint (Q2,Res3(1,5),Res3coef(1,5),npoints,s,mean3)
            call splint (Q2,Res3(1,6),Res3coef(1,6),npoints,s,error3)
          case (4)
            call splint (Q2,Res4(1,1),Res4coef(1,1),npoints,s,mean1)
            call splint (Q2,Res4(1,2),Res4coef(1,2),npoints,s,error1)
            call splint (Q2,Res4(1,3),Res4coef(1,3),npoints,s,mean2)
            call splint (Q2,Res4(1,4),Res4coef(1,4),npoints,s,error2)
            call splint (Q2,Res4(1,5),Res4coef(1,5),npoints,s,mean3)
            call splint (Q2,Res4(1,6),Res4coef(1,6),npoints,s,error3)
          case (5)
            call splint (Q2,Res5(1,1),Res5coef(1,1),npoints,s,mean1)
            call splint (Q2,Res5(1,2),Res5coef(1,2),npoints,s,error1)
            call splint (Q2,Res5(1,3),Res5coef(1,3),npoints,s,mean2)
            call splint (Q2,Res5(1,4),Res5coef(1,4),npoints,s,error2)
            call splint (Q2,Res5(1,5),Res5coef(1,5),npoints,s,mean3)
            call splint (Q2,Res5(1,6),Res5coef(1,6),npoints,s,error3)
          case (6)
            call splint (Q2,Res6(1,1),Res6coef(1,1),npoints,s,mean1)
            call splint (Q2,Res6(1,2),Res6coef(1,2),npoints,s,error1)
            call splint (Q2,Res6(1,3),Res6coef(1,3),npoints,s,mean2)
            call splint (Q2,Res6(1,4),Res6coef(1,4),npoints,s,error2)
            call splint (Q2,Res6(1,5),Res6coef(1,5),npoints,s,mean3)
            call splint (Q2,Res6(1,6),Res6coef(1,6),npoints,s,error3)
          case (7)
            call splint (Q2,Res7(1,1),Res7coef(1,1),npoints,s,mean1)
            call splint (Q2,Res7(1,2),Res7coef(1,2),npoints,s,error1)
            call splint (Q2,Res7(1,3),Res7coef(1,3),npoints,s,mean2)
            call splint (Q2,Res7(1,4),Res7coef(1,4),npoints,s,error2)
            call splint (Q2,Res7(1,5),Res7coef(1,5),npoints,s,mean3)
            call splint (Q2,Res7(1,6),Res7coef(1,6),npoints,s,error3)
          case (8)
            call splint (Q2,Res8(1,1),Res8coef(1,1),npoints,s,mean1)
            call splint (Q2,Res8(1,2),Res8coef(1,2),npoints,s,error1)
            call splint (Q2,Res8(1,3),Res8coef(1,3),npoints,s,mean2)
            call splint (Q2,Res8(1,4),Res8coef(1,4),npoints,s,error2)
            call splint (Q2,Res8(1,5),Res8coef(1,5),npoints,s,mean3)
            call splint (Q2,Res8(1,6),Res8coef(1,6),npoints,s,error3)
          case (9)
            call splint (Q2,Res9(1,1),Res9coef(1,1),npoints,s,mean1)
            call splint (Q2,Res9(1,2),Res9coef(1,2),npoints,s,error1)
            call splint (Q2,Res9(1,3),Res9coef(1,3),npoints,s,mean2)
            call splint (Q2,Res9(1,4),Res9coef(1,4),npoints,s,error2)
            call splint (Q2,Res9(1,5),Res9coef(1,5),npoints,s,mean3)
            call splint (Q2,Res9(1,6),Res9coef(1,6),npoints,s,error3)
          case (10)
            call splint (Q2,Res10(1,1),Res10coef(1,1),npoints,s,mean1)
            call splint (Q2,Res10(1,2),Res10coef(1,2),npoints,s,error1)
            call splint (Q2,Res10(1,3),Res10coef(1,3),npoints,s,mean2)
            call splint (Q2,Res10(1,4),Res10coef(1,4),npoints,s,error2)
            call splint (Q2,Res10(1,5),Res10coef(1,5),npoints,s,mean3)
            call splint (Q2,Res10(1,6),Res10coef(1,6),npoints,s,error3)
          case (11)
            call splint (Q2,Res11(1,1),Res11coef(1,1),npoints,s,mean1)
            call splint (Q2,Res11(1,2),Res11coef(1,2),npoints,s,error1)
            call splint (Q2,Res11(1,3),Res11coef(1,3),npoints,s,mean2)
            call splint (Q2,Res11(1,4),Res11coef(1,4),npoints,s,error2)
            call splint (Q2,Res11(1,5),Res11coef(1,5),npoints,s,mean3)
            call splint (Q2,Res11(1,6),Res11coef(1,6),npoints,s,error3)
          case (12)
            call splint (Q2,Res12(1,1),Res12coef(1,1),npoints,s,mean1)
            call splint (Q2,Res12(1,2),Res12coef(1,2),npoints,s,error1)
            call splint (Q2,Res12(1,3),Res12coef(1,3),npoints,s,mean2)
            call splint (Q2,Res12(1,4),Res12coef(1,4),npoints,s,error2)
            call splint (Q2,Res12(1,5),Res12coef(1,5),npoints,s,mean3)
            call splint (Q2,Res12(1,6),Res12coef(1,6),npoints,s,error3)
          case (13)
            call splint (Q2,Res13(1,1),Res13coef(1,1),npoints,s,mean1)
            call splint (Q2,Res13(1,2),Res13coef(1,2),npoints,s,error1)
            call splint (Q2,Res13(1,3),Res13coef(1,3),npoints,s,mean2)
            call splint (Q2,Res13(1,4),Res13coef(1,4),npoints,s,error2)
            call splint (Q2,Res13(1,5),Res13coef(1,5),npoints,s,mean3)
            call splint (Q2,Res13(1,6),Res13coef(1,6),npoints,s,error3)
          case (14)
            call splint (Q2,Res14(1,1),Res14coef(1,1),npoints,s,mean1)
            call splint (Q2,Res14(1,2),Res14coef(1,2),npoints,s,error1)
            call splint (Q2,Res14(1,3),Res14coef(1,3),npoints,s,mean2)
            call splint (Q2,Res14(1,4),Res14coef(1,4),npoints,s,error2)
            call splint (Q2,Res14(1,5),Res14coef(1,5),npoints,s,mean3)
            call splint (Q2,Res14(1,6),Res14coef(1,6),npoints,s,error3)
          case (15)
            call splint (Q2,Res15(1,1),Res15coef(1,1),npoints,s,mean1)
            call splint (Q2,Res15(1,2),Res15coef(1,2),npoints,s,error1)
            call splint (Q2,Res15(1,3),Res15coef(1,3),npoints,s,mean2)
            call splint (Q2,Res15(1,4),Res15coef(1,4),npoints,s,error2)
            call splint (Q2,Res15(1,5),Res15coef(1,5),npoints,s,mean3)
            call splint (Q2,Res15(1,6),Res15coef(1,6),npoints,s,error3)
          case (16)
            call splint (Q2,Res16(1,1),Res16coef(1,1),npoints,s,mean1)
            call splint (Q2,Res16(1,2),Res16coef(1,2),npoints,s,error1)
            call splint (Q2,Res16(1,3),Res16coef(1,3),npoints,s,mean2)
            call splint (Q2,Res16(1,4),Res16coef(1,4),npoints,s,error2)
            call splint (Q2,Res16(1,5),Res16coef(1,5),npoints,s,mean3)
            call splint (Q2,Res16(1,6),Res16coef(1,6),npoints,s,error3)
          case (17)
            call splint (Q2,Res17(1,1),Res17coef(1,1),npoints,s,mean1)
            call splint (Q2,Res17(1,2),Res17coef(1,2),npoints,s,error1)
            call splint (Q2,Res17(1,3),Res17coef(1,3),npoints,s,mean2)
            call splint (Q2,Res17(1,4),Res17coef(1,4),npoints,s,error2)
            call splint (Q2,Res17(1,5),Res17coef(1,5),npoints,s,mean3)
            call splint (Q2,Res17(1,6),Res17coef(1,6),npoints,s,error3)
          case (18)
            call splint (Q2,Res18(1,1),Res18coef(1,1),npoints,s,mean1)
            call splint (Q2,Res18(1,2),Res18coef(1,2),npoints,s,error1)
            call splint (Q2,Res18(1,3),Res18coef(1,3),npoints,s,mean2)
            call splint (Q2,Res18(1,4),Res18coef(1,4),npoints,s,error2)
            call splint (Q2,Res18(1,5),Res18coef(1,5),npoints,s,mean3)
            call splint (Q2,Res18(1,6),Res18coef(1,6),npoints,s,error3)
          case (19)
            call splint (Q2,Res19(1,1),Res19coef(1,1),npoints,s,mean1)
            call splint (Q2,Res19(1,2),Res19coef(1,2),npoints,s,error1)
            call splint (Q2,Res19(1,3),Res19coef(1,3),npoints,s,mean2)
            call splint (Q2,Res19(1,4),Res19coef(1,4),npoints,s,error2)
            call splint (Q2,Res19(1,5),Res19coef(1,5),npoints,s,mean3)
            call splint (Q2,Res19(1,6),Res19coef(1,6),npoints,s,error3)
          case default
            write (*,*)'Please select a resonance in the 1-19 range'
        end select
        a12 = random_gaussian_number ( mean1, error1, x(1), x(4) )
        a32 = random_gaussian_number ( mean2, error2, x(2), x(5) )
        s12 = random_gaussian_number ( mean3, error3, x(3), x(6) )
        return
      end subroutine mcamplitudes
!     ----------------------------------------
      subroutine readcouplings
        use numbers
        use couplingsandinterpolation
        implicit none
        ! internal variables
        character (25), dimension (nres) :: nameresfile
        integer :: k, i, j, nfile
        real (8) :: rvoid, cero
        nameresfile (1)  = namefolder // 'Res1.txt'
        nameresfile (2)  = namefolder // 'Res2.txt'
        nameresfile (3)  = namefolder // 'Res3.txt'
        nameresfile (4)  = namefolder // 'Res4.txt'
        nameresfile (5)  = namefolder // 'Res5.txt'
        nameresfile (6)  = namefolder // 'Res6.txt'
        nameresfile (7)  = namefolder // 'Res7.txt'
        nameresfile (8)  = namefolder // 'Res8.txt'
        nameresfile (9)  = namefolder // 'Res9.txt'
        nameresfile (10) = namefolder // 'Res10.txt'
        nameresfile (11) = namefolder // 'Res11.txt'
        nameresfile (12) = namefolder // 'Res12.txt'
        nameresfile (13) = namefolder // 'Res13.txt'
        nameresfile (14) = namefolder // 'Res14.txt'
        nameresfile (15) = namefolder // 'Res15.txt'
        nameresfile (16) = namefolder // 'Res16.txt'
        nameresfile (17) = namefolder // 'Res17.txt'
        nameresfile (18) = namefolder // 'Res18.txt'
        nameresfile (19) = namefolder // 'Res19.txt'
        nfile = 100
        do k = 1, nres
          open (nfile+k,file = nameresfile(k), status='old')
        enddo
        nfile = 100
        do i = 1, sizearray
          read (nfile+1,*)  Q2(i), ( Res1(i,j),j=1,6 )
          read (nfile+2,*)  rvoid, ( Res2(i,j),j=1,6 )
          read (nfile+3,*)  rvoid, ( Res3(i,j),j=1,6 )
          read (nfile+4,*)  rvoid, ( Res4(i,j),j=1,6 )
          read (nfile+5,*)  rvoid, ( Res5(i,j),j=1,6 )
          read (nfile+6,*)  rvoid, ( Res6(i,j),j=1,6 )
          read (nfile+7,*)  rvoid, ( Res7(i,j),j=1,6 )
          read (nfile+8,*)  rvoid, ( Res8(i,j),j=1,6 )
          read (nfile+9,*)  rvoid, ( Res9(i,j),j=1,6 )
          read (nfile+10,*) rvoid, ( Res10(i,j),j=1,6)
          read (nfile+11,*) rvoid, ( Res11(i,j),j=1,6)
          read (nfile+12,*) rvoid, ( Res12(i,j),j=1,6)
          read (nfile+13,*) rvoid, ( Res13(i,j),j=1,6)
          read (nfile+14,*) rvoid, ( Res14(i,j),j=1,6)
          read (nfile+15,*) rvoid, ( Res15(i,j),j=1,6)
          read (nfile+16,*) rvoid, ( Res16(i,j),j=1,6)
          read (nfile+17,*) rvoid, ( Res17(i,j),j=1,6)
          read (nfile+18,*) rvoid, ( Res18(i,j),j=1,6)
          read (nfile+19,*) rvoid, ( Res19(i,j),j=1,6)
        enddo
        npoints = sizearray
        cero = zero
!       interpolation
        do j = 1, 6
          call spline (Q2,Res1(1,j),npoints,cero,cero,Res1coef(1,j))
          call spline (Q2,Res2(1,j),npoints,cero,cero,Res2coef(1,j))
          call spline (Q2,Res3(1,j),npoints,cero,cero,Res3coef(1,j))
          call spline (Q2,Res4(1,j),npoints,cero,cero,Res4coef(1,j))
          call spline (Q2,Res5(1,j),npoints,cero,cero,Res5coef(1,j))
          call spline (Q2,Res6(1,j),npoints,cero,cero,Res6coef(1,j))
          call spline (Q2,Res7(1,j),npoints,cero,cero,Res7coef(1,j))
          call spline (Q2,Res8(1,j),npoints,cero,cero,Res8coef(1,j))
          call spline (Q2,Res9(1,j),npoints,cero,cero,Res9coef(1,j))
          call spline (Q2,Res10(1,j),npoints,cero,cero,Res10coef(1,j))
          call spline (Q2,Res11(1,j),npoints,cero,cero,Res11coef(1,j))
          call spline (Q2,Res12(1,j),npoints,cero,cero,Res12coef(1,j))
          call spline (Q2,Res13(1,j),npoints,cero,cero,Res13coef(1,j))
          call spline (Q2,Res14(1,j),npoints,cero,cero,Res14coef(1,j))
          call spline (Q2,Res15(1,j),npoints,cero,cero,Res15coef(1,j))
          call spline (Q2,Res16(1,j),npoints,cero,cero,Res16coef(1,j))
          call spline (Q2,Res17(1,j),npoints,cero,cero,Res17coef(1,j))
          call spline (Q2,Res18(1,j),npoints,cero,cero,Res18coef(1,j))
          call spline (Q2,Res19(1,j),npoints,cero,cero,Res19coef(1,j))
        enddo
        do i = 1, nres
          close (nfile+i)
        enddo
        return
      end subroutine readcouplings
!     -------------------------------------------------------------------------
!     Random numbers routines
!     -------------------------------------------------------------------------
      subroutine init_random_seed2 ()
        use iso_fortran_env, only: int64
        implicit none
        integer, allocatable :: seed(:)
        integer :: i, n, un, istat, dt(8), pid
        integer(int64) :: t
!     code
        call random_seed (size = n)
        allocate(seed(n))
!     First try if the OS provides a random number generator
        open(newunit=un, file="/dev/urandom", access="stream",
     +     form="unformatted", action="read",
     +     status="old", iostat=istat)
        if (istat == 0) then
          if (this_image()==1) read  (un) seed
!print *,"OS provides random number generator"
          close (un)
        else
          if (this_image()==1)
     +          print *,"OS does not provide random number generator"
!     Fallback to XOR:ing the current time and pid. The PID is
!     useful in case one launches multiple instances of the same
!     program in parallel.
          call system_clock(t)
          if (t == 0) then
              call date_and_time(values=dt)
              t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000
     +           + dt(2) * 31_int64 * 24 * 60 * 60 * 1000
     +           + dt(3) * 24_int64 * 60 * 60 * 1000
     +           + dt(5) * 60 * 60 * 1000
     +           + dt(6) * 60 * 1000 + dt(7) * 1000
     +           + dt(8)
          endif
          pid = getpid()
          t   = ieor(t, int(pid, kind(t)))
          do i = 1, n
              seed(i) = lcg(t)
            end do
          end if
          call random_seed (put=seed)
        contains
        function lcg (s)
          integer :: lcg
          integer(int64) :: s
          if (s == 0) then
            s = 104729
          else
            s = mod(s, 4294967296_int64)
          end if
          s = mod(s * 279470273_int64, 4294967291_int64)
          lcg = int(mod(s, int(huge(0), int64)), kind(0))
        end function lcg
      end subroutine init_random_seed2
!     -------------------------------------------------------------------------
      real(8) function random_gaussian_number(mean,sigma,x,y)
        implicit none
        real(8) :: mean,sigma,x,y,z,rho,theta
        real(8), parameter :: pi = 4.d0*datan(1.d0)
        rho = -log(x)
        theta = 2*pi*y
        z = sqrt(2*rho)*cos(theta)
        random_gaussian_number = sigma*z + mean
        return
      end function random_gaussian_number
!     ------------------------------------------------------------------
!     interpolating routines
!     ------------------------------------------------------------------
      SUBROUTINE SPLINE (X,Y,N,YP1,YPN,Y2)
        use numbers
        IMPLICIT NONE
        INTEGER, PARAMETER :: NPUNT1 = 701 ! npoints
        INTEGER   :: I,N,K
        REAL (8) :: X(NPUNT1),Y(NPUNT1),Y2(NPUNT1),U(NPUNT1)
        REAL (8) :: YP1,YPN,SIG,P,QN,UN
        IF (YP1.GT..99E30) THEN
          Y2(1) = zero
          U(1)  = zero
        ELSE
          Y2(1) = -one/two
          U(1)  = (three/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
        ENDIF
        DO 11 I = 2,N-1
          SIG   = (X(I)-X(I-1))/(X(I+1)-X(I-1))
          P     = SIG*Y2(I-1)+two
          Y2(I) = (SIG-one)/P
          U(I)  = (six*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     +      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
  11    CONTINUE
        IF (YPN.GT..99d30) THEN
          QN = zero
          UN = zero
        ELSE
          QN = one/two
          UN = (three/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
        ENDIF
        Y2(N) = (UN-QN*U(N-1))/(QN*Y2(N-1)+one)
        DO 12 K = N-1,1,-1
          Y2(K) = Y2(K)*Y2(K+1)+U(K)
  12    CONTINUE
        RETURN
      END
!     ------------------------------------------------------------------
      SUBROUTINE SPLINT (XA,YA,Y2A,N,X,Y)
        use numbers
        IMPLICIT NONE
        INTEGER, PARAMETER :: NPUNT1 = 701 ! npoints
        INTEGER  :: N,KLO,KHI,K
        REAL (8) :: XA(NPUNT1),YA(NPUNT1),Y2A(NPUNT1)
        REAL (8) :: A,B,Y,H,X
        KLO = 1
        KHI = N
 1      IF (KHI-KLO.GT.1) THEN
          K = ( KHI+KLO )/2
          IF (XA(K).GT.X) THEN
            KHI = K
          ELSE
            KLO = K
          ENDIF
          GOTO 1
        ENDIF
        H = XA(KHI) - XA(KLO)
        IF (H.EQ.zero) THEN
          PRINT*,KHI,KLO,N,X,XA(KHI),XA(KLO),XA(1)
          WRITE (*,*)'Bad XA input.'
          READ (*,*)
          DO 100 K = 1,KHI
            PRINT*,K,XA(K)
100         CONTINUE
        ENDIF
        A = ( XA(KHI)-X )/H
        B = ( X-XA(KLO) )/H
        Y = A*YA(KLO) + B*YA(KHI) +
     +   ( ( A**3-A )*Y2A(KLO)+( B**3-B )*Y2A(KHI) )*(H**2)/six
        RETURN
      END
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
