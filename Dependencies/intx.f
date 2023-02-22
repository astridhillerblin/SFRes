!=====================================================================!
! This file includes ready-made integration routines used to calculate
! truncated moments in resonance_files.f
!=====================================================================!
      SUBROUTINE DSG20R(A,B,N,X,NP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Y(10),X(2000)
      DATA Y/.9931285991d0,.9639719272d0,.9122344282d0,.8391169718d0,
     F .7463319064d0,.6360536807d0,.5108670019d0,.3737060887d0,
     F .2277858511d0,.0765265211d0/
      NP=20*N
      XN= N
      DINT=B - A
      DINT = DINT/XN
      DELT=DINT*0.5D0
      ORIG=A-DELT
      I1=-20
      DO 1 I=1,N
      ORIG=ORIG+DINT
      DORIG=ORIG+ORIG
      I1=I1+20
      I2=I1+21
      DO 2 J=1,10
      J1=I1+J
      J2=I2-J
      X(J1)=ORIG-DELT*Y(J)
 2    X(J2)=DORIG-X(J1)
 1    CONTINUE
      RETURN
      END
*********************************************************************
*********************************************************************
*********************************************************************
      SUBROUTINE DRG20R(A,B,N,CF,CRES)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION W(10),CF(2000)
      DATA W/.0176140071d0,.0406014298d0,.0626720483d0,.0832767415d0,
     F .1019301198d0,.1181945319d0,.1316886384d0,.1420961093d0,
     f .1491729864d0,.1527533871d0/
      cr=0.d0
      I1=-20
      DO 1 I=1,N
      I1=I1+20
      I2=I1+21
      DO 2 J=1,10
      J1=I1+J
      J2=I2-J
 2    CR=CR+W(J)*(CF(J1)+CF(J2))
 1    CONTINUE
      CRES=CR*0.5D0*(B-A)/DBLE(N)
      RETURN
      END    
