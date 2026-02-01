      !      PRINT *,'DUAMEL'
      !      PRINT *, UN1,UT,DT,Q(1),NU1,N,M,K,NTAU
      !      PRINT *,'DUAMEL?'
 
      ! SUBROUTINE DUAMEL(Q,UN1,UT,DT,N,MM,K,NTAU,QB)
      ! AWW adding documentation for input variables
      ! Q      : TCI (total channel inflow) vector ... .ie unrouted
      ! U1     : unit hydrograph vector
      ! UN1    : unit hydrograph shape parameter (for gamma dist)
      ! UT     : unit hydrograph scale parameter (for gamma dist)
      ! DT     : timestep of the UH function (in days or fractions thereof)
      ! N      : sim_length + uh_length ...ie length of U1
      ! M      : max UH length?
      ! QB     : routed flow vector
      ! K      : 
      ! NTAU   : 

      ! AWW adding documentation for local variables
      ! SP     : 

C
C===========================================================
C
C     THIS SUBROUTINE PERFORM UNIT HYDROGRAPH ROUTING
C
      SUBROUTINE DUAMEL_SYNC_UH(Q,U1,DT,N,MM,NTAU,QB)
      IMPLICIT REAL (A-H,O-Z)
      INTEGER A,B,M
      
      REAL, DIMENSION(N-MM), INTENT(IN) :: Q
      REAL, INTENT(IN) :: DT
      INTEGER, INTENT(IN) :: N, MM, NTAU
      REAL, DIMENSION(N), INTENT(OUT) :: QB

      REAL, DIMENSION(MM), INTENT(IN) :: U1

!      write(*,*) DT,N,MM,NTAU
!      U1 = 0

      M = MM

!      do L=1,M
!        write(*,*)'U1',L,U1(L)
!      end do 
!    6 CONTINUE
      IOC=N+NTAU
      IF(N.GT.M)GO TO 10
      DO 2 I=1,IOC
      QB(I)=0.
      A=1
      IF(I.GT.M)A=I-M+1
      B=I
      IF(I.GT.N)B=N
      DO 3 J=A,B
      IOR=I-J+1
!      if(I.lt.10)write(*,*)QB(I)+U1(J)*Q(IOR)
      QB(I)=QB(I)+Q(J)*U1(IOR)
    3 CONTINUE
    2 CONTINUE
      GO TO 11
   10 DO 4 I=1,IOC
      QB(I)=0.
      A=1
      IF(I.GT.N)A=I-N+1
      B=I
      IF(I.GT.M)B=M
      DO 5 J=A,B
      IOR=I-J+1
      !if(I.lt.10)write(*,*)QB(I)+U1(J)*Q(IOR)
      QB(I)=QB(I)+U1(J)*Q(IOR)
 5    CONTINUE
 4    CONTINUE
 11   RETURN
      END
