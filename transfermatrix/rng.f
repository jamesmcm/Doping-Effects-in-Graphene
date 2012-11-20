c$$$      PROGRAM RAND2
c$$$      IMPLICIT NONE
c$$$      REAL GETRAND
c$$$      REAL R
c$$$
c$$$      CALL seedtime()
c$$$
c$$$      R=GETRAND()
c$$$      WRITE(*,310) R
c$$$
c$$$ 310  FORMAT(F20.10)
c$$$
c$$$      STOP
c$$$      END
c$$$NOTE PROGRAM MUST BE COMPILED WITH TESTTIME.C OBJECT
c$$$gcc -c gettime.c
c$$$gfortran -c rng.f
c$$$gfortran -o program <allobjectfileshere.o>

      REAL FUNCTION GETRAND(SEED)
      IMPLICIT NONE
      INTEGER*4 SEED
      REAL ran2

C     SEED must be the same as the seed. To get a new sequence recall safeseed with the new seed, then repeatedly call GETRAND with num of the new seed.
      IF (SEED .EQ. 0) THEN
         WRITE (*,*) "Warning: Seed was set to zero."
      ELSE IF (SEED .GT. 0) THEN
         SEED=-1*SEED
      END IF
   
      GETRAND=ran2(SEED)
      RETURN
      END

      SUBROUTINE SEEDTIME()
      IMPLICIT NONE
      INTEGER TIME
      INTEGER gettime1
      REAL ran2
      REAL Z

      TIME=-1*gettime1()
      Z=ran2(TIME)

      RETURN
      END

      SUBROUTINE SAFESEED(SEED)
      IMPLICIT NONE
      INTEGER SEED
      REAL ran2
      REAL Z
      
      IF (SEED .EQ. 0) THEN
         WRITE (*,*) "Warning: Seed was set to zero."
      ELSE IF (SEED .GT. 0) THEN
         SEED=-1*SEED
      END IF
      
      Z=ran2(SEED)
         
      RETURN
      END
      
      REAL FUNCTION ran2(idum)
      IMPLICIT NONE
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,
     +     AM=1./IM1,IMM1=IM1-1, IA1=40014,IA2=40692,IQ1=53668,
     +     IQ2=52774,IR1=12211, IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,
     +     EPS=1.2e-7,RNMX=1.-EPS)

      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/


      if (idum.le.0) then 
         idum=max(-idum,1)
         idum2=idum
         do j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
         enddo
         iy=iv(1)
      endif
 
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1

      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2 

      if (idum2.lt.0) idum2=idum2+IM2 
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
