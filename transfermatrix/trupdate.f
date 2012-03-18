c$$   NEW SUBROUTINE TO GENERATE A,B,C,D
      SUBROUTINE GENABCD(LIMX, MULT, A, B, C, D,U)
      IMPLICIT NONE
      INTEGER LIMX
      DOUBLE COMPLEX UNITY, ZERO,NEG,HALF,NHALF, 
     + A(LIMX, LIMX), B(LIMX, LIMX),C(LIMX, LIMX), D(LIMX, LIMX),


     + MULT(2*LIMX, 2*LIMX),
     
C    M(A,B,C,D) COMPONENTS OF MATRIX, U: PHASE MATRIX
     + MA(LIMX,LIMX),MB(LIMX,LIMX),MC(LIMX,LIMX),MD(lIMX,LIMX),
     + U(LIMX,LIMX),TEMPX(LIMX,LIMX),TEMPY(LIMX,LIMX),
     + LAMONE(LIMX,LIMX),LAMTWO(LIMX,LIMX)
     
      UNITY = 1.0
      ZERO  = 0.0
      NEG   =-1.0
      HALF  = 0.5
      NHALF = -0.5

C     
      MA = MULT(1:LIMX, 1:LIMX)
      MB = MULT((LIMX+1):2*LIMX, 1:LIMX)
      MC = MULT(1:LIMX, LIMX+1:2*LIMX)
      MD = MULT((LIMX+1):2*LIMX, (LIMX+1):2*LIMX)

c     OLD INTITISE
c      MA = MULT(1:LIMX, 1:LIMX)
c      MB = MULT((LIMX+1):2*LIMX, 1:LIMX)
c      MC = MULT(1:LIMX, LIMX+1:2*LIMX)
c      MD = MULT((LIMX+1):2*LIMX, (LIMX+1):2*LIMX)

c$$$ M(A,B,C,D) INITALISE TO ZERO'S
      CALL SQZERO(A,LIMX)
      CALL SQZERO(B,LIMX)
      CALL SQZERO(C,LIMX)
      CALL SQZERO(D,LIMX)

C     U NOW INTIALISED IN FILLU
C      CALL SQUNITZ(U, ZI, LIMX)


      CALL SQDOT(TEMPX,MB,U,LIMX)
      CALL SQDOT(LAMONE,MD,U,LIMX)

      CALL ZGEMM('C','N',LIMX,LIMX,LIMX,UNITY,U,LIMX,MC,LIMX,UNITY,
     + TEMPY,LIMX)
      CALL ZGEMM('C','N',LIMX,LIMX,LIMX,UNITY,U,LIMX,LAMONE,LIMX,UNITY,
     + LAMTWO,LIMX)

c$$$ CREATED TEMP MATICIES -> GENERATE A
      CALL ZAXPY(LIMX*LIMX,HALF,MA,1,A,1)
      CALL ZAXPY(LIMX*LIMX,HALF,TEMPX,1,A,1)
      CALL ZAXPY(LIMX*LIMX,HALF,TEMPY,1,A,1)
      CALL ZAXPY(LIMX*LIMX,HALF,LAMTWO,1,A,1)
c$$$ GENERATE B
      CALL ZAXPY(LIMX*LIMX,HALF,MA,1,B,1)
      CALL ZAXPY(LIMX*LIMX,NHALF,TEMPX,1,B,1)
      CALL ZAXPY(LIMX*LIMX,HALF,TEMPY,1,B,1)
      CALL ZAXPY(LIMX*LIMX,NHALF,LAMTWO,1,B,1)
c$$$  GENERTATE C
      CALL ZAXPY(LIMX*LIMX,HALF,MA,1,C,1)
      CALL ZAXPY(LIMX*LIMX,HALF,TEMPX,1,C,1)
      CALL ZAXPY(LIMX*LIMX,NHALF,TEMPY,1,C,1)
      CALL ZAXPY(LIMX*LIMX,NHALF,LAMTWO,1,C,1)      
c$$$  GENERTATE D
      CALL ZAXPY(LIMX*LIMX,HALF,MA,1,D,1)
      CALL ZAXPY(LIMX*LIMX,NHALF,TEMPX,1,D,1)
      CALL ZAXPY(LIMX*LIMX,NHALF,TEMPY,1,D,1)
      CALL ZAXPY(LIMX*LIMX,HALF,LAMTWO,1,D,1)      


C$$$ I HAVE VERIFIED THAT AD-BC=1 (IDENTITY MATRIX) AS EXPECTED
      RETURN
      END
 
C$$$ ROUTINE TO GENERATE T AND R
      SUBROUTINE GENTANDRINC(LIMX, TINC,RINC,TTILDEINC,RTILDEINC,A,B,
     + C,D)
      IMPLICIT NONE
      INTEGER LIMX
      DOUBLE COMPLEX UNITY, ZERO
      DOUBLE COMPLEX A(LIMX, LIMX), B(LIMX, LIMX),
     +               C(LIMX, LIMX), D(LIMX, LIMX),
     + TINC(LIMX, LIMX), TTILDEINC(LIMX, LIMX),
     + RINC(LIMX, LIMX), RTILDEINC(LIMX, LIMX)
c     UNITMATRIX is no longer needed
c     UNITMATRIX(LIMX, LIMX)

      UNITY = 1.0
      ZERO = 0.0
C$$$ T~ = D^-1
      CALL SQCOPY (D, TTILDEINC, LIMX)
      CALL INVERTMATRIX(TTILDEINC, LIMX)
C$$$ R~ = BD^-1
      CALL SQDOT (RTILDEINC, B, TTILDEINC, LIMX)
C$$$ R = -D^-1 C
      CALL SQDOTAX (RINC, -UNITY, TTILDEINC, C, LIMX)
C$$$ R=-R
C$$$ T=(A-)? BD^-1 C
C$$$ Thus, We can simply reuse R here: T = A + B * R
C$$$ Also, removed a slow call to zgemm, replacing it with zcopy
      CALL SQCOPY(A, TINC, LIMX)
      CALL SQDOTUPD(UNITY, TINC, UNITY, B, RINC, LIMX)
      
      RETURN
      END

      SUBROUTINE UPDATETANDR(TINC, TTILDEINC, R, RTILDEINC, T, TTILDE,
     + RTILDE, LIMX, RINC)
      IMPLICIT NONE
      INTEGER LIMX
      DOUBLE COMPLEX T(LIMX, LIMX), TTILDE(LIMX, LIMX),
     +               R(LIMX, LIMX), RTILDE(LIMX, LIMX),
     +               TINC(LIMX, LIMX), TTILDEINC(LIMX, LIMX),
     + RINC(LIMX, LIMX), RTILDEINC(LIMX, LIMX)

C$$$ TEMPORARY LOCAL VARIABLES
      DOUBLE COMPLEX T1TEMP(LIMX, LIMX), TTILDE1TEMP(LIMX, LIMX),
     +               R1TEMP(LIMX, LIMX), RTILDE1TEMP(LIMX, LIMX)
      DOUBLE COMPLEX UNITY, ZERO
      DOUBLE COMPLEX BRACKET12(LIMX, LIMX), BRACKET21(LIMX, LIMX),
     +               TRTEMP(LIMX, LIMX), TRTEMP2(LIMX, LIMX),
     + UNITMATRIX(LIMX, LIMX), ALLZERO(LIMX, LIMX)
      UNITY = 1.0
      ZERO = 0.0
      CALL SQCOPY (T,      T1TEMP,      LIMX)
      CALL SQCOPY (TTILDE, TTILDE1TEMP, LIMX)
      CALL SQCOPY (R,      R1TEMP,      LIMX)
      CALL SQCOPY (RTILDE, RTILDE1TEMP, LIMX)
      CALL SQZERO (ALLZERO,    LIMX)
      CALL SQUNIT (UNITMATRIX, LIMX)

C$$$ BRACKET12 = (1 - RTILDE1.R2)^-1
      CALL SQCOPY (UNITMATRIX, BRACKET12, LIMX)
      CALL SQDOTUPD (UNITY, BRACKET12, -UNITY, RTILDE1TEMP, RINC, LIMX)
      CALL INVERTMATRIX(BRACKET12, LIMX)

C$$$ BRACKET21 = (1 - R2.RTILDE1)^-1
      CALL SQCOPY (UNITMATRIX, BRACKET21, LIMX)
      CALL SQDOTUPD (UNITY, BRACKET21, -UNITY, RINC, RTILDE1TEMP, LIMX)
      CALL INVERTMATRIX(BRACKET21, LIMX)


C$$$ T = T2.BRACKET12.T1
C$$$ Note: T2*Br12 can be reused later.
      CALL SQDOT (TRTEMP, TINC, BRACKET12, LIMX)
      CALL SQDOT (T, TRTEMP, T1TEMP, LIMX)

C$$$ RTILDE = RTILDE2 + T2.BRACKET12.RTILDE1.TTILDE2
      CALL SQCOPY (RTILDEINC, RTILDE, LIMX)
      CALL SQDOT (TRTEMP2, TRTEMP, RTILDE1TEMP, LIMX)

      CALL SQDOTUPD (UNITY, RTILDE, UNITY, TRTEMP2, TTILDEINC, LIMX)

C$$$ TTILDE = TTILDE1.BRACKET21.TTILDE2
C$$$ Note: T1~ * Br21 can be reused later. I reorganized the code
C$$$ to make it possible
      CALL SQDOT (TRTEMP, TTILDE1TEMP, BRACKET21, LIMX)
      CALL SQDOT (TTILDE, TRTEMP, TTILDEINC, LIMX)

C$$$ R = R1 + TTILDE1.BRACKET21.R2.T1. Note: TRTEMP from above is reused here
      CALL SQCOPY   (R1TEMP,  R, LIMX)
      CALL SQDOT    (TRTEMP2, TRTEMP, RINC, LIMX)
      CALL SQDOTUPD (UNITY, R, UNITY, TRTEMP2, T1TEMP, LIMX)
      RETURN
      END
