c$$   NEW SUBROUTINE TO GENERATE A,B,C,D
      SUBROUTINE GENABCD(LIMX, MULT, A, B, C, D,U)
      IMPLICIT NONE
      INTEGER LIMX
C    M(A,B,C,D) COMPONENTS OF MATRIX, U: PHASE MATRIX
      DOUBLE COMPLEX  MULT(2*LIMX, 2*LIMX)
      DOUBLE COMPLEX  U(LIMX, LIMX)
      DOUBLE COMPLEX  A(LIMX, LIMX), B(LIMX, LIMX),
     +                C(LIMX, LIMX), D(LIMX, LIMX)
      DOUBLE COMPLEX  MA(LIMX, LIMX), MB(LIMX, LIMX),
     +                MC(LIMX, LIMX), MD(lIMX, LIMX)
      DOUBLE COMPLEX  BU (LIMX, LIMX),  UC (LIMX, LIMX),
     +                DU (LIMX, LIMX),  UDU(LIMX,LIMX),
     +                TEMPBC  (LIMX, LIMX)
      DOUBLE COMPLEX UNITY/1.0/, ZERO/0.0/,   NEG/-1.0/,
     +               HALF/0.5/,  NHALF/-0.5/
     
C     Surely this negates any purpose in rewriting the code?
      MA = MULT(1:LIMX, 1:LIMX)
      MB = MULT(1:LIMX, (LIMX+1):2*LIMX)
      MC = MULT((LIMX+1):2*LIMX, 1:LIMX)
      MD = MULT((LIMX+1):2*LIMX, (LIMX+1):2*LIMX)

c$$$ M(A,B,C,D) INITALISE TO ZERO'S
      CALL SQZERO(A, LIMX)
      CALL SQZERO(B, LIMX)
      CALL SQZERO(C, LIMX)
      CALL SQZERO(D, LIMX)
      CALL SQZERO(TEMPBC, LIMX)

c     BU := Mb*U
      CALL SQDOT(BU, MB, U, LIMX)
c     DU := Md*U
      CALL SQDOT(DU, MD, U, LIMX)
c     UC := Udagger * Mc
      CALL ZGEMM('C', 'N', LIMX, LIMX, LIMX,
     +     UNITY, U,  LIMX,
     +            MC, LIMX,
     +      ZERO, UC, LIMX)
c     UDU := Udagger*DU = Udagger * Md * U     
      CALL ZGEMM('C', 'N', LIMX, LIMX, LIMX,
     +     UNITY, U, LIMX,
     +           DU, LIMX,
     +     ZERO, UDU, LIMX)

c$$$ CREATED TEMP MATICIES -> GENERATE A
c     Micro-optimization: 16 calls to ZAXPY can be replaced by
c     10 calls to ZAXPY + 3 call to ZCOPY -- AVS
c             
c     A := 1/2 Ma
      CALL SQUPDAXPY (A, HALF, MA, LIMX)
c     B := A = 1/2 Ma     
      CALL SQCOPY    (A, B, LIMX)
c     A := 1/2 (Ma + U*Md*U)
      CALL SQUPDAXPY (A, HALF, UDU, LIMX)
c     D := A = 1/2(Ma + U*Md*U)
      CALL SQCOPY    (A, D, LIMX)
c     TEMPBC := 1/2 (Mb*u + u*Mc)      
      CALL SQUPDAXPY (TEMPBC, HALF, BU, LIMX)      
      CALL SQUPDAXPY (TEMPBC, HALF, UC, LIMX)
c     A := 1/2 (Ma + U*Md*U + Mb*u + u * Mc)
      CALL SQUPDAXPY (A, UNITY, TEMPBC, LIMX)
c     D := 1/2(Ma + U*Md*U - Mb*u - u*Mc) 
      CALL SQUPDAXPY (D, NEG,   TEMPBC, LIMX)
c     C:= B := 1/2 (Ma - U*Md*U)
      CALL SQUPDAXPY (B, NHALF, UDU, LIMX)
      CALL SQCOPY    (B, C, LIMX)
c     TEMPBC := TEMPBC - BU = 1/2(-Mb*U + U*Mc)      
      CALL SQUPDAXPY (TEMPBC, NEG, BU, LIMX)
c     B  := 1/2(Ma - U*Md*U - Mb*U + U*Mc)
      CALL SQUPDAXPY (B, UNITY, TEMPBC, LIMX)
c     C  := 1/2(Ma - U*Md*U + Mb*U - U*Mc)
      CALL SQUPDAXPY (C, NEG,   TEMPBC, LIMX)
      
      RETURN
      END
 
C$$$ ROUTINE TO GENERATE T AND R
      SUBROUTINE GENTANDRINC(LIMX, TINC, RINC, TTILDEINC, RTILDEINC,
     +                       A, B, C, D)
      IMPLICIT NONE
      INTEGER LIMX
      DOUBLE COMPLEX UNITY/1.0/
      DOUBLE COMPLEX A(LIMX, LIMX), B(LIMX, LIMX),
     +               C(LIMX, LIMX), D(LIMX, LIMX),
     + TINC(LIMX, LIMX), TTILDEINC(LIMX, LIMX),
     + RINC(LIMX, LIMX), RTILDEINC(LIMX, LIMX)

C$$$ T~ = D^-1
      CALL SQCOPY (D, TTILDEINC, LIMX)
      CALL INVERTMATRIX(TTILDEINC, LIMX)
C$$$ R~ = BD^-1
      CALL SQDOT (RTILDEINC, B, TTILDEINC, LIMX)
C$$$ R = -D^-1 C
      CALL SQDOTAX (RINC, -UNITY, TTILDEINC, C, LIMX)
C$$$ T=(A-)? BD^-1 C
C$$$ Thus, We can simply reuse R here: T = A + B * R
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
     +               RINC(LIMX, LIMX), RTILDEINC(LIMX, LIMX)

C$$$ TEMPORARY LOCAL VARIABLES
      DOUBLE COMPLEX T1TEMP(LIMX, LIMX), TTILDE1TEMP(LIMX, LIMX),
     +               R1TEMP(LIMX, LIMX), RTILDE1TEMP(LIMX, LIMX)
      DOUBLE COMPLEX BRACKET(LIMX, LIMX),
     +               TRTEMP(LIMX, LIMX)
      DOUBLE COMPLEX UNITY/1.0/
      
      CALL SQCOPY (T,      T1TEMP,      LIMX)
      CALL SQCOPY (TTILDE, TTILDE1TEMP, LIMX)
      CALL SQCOPY (R,      R1TEMP,      LIMX)
      CALL SQCOPY (RTILDE, RTILDE1TEMP, LIMX)
      
C$$$ BRACKET := (1 - RTILDE1.R2)^-1 . T1
      CALL SQUNIT   (BRACKET, LIMX)
      CALL SQDOTUPD (UNITY, BRACKET, -UNITY, RTILDE1TEMP, RINC, LIMX)
      CALL INVSOLVE (BRACKET, T1TEMP, LIMX)
      
C$$$ T = T2.BRACKET
      CALL SQDOT  (T, TINC, BRACKET, LIMX)
      
C$$$ R = R1 + TTILDE1.(1-R2.R1~)^-1.R2.T1
C$$$   = R1 + T1~.R_2.(1-R1~.R2)^-1.T1 = R1 + T1~.R2.BRACKET
      CALL SQDOT  (TRTEMP, RINC, BRACKET, LIMX)
      CALL SQCOPY (R1TEMP, R, LIMX)
      CALL SQDOTUPD (UNITY, R, UNITY, TTILDE1TEMP, TRTEMP, LIMX)

C$$$ BRACKET := (1 - R2.RTILDE1)^-1.T2~
      CALL SQUNIT (BRACKET, LIMX)
      CALL SQDOTUPD (UNITY, BRACKET, -UNITY, RINC, RTILDE1TEMP, LIMX)
      CALL INVSOLVE (BRACKET, TTILDEINC, LIMX)
      
C$$$  TTILDE = TTILDE1.BRACKET
      CALL SQDOT (TTILDE, TTILDE1TEMP, BRACKET, LIMX)
      
C$$$ RTILDE = RTILDE2 + T2.inv(1 - R1~.R2).RTILDE1.TTILDE2
C$$$        = R2~ + T2.R1~.inv(1 - R2.R1~).T2~
C$$$        = R2~ + T2.R1~.Bracket      
      CALL SQDOT  (TRTEMP,   RTILDE1TEMP, BRACKET, LIMX)
      CALL SQCOPY (RTILDEINC, RTILDE, LIMX)
      CALL SQDOTUPD (UNITY, RTILDE, UNITY, TINC, TRTEMP, LIMX)
      
      RETURN
      END
