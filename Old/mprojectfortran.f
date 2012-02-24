      PROGRAM GRAPHENE
      INTEGER ,PARAMETER :: LIMX=3, LIMY=3, HSIZE=LIMX*LIMY*LIMX*LIMY,
     +     LSIZE=LIMX*LIMY, WRAPY=0, WRAPX=1
      INTEGER*4 I/1/, J/1/,L/1/, M/1/, P, COLSTAG, ROWSTAG, ABSTAG

      REAL HAMIL(LSIZE, LSIZE)
      DATA HAMIL/HSIZE*0.0/

C$$$ 10   IF (I .LE.4 ) THEN
C$$$         HAMIL(I,I)=1
C$$$         I=I+1
C$$$         GOTO 10
C$$$         ENDIF
C$$$         I=1

C$$$  Let M be Y, L be X
      DO M = 1, LIMY
         DO L = 1, LIMX
            P=L+(LIMX*(M-1))
C$$$  BEGIN Y Link upwards
            IF (P+LIMX .LE. LIMX*LIMY) THEN
               HAMIL(P, P+LIMX)=1
            ENDIF
            IF (P+LIMX .GT. LIMX*LIMY) THEN
               HAMIL(P,(P+LIMX)-(LIMX*LIMY))=WRAPY
            ENDIF
c$$$  BEGIN Y Link downwards
            IF (P-LIMX .GE. 1) THEN
               HAMIL(P, P-LIMX)=1
            ENDIF
            IF (P-LIMX .LT. 1) THEN
               HAMIL(P, (LIMX*LIMY)+(P-LIMX))=WRAPY
            ENDIF
c$$$  BEGIN X Linking
c$$$      int col_stagger = (x % 2) ? 1 : -1; 
c$$$  FORTRAN counts from 1!
c$$$      int row_stagger = (y % 2) ? 1 : -1;
c$$$      int ab_stagger  = col_stagger * row_stagger;
            IF (MOD(L,2) .EQ. 0) THEN
               COLSTAG=-1
            ENDIF
            IF (MOD(L,2) .EQ. 1) THEN
               COLSTAG=1
            ENDIF
            IF (MOD(M,2) .EQ. 0) THEN
               ROWSTAG=-1
            ENDIF
            IF (MOD(M,2) .EQ. 1) THEN
               ROWSTAG=1
            ENDIF
            IF (ABSTAG .GT. 0) THEN
               IF (((P-1) + ABSTAG)/LIMX .EQ. M-1) THEN
                  HAMIL(P, P+ABSTAG)=1
               ENDIF
               IF (((P-1) + ABSTAG)/LIMX .GT. M-1) THEN
                  HAMIL(P, ((M-1)*LIMX)+((P+ABSTAG)-((M)*LIMX)))=WRAPX
               ENDIF
            ENDIF
            IF (ABSTAG .LT. 0) THEN
               IF (((P-1) + ABSTAG)/LIMX .EQ. M-1) THEN
                  HAMIL(P, P+ABSTAG)=1
               ENDIF
               IF (((P-1) + ABSTAG)/LIMX .LT. M-1) THEN
                  HAMIL(P, ((M-1)*LIMX)+((LIMX+ABSTAG))+1)=WRAPX
               ENDIF
            ENDIF


         ENDDO
      ENDDO

C$$$  PRINT HAMILTONIAN         
      DO J = 1, LIMX*LIMY
         WRITE (*,20) (HAMIL(J,I), I = 1, LIMX*LIMY)
      END DO	
      
 20   FORMAT (9F2.0)
      END
      

c$$$#define RESTRICT(x, L) (((x) >= 0) ? ((x) % (L)) : ((x) + (L)))  // Restrict x to the interval [0:L-1]
c$$$                                                                 // This is needed for wrapping
c$$$                                                                 // Note: does not work when x < -L
c$$$#define SITE(x, y, L, M) (RESTRICT((x), (L)) + RESTRICT ((y), (M)) * (L) )  // Site index 
c$$$int fillH(double H[L*M][L*M]){
c$$$  int x=0;
c$$$  int y=0;
c$$$  for (y = 0; y < M; y ++) {
c$$$    for (x = 0; x < L; x++) {
c$$$      int current     = SITE (x, y, L, M); 
c$$$      int col_stagger = (x % 2) ? 1 : -1; 
c$$$      int row_stagger = (y % 2) ? 1 : -1;
c$$$      int ab_stagger  = col_stagger * row_stagger; // +1 on A, -1 on B sublattice 
c$$$      int h_neigh     = SITE (x + ab_stagger, y, L, M); // horizontal neighbour,
c$$$                                                        // left or right, depending on AB
c$$$      int up_neigh    = SITE (x, y + 1, L, M);          // A neighbour upstairs
c$$$      int down_neigh  = SITE (x, y - 1, L, M);          // A neighbour downstairs
c$$$      
c$$$      if ((y < M - 1) || WRAPY ) { //link above
c$$$          H[current][up_neigh] = 1;
c$$$      }
c$$$
c$$$      if ((y > 0) || WRAPY) {      //link below
c$$$          H[current][down_neigh] = 1;
c$$$      }
c$$$
c$$$      if (((x > 0) && (x <  L - 1)) || WRAPX) 
c$$$           H[current][h_neigh] = 1; 
c$$$      }
c$$$     
c$$$  }
c$$$  return 0;
c$$$  //array is always Hermitian or symmetric so should only need 
c$$$  //to run on half the values - fix this!
c$$$
c$$$}
