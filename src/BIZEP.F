C.... =================================================================
C.... PROGRAM FOR SOLUTION OF TWO DIMENSIONAL
C.... NON-LINEAR MATERIAL PROBLEMS BY THE
C.... BOUNDARY ELEMENT METHOD
C.... J. C. DE FARIA TELLES
C.... =================================================================
C.... NE - Number of elements
C.... NC - Number of cells
C.... NN - Number of nodes
C.... NP - Number of internal points
C.... IPL - Index for type of problem, 1-plane stress/2-plane strain
C.... E - Youngs modulus
C.... ET - Hardening slope for computation of H' (from Eq. 6.11/6.12)
C.... HL - H' for Taylor series approximation of stress
C.... SY - ???
C.... CC - Cohesion
C.... PHI - Internal Friction angle [in deg]
C.... PO - Poissons Ratio
C.... NN2 - Number of nodes times 2, NN2 = 2 * NN
C.... NT - Internal + boundary nodes, NT = NN + NP
C.... NT3 - ???
C.... C1 to C11 - Constants needed by fundamental solution
C.... PER - Load increment in %
C.... SM - ???
C.... C12 - ???
C.... IFIP - Index for boundary type, 0-zero traction, 1-prescribed displacement
C.... ISYM - Indicates symmetry (0-no sym, 1-y sym, 2-x sym, 3-xy sym)
C.... XSYM - X coordinates of symmetry plane
C.... YSYM - Y coordinates of symmetry plane
C.... INFB - Infinite boundary index (1-infinite case, 0-bounded case)
C.... C13 to C18 - ???
C.... SNPHI - Sine of internal friction angle
C.... IYIED - Yield criteria 1-Tesca/2-Von Mises/3-Mohr-Coulomb/4-Drucker-Prager
C.... ITECH - Iteration methodology (Section 7.5 page 291)
C.... D - Kronecker delta
C.... XI - Acissas of Gaussian integration points
C.... W - Weights of Gausssian integration points
C.... IDUP - Double nodes vector (=0 no double node, =N N double nodes)
C.... INC - Boundary nodes incidencies
C.... ICC - Cells conectivity
C.... NFIP - Number of nodes with displacement
C.... NDFIP - Number of nodes with nonzero tractions prescribed
C.... =================================================================
      PROGRAM BIZEP

      COMMON /RW/ IRE,IWR
      COMMON /A/ X(73),Y(73),D(2,2),A(6,3),W(6,3)
      COMMON /A1/ ISYM(73),INC(33,2),IDUP(73),ICC(110,3),AR(110),C(33)
      COMMON /A2/ DI(219,219),DB(76,219)
      COMMON /A3/ IFIP(76),U(76),P(76),S(73,4),SR(219),DSR(219),DSE(219)
     *,XM(76),XN(219),EPE(73),DEPE(73),SR3(73),SEQ(73),ASEQ(73),DSR3(73)

C.... BIND INPUT AND OUTPUT FILES TO IRE AND IRW VARIABLES
      OPEN(UNIT=11, FILE="src/input.dat", STATUS='OLD')
      OPEN(UNIT=12, FILE="src/results.out")
      IRE=11
      IWR=12

C.... TOLERANCE FOR ??
      TOL=0.001
      J=0

C     READ INPUT FILE
      CALL INPUT(NE,NC,NN,NP,IPL,E,ET,SY,PO,NN2,NT,NT3,C1,C2,
     *C3,C4,C5,C6,C7,C8,C9,C10,C11,PER,SM,C12,HL,IDSYM,XSYM,YSYM,INFB,
     *C13,C14,C15,C16,C17,C18,SNPHI,IYIED,ITECH)

C.... COMPUTE MATRICES
      CALL MATRX(NE,NC,NN,NN2,NT,NT3,C1,C2,C3,C4,C5,C6,C7,C8
     *,C9,C10,C11,PO,IDSYM,XSYM,YSYM,INFB,C15,C16,C17,C18)

C.... SOLVE ELASTIC PROBLEM AND CALCULATE LOAD AT FIRST YIELD
      CALL ELAST(C11,XLAMB,SY,NT3,NT,PER,NN2,FAC,IPRT,IYIED,SNPHI)
      GO TO 3

C.... INCREMENTAL PROCESS
    1 IF(XLAMB.GE.1.)GO TO 4

C.... COMPUTATION OF LOAD FACTOR
      XLAMB=XLAMB+FAC

C.... ITERATIONS
      DO 2 J=1,200
C.... SOLVE ELASTOPLASTIC PROBLEM AND COMPUTE RESIDUAL STRESS
C.... INCREMENTS
      IF(ITECH.EQ.0)CALL SOLVE(TOL,IC,SY,NT3,XLAMB,NT,NN2,HL,SNPHI,C11,F
     *AC,J,IYIED,C14,IPL,SM)
      IF(ITECH.NE.0)CALL SOLVE1(TOL,IC,SY,NT3,XLAMB,NT,NN2,HL,SNPHI,C11,
     *FAC,J,IYIED,C14,IPL,SM)
C.... VERIFY CONVERGENCE
      IF(IC.EQ.0)GO TO 3
    2 CONTINUE

C.... OUTPUT RESULTS
    3 CALL OUTPT(NN,NT,XLAMB,IPRT,J,SM,C13)
      GO TO 1
    4 STOP
      END PROGRAM BIZEP

      SUBROUTINE INPUT(NE,NC,NN,NP,IPL,E,ET,SY,PO,NN2,NT,NT3,C1,C2,
     *C3,C4,C5,C6,C7,C8,C9,C10,C11,PER,G,C12,HL,IDSYM,XSYM,YSYM,INFB,
     *C13,C14,C15,C16,C17,C18,SNPHI,IYIED,ITECH)
      COMMON /RW/ IRE,IWR
      COMMON /A/ X(73),Y(73),D(2,2),A(6,3),W(6,3)
      COMMON /A1/ ISYM(73),INC(33,2),IDUP(73),ICC(110,3),AR(110),C(33)
      COMMON /A3/ IFIP(76),U(76),P(76),S(73,4),SR(219),DSR(219),DSE(219)
     *,XM(76),XN(219),EPE(73),DEPE(73),SR3(73),SEQ(73),ASEQ(73),DSR3(73)
      CHARACTER(LEN=25) :: Title

      WRITE(IWR,1)
    1 FORMAT(//////////,24X,'* * *  B O U N D A R Y  E L E M E N T  M E
     *T H O D  A P P L I E D  T O  * * *',//,24X,'* * *  P L A N E  E L
     *A S T O P L A S T I C  P R O B L E M S            * * *',///)

C.... TITLE OF PROBLEM
      READ(IRE, 2) Title
      WRITE(IWR, *) 'STUDY CASE: ', Title
    2 FORMAT(A)

C.... GENERAL INFORMATION ABOUT THE PROBLEM
      READ(IRE, *) INFB, ITECH, NE, NC, NN, NP, IPL, IDSYM

C.... LOAD MATERIAL PROPERTIES AND YEILD CRITERIA
      READ(IRE, *) IYIED, E, ET, SY, PO, CC, PHI

      SNPHI = SIN(PHI * 0.017453292)
      GO TO (71, 71, 72, 73, 80), IYIED
   72 SY = CC * COS(PHI * 0.017453292)
      GO TO 71
   73 SY=6. * CC * COS(PHI  *0.017453292)/(1.73205080757 * (3. - SNPHI))
      GO TO 71
   80 SY=1.73205080757*CC*COS(PHI*0.017453292)/SQRT(3.+SNPHI**2)

   71 IF(INFB.EQ.0) GO TO 60

      WRITE(IWR,61)
   61 FORMAT(//,13X,'* INFINITE BOUNDARY *')

   60 IF(ITECH.EQ.0) WRITE(IWR,62)
   62 FORMAT(//,13X,'* STANDARD IN. STR. *')

      IF(ITECH.NE.0)WRITE(IWR,63)
   63 FORMAT(//,13X,'* MODIFIED IN. STR. *')

      WRITE(IWR, 4) NE, NC, NN, NP, IPL, IDSYM
    4 FORMAT(//,15X,'NO. ELEMENTS =',I5,//,15X,'NO.    CELLS =',I5,//,15
     *X,'NO.    NODES =',I5,//,15X,'NO.   POINTS =',I5,//,15X,'PROBL.  T
     *YPE =',I5,//,15X,'SYMME.  TYPE =',I5)

      WRITE(IWR,74) IYIED, E, ET, SY, PO, CC, PHI
   74 FORMAT(//,15X,'MATERIAL PROPERTIES',//,15X,'Y. CRIT.=',I10,
     *//,15X,'E       =',F10.0,//,15X,'ET      =',F10.0,//,15X,
     *'FY      =',F15.5,//,15X,'POISSON =',F15.5,//,15X,'C       =',
     *F15.5,//,15X,'PHI     =',F15.5,///,30X,'COORDINATES OF BOUNDARY NO
     *DES',//,12X,'NODE',14X,'X',15X,'Y',/)

      NN2 = NN * 2
      NT = NN + NP
      NT3 = 3 * NT
C.... NODES AND POINTS COORDINATES
      DO 5 I = 1, NT
        READ(IRE, *) K, X(K), Y(K), IDUP(K), ISYM(K)
C....   VERIFY IF DOUBLE NODE
        IF(IDUP(K).EQ.0) GO TO 5
        J = IDUP(K)
        IDUP(J) = K
        X(K) = X(J)
        Y(K) = Y(J)
    5 CONTINUE

      WRITE(IWR,7) (K, X(K), Y(K), K=1, NN)
    7 FORMAT(10X, I5, 5X, F15.4, 1X, F15.4)

      IF(NP.EQ.0) GO TO 9

      WRITE(IWR, 8)
    8 FORMAT(//,30X,'COORDINATES OF INTERNAL POINTS',//,11X,'POINT',14X,
     *'X',15X,'Y',/)

      K = NN + 1

      WRITE(IWR, 7)(J, X(J), Y(J), J=K, NT)
C.... NODES AND POINTS AT SYMMETRY LINES
    9 IF(IDSYM.EQ.0) GO TO 49

      WRITE(IWR, 42)
   42 FORMAT(//,30X,'BOUNDARY NODES AND INTERNAL POINTS AT SYMMETRY LINE
     *(S)',//,12X,'L. X',12X,'L. Y',/)

      DO 43 K=1, NT
        IF(ISYM(K).EQ.0) GO TO 43
        IZZ = ISYM(K)
        GO TO (44, 45, 46), IZZ
   44   YSYM = Y(K)
        WRITE(IWR, 47) K
   47   FORMAT(10X, I5)
        GO TO 43
   45   XSYM = X(K)
        WRITE(IWR, 48) K
   48   FORMAT(26X, I5)
        GO TO 43
   46   WRITE(IWR, 50) K, K
   50   FORMAT(10X, I5, 11X, I5)
   43 CONTINUE

C.... ELEMENT INCIDENCES
   49 WRITE(IWR, 10)
   10 FORMAT(//,30X,'ELEMENT INCIDENCES',//,13X,'EL',13X,'N. 1',12X,'N.
     *2',14X,'L',/)
      DO 11 I=1, NE
      READ(IRE, *) K, INC(K,1), INC(K,2)
      II = INC(K, 1)
      IF = INC(K, 2)
C.... COMPUTE LENGHT OF ELEMENT
   11 C(K) = SQRT((X(IF) - X(II))**2 + (Y(IF) - Y(II))**2)
      WRITE(IWR,13)(I, INC(I, 1), INC(I, 2), C(I), I=1, NE)
   13 FORMAT(10X, I5, 11X, I5, 11X, I5, 5X, F15.4)
      IF(NC.EQ.0) GO TO 14
C.... CELLS CONNECTIVITY
      WRITE(IWR, 17)
   17 FORMAT(//,30X,'CELLS CONNECTIVITY',//,12X,'CELL',12X,'P. 1',12X,'P
     *. 2',12X,'P. 3',14X,'A',/)
      DO 15 I=1, NC
      READ(IRE, *) K, ICC(K, 1), ICC(K, 2), ICC(K, 3)
      I1 = ICC(K, 1)
      I2 = ICC(K, 2)
      I3 = ICC(K, 3)
   15 AR(K)=0.5*ABS((Y(I2)-Y(I3))*(X(I1)-X(I3))-(X(I2)-X(I3))*(Y(I1)-Y(I
     *3)))
      WRITE(IWR,18)(I,ICC(I,1),ICC(I,2),ICC(I,3),AR(I),I=1,NC)
   18 FORMAT(10X,I5,11X,I5,11X,I5,11X,I5,5X,F15.4)
C.... CONSTANTS
C.... COMPUTE SHEAR MODULUS
   14 G=E/(2.*(1.+PO))
      C12=PO
C.... COMPUTE H' HARDENING FACTOR
      HL=E*ET/(E-ET)
      C13=PO/(1.+PO)
      C11=PO
C.... If IPL = 1 (plane-stress) them assertion is false
      IF(IPL-1) 40, 40, 41
   40 PO=C13
      C11=0.
   41 C2=3.-4.*PO
      C18=1./(1.-PO)
      C3=C18/12.56637062
      C4=1.-2.*PO
      C14=PO/C4
      C6=2.*C3*G
      C7=1.-4.*PO
      C1=C3/(2.*G)
      C5=C1/2.
      C8=2.*G*C18
      C9=PO*C18
      C10=(2.-PO)*C18
      C15=(5.-4.*PO)*C18/8.
      C16=C7*C18/8.
      C17=C2*C18/4.
C.... BOUNDARY VALUES PRESCRIBED
C.... Zero P am XM matrix,  XM to be right side vector, P vector of loads
      DO 19 I=1, NN2
      P(I) = 0.
      XM(I) = 0.
   19 IFIP(I) = 0
      READ(IRE, *) NFIP, NDFIP, PER
      WRITE(IWR, 21) NFIP, NDFIP, PER
   21 FORMAT(//,15X,'NO. DISPL. PRESC. =',I5,//,15X,'NO. TRACT. PRESC. =
     *',I5,//,15X,'LOAD INCREMENT    =',F10.4,'%',///,15X,'DISPLACEMENTS
     *',//,12X,'NODE',14X,'U',15X,'V',/)
      PER=PER/100.

      IF(NFIP.EQ.0) GO TO 22

      DO 23 I=1, NFIP
        READ(IRE, *) K, P(2*K-1), P(2*K), IFIP(2*K-1), IFIP(2*K)
        IND = IFIP(2*K-1) + 2*IFIP(2*K)
        GO TO (25, 26, 27), IND
   25   WRITE(IWR,28) K, P(2*K-1)
   28   FORMAT(10X, I5, 5X, F15.4)
        GO TO 23
   26   WRITE(IWR, 29) K, P(2*K)
   29   FORMAT(10X, I5, 21X, F15.4)
        GO TO 23
   27   WRITE(IWR, 30) K, P(2*K-1), P(2*K)
   30   FORMAT(10X, I5, 5X, F15.4, 1X, F15.4)
   23 CONTINUE

   22 IF(NDFIP.EQ.0) GO TO 31
      WRITE(IWR, 34)
   34 FORMAT(//,15X,'TRACTIONS',//,12X,'NODE',13X,'TX',14X,'TY',/)
      DO 32 I=1, NDFIP
      READ(IRE, *) K, P(2*K - 1), P(2*K)
   32 WRITE(IWR, 30) K, P(2*K-1), P(2*K)
C.... INTEGRATION POINTS
   31 A(1,3)=-0.932469514203152
      A(2,3)=-0.661209386466265
      A(3,3)=-0.238619186083197
      A(4,3)=-A(3,3)
      A(5,3)=-A(2,3)
      A(6,3)=-A(1,3)
      W(1,3)=0.171324492379170
      W(2,3)=0.360761573048139
      W(3,3)=0.467913934572691
      W(4,3)=W(3,3)
      W(5,3)=W(2,3)
      W(6,3)=W(1,3)
      A(1,2)=-0.861136311594053
      A(2,2)=-0.339981043584856
      A(3,2)=-A(2,2)
      A(4,2)=-A(1,2)
      W(1,2)=0.347854845137454
      W(2,2)=0.652145154862546
      W(3,2)=W(2,2)
      W(4,2)=W(1,2)
      A(1,1)=-0.577350269189626
      A(2,1)=-A(1,1)
      W(1,1)=1.
      W(2,1)=1.
      RETURN
      END

      SUBROUTINE MATRX(NE,NC,NN,NN2,NT,NT3,C1,C2,C3,C4,C5,C6,C7,C8,
     *C9,C10,C11,PO,IDSYM,XSYM,YSYM,INFB,C15,C16,C17,C18)
      COMMON /A/ X(73),Y(73),D(2,2),A(6,3),W(6,3)
      COMMON /A1/ ISYM(73),INC(33,2),IDUP(73),ICC(110,3),AR(110),C(33)
      COMMON /A2/ DI(219,219),DB(76,219)
      COMMON /A3/ IFIP(76),U(76),P(76),S(73,4),SR(219),DSR(219),DSE(219)
     *,XM(76),XN(219),EPE(73),DEPE(73),SR3(73),SEQ(73),ASEQ(73),DSR3(73)
      COMMON /A4/ H(3,4),G(3,4),DM(3,9)
C.... KRONECKER DELTA
      D(1,1)=1.
      D(2,2)=1.
      D(1,2)=0.
      D(2,1)=0.
C.... CLEAR ARRAYS
      DO 1 I=1,NT3
      XN(I)=0.
      DO 1 J=1,NT3
    1 DI(I,J)=0.
      DO 4 I=1,NN2
      DO 4 J=1,NT3
    4 DB(I,J)=0.
      IFA=1
      NIF=1
      IF(IDSYM.EQ.1) IFA = 2
      IF(IDSYM.NE.2)GO TO 60
      IFA=3
      NIF=2
   60 IF(IDSYM.EQ.3) IFA = 4
C.... TEST FOR INFINITE BOUNDARY
      IF(INFB.EQ.0) GO TO 90
      DO 91 I=1,NN2
      IF(IFIP(I).NE.0) GO TO 92
      DI(I,I)=1.
      GO TO 91
   92 XM(I)=-P(I)
   91 CONTINUE
C.... IBA=1 - SING. POINT AT BOUNDARY NODES
C.... IBA=2 - SING. POINT AT INTERNAL POINTS
   90 DO 2 IBA=1,2
      GO TO (80,81),IBA
   80 IIH=1
      IFH=NN
      IP=2
      GO TO 82
   81 IIH=NN+1
      IFH=NT
      IP=3
      ICOD=4
C.... FIND VECTOR XM AND MATRIX DB MODIFIED
      CALL GAUSP(NN2,NT3,C8)
C.... CLEAR THE PART OF MATRIX DI WHICH WAS USED FOR MATRIX
C.... HB REORDERED
      DO 93 I=1,NN2
      DO 93 J=1,NN2
   93 DI(I,J)=0.
   82 DO 2 ISY=1,IFA,NIF
      IF(IBA.EQ.1)GO TO 83
      IIS=2
      IFS=2
      GO TO 70
   83 GO TO (70,71,71,73),ISY
   71 IIS=4-ISY
      IFS=IIS
      GO TO 70
   73 IIS=1
      IFS=2
   70 DO 2 I=IIH,IFH
      XS=X(I)
      YS=Y(I)
      IF(ISY.EQ.2.OR.ISY.EQ.4)YS=2.*YSYM-YS
      IF(ISY.GE.3)XS=2.*XSYM-XS
C.... COMPUTE MATRIX HB REORDERED IN THE ARRAY OF MATRIX DI AND VECTOR
C.... XM AS INDEPENDENT TERM FOR BOUNDARY CONSTRAINT EQ., AS WELL AS
C.... LOWER PART OF BOTH, MATRIX DI MODIFIED AND VECTOR XN
      DO 10 J=1,NE
      II=INC(J,1)
      IF=INC(J,2)
      IF(I.GT.NN)GO TO 6
      ICOD=1
      IF(ISY.NE.1.AND.ISYM(I).NE.(ISY-1))GO TO 6
      IF(I.EQ.II.OR.I.EQ.IDUP(II))ICOD=2
      IF(I.EQ.IF.OR.I.EQ.IDUP(IF))ICOD=3
    6 CALL FUNC(ICOD,J,C(J),C1,C2,C3,C4,C5,C6,C7,PO,II,IF,IP,XS,YS,ISY,I
     *IS,IFS)
      DO 10 K=1,IP
      JJ=IP*(I-1)+K
      M=0
      DO 10 NX=1,2
      DO 10 NV=1,2
      M=M+1
      IC=2*INC(J,NX)+NV-2
      IF(I.GT.NN)GO TO 5
      IF(IFIP(IC).NE.0)GO TO 67
      DI(JJ,IC)=DI(JJ,IC)+H(K,M)
      XM(JJ)=XM(JJ)+G(K,M)*P(IC)
      GO TO 68
   67 DI(JJ,IC)=DI(JJ,IC)-G(K,M)*C8
      XM(JJ)=XM(JJ)-H(K,M)*P(IC)
   68 GO TO (61,62,63,64),ISY
   62 IF(NV-2)61,64,61
   63 IF(NV-1)61,64,61
   64 H(K,M)=-H(K,M)
   61 IF(IFIP(JJ+NV-K).NE.0)GO TO 69
      DI(JJ,JJ+NV-K)=DI(JJ,JJ+NV-K)-H(K,M)
      GO TO 10
   69 XM(JJ)=XM(JJ)+H(K,M)*P(JJ+NV-K)
      GO TO 10
    5 IF(IFIP(IC).NE.0)GO TO 75
      DO 84 JX=1,NT3
   84 DI(JJ,JX)=DI(JJ,JX)-H(K,M)*DB(IC,JX)
      XN(JJ)=XN(JJ)+G(K,M)*P(IC)-H(K,M)*XM(IC)
      GO TO 10
   75 DO 85 JX=1,NT3
   85 DI(JJ,JX)=DI(JJ,JX)+G(K,M)*DB(IC,JX)
      XN(JJ)=XN(JJ)-H(K,M)*P(IC)+G(K,M)*XM(IC)
   10 CONTINUE
      IF(NC.EQ.0)GO TO 2
C.... COMPUTE MATRIX DB AND LOWER PART OF MATRIX DI
      DO 17 J=1,NC
      ICODC=0
      I1=ICC(J,1)
      I2=ICC(J,2)
      I3=ICC(J,3)
      IF(ISY.EQ.1)GO TO 65
      IF(ISYM(I).NE.3.AND.ISYM(I).NE.(ISY-1))GO TO 25
   65 IF(I1.NE.I.AND.I.NE.IDUP(I1))GO TO 22
      IIC=I2
      IFC=I3
      ICODC=1
      GO TO 25
   22 IF(I2.NE.I.AND.I.NE.IDUP(I2))GO TO 23
      IIC=I3
      IFC=I1
      ICODC=2
      GO TO 25
   23 IF(I3.NE.I.AND.I.NE.IDUP(I3))GO TO 25
      IIC=I1
      IFC=I2
      ICODC=3
   25 CALL FANC(ICOD,ICODC,J,AR(J),C1,C2,C3,C4,C5,C6,C7,PO,I1,I2,I3,IIC,
     *IFC,IP,C11,XS,YS,ISY,IIS,IFS)
      DO 17 K=1,IP
      JJ=IP*(I-1)+K
      M=0
      DO 17 NX=1,3
      DO 17 NV=1,3
      M=M+1
      IC=3*ICC(J,NX)+NV-3
      IF(I.GT.NN)GO TO 7
      DB(JJ,IC)=DB(JJ,IC)+DM(K,M)
      GO TO 17
    7 DI(JJ,IC)=DI(JJ,IC)+DM(K,M)
   17 CONTINUE
      IF(ISY.GT.1.OR.IBA.EQ.1)GO TO 2
C.... INDEPENDENT TERMS IN MATRIX DI
      DI(3*I-2,3*I-2)=DI(3*I-2,3*I-2)+C15
      DI(3*I-2,3*I)=DI(3*I-2,3*I)-C16
      DI(3*I-1,3*I-1)=DI(3*I-1,3*I-1)+C17
      DI(3*I,3*I-2)=DI(3*I,3*I-2)-C16
      DI(3*I,3*I)=DI(3*I,3*I)+C15
    2 CONTINUE
C.... COMPUTE UPPER PART OF BOTH, MATRIX DI MODIFIED AND
C.... VECTOR XN
      DO 30 I=1,NE
      II=INC(I,1)
      IF=INC(I,2)
      CC1=(Y(IF)-Y(II))/C(I)
      CC2=(X(II)-X(IF))/C(I)
      CALL FENC(C8,C9,C10,C11,CC1,CC2,I,C18,PO)
      DO 31 JP=1,2
      IIF=INC(I,JP)
      XFAC=2.
      IF(IDUP(IIF).NE.0.OR.ISYM(IIF).NE.0)XFAC=1.
      DO 31 IR=1,3
      M=0
      JJ=3*IIF+IR-3
      DO 32 IP=1,2
      IO=2*IIF+IP-2
      IF(IFIP(IO).NE.0)GO TO 76
      XN(JJ)=XN(JJ)+G(IR,IP)*P(IO)/XFAC
      GO TO 78
   76 XN(JJ)=XN(JJ)+G(IR,IP)*XM(IO)/XFAC
      DO 86 JX=1,NT3
   86 DI(JJ,JX)=DI(JJ,JX)+G(IR,IP)*DB(IO,JX)/XFAC
   78 DO 32 JR=1,2
      M=M+1
      IO=2*INC(I,IP)+JR-2
      IF(IFIP(IO).NE.0)GO TO 79
      XN(JJ)=XN(JJ)-H(IR,M)*XM(IO)/XFAC
      DO 87 JX=1,NT3
   87 DI(JJ,JX)=DI(JJ,JX)-H(IR,M)*DB(IO,JX)/XFAC
      GO TO 32
   79 XN(JJ)=XN(JJ)-H(IR,M)*P(IO)/XFAC
   32 CONTINUE
      DO 31 JR=1,3
      IO=3*IIF+JR-3
   31 DI(JJ,IO)=DI(JJ,IO)+DM(IR,JR)/XFAC
   30 CONTINUE
      RETURN
      END
      SUBROUTINE FUNC(ICOD,JA,C,C1,C2,C3,C4,C5,C6,C7,PO,II,IF,IP,XS,YS,I
     *SY,IIS,IFS)
C.... INTEGRALS OVER ALL BOUNDARY ELEMENTS
      COMMON /A/ X(73),Y(73),D(2,2),A(6,3),W(6,3)
      COMMON /A4/ H(3,4),G(3,4),DM(3,9)
      DIMENSION DXY(2),BN(2),B(2),DR(2),U(2,2),P(2,2),V(2,2,2),PL(2,2,2)
      DO 5 KK=1,IP
      DO 5 L=1,4
      G(KK,L)=0.
    5 H(KK,L)=0.
      DXY(1)=X(IF)-X(II)
      DXY(2)=Y(IF)-Y(II)
      GO TO (1,2,2,1),ICOD
    1 BN(1)=DXY(2)/C
      BN(2)=-DXY(1)/C
C.... SELECT NO. INTEGRATION POINTS
      SEL=0.5*SQRT((2.*XS-X(II)-X(IF))**2+(2.*YS-Y(II)-Y(IF))**2)/C
      NPI=4
      IF(SEL.LE.1.5)NPI=6
      IF(SEL.GT.5.5)NPI=2
      INP=NPI/2
C.... COMPUTE MATRICES H AND G NUMERICALLY
      DO 50 KK=1,NPI
      XMXI=0.5*(1.+A(KK,INP))*DXY(1)+X(II)-XS
      YMYI=0.5*(1.+A(KK,INP))*DXY(2)+Y(II)-YS
      R=SQRT(XMXI**2+YMYI**2)
      B(1)=-0.25*(A(KK,INP)-1.)*C
      B(2)=0.25*(A(KK,INP)+1.)*C
      DR(1)=XMXI/R
      DR(2)=YMYI/R
      DRDN=DR(1)*BN(1)+DR(2)*BN(2)
      IF(ICOD.EQ.4)GO TO 10
C.... COMPUTE MATRICES H AND G (BOUNDARY CONSTRAINT EQUATION)
      DO 6 I=1,2
      DO 6 J=1,2
      U(I,J)=-C1*(C2*ALOG(R)*D(I,J)-DR(I)*DR(J))
    6 P(I,J)=-C3*((C4*D(I,J)+2.*DR(I)*DR(J))*DRDN+C4*(DR(J)*BN(I)-DR(I)*
     *BN(J)))/R
      DO 7 LA=1,2
      IC=0
      DO 7 LL=1,2
      DO 7 JJ=1,2
      IC=IC+1
      G(LA,IC)=G(LA,IC)+U(LA,JJ)*B(LL)*W(KK,INP)
    7 H(LA,IC)=H(LA,IC)+P(LA,JJ)*B(LL)*W(KK,INP)
      GO TO 50
C.... COMPUTE MATRICES H AND G (STRESSES AT INTERNAL POINTS)
   10 DO 11 I=1,2
      DO 11 J=I,2
      DO 11 K=1,2
      V(I,J,K)=C3*(C4*(DR(J)*D(K,I)+DR(I)*D(K,J)-DR(K)*D(I,J))+2.*DR(I)*
     *DR(J)*DR(K))/R
      B1=2.*DRDN*(C4*DR(K)*D(I,J)+PO*(DR(J)*D(I,K)+DR(I)*D(J,K))-4.*DR(I
     *)*DR(J)*DR(K))
      B2=2.*PO*(BN(I)*DR(J)*DR(K)+BN(J)*DR(I)*DR(K))
      B3=C4*(2.*BN(K)*DR(I)*DR(J)+BN(J)*D(I,K)+BN(I)*D(J,K))
   11 PL(I,J,K)=C6*(B1+B2+B3-C7*BN(K)*D(I,J))/R**2
      IL=0
      DO 12 I=1,2
      DO 12 J=I,2
      IL=IL+1
      IC=0
      DO 12 IAA=1,2
      DO 12 JAA=1,2
      IC=IC+1
      G(IL,IC)=G(IL,IC)+B(IAA)*V(I,J,JAA)*W(KK,INP)
   12 H(IL,IC)=H(IL,IC)+B(IAA)*PL(I,J,JAA)*W(KK,INP)
   50 CONTINUE
      GO TO 18
C.... COMPUTE MATRICES H AND G ANALYTICALLY (BOUNDARY CONSTRAINT EQ.)
    2 AL=C5*C2*C
      AA=AL*(0.5-ALOG(C))
      DO 15 I=1,2
      DO 15 J=1,4
      IT=(J/2)*2+2-J
      G(I,J)=C5*DXY(I)*DXY(IT)/C
      IF(IT.EQ.I)G(I,J)=G(I,J)+AA
   15 CONTINUE
      IAA=-2
      IF(ICOD.EQ.3)IAA=0
      G(1,3+IAA)=G(1,3+IAA)+AL
      G(2,4+IAA)=G(2,4+IAA)+AL
      H(1,2-IAA)=C3*C4*(1.+IAA)
      H(2,1-IAA)=-H(1,2-IAA)
C.... SYMMETRY TEST
   18 IF(ISY.EQ.1.OR.(ISY+ICOD).EQ.8)GO TO 8
      DO 24 I=IIS,IFS
      DO 24 J=1,4
      H(I,J)=-H(I,J)
   24 G(I,J)=-G(I,J)
    8 RETURN
      END
      SUBROUTINE FANC(ICOD,ICODC,JA,AR,C1,C2,C3,C4,C5,C6,C7,PO,I1,I2,I3,
     *IIC,IFC,IP,C11,XS,YS,ISY,IIS,IFS)
C.... INTEGRALS OVER ALL INTERNAL CELLS
      COMMON /A/ X(73),Y(73),D(2,2),A(6,3),W(6,3)
      COMMON /A4/ H(3,4),G(3,4),DM(3,9)
      DIMENSION B(3),DR(2),DX(2),DY(2),EI(2,2,2),EII(2,2,2,2),AA(3),BB(3
     *),Z(3),CSI(3),FI(3)
      DO 1 K=1,IP
      DO 1 L=1,9
    1 DM(K,L)=0.
      AA(1)=X(I3)-X(I2)
      BB(1)=Y(I2)-Y(I3)
      AA(2)=X(I1)-X(I3)
      BB(2)=Y(I3)-Y(I1)
      AA(3)=X(I2)-X(I1)
      BB(3)=Y(I1)-Y(I2)
      IF(ICODC.EQ.0)GO TO 30
      FIMIN=ATAN2((Y(IIC)-YS),(X(IIC)-XS))
      FIMAX=ATAN2((Y(IFC)-YS),(X(IFC)-XS))
      IF(FIMAX.LT.FIMIN)FIMAX=FIMAX+6.283185308
      ITL=1
      GO TO 16
   30 CSI(1)=(X(I2)*Y(I3)-X(I3)*Y(I2)+BB(1)*XS+AA(1)*YS)/(2.*AR)
      CSI(2)=(X(I3)*Y(I1)-X(I1)*Y(I3)+BB(2)*XS+AA(2)*YS)/(2.*AR)
      CSI(3)=(X(I1)*Y(I2)-X(I2)*Y(I1)+BB(3)*XS+AA(3)*YS)/(2.*AR)
      ICO=0
      FI(1)=ATAN2((Y(I1)-YS),(X(I1)-XS))
      FI(2)=ATAN2((Y(I2)-YS),(X(I2)-XS))
      FI(3)=ATAN2((Y(I3)-YS),(X(I3)-XS))
   35 FIMAX=-6.283185308
      FIMIN=-FIMAX
      ICO=ICO+1
      DO 32 IF=1,3
      IF(FI(IF).GE.FIMIN)GO TO 33
      FIMIN=FI(IF)
      IMI=IF
   33 IF(FI(IF).LE.FIMAX)GO TO 32
      FIMAX=FI(IF)
      IMA=IF
   32 CONTINUE
      IF((FIMAX-FIMIN).LE.3.1415926535)GO TO 34
      FI(IMI)=FI(IMI)+6.283185308
      IF(ICO.LE.2)GO TO 35
   34 IME=6-IMA-IMI
      ITL=2
   16 DO 38 IT=1,ITL
      IF(ITL.EQ.1)GO TO 36
      IF(IT.EQ.2)GO TO 37
      IR=IMA
      FIMAX=FI(IME)
      GO TO 36
   37 IR=IMI
      FIMAX=FI(IMA)
      FIMIN=FI(IME)
   36 DFI=(FIMAX-FIMIN)/2.
      IF(DFI.LE.0.0001)GO TO 38
C.... SELECT NO. INTEGRATION POINTS
      NPI=4
      IF(DFI.GT.0.5672)NPI=6
      IF(DFI.LT.0.1309)NPI=2
      INP=NPI/2
      DO 2 KK=1,NPI
      FII=0.5*(A(KK,INP)*(FIMAX-FIMIN)+FIMAX+FIMIN)
      DR(1)=COS(FII)
      DR(2)=SIN(FII)
      DO 21 KX=1,3
   21 Z(KX)=(BB(KX)*DR(1)+AA(KX)*DR(2))/(2.*AR)
      IF(ICODC.EQ.0)GO TO 10
      IF(ICOD.EQ.4)GO TO 18
      DO 19 KX=1,3
      B(KX)=-0.5*DFI/Z(ICODC)
      IF(KX.EQ.ICODC)GO TO 19
      B(KX)=-B(KX)*Z(KX)/Z(ICODC)
   19 CONTINUE
      GO TO 20
   18 DO 22 KX=1,3
      IF(KX.EQ.ICODC)GO TO 22
      B(KX)=-DFI*Z(KX)/Z(ICODC)
   22 CONTINUE
      B(ICODC)=-DFI*ALOG(-Z(ICODC))
      GO TO 8
   10 R1=-ABS(CSI(IME))/Z(IME)
      R2=-ABS(CSI(IR))/Z(IR)
      IF(ICOD.EQ.4)GO TO 39
      DO 40 KX=1,3
   40 B(KX)=DFI*(CSI(KX)*(R1+R2)+0.5*Z(KX)*(SIGN(R1**2,R1)+SIGN(R2**2,R2
     *)))
      GO TO 20
   39 DO 41 KX=1,3
   41 B(KX)=DFI*(CSI(KX)*SIGN(1.,R1)*ALOG(-R1/R2)+Z(KX)*(R1+R2))
      GO TO 8
C.... COMPUTE MATRIX DM (BOUNDARY CONSTRAINT EQUATION)
   20 DO 4 J=1,2
      DO 4 K=J,2
      DO 4 I=1,2
    4 EI(J,K,I)=-C1*(2.-D(J,K))*(C4*(D(I,J)*DR(K)+D(K,I)*DR(J))-D(J,K)*D
     *R(I)+2.*DR(I)*DR(J)*DR(K))
      IC=0
      DO 5 IAR=1,3
      DO 5 I=1,2
      DO 5 J=I,2
      IC=IC+1
      DO 5 IL=1,2
    5 DM(IL,IC)=DM(IL,IC)+B(IAR)*EI(I,J,IL)*W(KK,INP)
      GO TO 2
C.... COMPUTE MATRIX DM (STRESSES AT INTERNAL POINTS)
    8 DO 6 I=1,2
      DO 6 J=I,2
      DO 6 K=1,2
      DO 6 L=K,2
      B1=C4*(D(I,K)*D(L,J)+D(J,K)*D(L,I)-D(I,J)*D(K,L)+2.*D(I,J)*DR(K)*D
     *R(L))
      B2=2.*PO*(DR(K)*DR(J)*D(I,L)+DR(L)*DR(J)*D(K,I)+DR(I)*DR(K)*D(J,L)
     *+DR(L)*DR(I)*D(J,K))
      B3=2.*D(K,L)*DR(I)*DR(J)-8.*DR(I)*DR(J)*DR(K)*DR(L)
    6 EII(I,J,K,L)=C3*(B1+B2+B3)*(2.-D(K,L))
      DO 7 IAR=1,3
      IL=0
      DO 7 I=1,2
      DO 7 J=I,2
      IC=(IAR-1)*3
      IL=IL+1
      DO 7 K=1,2
      DO 7 L=K,2
      IC=IC+1
    7 DM(IL,IC)=DM(IL,IC)+B(IAR)*EII(I,J,K,L)*W(KK,INP)
    2 CONTINUE
   38 CONTINUE
C.... SYMMETRY TEST
      IF(ISY.EQ.1.OR.(ISY+ICOD).EQ.8)GO TO 50
      DO 55 I=IIS,IFS
      DO 55 J=1,9
   55 DM(I,J)=-DM(I,J)
   50 RETURN
      END
      SUBROUTINE FENC(C8,C9,C10,C11,CC1,CC2,I,C18,PO)
C.... EXPRESSIONS FOR STRESSES AT BOUNDARY NODES
      COMMON /A1/ ISYM(73),INC(33,2),IDUP(73),ICC(110,3),AR(110),C(33)
      COMMON /A4/ H(3,4),G(3,4),DM(3,9)
C.... MATRIX H
      CO=-C8/C(I)
      H(1,1)=CO*CC2**3
      H(1,2)=-CO*CC1*CC2**2
      H(1,3)=-H(1,1)
      H(1,4)=-H(1,2)
      H(2,1)=H(1,2)
      H(2,2)=CO*CC2*CC1**2
      H(2,3)=-H(1,2)
      H(2,4)=-H(2,2)
      H(3,1)=H(2,2)
      H(3,2)=-CO*CC1**3
      H(3,3)=-H(2,2)
      H(3,4)=-H(3,2)
C.... MATRIX G
      G(1,1)=CC1**3+C10*CC1*CC2**2
      G(1,2)=-CC2*CC1**2+C9*CC2**3
      G(2,1)=CC2**3-C9*CC2*CC1**2
      G(2,2)=CC1**3-C9*CC1*CC2**2
      G(3,1)=-CC1*CC2**2+C9*CC1**3
      G(3,2)=CC2**3+C10*CC2*CC1**2
C.... MATRIX DM
      DM(1,1)=C18*(PO-CC2**2)*CC2**2+1.
      DM(1,2)=2.*C18*CC1*CC2**3
      DM(1,3)=C18*(PO-CC1**2)*CC2**2
      DM(2,1)=-C18*CC1*CC2*(PO-CC2**2)
      DM(2,2)=-2.*C18*CC1**2*CC2**2+1.
      DM(2,3)=-C18*CC1*CC2*(PO-CC1**2)
      DM(3,1)=C18*(PO-CC2**2)*CC1**2
      DM(3,2)=2.*C18*CC2*CC1**3
      DM(3,3)=C18*(PO-CC1**2)*CC1**2+1.
      RETURN
      END
      SUBROUTINE GAUSP(N,NT3,C8)
C.... THIS SUBROUTINE CALCULATES VECTOR XM AND MODIFIES MATRIX DB
C.... USING GAUSS ELIMINATION PROCESS
      COMMON /RW/ IRE,IWR
      COMMON /A2/ A(219,219),DB(76,219)
      COMMON /A3/ IFIP(76),U(76),P(76),S(73,4),SR(219),DSR(219),DSE(219)
     *,B(76),XN(219),EPE(73),DEPE(73),SR3(73),SEQ(73),ASEQ(73),DSR3(73)
      N1=N-1
      DO 100 K=1,N1
      K1=K+1
      C=A(K,K)
      IF(ABS(C))1,1,3
    1 DO 7 J=K1,N
      IF(ABS(A(J,K)))7,7,5
    5 DO 6 L=K,N
      C=A(K,L)
      A(K,L)=A(J,L)
    6 A(J,L)=C
      C=B(K)
      B(K)=B(J)
      B(J)=C
      DO 50 IM=1,NT3
      C=DB(K,IM)
      DB(K,IM)=DB(J,IM)
   50 DB(J,IM)=C
      C=A(K,K)
      GO TO 3
    7 CONTINUE
    8 WRITE(IWR,2)K
    2 FORMAT(/////,20X,'* * * SINGULARITY IN ROW',I5,' * * *')
      GO TO 300
    3 C=A(K,K)
      DO 4 J=K1,N
    4 A(K,J)=A(K,J)/C
      B(K)=B(K)/C
      DO 51 IM=1,NT3
   51 DB(K,IM)=DB(K,IM)/C
      DO 10 I=K1,N
      C=A(I,K)
      DO 9 J=K1,N
    9 A(I,J)=A(I,J)-C*A(K,J)
      B(I)=B(I)-C*B(K)
      DO 10 IM=1,NT3
   10 DB(I,IM)=DB(I,IM)-C*DB(K,IM)
  100 CONTINUE
      IF(ABS(A(N,N)))8,8,101
  101 B(N)=B(N)/A(N,N)
      DO 52 IM=1,NT3
   52 DB(N,IM)=DB(N,IM)/A(N,N)
      DO 200 L=1,N1
      K=N-L
      K1=K+1
      DO 200 J=K1,N
      B(K)=B(K)-A(K,J)*B(J)
      DO 200 IM=1,NT3
  200 DB(K,IM)=DB(K,IM)-A(K,J)*DB(J,IM)
      DO 70 I=1,N
      IF(IFIP(I).EQ.0)GO TO 70
      B(I)=B(I)*C8
      DO 71 J=1,NT3
   71 DB(I,J)=DB(I,J)*C8
   70 CONTINUE
  300 RETURN
      END
      SUBROUTINE ELAST(C11,XLAMB,SY,NT3,NT,PER,NN2,FAC,IPRT,IYIED,SNPHI)
      COMMON /A3/ IFIP(76),U(76),P(76),S(73,4),SR(219),DSR(219),DSE(219)
     *,XM(76),XN(219),EPE(73),DEPE(73),SR3(73),SEQ(73),ASEQ(73),DSR3(73)
      COMMON /A5/ DEVS(4),A(4),D(4)
      DIMENSION STR(4)
      DO 1 I=1,NT
      EPE(I)=0.
    1 SR3(I)=0.
      DO 2 I=1,NT3
      SR(I)=0.
    2 DSR(I)=0.
C.... FIND MOST HIGLY STRESSED NODE OR POINT
      SMAX=0.
      DO 3 I=1,NT
      STR(1)=XN(3*I-2)
      STR(2)=XN(3*I-1)
      STR(3)=XN(3*I)
      STR(4)=C11*(STR(1)+STR(3))
C.... CALCULATE EFFECTIVE STRESS
      CALL EFFST(SE,IYIED,STR,STEFF,THETA,VARJ2,SINT3,SNPHI)
      SEQ(I)=SE
      IF(SE.GT.SMAX)SMAX=SE
    3 CONTINUE
      XLAMB=SY/SMAX
      FAC=XLAMB*PER
      IPRT=1
      IF(XLAMB.GT.1.)IPRT=2
C.... RESULTS FOR LOAD AT FIRST YIELD
      DO 4 I=1,NN2
    4 U(I)=XM(I)*XLAMB
      DO 5 I=1,NT
      SEQ(I)=SEQ(I)*XLAMB
      S(I,1)=XN(3*I-2)*XLAMB
      S(I,2)=XN(3*I-1)*XLAMB
      S(I,3)=XN(3*I)*XLAMB
    5 S(I,4)=C11*(S(I,1)+S(I,3))
      RETURN
      END
      SUBROUTINE SOLVE(TOL,IC,SY,NT3,XLAMB,NT,NN2,HL,SNPHI,C11,FAC,ITR,I
     *YIED,C14,IPL,G)
      COMMON /A2/ DI(219,219),DB(76,219)
      COMMON /A3/ IFIP(76),U(76),P(76),S(73,4),SR(219),DSR(219),DSE(219)
     *,XM(76),XN(219),EPE(73),DEPE(73),SR3(73),SEQ(73),ASEQ(73),DSR3(73)
      COMMON /A5/ DEVS(4),A(4),D(4)
      DIMENSION DSTE(4),STOT(4),STR(4)
      IC=0
C.... CALCULATE ELASTIC STRESS INCREMENT
      DO 1 I=1,NT3
      DSE(I)=0.
      IF(ITR.EQ.1)DSE(I)=XN(I)*FAC
    1 CONTINUE
      DO 21 J=1,NT3
      IF(DSR(J).EQ.0.)GO TO 21
      DO 20 I=1,NT3
   20 DSE(I)=DSE(I)+DI(I,J)*DSR(J)
   21 CONTINUE
C.... LOOP OVER EACH NODE OR POINT
      DO 2 I=1,NT
      IF(ITR.EQ.1)DEPE(I)=0.
      PREYS=SY+EPE(I)*HL
      DSTE(1)=DSE(3*I-2)
      DSTE(2)=DSE(3*I-1)
      DSTE(3)=DSE(3*I)
      DSTE(4)=C11*(DSTE(1)+DSTE(3))
      DO 7 J=1,4
    7 STR(J)=S(I,J)+DSTE(J)
C.... EVALUATE EFFECTIVE STRESS
      CALL EFFST(SE,IYIED,STR,STEFF,THETA,VARJ2,SINT3,SNPHI)
C.... CHECK FOR YIELDING ON THE PREVIOUS ITERATIONS
      IF((SEQ(I)-PREYS).GE.0.)GO TO 3
C.... CHECK FOR YIELDING DURING THIS ITERATION
      ESCUR=SE-PREYS
      IF(ESCUR.LE.0.)GO TO 4
C.... COMPUTE FACTOR R
      RFAC=ESCUR/(SE-SEQ(I))
      GO TO 5
    3 ESCUR=SE-SEQ(I)
      IF(ESCUR.LE.0.)GO TO 4
      RFAC=1.
C.... EVALUATE NUMBER OF STRESS REDUCTION STEPS
    5 MSTEP=ESCUR*8./SY+1.
      ASTEP=MSTEP
      REDUC=1.-RFAC
      DO 6 J=1,4
      STOT(J)=S(I,J)+REDUC*DSTE(J)
    6 DSTE(J)=RFAC*DSTE(J)/ASTEP
      EPCON=0.
C.... LOOP OVER EACH STRESS REDUCTION STEP
      DO 8 IST=1,MSTEP
C.... CALCULATE VECTORS A AND D
      CALL EFFST(SE,IYIED,STOT,STEFF,THETA,VARJ2,SINT3,SNPHI)
      CALL FLOW(ABETA,IYIED,STEFF,THETA,VARJ2,SINT3,HL,SNPHI,C14,IPL,G)
C.... COMPUTE PLASTIC MULTIPLIER
      AGASH=0.
      DO 9 J=1,4
    9 AGASH=AGASH+A(J)*DSTE(J)
      DLAMB=AGASH*ABETA
      IF(DLAMB.LT.0.)DLAMB=0.
C.... COMPUTE ELASTOPLASTIC STRESSES
      BGASH=0.
      DO 10 J=1,4
      BGASH=BGASH+A(J)*STOT(J)
   10 STOT(J)=STOT(J)+DSTE(J)-DLAMB*D(J)
C.... CALCULATE EQUIVALENT PLASTIC STRAIN INCREMENT
      DELEP=DLAMB*BGASH/SE
      EPCON=EPCON+DELEP
C.... UPDATE EQUIVALENT PLASTIC STRAIN
      EPE(I)=EPE(I)+DELEP
    8 CONTINUE
C.... COMPUTE EFFECTIVE STRESS
      CALL EFFST(SE,IYIED,STOT,STEFF,THETA,VARJ2,SINT3,SNPHI)
C.... COMPUTE EQUIVALENT YIELD STRESS
      CURYS=SY+EPE(I)*HL
C.... SCALE DOWN STRESSES TO YIELD SURFACE
      BRING=1.
      IF(SE.GT.CURYS)BRING=CURYS/SE
      DO 11 J=1,4
   11 S(I,J)=BRING*STOT(J)
      SEQ(I)=BRING*SE
C.... INDEX OF CONVERGENCE
      IF(EPCON.GT.DEPE(I)*TOL)IC=IC+1
      DEPE(I)=DEPE(I)+EPCON
C.... COMPUTE RESIDUAL STRESSES
      DO 14 J=1,3
      JJ=3-J
      DSR(3*I-JJ)=STR(J)-S(I,J)
   14 SR(3*I-JJ)=SR(3*I-JJ)+DSR(3*I-JJ)
      SR3(I)=SR3(I)+STR(4)-S(I,4)
      GO TO 2
C.... COMPUTE STRESSES FOR ELASTIC NODES OR POINTS
    4 DO 13 J=1,4
      JJ=3-J
      IF(J.LT.4)DSR(3*I-JJ)=0.
   13 S(I,J)=STR(J)
      SEQ(I)=SE
    2 CONTINUE
      IF(IC.NE.0)GO TO 12
C.... COMPUTE BOUNDARY UNKNOWNS
      DO 15 I=1,NN2
   15 U(I)=XM(I)*XLAMB
      DO 22 J=1,NT3
      IF(SR(J).EQ.0.)GO TO 22
      DO 23 I=1,NN2
   23 U(I)=U(I)+DB(I,J)*SR(J)
   22 CONTINUE
   12 RETURN
      END
      SUBROUTINE SOLVE1(TOL,IC,SY,NT3,XLAMB,NT,NN2,HL,SNPHI,C11,FAC,ITR,
     *IYIED,C14,IPL,G)
      COMMON /A2/ DI(219,219),DB(76,219)
      COMMON /A3/ IFIP(76),U(76),P(76),S(73,4),SR(219),DSR(219),DSE(219)
     *,XM(76),XN(219),EPE(73),DEPE(73),SR3(73),SEQ(73),ASEQ(73),DSR3(73)
      COMMON /A5/ DEVS(4),A(4),D(4)
      DIMENSION DSTE(4),STOT(4),STR(4)
      IC=0
C.... CALCULATE ELASTIC STRESS
      DO 1 I=1,NT3
    1 DSE(I)=XLAMB*XN(I)
      DO 20 J=1,NT3
      EPLAS=SR(J)+DSR(J)
      IF(EPLAS.EQ.0.)GO TO 20
      DO 21 I=1,NT3
   21 DSE(I)=DSE(I)+DI(I,J)*EPLAS
   20 CONTINUE
C.... LOOP OVER EACH NODE OR POINT
      DO 2 I=1,NT
      PREYS=SY+EPE(I)*HL
      EPCON=0.
C.... COMPUTE ELASTIC STRESS INCREMENT
      DSTE(1)=DSE(3*I-2)-S(I,1)-SR(3*I-2)
      DSTE(2)=DSE(3*I-1)-S(I,2)-SR(3*I-1)
      DSTE(3)=DSE(3*I)-S(I,3)-SR(3*I)
      DSTE(4)=C11*(DSTE(1)+DSTE(3))
      DO 7 J=1,4
    7 STR(J)=S(I,J)+DSTE(J)
C.... EVALUATE EFFECTIVE STRESS
      CALL EFFST(SE,IYIED,STR,STEFF,THETA,VARJ2,SINT3,SNPHI)
C.... CHECK FOR YIELDING ON THE PREVIOUS INCREMENTS
      IF((SEQ(I)-PREYS).GE.0.)GO TO 3
C.... CHECK FOR YIELDING DURING THIS ITERATION
      ESCUR=SE-PREYS
      IF(ESCUR.LE.0.)GO TO 4
C.... COMPUTE FACTOR R
      RFAC=ESCUR/(SE-SEQ(I))
      GO TO 5
    3 ESCUR=SE-SEQ(I)
      IF(ESCUR.LE.0.)GO TO 4
      RFAC=1.
C.... EVALUATE NUMBER OF STRESS REDUCTION STEPS
    5 MSTEP=ESCUR*8./SY+1.
      ASTEP=MSTEP
      REDUC=1.-RFAC
      DO 6 J=1,4
      STOT(J)=S(I,J)+REDUC*DSTE(J)
    6 DSTE(J)=RFAC*DSTE(J)/ASTEP
C.... LOOP OVER EACH STRESS REDUCTION STEP
      DO 8 IST=1,MSTEP
C.... CALCULATE VECTORS A AND D
      CALL EFFST(SE,IYIED,STOT,STEFF,THETA,VARJ2,SINT3,SNPHI)
      CALL FLOW(ABETA,IYIED,STEFF,THETA,VARJ2,SINT3,HL,SNPHI,C14,IPL,G)
C.... COMPUTE PLASTIC MULTIPLIER
      AGASH=0.
      DO 9 J=1,4
    9 AGASH=AGASH+A(J)*DSTE(J)
      DLAMB=AGASH*ABETA
      IF(DLAMB.LT.0.)DLAMB=0.
C.... COMPUTE ELASTOPLASTIC STRESSES
      BGASH=0.
      DO 10 J=1,4
      BGASH=BGASH+A(J)*STOT(J)
   10 STOT(J)=STOT(J)+DSTE(J)-DLAMB*D(J)
C.... CALCULATE EQUIVALENT PLASTIC STRAIN INCREMENT
      EPCON=EPCON+DLAMB*BGASH/SE
    8 CONTINUE
C.... COMPUTE EFFECTIVE STRESS
      CALL EFFST(SE,IYIED,STOT,STEFF,THETA,VARJ2,SINT3,SNPHI)
C.... COMPUTE EQUIVALENT YIELD STRESS
      CURYS=SY+(EPE(I)+EPCON)*HL
C.... SCALE DOWN STRESSES TO YIELD SURFACE
      BRING=1.
      IF(SE.GT.CURYS)BRING=CURYS/SE
      DO 11 J=1,4
   11 STOT(J)=BRING*STOT(J)
      ASEQ(I)=BRING*SE
C.... COMPUTE RESIDUAL STRESS INCREMENTS
      DO 14 J=1,3
      JJ=3-J
   14 DSR(3*I-JJ)=STR(J)-STOT(J)
      DSR3(I)=STR(4)-STOT(4)
      GO TO 22
C.... ELASTIC NODES OR POINTS
    4 DO 13 J=1,3
      JJ=3-J
   13 DSR(3*I-JJ)=0.
      DSR3(I)=0.
      ASEQ(I)=SE
C.... INDEX OF CONVERGENCE
   22 IF(ABS(EPCON-DEPE(I)).GT.DEPE(I)*TOL)IC=IC+1
      DEPE(I)=EPCON
    2 CONTINUE
      IF(IC.NE.0)GO TO 12
      DO 30 I=1,NT
C.... UPDATE EQUIVALENT PLASTIC STRAIN
      EPE(I)=EPE(I)+DEPE(I)
C.... UPDATE EFFECTIVE STRESS
      SEQ(I)=ASEQ(I)
C.... COMPUTE RESIDUAL STRESSES AND TRUE STRESSES
      DO 31 J=1,3
      JJ=3-J
      SR(3*I-JJ)=SR(3*I-JJ)+DSR(3*I-JJ)
   31 S(I,J)=DSE(3*I-JJ)-SR(3*I-JJ)
      SR3(I)=SR3(I)+DSR3(I)
   30 S(I,4)=C11*(DSE(3*I-2)+DSE(3*I))-SR3(I)
C.... COMPUTE BOUNDARY UNKNOWNS
      DO 15 I=1,NN2
   15 U(I)=XM(I)*XLAMB
      DO 23 J=1,NT3
      IF(SR(J).EQ.0.)GO TO 23
      DO 24 I=1,NN2
   24 U(I)=U(I)+DB(I,J)*SR(J)
   23 CONTINUE
   12 RETURN
      END
      SUBROUTINE EFFST(SE,IYIED,SA,STEFF,THETA,VARJ2,SINT3,SNPHI)
C.... THIS SUBROUTINE CALCULATES THE DEVIATORIC STRESSES,THE EFFECTIVE
C.... STRESS AND THE STRESS INVARIANTS
      COMMON /A5/ DEVS(4),A(4),D(4)
      DIMENSION SA(4)
      ROOT3=1.73205080757
      SMEAN=(SA(1)+SA(3)+SA(4))/3.
      DEVS(1)=SA(1)-SMEAN
      DEVS(2)=SA(2)
      DEVS(3)=SA(3)-SMEAN
      DEVS(4)=SA(4)-SMEAN
      VARJ2=DEVS(2)*DEVS(2)+0.5*(DEVS(1)*DEVS(1)+DEVS(3)*DEVS(3)+DEVS(4)
     **DEVS(4))
      VARJ3=DEVS(4)*(DEVS(4)*DEVS(4)-VARJ2)
      STEFF=SQRT(VARJ2)
      SINT3=-3.*ROOT3*VARJ3/(2.*VARJ2*STEFF)
      IF(ABS(SINT3).GT.1.)SINT3=SIGN(1.,SINT3)
      THETA=ASIN(SINT3)/3.
      GO TO (1,2,3,4,5), IYIED
C.... TRESCA
    1 SE=2.*COS(THETA)*STEFF
      RETURN
C.... VON MISES
    2 SE=ROOT3*STEFF
      RETURN
C.... MOHR-COULOMB
    3 SE=SMEAN*SNPHI+STEFF*(COS(THETA)-SIN(THETA)*SNPHI/ROOT3)
      RETURN
C.... DRUCKER-PRAGER
    4 SE=6.*SMEAN*SNPHI/(ROOT3*(3.-SNPHI))+STEFF
      RETURN
    5 SE=3.*SMEAN*SNPHI/(ROOT3*SQRT(3.+SNPHI**2))+STEFF
      RETURN
      END
      SUBROUTINE FLOW(ABETA,IYIED,STEFF,THETA,VARJ2,SINT3,HL,SNPHI,C14,I
     *PL,G)
C.... THIS SUBROUTINE CALCULATES VECTORS A AND D
      COMMON /A5/ DEVS(4),A(4),D(4)
      DIMENSION A1(4),A2(4),A3(4)
      ROOT3=1.73205080757
C.... VECTOR A1
      A1(1)=1.
      A1(2)=0.
      A1(3)=1.
      A1(4)=1.
C.... VECTOR A2
      A2(1)=DEVS(1)/(2.*STEFF)
      A2(2)=DEVS(2)/STEFF
      A2(3)=DEVS(3)/(2.*STEFF)
      A2(4)=DEVS(4)/(2.*STEFF)
C.... VECTOR A3
      A3(1)=DEVS(3)*DEVS(4)+VARJ2/3.
      A3(2)=-2.*DEVS(2)*DEVS(4)
      A3(3)=DEVS(1)*DEVS(4)+VARJ2/3.
      A3(4)=DEVS(1)*DEVS(3)-DEVS(2)*DEVS(2)+VARJ2/3.
      GO TO (1,2,3,4,5),IYIED
C.... TRESCA
    1 CONS1=0.
      IF(ABS(THETA)*57.29577951308.LT.29.)GO TO 20
      CONS2=ROOT3
      CONS3=0.
      GO TO 40
   20 SINTH=SIN(THETA)
      CONS2=2.*(COS(THETA)+SINTH*TAN(3.*THETA))
      CONS3=ROOT3*SINTH/(VARJ2*COS(3.*THETA))
      GO TO 40
C.... VON MISES
    2 CONS1=0.
      CONS2=ROOT3
      CONS3=0.
      GO TO 40
C.... MOHR-COULOMB
    3 CONS1=SNPHI/3.
      IF(ABS(THETA)*57.29577951308.LT.29.)GO TO 30
      PLUMI=1.
      IF(THETA.GT.0.)PLUMI=-1.
      CONS2=0.5*(ROOT3+PLUMI*SNPHI/ROOT3)
      CONS3=0.
      GO TO 40
   30 COSTH=COS(THETA)
      TANTH=TAN(THETA)
      TANT3=TAN(3.*THETA)
      CONS2=COSTH*((1.+TANTH*TANT3)+SNPHI*(TANT3-TANTH)/ROOT3)
      CONS3=(ROOT3*SIN(THETA)+SNPHI*COSTH)/(2.*VARJ2*COS(3.*THETA))
      GO TO 40
C.... DRUCKER-PRAGER
    4 CONS1=2.*SNPHI/(ROOT3*(3.-SNPHI))
      CONS2=1.
      CONS3=0.
      GO TO 40
    5 CONS1=SNPHI/(ROOT3*SQRT(3.+SNPHI**2))
      CONS2=1.
      CONS3=0.
C.... VECTOR A
   40 DO 50 J=1,4
   50 A(J)=CONS1*A1(J)+CONS2*A2(J)+CONS3*A3(J)
C.... VECTOR D
      XML=C14*(A(1)+A(3))
      IF(IPL.EQ.2)XML=XML+C14*A(4)
      D(1)=2.*G*(A(1)+XML)
      D(2)=G*A(2)
      D(3)=2.*G*(A(3)+XML)
      D(4)=2.*G*(A(4)+XML)
      IF(IPL.EQ.1)D(4)=0.
C.... COMPUTE PART OF PLASTIC MULTIPLIER
      DENOM=HL
      DO 80 J=1,4
   80 DENOM=DENOM+A(J)*D(J)
      ABETA=1./DENOM
      RETURN
      END
      SUBROUTINE OUTPT(NN,NT,XLAMB,IPRT,ITR,G,C13)
C.... OUTPUT RESULTS
      COMMON /RW/ IRE,IWR
      COMMON /A3/ IFIP(76),U(76),P(76),S(73,4),SR(219),DSR(219),DSE(219)
     *,XM(76),XN(219),EPE(73),DEPE(73),SR3(73),SEQ(73),ASEQ(73),DSR3(73)
      DIMENSION UA(2),PA(2),EP(3)
      GO TO (1,2,3),IPRT
    2 WRITE(IWR,4)
    4 FORMAT(///,30X,'* * * NO PLASTIC STRAIN OCCURS IN THIS ANALYSIS *
     ** *')
    1 WRITE(IWR,5)
    5 FORMAT(///,30X,'L O A D  A T  F I R S T  Y I E L D')
      IPRT=3
    3 WRITE(IWR,6)XLAMB,ITR
    6 FORMAT(///,30X,'* LOAD FACTOR =',F12.7,' * NO. ITERATIONS =',I5,'
     **',///,15X,'BOUNDARY DISPLACEMENTS AND TRACTIONS',//,12X,'NODE',14
     *X,'U',15X,'V',14X,'TX',14X,'TY',/)
C.... BOUNDARY DISPLACEMENTS AND TRACTIONS
      DO 7 I=1,NN
      DO 10 J=1,2
      JA=2*I+J-2
      IF(IFIP(JA)-1)8,9,9
    8 UA(J)=U(JA)
      PA(J)=P(JA)*XLAMB
      GO TO 10
    9 UA(J)=P(JA)*XLAMB
      PA(J)=U(JA)
   10 CONTINUE
    7 WRITE(IWR,11)I,UA(1),UA(2),PA(1),PA(2)
   11 FORMAT(10X,I5,5X,F15.4,1X,F15.4,1X,F15.4,1X,F15.4)
C.... STRESSES AND PLASTIC STRAINS AT NODES AND POINTS
      WRITE(IWR,12)
   12 FORMAT(//,15X,'STRESSES AND PLASTIC STRAINS AT NODES AND INTERNAL
     * POINTS',//,11X,'NO/PT',13X,'SX',14X,'SXY',13X,'SY',14X,'EPX',12X
     *,'EPXY',13X,'EPY',/)
      DO 14 I=1,NT
      EP(1)=(SR(3*I-2)-C13*(SR(3*I-2)+SR(3*I)+SR3(I)))/(2.*G)
      EP(2)=SR(3*I-1)/(2.*G)
      EP(3)=(SR(3*I)-C13*(SR(3*I-2)+SR(3*I)+SR3(I)))/(2.*G)
   14 WRITE(IWR,13)I,S(I,1),S(I,2),S(I,3),EP(1),EP(2),EP(3)
   13 FORMAT(10X,I5,5X,F15.4,1X,F15.4,1X,F15.4,9X,E14.7,2X,E14.7,2X,E14.
     *7)
      RETURN
      END
