      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(*),DPRED(*),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
C
      DOUBLE PRECISION STRESS,STATEV,DDSDDE,SSE,SPD,SCD
      DOUBLE PRECISION RPL,DDSDDT,DRPLDE,DRPLDT
      DOUBLE PRECISION STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP
      DOUBLE PRECISION PREDEF,DPRED,PROPS,COORDS,DROT,PNEWDT
      DOUBLE PRECISION CELENT,DFGRD0,DFGRD1
C
      DOUBLE PRECISION ZERO,ONE,TWO,THREE,HALF
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0,
     1           THREE=3.D0, HALF=0.5D0)
C
      DOUBLE PRECISION KF_CIRC,KF_LONG,MU0,KAPPA_B
      DOUBLE PRECISION F(3,3),BMAT(3,3),DELTA(3,3)
      DOUBLE PRECISION SIGMA_PAS(3,3),SIGMA_ANISO_C(3,3)
      DOUBLE PRECISION SIGMA_ANISO_L(3,3),SIGMA_TOT(3,3)
      DOUBLE PRECISION C_PAS(3,3,3,3),C_ANISO_C(3,3,3,3)
      DOUBLE PRECISION C_ANISO_L(3,3,3,3),C_TOT(3,3,3,3)
      DOUBLE PRECISION JAC,LOGJ,P_SCALAR
      DOUBLE PRECISION ECIRC_REF(3),ELONG_REF(3)
      DOUBLE PRECISION ACIRC_CUR(3),ALONG_CUR(3)
      DOUBLE PRECISION ECIRC_CUR(3),ELONG_CUR(3)
      DOUBLE PRECISION LAMBDA_C,LAMBDA_L,NORM_AC,NORM_AL
      DOUBLE PRECISION I4C,I4L,BR_C,BR_L
      INTEGER I1,I2,I3,I4
      DOUBLE PRECISION DET3X3
C
C     PROPS:
C     1 = KF_CIRC
C     2 = KF_LONG
C     3 = MU0
C
      KF_CIRC = PROPS(1)
      KF_LONG = PROPS(2)
      MU0     = PROPS(3)
      KAPPA_B = 100.D0*MU0
C
      IF (NDI .NE. 3 .OR. NSHR .NE. 3 .OR. NTENS .NE. 6) THEN
         DO I1=1,NTENS
            STRESS(I1)=ZERO
            DDSDDT(I1)=ZERO
            DRPLDE(I1)=ZERO
            DO I2=1,NTENS
               DDSDDE(I1,I2)=ZERO
            END DO
         END DO
         RPL = ZERO
         DRPLDT = ZERO
         RETURN
      END IF
C
      DO I1=1,3
         DO I2=1,3
            F(I1,I2)     = DFGRD1(I1,I2)
            BMAT(I1,I2)  = ZERO
            DELTA(I1,I2) = ZERO
         END DO
         DELTA(I1,I1)=ONE
      END DO
C
      JAC = DET3X3(F)
      IF (JAC .LE. 1.D-12) JAC = 1.D-12
      LOGJ = DLOG(JAC)
C
      DO I1=1,3
         DO I2=1,3
            DO I3=1,3
               BMAT(I1,I2)=BMAT(I1,I2)+F(I1,I3)*F(I2,I3)
            END DO
         END DO
      END DO
C
C     Isotropic compressible neo-Hookean baseline
      P_SCALAR = KAPPA_B*LOGJ - MU0
      DO I1=1,3
         DO I2=1,3
            SIGMA_PAS(I1,I2) = (MU0/JAC)*BMAT(I1,I2)
         END DO
         SIGMA_PAS(I1,I1) = SIGMA_PAS(I1,I1) + P_SCALAR/JAC
      END DO
C
C     Reference directions:
C     e2 = circumferential, e3 = longitudinal
      ECIRC_REF(1)=ZERO
      ECIRC_REF(2)=ONE
      ECIRC_REF(3)=ZERO
      ELONG_REF(1)=ZERO
      ELONG_REF(2)=ZERO
      ELONG_REF(3)=ONE
C
      DO I1=1,3
         ACIRC_CUR(I1)=F(I1,1)*ECIRC_REF(1)
     1              + F(I1,2)*ECIRC_REF(2)
     2              + F(I1,3)*ECIRC_REF(3)
         ALONG_CUR(I1)=F(I1,1)*ELONG_REF(1)
     1              + F(I1,2)*ELONG_REF(2)
     2              + F(I1,3)*ELONG_REF(3)
      END DO
C
      NORM_AC = DSQRT(ACIRC_CUR(1)*ACIRC_CUR(1)
     1               + ACIRC_CUR(2)*ACIRC_CUR(2)
     2               + ACIRC_CUR(3)*ACIRC_CUR(3))
      NORM_AL = DSQRT(ALONG_CUR(1)*ALONG_CUR(1)
     1               + ALONG_CUR(2)*ALONG_CUR(2)
     2               + ALONG_CUR(3)*ALONG_CUR(3))
C
      IF (NORM_AC .LE. 1.D-12) NORM_AC = 1.D-12
      IF (NORM_AL .LE. 1.D-12) NORM_AL = 1.D-12
C
      LAMBDA_C = NORM_AC
      LAMBDA_L = NORM_AL
      I4C = LAMBDA_C*LAMBDA_C
      I4L = LAMBDA_L*LAMBDA_L
C
      DO I1=1,3
         ECIRC_CUR(I1)=ACIRC_CUR(I1)/LAMBDA_C
         ELONG_CUR(I1)=ALONG_CUR(I1)/LAMBDA_L
      END DO
C
      BR_C = I4C - ONE
      BR_L = I4L - ONE
      IF (BR_C .LT. ZERO) BR_C = ZERO
      IF (BR_L .LT. ZERO) BR_L = ZERO
C
      DO I1=1,3
         DO I2=1,3
            SIGMA_ANISO_C(I1,I2) =
     1         (THREE*KF_CIRC/JAC)*I4C*BR_C*BR_C*
     2         ECIRC_CUR(I1)*ECIRC_CUR(I2)
C
            SIGMA_ANISO_L(I1,I2) =
     1         (THREE*KF_LONG/JAC)*I4L*BR_L*BR_L*
     2         ELONG_CUR(I1)*ELONG_CUR(I2)
C
            SIGMA_TOT(I1,I2)=SIGMA_PAS(I1,I2)
     1                      + SIGMA_ANISO_C(I1,I2)
     2                      + SIGMA_ANISO_L(I1,I2)
         END DO
      END DO
C
      DO I1=1,3
         DO I2=1,3
            DO I3=1,3
               DO I4=1,3
                  C_PAS(I1,I2,I3,I4) =
     1              (KAPPA_B/JAC)*DELTA(I1,I2)*DELTA(I3,I4)
     2            - (P_SCALAR/JAC)*
     3              ( DELTA(I1,I3)*DELTA(I2,I4)
     4              + DELTA(I1,I4)*DELTA(I2,I3) )
C
                  C_PAS(I1,I2,I3,I4) = C_PAS(I1,I2,I3,I4)
     1            + HALF*( SIGMA_PAS(I1,I3)*DELTA(I2,I4)
     2            +        SIGMA_PAS(I1,I4)*DELTA(I2,I3)
     3            +        SIGMA_PAS(I2,I3)*DELTA(I1,I4)
     4            +        SIGMA_PAS(I2,I4)*DELTA(I1,I3) )
C
                  C_ANISO_C(I1,I2,I3,I4) =
     1              (12.D0*KF_CIRC*I4C*I4C*BR_C/JAC)
     2            * ECIRC_CUR(I1)*ECIRC_CUR(I2)
     3            * ECIRC_CUR(I3)*ECIRC_CUR(I4)
C
                  C_ANISO_L(I1,I2,I3,I4) =
     1              (12.D0*KF_LONG*I4L*I4L*BR_L/JAC)
     2            * ELONG_CUR(I1)*ELONG_CUR(I2)
     3            * ELONG_CUR(I3)*ELONG_CUR(I4)
C
                  C_TOT(I1,I2,I3,I4)=C_PAS(I1,I2,I3,I4)
     1                                 + C_ANISO_C(I1,I2,I3,I4)
     2                                 + C_ANISO_L(I1,I2,I3,I4)
               END DO
            END DO
         END DO
      END DO

      DO I1=1,3
         DO I2=1,3
            DO I3=1,3
               DO I4=1,3
                  C_TOT(I1,I2,I3,I4) = C_TOT(I1,I2,I3,I4)
     1            + HALF*(
     2              (SIGMA_ANISO_C(I1,I3)+SIGMA_ANISO_L(I1,I3))
     3              *DELTA(I2,I4)
     4            + (SIGMA_ANISO_C(I1,I4)+SIGMA_ANISO_L(I1,I4))
     5              *DELTA(I2,I3)
     6            + (SIGMA_ANISO_C(I2,I3)+SIGMA_ANISO_L(I2,I3))
     7              *DELTA(I1,I4)
     8            + (SIGMA_ANISO_C(I2,I4)+SIGMA_ANISO_L(I2,I4))
     9              *DELTA(I1,I3) )
               END DO
            END DO
         END DO
      END DO
C
      CALL VOIGT_MAP_3D(C_TOT,DDSDDE)
C
      STRESS(1)=SIGMA_TOT(1,1)
      STRESS(2)=SIGMA_TOT(2,2)
      STRESS(3)=SIGMA_TOT(3,3)
      STRESS(4)=SIGMA_TOT(1,2)
      STRESS(5)=SIGMA_TOT(1,3)
      STRESS(6)=SIGMA_TOT(2,3)
C
      SSE = ZERO
      SPD = ZERO
      SCD = ZERO
      RPL = ZERO
      DRPLDT = ZERO
      DO I1=1,NTENS
         DDSDDT(I1)=ZERO
         DRPLDE(I1)=ZERO
      END DO
C
      IF (NSTATV .GE. 1) STATEV(1)=I4C
      IF (NSTATV .GE. 2) STATEV(2)=I4L
C
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION DET3X3(A)
      INCLUDE 'ABA_PARAM.INC'
      DOUBLE PRECISION A(3,3)
      DET3X3 = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
     1       -A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
     2       +A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      RETURN
      END
C
      SUBROUTINE VOIGT_MAP_3D(C4,DDSDDE)
      INCLUDE 'ABA_PARAM.INC'
      DOUBLE PRECISION C4(3,3,3,3),DDSDDE(6,6)
      INTEGER IV(6),JV(6),K1,K2,A1,A2,A3,A4
      IV(1)=1
      JV(1)=1
      IV(2)=2
      JV(2)=2
      IV(3)=3
      JV(3)=3
      IV(4)=1
      JV(4)=2
      IV(5)=1
      JV(5)=3
      IV(6)=2
      JV(6)=3
      DO K1=1,6
         DO K2=1,6
            A1 = IV(K1)
            A2 = JV(K1)
            A3 = IV(K2)
            A4 = JV(K2)
            DDSDDE(K1,K2) = C4(A1,A2,A3,A4)
         END DO
      END DO
      RETURN
      END
