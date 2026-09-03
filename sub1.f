C ============================================================
C Model 1 UMATHT: fibrosis scalar transport
C
C The Abaqus temperature degree of freedom is used as the
C fibrosis scalar m.  This routine supplies the storage term,
C diffusive flux, and diffusion tangent for the thermal analogy.
C ============================================================
      SUBROUTINE UMATHT(U,DUDT,DUDG,FLUX,DFDT,DFDG,
     1 STATEV,TEMP,DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED,
     2 CMNAME,NTGRD,NSTATV,PROPS,NPROPS,COORDS,PNEWDT,
     3 NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      IMPLICIT NONE

      INTEGER NTGRD,NSTATV,NPROPS,NOEL,NPT,LAYER,KSPT,KSTEP,KINC

      CHARACTER*80 CMNAME
      DOUBLE PRECISION U,DUDT,DUDG(NTGRD),FLUX(NTGRD),DFDT(NTGRD),
     1 DFDG(NTGRD,NTGRD)
      DOUBLE PRECISION STATEV(NSTATV),TEMP,DTEMP,DTEMDX(NTGRD),
     1 TIME(2),DTIME
      DOUBLE PRECISION PREDEF(*),DPRED(*),PROPS(NPROPS),
     1 COORDS(3),PNEWDT

      DOUBLE PRECISION M,C_M,D_M
      DOUBLE PRECISION ZERO,ONE
      INTEGER I1,I2

      PARAMETER (ZERO=0.D0, ONE=1.D0)

C PROPS(1) = c_m, temporal scaling parameter for m.
C PROPS(2) = D_m, isotropic fibrosis diffusivity.

C     Symbol names follow the manuscript equations where possible.


      C_M = PROPS(1)
      D_M = PROPS(2)


C Current scalar field value used by Abaqus in this increment.
C TEMP is the previous value and DTEMP is the current increment.
      M = TEMP + DTEMP



C Storage term: U = c_m m, with tangent dU/dm = c_m.
      U    = C_M*M
      DUDT = C_M

      DO I1=1,NTGRD
         DUDG(I1) = ZERO
      END DO

      DO I1=1,NTGRD

C Diffusive flux: q = -D_m grad(m).
         FLUX(I1) = -D_M*DTEMDX(I1)
         DFDT(I1) = ZERO
         DO I2=1,NTGRD
            DFDG(I1,I2) = ZERO
         END DO
         DFDG(I1,I1) = -D_M
      END DO

      RETURN
      END




C ============================================================
C Model 1 UMAT: fibrosis-driven growth and passive mechanics
C
C This routine supplies the logistic source term, updates the
C radial growth multiplier theta_h, computes the elastic response
C using F = F_e F_h, and returns stress/tangent information.
C ============================================================
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

      IMPLICIT NONE

      INTEGER NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER,KSPT,KINC

      CHARACTER*80 CMNAME
      DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),SSE,SPD,SCD
      DOUBLE PRECISION RPL,DDSDDT(NTENS),DRPLDE(NTENS),DRPLDT
      DOUBLE PRECISION STRAN(NTENS),DSTRAN(NTENS),TIME(2),
     1 DTIME,TEMP,DTEMP
      DOUBLE PRECISION PREDEF(*),DPRED(*),PROPS(NPROPS),
     1 COORDS(3),DROT(3,3),PNEWDT
      DOUBLE PRECISION CELENT,DFGRD0(3,3),DFGRD1(3,3)
      INTEGER JSTEP(4)

      DOUBLE PRECISION ZERO,ONE,TWO,THREE,FOUR,HALF
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0,
     1           THREE=3.D0, FOUR=4.D0, HALF=0.5D0)

C Main material parameters for Model 1:
C   PROPS(1) = mu_0, healthy shear modulus
C   PROPS(2) = k_c, circumferential fiber stiffness
C   PROPS(3) = k_l, longitudinal fiber stiffness
C   PROPS(4) = gamma_g, growth rate parameter
C   PROPS(5) = kappa_fac, bulk-to-shear modulus factor
C   PROPS(6) = s_0, logistic source rate
C
C State variables used by Model 1:
C   STATEV(1) = theta_h, radial growth multiplier

      DOUBLE PRECISION MU_0,K_C,K_L
      DOUBLE PRECISION GAMMA_G, KAPPA_FAC
      DOUBLE PRECISION MU,KAPPA
      DOUBLE PRECISION M, GAMMA_G_M, GAMMA_G_DT
      DOUBLE PRECISION THETA_H_OLD, THETA_H_NEW

      DOUBLE PRECISION F_TOTAL(3,3),F_E(3,3)
      DOUBLE PRECISION F_H_INV(3,3),F_H(3,3)
      DOUBLE PRECISION B_E(3,3)
      DOUBLE PRECISION IDENTITY(3,3),DELTA(3,3)
      DOUBLE PRECISION N_R(3),THETA_H_MINUS_ONE

      DOUBLE PRECISION SIGMA_ISO(3,3),SIGMA(3,3)
      DOUBLE PRECISION SIGMA_C(3,3),SIGMA_L(3,3)

      DOUBLE PRECISION C_ISO(3,3,3,3),C_TOTAL(3,3,3,3)
      DOUBLE PRECISION C_C(3,3,3,3)
      DOUBLE PRECISION C_L(3,3,3,3)
      DOUBLE PRECISION C_GEOM(3,3,3,3)

      DOUBLE PRECISION D_SIGMA_D_M(3,3),N_R_CUR(3)
      DOUBLE PRECISION J_E, LOG_J_E
      DOUBLE PRECISION PRESSURE_TERM, D_PRESSURE_D_M

      DOUBLE PRECISION N_C(3),N_L(3)
      DOUBLE PRECISION F_E_N_C(3),F_E_N_L(3)
      DOUBLE PRECISION N_C_CUR(3),N_L_CUR(3)

      DOUBLE PRECISION LAMBDA_C_E,LAMBDA_L_E
      DOUBLE PRECISION I4C_E,I4L_E
      DOUBLE PRECISION BRACKET_C,BRACKET_L

      DOUBLE PRECISION S0
      DOUBLE PRECISION LAMBDA_R_E

      INTEGER I1,I2,I3,I4
      DOUBLE PRECISION DET3X3


C This constitutive implementation is written for 3D continuum elements
C with six Voigt stress/tangent components.
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

      RPL    = ZERO
      DRPLDT = ZERO
      DO I1=1,NTENS
         DDSDDT(I1)=ZERO
         DRPLDE(I1)=ZERO
      END DO


C Read material parameters from the Abaqus input file in manuscript order.
      MU_0      = PROPS(1)
      K_C       = PROPS(2)
      K_L       = PROPS(3)
      GAMMA_G   = PROPS(4)
      KAPPA_FAC = PROPS(5)
      S0        = PROPS(6)

      MU    = MU_0
      KAPPA = KAPPA_FAC*MU_0


C Current fibrosis value m at the integration point.
      M         = TEMP + DTEMP
      GAMMA_G_M = GAMMA_G*M


C Logistic reaction source S = s_0 m (1 - m) and its tangent.
      RPL    = S0*M*(ONE - M)
      DRPLDT = S0*(ONE - TWO*M)


C Retrieve and update the radial growth multiplier theta_h.
C The exponential form is the exact update for constant m over DTIME.
      THETA_H_OLD = ONE
      IF (NSTATV .GE. 1) THETA_H_OLD = STATEV(1)

      THETA_H_NEW = THETA_H_OLD*DEXP(GAMMA_G_M*DTIME)

      DO I1=1,3
         DO I2=1,3
            F_TOTAL(I1,I2)    = DFGRD1(I1,I2)
            F_E(I1,I2)  = ZERO
            IDENTITY(I1,I2)        = ZERO
            DELTA(I1,I2)           = ZERO
            B_E(I1,I2)   = ZERO
            F_H_INV(I1,I2)       = ZERO
            F_H(I1,I2)    = ZERO
            SIGMA_ISO(I1,I2)       = ZERO
            SIGMA(I1,I2)     = ZERO
            SIGMA_C(I1,I2) = ZERO
            SIGMA_L(I1,I2) = ZERO
            D_SIGMA_D_M(I1,I2) = ZERO
         END DO
         IDENTITY(I1,I1) = ONE
         DELTA(I1,I1)    = ONE
      END DO


C Local material basis: local-1 is radial, local-2 is circumferential,
C and local-3 is longitudinal under the cylindrical orientation.
      N_R(1) = ONE
      N_R(2) = ZERO
      N_R(3) = ZERO
      THETA_H_MINUS_ONE = THETA_H_NEW - ONE


C Growth tensor F_h = I + (theta_h - 1) N_r x N_r.
      DO I1=1,3
         DO I2=1,3
            F_H(I1,I2) = IDENTITY(I1,I2)
     1                           + THETA_H_MINUS_ONE*N_R(I1)*N_R(I2)
         END DO
      END DO


C Elastic deformation gradient F_e = F F_h^{-1}.
      CALL M3INV(F_H,F_H_INV)

      CALL M3MULT(F_TOTAL,F_H_INV,F_E)


C Elastic Jacobian and elastic left Cauchy-Green tensor b_e.
      J_E = DET3X3(F_E)
      LOG_J_E = DLOG(J_E)

      DO I1=1,3
         DO I2=1,3
            DO I3=1,3
               B_E(I1,I2) = B_E(I1,I2)
     1                             + F_E(I1,I3)
     2                             * F_E(I2,I3)
            END DO
         END DO
      END DO


C Reference fiber families: circumferential N_c and longitudinal N_l.
      N_C(1) = ZERO
      N_C(2) = ONE
      N_C(3) = ZERO

      N_L(1) = ZERO
      N_L(2) = ZERO
      N_L(3) = ONE

      DO I1=1,3
         F_E_N_C(I1) = F_E(I1,1)*N_C(1)
     1                    + F_E(I1,2)*N_C(2)
     2                    + F_E(I1,3)*N_C(3)

         F_E_N_L(I1) = F_E(I1,1)*N_L(1)
     1                    + F_E(I1,2)*N_L(2)
     2                    + F_E(I1,3)*N_L(3)
      END DO


C Elastic fiber stretches and invariants.
      LAMBDA_C_E = DSQRT(F_E_N_C(1)*F_E_N_C(1)
     1                  + F_E_N_C(2)*F_E_N_C(2)
     2                  + F_E_N_C(3)*F_E_N_C(3))

      LAMBDA_L_E = DSQRT(F_E_N_L(1)*F_E_N_L(1)
     1                  + F_E_N_L(2)*F_E_N_L(2)
     2                  + F_E_N_L(3)*F_E_N_L(3))


      I4C_E = LAMBDA_C_E*LAMBDA_C_E
      I4L_E = LAMBDA_L_E*LAMBDA_L_E

      DO I1=1,3
         N_C_CUR(I1) = F_E_N_C(I1)/LAMBDA_C_E
         N_L_CUR(I1) = F_E_N_L(I1)/LAMBDA_L_E
      END DO


C Compressible neo-Hookean isotropic Cauchy stress.
      PRESSURE_TERM = KAPPA*LOG_J_E - MU

      DO I1=1,3
         DO I2=1,3
            SIGMA_ISO(I1,I2) =
     1         (MU/J_E)*B_E(I1,I2)
         END DO
         SIGMA_ISO(I1,I1) = SIGMA_ISO(I1,I1)
     1                   + (PRESSURE_TERM/J_E)
      END DO


C Tension-only fiber response: negative brackets are clipped to zero.
      BRACKET_C = I4C_E - ONE
      BRACKET_L = I4L_E - ONE

      IF (BRACKET_C .LT. ZERO) BRACKET_C = ZERO
      IF (BRACKET_L .LT. ZERO) BRACKET_L = ZERO

      DO I1=1,3
         DO I2=1,3
            SIGMA_C(I1,I2) =
     1         (THREE*K_C/J_E)
     2         * I4C_E*BRACKET_C*BRACKET_C
     3         * N_C_CUR(I1)*N_C_CUR(I2)

            SIGMA_L(I1,I2) =
     1         (THREE*K_L/J_E)
     2         * I4L_E*BRACKET_L*BRACKET_L
     3         * N_L_CUR(I1)*N_L_CUR(I2)

            SIGMA(I1,I2) = SIGMA_ISO(I1,I2)
     1                       + SIGMA_C(I1,I2)
     2                       + SIGMA_L(I1,I2)
         END DO
      END DO


C Assemble material and geometric spatial tangent contributions.
      DO I1=1,3
         DO I2=1,3
            DO I3=1,3
               DO I4=1,3

                  C_ISO(I1,I2,I3,I4) =
     1              (KAPPA/J_E)
     2              * DELTA(I1,I2)*DELTA(I3,I4)
     3            - (PRESSURE_TERM/J_E)
     4              * ( DELTA(I1,I3)*DELTA(I2,I4)
     5                + DELTA(I1,I4)*DELTA(I2,I3) )

                  C_C(I1,I2,I3,I4) =
     1              (12.D0*K_C*I4C_E*I4C_E
     2              *BRACKET_C/J_E)
     3              * N_C_CUR(I1)*N_C_CUR(I2)
     4              * N_C_CUR(I3)*N_C_CUR(I4)

                  C_L(I1,I2,I3,I4) =
     1              (12.D0*K_L*I4L_E*I4L_E
     2              *BRACKET_L/J_E)
     3              * N_L_CUR(I1)*N_L_CUR(I2)
     4              * N_L_CUR(I3)*N_L_CUR(I4)

                  C_GEOM(I1,I2,I3,I4) =
     1              HALF*
     2              ( DELTA(I1,I3)*SIGMA(I2,I4)
     3              + DELTA(I1,I4)*SIGMA(I2,I3)
     4              + DELTA(I2,I3)*SIGMA(I1,I4)
     5              + DELTA(I2,I4)*SIGMA(I1,I3) )

                  C_TOTAL(I1,I2,I3,I4) =
     1               C_ISO(I1,I2,I3,I4)
     2             + C_C(I1,I2,I3,I4)
     3             + C_L(I1,I2,I3,I4)
     4             + C_GEOM(I1,I2,I3,I4)

               END DO
            END DO
         END DO
      END DO


C Algorithmic stress sensitivity to the fibrosis scalar m.
      GAMMA_G_DT = GAMMA_G*DTIME

      LAMBDA_R_E = DSQRT(F_E(1,1)*F_E(1,1)
     1              + F_E(2,1)*F_E(2,1)
     2              + F_E(3,1)*F_E(3,1))
      DO I1=1,3
         N_R_CUR(I1) = F_E(I1,1)/LAMBDA_R_E
      END DO

      D_PRESSURE_D_M =
     1   (GAMMA_G_DT/J_E)
     2 * (KAPPA*LOG_J_E - MU - KAPPA)

      DO I1=1,3
         DO I2=1,3
            D_SIGMA_D_M(I1,I2) =
     1         (MU*GAMMA_G_DT/J_E)
     2         * ( B_E(I1,I2)
     3           - TWO*LAMBDA_R_E*LAMBDA_R_E*N_R_CUR(I1)*N_R_CUR(I2) )
         END DO
         D_SIGMA_D_M(I1,I1) = D_SIGMA_D_M(I1,I1)
     1                        + D_PRESSURE_D_M
      END DO

      DO I1=1,3
         DO I2=1,3
            D_SIGMA_D_M(I1,I2) = D_SIGMA_D_M(I1,I2)
     1                            + GAMMA_G_DT*SIGMA_C(I1,I2)
     2                            + GAMMA_G_DT*SIGMA_L(I1,I2)
         END DO
      END DO


C Map full tensor quantities to Abaqus Voigt ordering.
      CALL VOIGT_MAP_3D(C_TOTAL,DDSDDE)

      STRESS(1) = SIGMA(1,1)
      STRESS(2) = SIGMA(2,2)
      STRESS(3) = SIGMA(3,3)
      STRESS(4) = SIGMA(1,2)
      STRESS(5) = SIGMA(1,3)
      STRESS(6) = SIGMA(2,3)

      DDSDDT(1) = D_SIGMA_D_M(1,1)
      DDSDDT(2) = D_SIGMA_D_M(2,2)
      DDSDDT(3) = D_SIGMA_D_M(3,3)
      DDSDDT(4) = D_SIGMA_D_M(1,2)
      DDSDDT(5) = D_SIGMA_D_M(1,3)
      DDSDDT(6) = D_SIGMA_D_M(2,3)

      SSE = ZERO
      SPD = ZERO
      SCD = ZERO


C Store theta_h for the next increment.
      IF (NSTATV .GE. 1) STATEV(1) = THETA_H_NEW

      RETURN
      END




C ============================================================
C Linear algebra helper routines
C ============================================================
      DOUBLE PRECISION FUNCTION DET3X3(A)
C DET3X3: returns the determinant of a 3x3 matrix.

      INCLUDE 'ABA_PARAM.INC'

      DOUBLE PRECISION A(3,3)

      DET3X3 = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
     1       -A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
     2       +A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

      RETURN
      END


      SUBROUTINE MCOFAC(A,C)
C MCOFAC: returns the cofactor matrix of a 3x3 matrix.

      INCLUDE 'ABA_PARAM.INC'

      DOUBLE PRECISION A(3,3),C(3,3)

      C(1,1)= A(2,2)*A(3,3)-A(2,3)*A(3,2)
      C(1,2)=-(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      C(1,3)= A(2,1)*A(3,2)-A(2,2)*A(3,1)
      C(2,1)=-(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      C(2,2)= A(1,1)*A(3,3)-A(1,3)*A(3,1)
      C(2,3)=-(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      C(3,1)= A(1,2)*A(2,3)-A(1,3)*A(2,2)
      C(3,2)=-(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      C(3,3)= A(1,1)*A(2,2)-A(1,2)*A(2,1)

      RETURN
      END


      SUBROUTINE MTRANS(A,AT)
C MTRANS: returns the transpose of a 3x3 matrix.

      INCLUDE 'ABA_PARAM.INC'

      DOUBLE PRECISION A(3,3),AT(3,3)
      INTEGER I,J

      DO I=1,3
         DO J=1,3
            AT(J,I)=A(I,J)
         END DO
      END DO

      RETURN
      END


      SUBROUTINE M3INV(A,AINV)
C M3INV: returns the inverse of a 3x3 matrix via its adjugate.

      INCLUDE 'ABA_PARAM.INC'

      DOUBLE PRECISION A(3,3),AINV(3,3)
      DOUBLE PRECISION DET,COF(3,3),ADJ(3,3)
      INTEGER I,J
      DOUBLE PRECISION DET3X3

      DET = DET3X3(A)

      CALL MCOFAC(A,COF)
      CALL MTRANS(COF,ADJ)

      DO I=1,3
         DO J=1,3
            AINV(I,J)=ADJ(I,J)/DET
         END DO
      END DO

      RETURN
      END


      SUBROUTINE M3MULT(A,B,C)
C M3MULT: returns the matrix product C = A*B of two 3x3 matrices.

      INCLUDE 'ABA_PARAM.INC'

      DOUBLE PRECISION A(3,3),B(3,3),C(3,3)
      INTEGER I,J,K

      DO I=1,3
         DO J=1,3
            C(I,J)=0.D0
            DO K=1,3
               C(I,J)=C(I,J)+A(I,K)*B(K,J)
            END DO
         END DO
      END DO

      RETURN
      END


      SUBROUTINE VOIGT_MAP_3D(C4,DDSDDE)
C VOIGT_MAP_3D: maps a 4th-order tensor to a 6x6 Abaqus Voigt matrix.

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