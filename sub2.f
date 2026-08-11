C ============================================================
C Model 2 UMATHT: monodomain electrical transport
C
C The Abaqus temperature degree of freedom is used as the
C normalized ICC potential phi.  This routine supplies the
C storage term, electrical flux, and diffusion tangent.
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

      DOUBLE PRECISION PHI,C_PHI,D_0,D_EFF,M,ETA_D
C
C PHI is the normalized electrical potential.
C M is the imported fibrosis field from Model 1.
C C_PHI, D_0, and ETA_D define c_phi and d_eff.
C     Symbol names follow the manuscript equations where possible.
      INTEGER I1,I2



C Current electrical potential phi and imported fibrosis scalar m.
      PHI     = TEMP + DTEMP
      M       = PREDEF(1)
      C_PHI   = PROPS(1)
      D_0      = PROPS(2)
      ETA_D = PROPS(3)


C Effective diffusivity: d_eff = d_0 (1 - eta_D m).
      D_EFF = D_0*(1.D0 - ETA_D*M)


C Storage term: U = c_phi phi, with tangent dU/dphi = c_phi.
      U    = C_PHI*PHI
      DUDT = C_PHI

      DO I1=1,NTGRD
         DUDG(I1) = 0.D0

C Electrical flux: q = -d_eff grad(phi).
         FLUX(I1) = -D_EFF*DTEMDX(I1)
         DFDT(I1) = 0.D0
         DO I2=1,NTGRD
            DFDG(I1,I2) = 0.D0
         END DO
         DFDG(I1,I1) = -D_EFF
      END DO

      RETURN
      END




C ============================================================
C Model 2 UMAT: electromechanics and active smooth-muscle stress
C
C This routine updates the FHN recovery variable w, computes
C active tension, evaluates passive and active stresses, and
C returns mechanical/electromechanical tangents.
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

C Main Model 2 parameters are read from PROPS in the same order
C used in the Abaqus input deck.  The names match the manuscript
C symbols where Fortran naming permits.
C
C State variables used by Model 2:
C   STATEV(1) = w, FHN recovery variable
C   STATEV(2) = T_act, active smooth-muscle tension
C   STATEV(3) = phi_smc, dimensional SMC membrane potential
C   STATEV(4) = I4c, circumferential fiber invariant
C   STATEV(5) = I4l, longitudinal fiber invariant
C   STATEV(6) = m, imported fibrosis scalar
C   STATEV(7:9)   = N_c components
C   STATEV(10:12) = e_r components

      DOUBLE PRECISION PHI,M
      DOUBLE PRECISION W_OLD,W_NEW,D_W_D_PHI,D_W_D_EPSILON
      DOUBLE PRECISION T_ACT_OLD,T_ACT_NEW
      DOUBLE PRECISION PHI_SMC,A_GATE,T_ACT_TARGET,T_MAX,
     1 D_A_GATE_D_PHI,D_T_ACT_D_PHI
      DOUBLE PRECISION XI
      DOUBLE PRECISION PHI_R,PHI_OFF,A_SW,PHI_TH,DELTA_S
      DOUBLE PRECISION TAU_T,T0,DT_TAU
      DOUBLE PRECISION KAPPA_FHN,A_0,A_THRESH,B_PARAM,C_PARAM
      DOUBLE PRECISION EPSILON,EPSILON_RAW,DENOM_W
      DOUBLE PRECISION EPSILON_0,EPSILON_Z,EPSILON_LAMBDA
      DOUBLE PRECISION ETA_MU,ETA_T,ETA_A
      DOUBLE PRECISION Z_COORD,D_EPSILON_D_LAMBDA_C,D_R_D_LAMBDA_C


      DOUBLE PRECISION MU_0,K_C,K_L
      DOUBLE PRECISION MU,KAPPA,KAPPA_FAC
      DOUBLE PRECISION F(3,3),F_E(3,3),F_H_INV(3,3)
      DOUBLE PRECISION B_LEFT_CG(3,3),IDENT(3,3)
      DOUBLE PRECISION SIGMA_PAS(3,3),SIGMA_ACT(3,3),SIGMA(3,3)
      DOUBLE PRECISION SIGMA_C(3,3),SIGMA_L(3,3)
      DOUBLE PRECISION C_PAS(3,3,3,3),C_ACT(3,3,3,3),C_TOTAL(3,3,3,3)
      DOUBLE PRECISION C_C(3,3,3,3),C_L(3,3,3,3)
      DOUBLE PRECISION C_GEOM(3,3,3,3)
      DOUBLE PRECISION J,LOG_J,PRESSURE_TERM

      DOUBLE PRECISION E_R(3),N_C(3),N_L(3)
      DOUBLE PRECISION F_N_C(3),F_N_L(3)
      DOUBLE PRECISION N_C_CUR(3),N_L_CUR(3)
      DOUBLE PRECISION E_Z(3),NORM_VALUE

      DOUBLE PRECISION LAMBDA_C,LAMBDA_L
      DOUBLE PRECISION I4C,I4L,BRACKET_C,BRACKET_L
      DOUBLE PRECISION P_C(3,3)

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


C Read the primary electrical field and the transferred triad fields.
C PREDEF(1) is m, PREDEF(2:4) is e_r, PREDEF(5:7) is e_theta,
C and PREDEF(8:10) is e_z.
      PHI = TEMP + DTEMP
      M   = PREDEF(1)

      E_R(1) = PREDEF(2)
      E_R(2) = PREDEF(3)
      E_R(3) = PREDEF(4)

      N_C(1) = PREDEF(5)
      N_C(2) = PREDEF(6)
      N_C(3) = PREDEF(7)

      E_Z(1) = PREDEF(8)
      E_Z(2) = PREDEF(9)
      E_Z(3) = PREDEF(10)



C Normalize imported basis vectors before using them in tensor operations.
      NORM_VALUE = DSQRT(E_Z(1)*E_Z(1)+E_Z(2)*E_Z(2)+E_Z(3)*E_Z(3))
      E_Z(1) = E_Z(1)/NORM_VALUE
      E_Z(2) = E_Z(2)/NORM_VALUE
      E_Z(3) = E_Z(3)/NORM_VALUE

      NORM_VALUE = DSQRT(E_R(1)*E_R(1)+E_R(2)*E_R(2)+E_R(3)*E_R(3))
      E_R(1) = E_R(1)/NORM_VALUE
      E_R(2) = E_R(2)/NORM_VALUE
      E_R(3) = E_R(3)/NORM_VALUE

      N_L(1) = E_Z(1)
      N_L(2) = E_Z(2)
      N_L(3) = E_Z(3)


C Retrieve history variables from the previous increment.
      W_OLD    = ZERO
      T_ACT_OLD = ZERO
      IF (NSTATV .GE. 1) W_OLD    = STATEV(1)
      IF (NSTATV .GE. 2) T_ACT_OLD = STATEV(2)



C Read FHN, active-stress, passive-mechanics, and disease parameters.
      KAPPA_FHN          = PROPS(1)
      A_0 = PROPS(2)
      B_PARAM          = PROPS(3)
      C_PARAM          = PROPS(4)
      PHI_R          = PROPS(5)
      PHI_OFF     = PROPS(6)
      A_SW          = PROPS(7)
      PHI_TH         = PROPS(8)
      DELTA_S        = PROPS(9)
      T0             = PROPS(10)
      TAU_T          = PROPS(11)
      MU_0            = PROPS(12)
      EPSILON_0      = PROPS(13)
      EPSILON_Z = PROPS(14)
      K_C        = PROPS(15)
      K_L        = PROPS(16)
      EPSILON_LAMBDA    = PROPS(17)
      ETA_MU         = PROPS(18)
      ETA_T       = PROPS(19)
      ETA_A          = PROPS(20)


C Disease-dependent constitutive parameters.
      MU      = MU_0*(ONE + ETA_MU*M)
C Bulk modulus is hardcoded at 100x the shear modulus to enforce
C near-incompressibility.
      KAPPA = 100.D0*MU
      T_MAX  = T0*(ONE - ETA_T*M)
      A_THRESH = A_0*(ONE - ETA_A*M)


C Model 2 has no additional growth, so F_h^{-1} is the identity
C and the elastic deformation gradient equals the total gradient.
      DO I1=1,3
         DO I2=1,3
            F(I1,I2)     = DFGRD1(I1,I2)
            F_E(I1,I2)    = ZERO
            IDENT(I1,I2) = ZERO
            B_LEFT_CG(I1,I2)  = ZERO
         END DO
         IDENT(I1,I1) = ONE
      END DO


      DO I1=1,3
         DO I2=1,3
            F_H_INV(I1,I2) = IDENT(I1,I2)
         END DO
      END DO


      DO I1=1,3
         DO I2=1,3
            DO I3=1,3
               F_E(I1,I2) = F_E(I1,I2) + F(I1,I3)*F_H_INV(I3,I2)
            END DO
         END DO
      END DO



C Jacobian and left Cauchy-Green tensor for passive mechanics.
      J = DET3X3(F_E)
      LOG_J = DLOG(J)

      DO I1=1,3
         DO I2=1,3
            DO I3=1,3
               B_LEFT_CG(I1,I2) = B_LEFT_CG(I1,I2)
     1                          + F_E(I1,I3)*F_E(I2,I3)
            END DO
         END DO
      END DO



C Push forward circumferential and longitudinal fiber directions.
      DO I1=1,3
         F_N_C(I1) = F_E(I1,1)*N_C(1)
     1                + F_E(I1,2)*N_C(2)
     2                + F_E(I1,3)*N_C(3)

         F_N_L(I1) = F_E(I1,1)*N_L(1)
     1                + F_E(I1,2)*N_L(2)
     2                + F_E(I1,3)*N_L(3)
      END DO

      LAMBDA_C = DSQRT(F_N_C(1)*F_N_C(1)
     1               + F_N_C(2)*F_N_C(2)
     2               + F_N_C(3)*F_N_C(3))
      LAMBDA_L = DSQRT(F_N_L(1)*F_N_L(1)
     1               + F_N_L(2)*F_N_L(2)
     2               + F_N_L(3)*F_N_L(3))

      I4C = LAMBDA_C*LAMBDA_C
      I4L = LAMBDA_L*LAMBDA_L

      DO I1=1,3
         N_C_CUR(I1) = F_N_C(I1)/LAMBDA_C
         N_L_CUR(I1) = F_N_L(I1)/LAMBDA_L
      END DO



C Stretch-dependent recovery rate with a small numerical floor.
      Z_COORD = COORDS(3)
      EPSILON_RAW = EPSILON_0 + EPSILON_Z*Z_COORD
     1       + EPSILON_LAMBDA*(LAMBDA_C - ONE)

      IF (EPSILON_RAW .GT. 1.D-12) THEN
         EPSILON = EPSILON_RAW
         D_EPSILON_D_LAMBDA_C = EPSILON_LAMBDA
      ELSE
         EPSILON = 1.D-12
         D_EPSILON_D_LAMBDA_C = ZERO
      END IF



C Backward-Euler update of the FHN recovery variable w.
C DENOM_W caches ONE + DTIME*EPSILON so it is formed once, not three times.
      DENOM_W = ONE + DTIME*EPSILON

      W_NEW = ( W_OLD + DTIME*EPSILON*B_PARAM*(PHI-C_PARAM) )
     1      / DENOM_W

      D_W_D_PHI = ( DTIME*EPSILON*B_PARAM ) / DENOM_W

      D_W_D_EPSILON = DTIME*(B_PARAM*(PHI-C_PARAM) - W_OLD)
     1       / ( DENOM_W*DENOM_W )



C Map normalized potential phi to dimensional SMC membrane potential.
      PHI_SMC = PHI_R + PHI_OFF + A_SW*PHI
      XI = (PHI_SMC - PHI_TH)/DELTA_S


C Logistic activation gate.  The large-|XI| branches avoid overflow
C while preserving the limiting values of the sigmoid.
      IF (XI .GT. 50.D0) THEN
         A_GATE = ONE
         D_A_GATE_D_PHI = ZERO
      ELSE IF (XI .LT. -50.D0) THEN
         A_GATE = ZERO
         D_A_GATE_D_PHI = ZERO
      ELSE
         A_GATE = ONE/(ONE + DEXP(-XI))
         D_A_GATE_D_PHI = A_GATE*(ONE-A_GATE)*(A_SW/DELTA_S)
      END IF


C First-order active tension update and its phi derivative.
C DT_TAU caches DTIME/TAU_T so the ratio is formed once, not four times.
      DT_TAU       = DTIME/TAU_T
      T_ACT_TARGET = T_MAX*A_GATE
      T_ACT_NEW = ( T_ACT_OLD + DT_TAU*T_ACT_TARGET )
     1         / ( ONE + DT_TAU )

      D_T_ACT_D_PHI = ( ONE/(ONE + DT_TAU) )
     1           * DT_TAU * T_MAX * D_A_GATE_D_PHI



C FHN source term and consistent scalar tangent.
      RPL = KAPPA_FHN*PHI*(PHI-A_THRESH)*(ONE-PHI) - W_NEW

      DRPLDT = KAPPA_FHN*(-THREE*PHI*PHI
     1         + TWO*(ONE+A_THRESH)*PHI - A_THRESH) - D_W_D_PHI


C Strain sensitivity of the source term from mechanoelectrical feedback.
      D_R_D_LAMBDA_C = -D_W_D_EPSILON*D_EPSILON_D_LAMBDA_C

      DRPLDE(1) = D_R_D_LAMBDA_C*(LAMBDA_C*N_C_CUR(1)*N_C_CUR(1))
      DRPLDE(2) = D_R_D_LAMBDA_C*(LAMBDA_C*N_C_CUR(2)*N_C_CUR(2))
      DRPLDE(3) = D_R_D_LAMBDA_C*(LAMBDA_C*N_C_CUR(3)*N_C_CUR(3))
      DRPLDE(4) = D_R_D_LAMBDA_C*(LAMBDA_C*N_C_CUR(1)*N_C_CUR(2))
      DRPLDE(5) = D_R_D_LAMBDA_C*(LAMBDA_C*N_C_CUR(1)*N_C_CUR(3))
      DRPLDE(6) = D_R_D_LAMBDA_C*(LAMBDA_C*N_C_CUR(2)*N_C_CUR(3))



C Passive stress: compressible neo-Hookean matrix plus tension-only fibers.
      PRESSURE_TERM = KAPPA*LOG_J - MU

      DO I1=1,3
         DO I2=1,3
            SIGMA_PAS(I1,I2) = (MU/J)*B_LEFT_CG(I1,I2)
         END DO
         SIGMA_PAS(I1,I1) = SIGMA_PAS(I1,I1) + (PRESSURE_TERM/J)
      END DO


      BRACKET_C = I4C - ONE
      BRACKET_L = I4L - ONE
      IF (BRACKET_C .LT. ZERO) BRACKET_C = ZERO
      IF (BRACKET_L .LT. ZERO) BRACKET_L = ZERO

      DO I1=1,3
         DO I2=1,3
            SIGMA_C(I1,I2) =
     1         (THREE*K_C/J)*I4C*BRACKET_C*BRACKET_C*
     2         N_C_CUR(I1)*N_C_CUR(I2)

            SIGMA_L(I1,I2) =
     1         (THREE*K_L/J)*I4L*BRACKET_L*BRACKET_L*
     2         N_L_CUR(I1)*N_L_CUR(I2)
         END DO
      END DO



C Active Cauchy stress acts along the current circumferential direction.
      DO I1=1,3
         DO I2=1,3
            SIGMA_ACT(I1,I2) = T_ACT_NEW*N_C_CUR(I1)*N_C_CUR(I2)

            SIGMA(I1,I2) = SIGMA_PAS(I1,I2)
     1                       + SIGMA_C(I1,I2)
     2                       + SIGMA_L(I1,I2)
     3                       + SIGMA_ACT(I1,I2)

            P_C(I1,I2) = IDENT(I1,I2)
     1                  - N_C_CUR(I1)*N_C_CUR(I2)
         END DO
      END DO



C Assemble passive, active, and geometric tangent contributions.
      DO I1=1,3
         DO I2=1,3
            DO I3=1,3
               DO I4=1,3

                  C_PAS(I1,I2,I3,I4) =
     1              (KAPPA/J)*IDENT(I1,I2)*IDENT(I3,I4)
     2            - (PRESSURE_TERM/J)*
     3              ( IDENT(I1,I3)*IDENT(I2,I4)
     4              + IDENT(I1,I4)*IDENT(I2,I3) )

                  C_C(I1,I2,I3,I4) =
     1              (12.D0*K_C*I4C*I4C*BRACKET_C/J)
     2            * N_C_CUR(I1)*N_C_CUR(I2)
     3            * N_C_CUR(I3)*N_C_CUR(I4)

                  C_L(I1,I2,I3,I4) =
     1              (12.D0*K_L*I4L*I4L*BRACKET_L/J)
     2            * N_L_CUR(I1)*N_L_CUR(I2)
     3            * N_L_CUR(I3)*N_L_CUR(I4)

                  C_ACT(I1,I2,I3,I4) =
     1              HALF*T_ACT_NEW*(
     2              P_C(I1,I3)*N_C_CUR(I2)*N_C_CUR(I4)
     3            + P_C(I1,I4)*N_C_CUR(I2)*N_C_CUR(I3)
     4            + P_C(I2,I3)*N_C_CUR(I1)*N_C_CUR(I4)
     5            + P_C(I2,I4)*N_C_CUR(I1)*N_C_CUR(I3) )

                  C_GEOM(I1,I2,I3,I4) =
     1              HALF*(
     2              IDENT(I1,I3)*SIGMA(I2,I4)
     3            + IDENT(I1,I4)*SIGMA(I2,I3)
     4            + IDENT(I2,I3)*SIGMA(I1,I4)
     5            + IDENT(I2,I4)*SIGMA(I1,I3) )

                  C_TOTAL(I1,I2,I3,I4) =
     1               C_PAS(I1,I2,I3,I4)
     2             + C_C(I1,I2,I3,I4)
     3             + C_L(I1,I2,I3,I4)
     4             + C_ACT(I1,I2,I3,I4)
     5             + C_GEOM(I1,I2,I3,I4)

               END DO
            END DO
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

      DDSDDT(1) = D_T_ACT_D_PHI*N_C_CUR(1)*N_C_CUR(1)
      DDSDDT(2) = D_T_ACT_D_PHI*N_C_CUR(2)*N_C_CUR(2)
      DDSDDT(3) = D_T_ACT_D_PHI*N_C_CUR(3)*N_C_CUR(3)
      DDSDDT(4) = D_T_ACT_D_PHI*N_C_CUR(1)*N_C_CUR(2)
      DDSDDT(5) = D_T_ACT_D_PHI*N_C_CUR(1)*N_C_CUR(3)
      DDSDDT(6) = D_T_ACT_D_PHI*N_C_CUR(2)*N_C_CUR(3)

      SSE = ZERO
      SPD = ZERO
      SCD = ZERO


C Store selected history and diagnostic quantities for output.
      IF (NSTATV .GE. 1) STATEV(1) = W_NEW
      IF (NSTATV .GE. 2) STATEV(2) = T_ACT_NEW
      IF (NSTATV .GE. 3) STATEV(3) = PHI_SMC
      IF (NSTATV .GE. 4) STATEV(4) = I4C
      IF (NSTATV .GE. 5) STATEV(5) = I4L
      IF (NSTATV .GE. 6) STATEV(6) = M

      IF (NSTATV .GE. 7) STATEV(7)  = N_C(1)
      IF (NSTATV .GE. 8) STATEV(8)  = N_C(2)
      IF (NSTATV .GE. 9) STATEV(9)  = N_C(3)
      IF (NSTATV .GE. 10) STATEV(10) = E_R(1)
      IF (NSTATV .GE. 11) STATEV(11) = E_R(2)
      IF (NSTATV .GE. 12) STATEV(12) = E_R(3)


      RETURN
      END




C ============================================================
C Linear algebra and Voigt mapping helper routines
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