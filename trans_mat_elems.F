* =====================================================================
*
*     Transition Matrix Elements for A_{Hg}, A_{Hq}, A_{qq,H}, A_{gq,H}
*           and A_{gg,H}
*
* =====================================================================
*
*     This file contains the exact functions up to NNLO and then some
*     approximations for the N3LO contributions based on available 
*     moments.
*
*     J. M. - 11/11/2020
*
* =====================================================================
*     NLO
* =====================================================================
*     A_{Hg}
! ---------------------------------------------------------------------
!
! ..This is the regular piece.
!
      FUNCTION AHg1A (X, eps)
      IMPLICIT REAL*8 (A-Z)
*
      TR = 1./2.
*       
      AHg1A = -2.*TR*(X**2. + (1.-X)**2.)*LOG(eps)
*
      RETURN 
      END                                           
!
* =====================================================================
*     NNLO
* =====================================================================
*     A_{Hg}
! ---------------------------------------------------------------------
!
! ..File: xa2hgp.f (provided by A.Vogt for G.Salam's HOPPET)
!
!
! ..Calculation of the alpha_s^2 heavy-quark singlet operator matrix
!    element (OME) a^2_Hg(x) in the MS(bar) scheme for mu_f = m_H via
!    a compact parametrization involving only logarithms.
!    The coupling constant is normalized as  a_s = alpha_s/(4*pi).
!
! ..This quantity, presented in Appendix B of M. Buza, Y. Matiounine,
!    J. Smith and W.L. van Neerven, Eur. Phys. J. C1 (1998) 301 (BSMN),
!    is required for the N_f matching of the NNLO parton densities.
!
!  ..The distributions (in the mathematical sense) are given as in eq.
!    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
!    The name-endings A, B, and C of the functions below correspond to
!    the kernel superscripts [2], [3], and [1] in that equation.
!
!  ..The relative accuracy of the OME, as well as of its convolutions
!    with gluon distributions, amounts to a few thousandth.
!
!  ..The user should also cite the original calculation by BSMN.
!
! =====================================================================
!
! ..This is the regular piece.
!
      FUNCTION AHG2A (Y)
      IMPLICIT REAL*8 (A-Z)
* 
      DL  = LOG(Y)
      DL1 = LOG(1.-Y)
*       
      AHg2A = - 24.89d0 / Y - 187.8d0 + 249.6d0 * Y
     &     - 146.8d0 * DL**2 * DL1
     &     - 1.556d0 * DL**3  - 3.292d0 * DL**2  - 93.68d0 * DL
     &     - 1.111d0 * DL1**3 - 0.400d0 * DL1**2 - 2.770d0 * DL1
*
      RETURN 
      END                                           
!
! ---------------------------------------------------------------------
!
! ..This is the 'local' piece, which has no counterpart in W. van
!    Neerven's program, as it does not exist in the exact expressions.
*       
      FUNCTION AHG2C (Y) 
      IMPLICIT REAL*8 (A-Z) 
*       
      AHg2C = - 0.006d0 
*       
      RETURN 
      END
*
! =====================================================================
*     A_{Hq}
! ---------------------------------------------------------------------
!
! ..This is the regular piece.
!
      FUNCTION AHQ2A(X)
      IMPLICIT REAL*8 (A-Z)
* 
      DL  = LOG(X) 
      DL1 = LOG(1.-X)
*       
      AHQ2A = -(896./81.)*(1./X) - 88.76*X + 59.07*X**2. + 3.41*DL 
     .+ 8.078*DL**2. + 16.42*DL*DL1**2. + 40.98*X**3. + 97.04*DL*DL1
*
      RETURN 
      END  
*
! =====================================================================
*     A_{qq,H}
! ---------------------------------------------------------------------
!
! ..This is the regular piece.
!
      FUNCTION AqqH2A(X)
      IMPLICIT REAL*8 (A-Z)
*
      CF = 4./3.
      TR = 1./2.
* 
      DL  = LOG(X) 
      DL1 = LOG(1.-X)
*       
      AqqH2A = CF*TR*((((2./3.)*DL**2. + (20./9.)*DL)*(1 + X**2.)/(1.-X)
     .+ (8./3.*(1.-X)*DL + 44./27. - (268./27.)*X)))
*
      RETURN 
      END  
*
! ---------------------------------------------------------------------
!
! ..This is the singular piece.
!
      FUNCTION AqqH2B(X)
      IMPLICIT REAL*8 (A-Z)
* 
      CF = 4./3.
      TR = 1./2.
      DL  = LOG(X) 
      DL1 = LOG(1.-X)
*       
      AqqH2B = CF*TR*((224./27.)*(1./(1.-X)))
*
      RETURN 
      END  
*
! ---------------------------------------------------------------------
!
! ..This is the residual piece.
!
      FUNCTION AqqH2C(X)
      IMPLICIT REAL*8 (A-Z)
* 
      CF = 4./3.
      TR = 1./2.
*
      DL  = LOG(X) 
      DL1 = LOG(1.-X)
*
      z2=1.6449
      z3=1.20205
*       
      AqqH2C = CF*TR*((224./27.)*DL1 + 73./18. + (40./9.)*z2 
     .- (8./3.)*z3)
*
      RETURN 
      END  
*
! ---------------------------------------------------------------------
!
! ..This is the regular piece.
!
      FUNCTION AqqH2fullA(Z, eps)
      IMPLICIT REAL*8 (A-Z)
*
      CF = 4./3.
      TR = 1./2.
* 
      ALZ  = LOG(Z) 
*       
      AqqH2fullA = CF*TR*((8./3.*log(z)*log(eps)+(2./3.*alz*alz
     .+20./9.*alz))*(1.+z*z)/(1.-z)-4./3.*(1.+z)*(log(eps))**2.
     .+(8./9.-88./9.*z)*log(eps)
     .+(8./3.*(1.-z)*alz+44./27.-268./27.*z))
*
      RETURN 
      END  
*
! ---------------------------------------------------------------------
!
! ..This is the singular piece.
!
      FUNCTION AqqH2fullB(Z, eps)
      IMPLICIT REAL*8 (A-Z)
* 
      CF = 4./3.
      TR = 1./2.
*       
      AqqH2fullB = CF*TR*(8./3.*(log(eps))**2.+80./9.*log(eps)
     .+224./27.)/(1.-Z)
*
      RETURN 
      END  
*
! ---------------------------------------------------------------------
!
! ..This is the residual piece.
!
      FUNCTION AqqH2fullC(Z, eps)
      IMPLICIT REAL*8 (A-Z)
* 
      CF = 4./3.
      TR = 1./2.
*
      AL  = LOG(Z) 
      AL1 = LOG(1.-Z)
*
      z2=1.6449
      z3=1.20205
*       
      AqqH2fullC = CF*TR*((8./3.*(log(eps))**2. + 80./9.*log(eps)
     .+ 224./27.)*AL1 + (2*(log(eps))**2. + (16./3.*z2 
     .+ 2./3.)*log(eps) + (73./18. + 40./9.*z2 - 8./3.*z3)))
*
      RETURN 
      END  
*
! =====================================================================
*     A_{gq,H}
! ---------------------------------------------------------------------
!
! ..This is the regular piece.
!
      FUNCTION AgqH2A(X)
      IMPLICIT REAL*8 (A-Z)
*
      CF = 4./3.
      TR = 1./2.
* 
      DL  = LOG(X) 
      DL1 = LOG(1.-X)
*       
      AgqH2A = CF*TR*((4./3.)*((2./X) - 2. + X)*DL1**2. + (8./9.)*((10./X) 
     .- 10. + 8.*X)*DL1 + ((448./X) - 448. + 344.*x)/27.)
*
      RETURN 
      END  
*
! =====================================================================
*     A_{gg,H}
! ---------------------------------------------------------------------
!
! ..This is the regular piece.
!
      FUNCTION AggH2A(X)
      IMPLICIT REAL*8 (A-Z)
*
      CF = 4./3.
      TR = 1./2.
      CA = 3.
* 
      DL  = LOG(X) 
      DL1 = LOG(1.-X)
*       
      AggH2A = CF*TR*((4./3.)*(1.+X)*DL**3. + (6. + 10.*X)*DL**2. + (32. 
     .+ 48.*X)*DL - (8./X) + 80. - 48.*X - 24.*X**2.) 
     .+ CA*TR*((4./3.)*(1. + X)*DL**2. + (52. + 88.*X)*DL/9. - (4./3.)*X*DL1 
     .+ (556./27.)*(1./X) - (628./27.) + (548./27.)*X - (700./27.)*X**2.)
*
      RETURN 
      END  
*
! ---------------------------------------------------------------------
!
! ..This is the singular piece.
!
      FUNCTION AggH2B(X)
      IMPLICIT REAL*8 (A-Z)
* 
      CA = 3.
      TR = 1./2.
*
      DL  = LOG(X) 
      DL1 = LOG(1.-X)
*       
      AggH2B = CA*TR*((224./27.)*(1./(1.-X)))
*
      RETURN 
      END  
*
! ---------------------------------------------------------------------
!
! ..This is the residual piece.
!
      FUNCTION AggH2C(X)
      IMPLICIT REAL*8 (A-Z)
* 
      CA = 3.
      TR = 1./2.
      CF = 4./3.
*
      DL  = LOG(X) 
      DL1 = LOG(1.-X)
*       
      AggH2C = CA*TR*((224./27.)*DL1 + 10./9.) - CF*TR*15.
*
      RETURN 
      END  
*
* =====================================================================
*     N3LO
* =====================================================================
*     A_{Hg}
! ---------------------------------------------------------------------
!
! ..This is the function which approximates AHg(3) based on the 
!   variational parameter a (xp(56)).
!   The method of simultaneous equations is used to estimate the 
!   coefficients of a fixed functional form using the hard coded Mellin
!   moments.
!
      RECURSIVE FUNCTION AHg3A(x, a) RESULT(AHg3A_res)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER n_dimension
      DOUBLE PRECISION N(5), C(5), N_exact(5), 
     &        N_a(5)
      DOUBLE PRECISION inv_A(5, 5)
      DATA PI,PI2/3.14159,9.8696/
*
      n_dimension=5
C 
C     Known moments minus those for the known exact small-x contribution (Vogt)
      N_exact(1) = 1996.50147366
      N_exact(2) = 232.27499707
      N_exact(3) = 7.44963895
      N_exact(4) = -77.93698013
      N_exact(5) = -121.36688921
C 
C     Moments for the varied function (1/x) multiplied by a 
      N_a(1) = 1.0D0
      N_a(2) = 0.3333333333333333D0
      N_a(3) = 0.20D0
      N_a(4) = 0.14285714285714285D0
      N_a(5) = 0.1111111111111111D0
*
      DO i=1, n_dimension
        N(i) = N_exact(i)-a*N_a(i)
      ENDDO
*
      inv_A = reshape((/-1.13879194D-01, -7.30104359D-01, -1.16019755D+02,
     &               1.47898113D+02, -2.48866805D+01, 2.60826752D+00,
     &               1.65606748D+01, 1.87421683D+03, -2.53062040D+03,  
     &               2.76071734D+02, -1.16231560D+01, -7.31022502D+01, 
     &               -6.56760263D+03, 9.21295249D+03, -8.50419108D+02, 
     &               1.73021125D+01, 1.07828576D+02, 8.21779251D+03, 
     &               -1.18144214D+04,  9.95444799D+02, -8.20269840D+00, 
     &               -5.06719525D+01, -3.41032145D+03, 4.98742337D+03, 
     &               -3.96388523D+02/), shape=(/5,5/))

C     Matrix multiplication of inv_A and N
*
      DO i=1,n_dimension
        C(i) = 0.0D0
        DO j=1,n_dimension
          C(i) = C(i) + inv_A(i, j)*N(j)
        ENDDO
      ENDDO
*
      AHg3A_res = C(1)*LOG(1.-X)**5.+C(2)*LOG(1.-X)**4.+C(3)*X 
     &          +C(4)*X**2.+C(5)*LOG(X)+a*(1./X)
     &          +AHg3exactA(X)
*
      if(x.gt.0.924D0) AHg3A_res = AHg3A(0.924D0, a)*(1.-X)
      ! AHg3A_res = 0.0D0
*
      RETURN 
      END 
*
! ---------------------------------------------------------------------
!
! ..This is the regular piece of the known exact contribution to AHg(3).
!    Also contained is the approximate expression which converts to MSbar.
!      
      FUNCTION AHg3exactA(Y) 
      IMPLICIT DOUBLE PRECISION (A-Z) 
      DATA PI,PI2/3.14159D0,9.8696D0/
*     
      AHg3exactA = ((41984.0D0/27.0D0)+160.0D0*(PI2/6.0D0)
     &             -224.0D0*1.20206D0)*LOG(Y)/Y
*
      RETURN
      END
*
! =====================================================================
*     A_{Hq}
! ---------------------------------------------------------------------
!
! ..This is the approximate parameterisation of the AHq(3) function.
!   Including the conversion function to MSbar.
!      
      RECURSIVE FUNCTION AHq3A(Y) RESULT(AHq3A_res)
      IMPLICIT DOUBLE PRECISION (A-Z) 
      DATA PI,PI2/3.14159,9.8696/
*
      AHq3A_res = ((1.-Y)**2)*(-158.3709570235676D0*log(1.-Y)**3.
     .+2151.0528216721095D0*log(1.-Y)**2.-14178.20613565087D0)
     .+0.7136042069623996D0*Y+15618.788026509601D0*Y*log(Y)**2.
     .+2117.261727957378D0*Y*log(Y)+688.396D0*log(Y)/Y
     .+(1.-Y)*3812.8990D0*(1/Y)+1.6D0*(log(Y)**5)-20.345678D0*(log(Y)**4)
     .+165.11455D0*(log(Y)**3)-604.63554D0*(log(Y)**2)+3524.9967D0*log(Y)
     .+(1.-Y)*(0.24691358D0*(log(1.-Y)**4)-4.44444444D0*(log(1.-Y)**3)
     .-2.28230742D0*(log(1.-Y)**2)-357.426943D0*log(1.-Y)+116.478169D0)

      ! if (y.gt.0.4D0) AHq3A_res = AHq3A(0.4D0)*(1.-Y)/0.6D0

      ! AHq3A_res = 0.0D0
*
      RETURN
      END
*
! =====================================================================
*     A_{qq,H}
! ---------------------------------------------------------------------
!
! ..This is the function which approximates AqqH(3) based on the 
!   variational parameter a (xp(?)).
!   The method of simultaneous equations is used to estimate the 
!   coefficients of a fixed functional form using the hard coded Mellin
!   moments.
!
      SUBROUTINE AqqH3AB(x, a, AqqH3A, AqqH3B)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER n_dimension
      DOUBLE PRECISION N(7), C(7), N_exact(7), N_a(7)
      DOUBLE PRECISION inv_A(7, 7)
      DATA PI,PI2/3.14159,9.8696/
*
      n_dimension=7
C 
C     Known moments minus those for the known exact small-x contribution (Vogt)
      N_exact(1) = -31.101
      N_exact(2) = -27.989
      N_exact(3) = -20.146
      N_exact(4) = -13.313
      N_exact(5) = -7.629
      N_exact(6) = -2.854
      N_exact(7) = 1.228
C 
C     Moments for the varied function (1/x) multiplied by a 
      N_a(1) = -0.375D0
      N_a(2) = -0.023437499998337614D0
      N_a(3) = -0.004629629629528962D0
      N_a(4) = -0.0014648437499999065D0
      N_a(5) = -0.0005999999999999996D0
      N_a(6) = -0.0002893518518518518D0
      N_a(7) = -0.00015618492294877136D0
*
      DO i=1, n_dimension
        N(i) = N_exact(i) - a*N_a(i)
      ENDDO
*
      inv_A = reshape((/2.16367103D+01, 2.20120178D+02, -7.28098793D-01,
     &                  1.08859149D+03, 3.16215872D+01, 7.57356352D+02,
     &                  -1.31967754D+02, -1.30805003D+03, -1.30363921D+04,
     &                  4.71551230D+01, -6.03760425D+04,-1.06557895D+03,
     &                  -4.38642686D+04, 5.99320146D+03, 1.36620524D+04, 
     &                  1.34037998D+05, -5.20816646D+02, 5.95121449D+05, 
     &                  8.47026017D+03, 4.43696299D+05, -5.19148509D+04, 
     &                  -5.23812588D+04, -5.07554870D+05,2.09319457D+03, 
     &                  -2.18873394D+06, -2.74098723D+04, -1.65937417D+06, 
     &                  1.75138916D+05, 9.13767640D+04, 8.76404430D+05, 
     &                  -3.80468317D+03, 3.69828317D+06, 4.24347906D+04, 
     &                  2.83722814D+06, -2.77933570D+05,-7.40399878D+04, 
     &                  -7.04036144D+05, 3.19808122D+03,-2.92103432D+06, 
     &                  -3.14155811D+04, -2.26097908D+06,2.09168184D+05, 
     &                  2.26683150D+04, 2.13961524D+05, -1.01233621D+03, 
     &                  8.75640617D+05, 8.95431572D+03, 6.82526398D+05, 
     &                  -6.03194325D+04/), shape=(/7,7/))
*
*     Matrix multiplication of inv_A and N
*
      DO i=1,n_dimension
        C(i) = 0.0D0
        DO j=1,n_dimension
          C(i) = C(i) + inv_A(i, j)*N(j)
        ENDDO
      ENDDO
*
      AqqH3A = C(1)*LOG(1.-X)**3. + C(2)*LOG(1.-X)**2.
     &          +C(4)*X + C(5)*LOG(X)**2. + C(6)*LOG(1.-X) + C(7) 
     &          + a*LOG(X)**3.
*
      AqqH3B = C(3)*(1./(1.-X))
*
      if (X.gt.0.4D0) AqqH3A = 0.0D0
      if (X.gt.0.4D0) AqqH3B = 0.0D0
*
      RETURN 
      END 
*
! =====================================================================
*     A_{gq,H}
! ---------------------------------------------------------------------
!
! ..This is the approximate function the AgqH(3). Including its 
!   conversion to MSbar.
!
      RECURSIVE FUNCTION AgqH3A(Y) RESULT(AgqH3A_res)
      IMPLICIT DOUBLE PRECISION (A-Z) 
      DATA PI,PI2/3.14159,9.8696/
*
      AgqH3A_res = -237.1720947621626*log(1.-Y)**3. - 201.496990873891*log(1.-Y)**2. 
     &        + 7247.6979315802555*log(1.-Y) + 39967.310031716734*Y**2. 
     &        - 22017.709314369586 - 28459.052011351927*log(Y) 
     &        - 14511.477130954034*log(Y)**2. + 341.543*(log(1./Y)/Y) 
     &        + 1814.73*(1./Y) - (580./243.)*log(1.-Y)**4. 
     &        - (17624./729.)*log(1.-Y)**3. - 135.699*log(1.-Y)**2.

      if (y.gt.0.7D0) AgqH3A_res = AgqH3A(0.7D0)*(1.-Y)/0.3D0
*
      RETURN
      END                                          
!
! =====================================================================
*     A_{gg,H}
! ---------------------------------------------------------------------
!
! ..This is the function which approximates the AggH(3) function based
!   based on the variational parameter a (xp(57)).
!   It uses simultaneous equations to estimate the coefficients of a 
!   fixed functional form using the hard coded moments.
!
      RECURSIVE FUNCTION AggH3A(X, a) RESULT(AggH3A_res)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER n_dimension
      DOUBLE PRECISION N(5), C(5), N_exact(5), 
     &        N_a(5)
      DOUBLE PRECISION inv_A(5, 5)
      DATA PI,PI2/3.14159,9.8696/
*
      n_dimension=5
*
*     Known moments minus the known exact small-x contributions
      N_exact(1) = -441.3444D0
      N_exact(2) = -96.32D0
      N_exact(3) = -13.33D0
      N_exact(4) = 31.1D0
      N_exact(5) = 60.84D0
* 
*     Moments for the varied function (log(x)/x) multiplied by a 
      N_a(1) = -0.9999999999999999D0
      N_a(2) = -0.11111111111076043D0
      N_a(3) = -0.040000000000137814D0
      N_a(4) = -0.020408163265306242D0
      N_a(5) = -0.01234567901234568D0
*
      DO i=1, n_dimension
        N(i) = N_exact(i) - a * N_a(i)
      ENDDO
*
!       inv_A = reshape((/1.05434437D-01, -3.76135977D+00, 2.08675271D+01, 
!      &                  -3.56461172D+01, 1.86043428D+01, 6.29487422D-01, 
!      &                  -2.19141532D+01,  1.19081094D+02, -1.99741697D+02, 
!      &                  1.02541857D+02, 2.85832310D+01, -6.90092962D+02, 
!      &                  2.82378438D+03, -3.88590569D+03, 1.72638239D+03, 
!      &                  8.06540899D+00, 1.98982084D+01, -4.65536179D+02, 
!      &                  9.44684628D+02, -5.08912706D+02, 7.21897161D+00, 
!      &                  -2.13313290D+02,1.00699129D+03, -1.50242240D+03, 
!      &                  7.03070750D+02/), shape=(/5,5/))

      inv_A = reshape((/1.43759679D+01, -3.46333748D+02, 1.61699126D+03,
     &                  -2.51174722D+03, 1.23814937D+03, 1.21751771D+02, 
     &                  -2.87863393D+03, 1.32164568D+04, -2.02263161D+04, 
     &                  9.83890667D+03, 2.91559485D+02, -5.88750691D+03,
     &                  2.44589984D+04, -3.49174315D+04, 1.61247555D+04,
     &                  -2.80668008D+01,3.50042899D+02, -1.18494800D+03,
     &                  1.50037935D+03,-6.38918620D+02, -3.82520602D+01,
     &                  1.94881955D+01, 2.01486391D+03, -5.00918023D+03,
     &                  3.06450021D+03/), shape=(/5,5/))
*
*     Matrix multiplication of inv_A and N
      DO i=1,n_dimension
        C(i) = 0
        DO j=1,n_dimension
          C(i) = C(i) + inv_A(j, i) * N(j)
        ENDDO
      ENDDO
*
!       AggH3A_res = (C(1) * LOG(1.-X)**4. + C(2) * LOG(1.-X)**3. 
!      &          + C(3) * LOG(X) + C(4) * LOG(X)**2. 
!      &          + C(5) * (1./x) + a * (LOG(X)/X))

      AggH3A_res = (C(1) * LOG(1.-X)**2. + C(2) * LOG(1.-X) 
     &          + C(3) * X**2. + C(4) * LOG(X) 
     &          + C(5) * X + a * (LOG(X)/X))
*
!      if (x.gt.0.4D0) AggH3A_res = AggH3A(0.4D0, a)*(1.-X)/0.6D0
      ! AggH3A_res = 0.0D0
*
      RETURN 
      END 