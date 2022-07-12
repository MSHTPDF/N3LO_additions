* =====================================================================
*
*     Splitting Functions for P_{qq}^{NS}, P_{qq}^{PS}, P_{qg}, P_{gq}
*           and P_{gg}
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
*     LO
* =====================================================================
*     P_{qq}^{NS}
* ---------------------------------------------------------------------
*
*
*
      FUNCTION Pqqns0A(X)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      CF = 4./3.
*
      Pqqns0A = - 2*CF*(1. + X)
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
*
      FUNCTION Pqqns0B(X)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      CF = 4./3.
*
      Pqqns0B = 2*CF*(2./(1.-X))
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
*
      FUNCTION Pqqns0C(X)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      CF = 4./3.
      AL1 = LOG(1.-X)
*
      Pqqns0C = 2.*CF*(3./2. + 2 * AL1)
*
      RETURN
      END
*   
* =====================================================================
*     P_{qg}
! ---------------------------------------------------------------------
*
*
*
      FUNCTION Pqg0A(X, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      TR = 1./2.
*
      Pqg0A = 4.*TR*FLAVOR*(X**2. + (1. - X)**2.)
*
      RETURN
      END
*
* =====================================================================
*     P_{gq}
! ---------------------------------------------------------------------
*
*
*
      FUNCTION Pgq0A(x)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      CF = 4./3.
*
      Pqg0A = 2.*CF*(1. + (1 - X)**2.) / X
*
      RETURN
      END
*
* =====================================================================
*     P_{gg}
! ---------------------------------------------------------------------
*
*
*
      FUNCTION Pgg0A(X)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      CA = 3.
*
      Pgg0A = 4.*CA*(1./X - 2. + X * (1. - X))
*
      RETURN
      END
*
! ---------------------------------------------------------------------
*
*
*
      FUNCTION Pgg0B(X)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      CA = 3.
*
      Pgg0B = 4.*CA*(1./(1. - X)) 
*
      RETURN
      END
*
! ---------------------------------------------------------------------
*
*
*
      FUNCTION Pgg0C(X, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      CA = 3.
      TR = 1./2.
      AL1 = LOG(1.-X)
*
      Pgg0C = 2.*(11.*CA - 4.*FLAVOR*TR)/6. + 4.*CA*AL1 
*
      RETURN
      END
*
* =====================================================================
*     NLO
* =====================================================================
*     P_{qq}^{NS}
! ---------------------------------------------------------------------
*
      FUNCTION PqqVns1A(X, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA PI,PI2/3.14159,9.8696/
*
      CF = 4./3.
      CA = 3.
      TR = 1./2.
      TF = TR*FLAVOR
*
      AL = LOG(X)
      AL1 = LOG(1.-X)
*      
      PQQ = 2./(1.-X) - 1.-X
      PQQA = -1.-X

      term1 = CF**2.*(-(2.*AL*AL1 + (3./2.)*AL)*PQQ - ((3./2.) 
     .+ (7./2.)*X)*AL - (1./2.)*(1.+X)*AL**2. - 5.*(1.-X))
*
      term2 = CF*CA*(((1./2.)*AL**2. + (11./6.)*AL)*PQQ + ((67./18.) 
     .- PI2/6.)*PQQA + (1.+X)*AL + (20./3.)*(1.-X))
*
      term3 = CF*TF*(-(2./3.)*AL*PQQ - (10./9.)*PQQA - (4./3.)*(1.-X))
*
      PqqVns1A = 4.*(term1 + term2 + term3)
*
      RETURN
      END
*
! ---------------------------------------------------------------------
*
      FUNCTION PqqVbarns1A(X, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA PI,PI2/3.14159,9.8696/
*
      CF = 4./3.
      CA = 3.
      TR = 1./2.
      TF = TR*FLAVOR
*
      AL = LOG(X)
      AL1 = LOG(1.-X)
*      
      PQQ_M = 2./(1.+X) - 1.+X
*
      PqqVbarns1A = 4.*CF*(CF - CA/2.)*(2.*PQQ_M*S2(X) + 2.*(1.+X)*AL 
     .+ 4.*(1.-X))
*
      RETURN
      END
*
! ---------------------------------------------------------------------
*
      FUNCTION PqqVns1B(X, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA PI,PI2/3.14159,9.8696/
*
      CF = 4./3.
      CA = 3.
      TR = 1./2.
      TF = TR*FLAVOR
*      
      PQQB = 2./(1.-X)
*
      PqqVns1B = 4.*(CF*CA*(67./18. - PI2/6.)*PQQB + CF*TF*(-10./9.)*PQQB)
*
      RETURN
      END
*
! ---------------------------------------------------------------------
*
      FUNCTION Pqqns1C(X, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA PI,PI2/3.14159,9.8696/

      CF = 4./3.
      CA = 3.
      TR = 1./2.
      TF = TR*FLAVOR
      Z3 = 1.2021
      
      AL1 = LOG(1.-X)

      term1 = 4.*(CF*CA*(67./18. - PI2/6.)*2.*AL1
     .+ CF*TF*(-10./9.)*2.*AL1)
*
      term2 = 4.*(CF**2.*((3./8.) - PI2/2. + 6.*Z3) 
     .+ CF*CA*(17./24. + (11./18.)*PI2 - 3.*Z3) - CF*TF*((1./6.) 
     .+ (2./9.)*PI2))
*
      Pqqns1C = term1 + term2

      RETURN
      END
*
* =====================================================================
*     P_{qq}^{PS}
* ---------------------------------------------------------------------
*
*
*
      FUNCTION Pqqps1A(X, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA PI,PI2/3.14159,9.8696/
*
      CF = 4./3.
      CA = 3.
      TR = 1./2.
      TF = TR*FLAVOR
*
      AL = LOG(X)
      AL1 = LOG(1.-X)
*
      Pqqps1A = 4.*(CF*TF*(-16./3. + (40./3.)*X + (10.*X + (16./3.)*X**2. 
     .+ 2.)*AL - (112./9.)*X**2. + (40./(9.*X)) - 2.*(1.+X)*AL**2. 
     .+ (4./3.)*(1.-X)))
*
      RETURN
      END
*
* =====================================================================
*     P_{qg}
! ---------------------------------------------------------------------
*
*
*
      FUNCTION Pqg1A(Z, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA PI,PI2/3.14159,9.8696/
*
      CF = 4./3.
      CA = 3.
      TR = 1./2.
      TF = TR*FLAVOR
*
      ALZ = LOG(Z)
      AL1 = LOG(1.-Z)
*
      Pqg1A = -4.*CF*TF*(-14.+29.*Z-20.*Z*Z-(3.-4.*Z+8.*Z*Z)*ALZ-
     .(1.-2.*Z+4.*Z*Z)*ALZ**2-8.*Z*(1.-Z)*AL1+2.*(1.-2.*Z+2.*Z*Z)*
     .(2.*ALZ*AL1+PI2/3.-AL1**2))-4.*CA*TF*(-2./9.*(20./Z-18.+225.*Z-
     .218.*Z*Z)-2./3.*(3.+24.*Z+44.*Z*Z)*ALZ+(3.+6.*Z+2.*Z*Z)*
     .ALZ**2+8.*Z*(1.-Z)*AL1+(1.-2.*Z+2.*Z*Z)*(2*AL1**2-PI2/3.)
     .-2.*(1.+2.*Z+2.*Z*Z)*S2(Z))
*
      RETURN
      END
*
* =====================================================================
*     P_{gq}
! ---------------------------------------------------------------------
*
*
*
      FUNCTION Pgq1A(X, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA PI,PI2/3.14159,9.8696/
*
      CF = 4./3.
      CA = 3.
      TR = 1./2.
      TF = TR*FLAVOR
*
      AL = LOG(X)
      AL1 = LOG(1.-X)
*
      Pgq1A = -4.*(4./3.*CF*TF*(2./3.*(5./X - 5. + 4.*X)+(2./X - 2. + X)*AL1)
     .+ (1./2.)*CF*CF*(5. + 7.*X-(4. + 7.*X)*AL + 2.*(6./X - 6. + 5*X)*AL1
     .+ (2. - X)*AL**2 + 2.*(2./X - 2. + X)*AL1**2.) + CA*CF*(-1./9.*(9./X + 19.
     .+ 37.*X + 44.*X*X)+1./3.*(36. + 15.*X + 8.*X*X)*AL-.5*(2./X + 6.
     .+ 3.*X)*AL**2. - 1./3.*(22./X - 22. + 17.*X)*AL1 + (2./X - 2. + X)*(
     .2.*AL*AL1 + PI2/6. - AL1**2.) + (2./X + 2. + X)*S2(X)))
*
      RETURN
      END
*
* =====================================================================
*     P_{gg}
! ---------------------------------------------------------------------
*
*
*
      FUNCTION Pgg1A(X, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA PI,PI2/3.14159,9.8696/
*
      CF = 4./3.
      CA = 3.
      TR = 1./2.
      TF = TR*FLAVOR
*
      AL = LOG(X)
      AL1 = LOG(1.-X)

      PGG = 1./(1.-X) + (1./X) - 2. + X*(1.-X)
      PGG_M = 1./(1.+X) - (1./X) - 2. - X*(1.+X)
      PGGA = (1./X) - 2. + X*(1.-X)
*
      term1 = CF*TF*(-16. + 8.*X + (20./3.)*X**2. + (4./3.)*(1./x) - (6. 
     .+ 10.*X)*AL - (2. + 2.*X)*AL**2.)
*
      term2 = CA*TF*(2. - 2.*X + (26./9.)*(X**2. - (1./X)) - (4./3.)*(1.+X)*AL 
     .- (20./9.)*PGGA)
*
      term3 = CA**2.*((27./2.)*(1.-X) + (67./9.)*(X**2. - (1./X)) - ((25./3.) 
     .- (11./3.)*X + (44./3.)*X**2.)*AL + 4.*(1. + X)*AL**2. + 2.*PGG_M*S2(x) 
     .+ (-4.*AL*AL1 + AL**2.)*PGG + ((67./9.) - PI2/3.)*PGGA)
*
      Pgg1A = 4.*(term1 + term2 + term3)
*
      RETURN
      END
*
! ---------------------------------------------------------------------
*
*
*
      FUNCTION Pgg1B(X, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA PI,PI2/3.14159,9.8696/
*
      CF = 4./3.
      CA = 3.
      TR = 1./2.
      TF = TR*FLAVOR
*
      AL = LOG(X)
      AL1 = LOG(1.-X)

      PGGB = 1./(1.-X)
*
      term2b = CA*TF*(- (20./9.)*PGGB)
*
      term3b = CA**2.*(((67./9.) - PI2/3.)*PGGB)
*
      Pgg1B = 4.*(term2b + term3b)
*
      RETURN
      END
*
! ---------------------------------------------------------------------
*
*
*
      FUNCTION Pgg1C(X, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA PI,PI2/3.14159,9.8696/
*
      CF = 4./3.
      CA = 3.
      TR = 1./2.
      TF = TR*FLAVOR
      Z3 = 1.2021
*
      AL = LOG(X)
      AL1 = LOG(1.-X)
*
      term2b = CA*TF*(- (20./9.)*AL1)
*
      term3b = CA**2.*(((67./9.) - PI2/3.)*AL1)
*
      Pgg1C = 4.*(term2b + term3b) + 4.*(CA**2.*((8./3.) + 3.*Z3) 
     .- CF*TF - (4./3.)*CA*TF)
*
      RETURN
      END
*
* =====================================================================
*     NNLO
* =====================================================================
*     P_{qq}^{NS}
! ---------------------------------------------------------------------
C ===================================================================

*
* ..File: xpns2p.f 
*
*                           __
* ..The parametrized 3-loop MS non-singlet splitting functions P^(2)
*    for the evolution of unpolarized partons densities, mu_r = mu_f.
*    The expansion parameter is alpha_s/(4 pi).
*
* ..The distributions (in the mathematical sense) are given as in eq.
*    (B.26) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*    The name-endings A, B, and C of the functions below correspond to
*    the kernel superscripts [2], [3], and [1] in that equation.
*
*  ..The relative accuracy of these parametrizations, as well as of
*    the convolution results, is better than one part in thousand.
*
* ..References: S. Moch, J. Vermaseren and A. Vogt,
*               hep-ph/0209100 = Nucl. Phys. B646 (2002) 181,
*               hep-ph/0403192 (submitted to Nucl. Phys. B)
*
* =====================================================================
*
*
* ..This is the regular piece of P2_NS+.  The rational coefficients are
*    exact, the rest has been fitted for x between 10^-6 and 1 - 10^-6.
*    The N_f^2 part is exact and was first determined in N-space by 
*    J.A. Gracey in Phys. Lett. B322 (1994) 141.
*
      FUNCTION P2NSPA (Y, NF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      INTEGER NF
*
      DL  = LOG (Y)
      DL1 = LOG (1.-Y)
      D81 = 1./81.D0
*
      P2NSPA =   1641.1 - 3135.* Y + 243.6 * Y**2 - 522.1 * Y**3
     ,            + 128.*D81 * DL**4 + 2400.*D81 * DL**3
     ,            + 294.9 * DL**2 + 1258.* DL
     ,            + 714.1 * DL1 + DL*DL1 * (563.9 + 256.8 * DL)
     ,        + NF * ( -197.0 + 381.1 * Y + 72.94 * Y**2 + 44.79 * Y**3
     ,            - 192.*D81 * DL**3  - 2608.*D81 * DL**2 - 152.6 * DL
     ,            - 5120.*D81 * DL1 - 56.66 * DL*DL1 - 1.497 * Y*DL**3 )
     ,        + NF**2 * ( 32.* Y*DL/(1.-Y) * (3.* DL + 10.) + 64.
     ,            + (48.* DL**2 + 352.* DL + 384.) * (1.-Y) ) * D81
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
* ..This is the regular piece of P2_NS-.  The rational coefficients are 
*    exact, the rest has been fitted for x between 10^-6 and 1 - 10^-6.
*    The N_f^2 part is exact (and identical to that of P2_NS+). 
*
      FUNCTION P2NSMA (Y, NF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      INTEGER NF
*
      DL  = LOG (Y)
      DL1 = LOG (1.-Y)
      D81 = 1./81.D0
*
      P2NSMA =   1860.2 - 3505.* Y + 297.0 * Y**2 - 433.2 * Y**3
     ,            + 116.*D81 * DL**4 + 2880.*D81 * DL**3 
     ,            + 399.2 * DL**2 + 1465.2 * DL
     ,            + 714.1 * DL1 + DL*DL1 * (684.0 + 251.2 * DL)
     ,        + NF * ( -216.62 + 406.5 * Y + 77.89 * Y**2 + 34.76 * Y**3
     ,            - 256.*D81 * DL**3  - 3216.*D81 * DL**2 - 172.69 * DL 
     ,            - 5120.*D81 * DL1 - 65.43 * DL*DL1 - 1.136 * Y*DL**3 )
     ,        + NF**2 * ( 32.* Y*DL/(1.-Y) * (3.* DL + 10.) + 64.
     ,            + (48.* DL**2 + 352.* DL + 384.) * (1.-Y) ) * D81
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
* ..This is the singular piece of both P2_NS+ and P2_NS-. It is exact 
*    up to the truncation of the irrational coefficients.
*
      FUNCTION P2NSB (Y, NF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      INTEGER NF
*
      P2NSB = ( 1174.898 - NF * 183.187 - NF**2 * 64./81.D0 ) / (1.-Y)
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece of P2_NS+. The coefficients of delta(1-x)
*    have been partly shifted relative to the exact (truncated) values.
*
      FUNCTION P2NSPC (Y, NF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      INTEGER NF
*
      DL1 = LOG (1.-Y)
*
      P2NSPC =       1174.898 * DL1 + 1295.624 - 0.24
     ,        - NF * ( 183.187 * DL1 + 173.938 - 0.011 )
     ,        + NF**2 * ( - 64./81.D0 * DL1 + 1.13067 )
*
      RETURN
      END
*
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece of P2_NS-. The coefficients of delta(1-x) 
*    have been partly shifted relative to the exact (truncated) values.
*
      FUNCTION P2NSMC (Y, NF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      INTEGER NF
*
      DL1 = LOG (1.-Y)
*
      P2NSMC =       1174.898 * DL1 + 1295.624 - 0.154
     ,        - NF * ( 183.187 * DL1 + 173.938  - 0.005 )
     ,        + NF**2 * ( - 64./81.D0 * DL1 + 1.13067 )
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
* ..This is P2_NSS, the difference of P2_NSV and P2_NS-.
*
      FUNCTION P2NSSA (Y, NF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

*
      INTEGER NF
*
      D27 = 1./27.D0
      DL  = LOG (Y)
      Y1  = 1.- Y
      DL1 = LOG (Y1)
*
      P2NSSA = Y1* ( 151.49 + 44.51 * Y - 43.12 * Y**2 + 4.820 * Y**3 )
     1          + 40.*D27 * DL**4 - 80.*D27 * DL**3 + 6.892 * DL**2 
     2          + 178.04 * DL + DL*DL1 * ( - 173.1 + 46.18 * DL )
     4          + Y1*DL1 * ( - 163.9 / Y - 7.208 * Y ) 
*
      P2NSSA  = NF * P2NSSA
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
*
*
      function p2nsa(y,ipm,nf)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      ipm2=ipm*ipm
      p2nsa = 0.d0              ! G.W. 07/05/2007
c$$$      if(ipm2.ne.1) go to 1
      if(ipm2.ne.1.and.ipm.ne.-2) go to 1 ! G.W. 04/07/2007
      if(ipm.eq.1) p2nsa=p2nspa(y,nf)
      if(ipm.eq.-1) p2nsa=p2nsma(y,nf)
c$$$      if(ipm.eq.-2) p2nsa=p2nsma(y,nf)+p2nssa(y,nf) ! G.W. 04/07/2007
      if(ipm.eq.-2) p2nsa=p2nsma(y,nf)+p2nssa(y,2) ! G.W. 05/07/2007
      return
    1 print 99,ipm
      stop
   99 format(1x,'ipm=',i3)
      end
*
* ---------------------------------------------------------------------
*
*
*
*
      function p2nsc(y,ipm,nf)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      ipm2=ipm*ipm
      p2nsc = 0.d0              ! G.W. 07/05/2007
c$$$      if(ipm2.ne.1) go to 1
      if(ipm2.ne.1.and.ipm.ne.-2) go to 1 ! G.W. 04/07/2007
      if(ipm.eq.1) p2nsc=p2nspc(y,nf)
      if(ipm.eq.-1) p2nsc=p2nsmc(y,nf)
      if(ipm.eq.-2) p2nsc=p2nsmc(y,nf) ! G.W. 04/07/2007
      return
    1 print 99,ipm
      stop
   99 format(1x,'ipm=',i3)
      end
*
* =====================================================================
*     P_{qq}^{PS}
* ---------------------------------------------------------------------
* =================================================================av==
*
* ..File: xpij2p.f 
*
*                           __
* ..The parametrized 3-loop MS singlet splitting functions P^(2) for 
*    the evolution of unpol. singlet parton densities at mu_r = mu_f.
*    The expansion parameter is alpha_s/(4 pi).
*
* ..The distributions (in the mathematical sense) are given as in eq.
*   (B.27) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*   The name-endings A, B, and C of the functions below correspond to 
*   the kernel superscripts [2], [3], and [1] in that equation.
*
* ..The relative accuracy of these parametrisations, as well as of
*    the convolution results, is better than one part in thousand.

* ..The coefficients of 1/(1-x)_+, (ln x)/x and 1/x are exact (up
*    to a truncation of irrational coefficients).  Furthermore all
*    coefficients written as fractions (e.g., 160./27.D0) are exact.
*    The other terms at x < 1 have fitted to the exact results for x 
*    between 10^-6 and 1 - 10^-6.  The coefficient of delta(1-x) of
*    P_gg^(2) have been slightly adjusted using the second moments.
*
* ..References: S. Moch, J. Vermaseren and A. Vogt,
*               hep-ph/0403192 (to appear in Nucl. Phys. B)
*               A. Vogt, S. Moch and J. Vermaseren,
*               hep-ph/0404111 (submitted to Nucl. Phys. B)
* 
* =====================================================================
*
*
* ..The (regular) pure-singlet splitting functions P_ps^(2).
*    P_qq^(2) is obtained by adding the non-singlet quantity P_NS^(2)+.
*    A parametrization of the latter is provided in the file  xpns2p.f.

      FUNCTION P2PSA (Y, NF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      INTEGER NF
*
      DL  = LOG (Y)
      DL1 = LOG (1.-Y)
*
      P2PS1 = - 3584./(27.D0*Y) * DL - 506.0/ Y + 160./27.D0 * DL**4
     ,         - 400./9.D0 * DL**3 + 131.4 * DL**2 - 661.6 * DL
     ,         - 5.926  * DL1**3 - 9.751 * DL1**2 - 72.11 * DL1
     ,         + 177.4 + 392.9 * Y - 101.4 * Y**2 - 57.04 * DL*DL1
      P2PS2 =   256./(81.*Y) + 32./27.D0 * DL**3 + 17.89 * DL**2
     ,         + 61.75 * DL + 1.778 * DL1**2 + 5.944 * DL1 + 100.1
     ,         - 125.2 * Y + 49.26 * Y**2 - 12.59 * Y**3 
     ,         - 1.889 * DL*DL1 
*
      P2PSA = (1.-Y) * ( P2PS1 + NF * P2PS2 )
*
      RETURN
      END
*
* =====================================================================
*     P_{qg}
! ---------------------------------------------------------------------
*
*
* ..The gluon->quark splitting functions P_qg^(2).
*
      FUNCTION P2QGA (Y, NF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

*
      INTEGER  NF
*
      DL  = LOG (Y)
      DL1 = LOG (1.-Y)
*
      P2QG1 = - 896./(3.D0*Y) * DL - 1268.3 / Y + 536./27.D0 * DL**4 
     ,         - 44./3.D0 * DL**3 + 881.5 * DL**2 + 424.9 * DL 
     ,         + 100./27.D0 * DL1**4 - 70./9.D0 * DL1**3 
     ,         - 120.5 * DL1**2 + 104.42 * DL1
     ,         + 2522. - 3316.* Y + 2126.* Y**2
     ,         + DL*DL1 * (1823. - 25.22 * DL) - 252.5 * Y*DL**3  
      P2QG2 =   1112./(243.D0*Y) - 16./9.D0 * DL**4 
     ,         - 376./27.D0 * DL**3 - 90.8 * DL**2 - 254.0 * DL   
     ,         + 20./27.D0 * DL1**3 + 200./27.D0 * DL1**2 - 5.496 * DL1
     ,         - 252.0  + 158.0 * Y + 145.4 * Y**2 - 139.28 * Y**3
     ,         - DL*DL1 * ( 53.09  + 80.616 * DL) - 98.07 * Y*DL**2
     ,         + 11.70 * Y*DL**3
* 
      P2QGA = ( P2QG1 + NF * P2QG2 )
*
      RETURN
      END
*
* =====================================================================
*     P_{gq}
! ---------------------------------------------------------------------

*
*
* ..The quark->gluon splitting functions P_gq^(2).  P2GQ2 is exact.
*
      FUNCTION P2GQA (Y, NF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

*
      INTEGER NF
*
      DL  = LOG (Y)
      DL1 = LOG (1.-Y)
*
      P2GQ0 =   1189.3 * DL/Y + 6163.1 / Y - 4288./81.D0 * DL**4
     ,         + 1568./9.D0 * DL**3 - 1794. * DL**2 + 4033. * DL
     ,         + 400./81.D0 * DL1**4 + 2200./27.D0 * DL1**3
     ,         + 606.3 * DL1**2 + 2193.* DL1 
     ,         - 4307. + 489.3 * Y + 1452.* Y**2 + 146.0 * Y**3
     ,         - 447.3 * DL**2*DL1 - 972.9 * Y*DL**2
      P2GQ1 =   71.082 * DL/Y  - 46.41 / Y + 128./27.D0 * DL**4
     ,         + 704/81.D0 * DL**3 + 20.39  * DL**2 + 174.8 * DL
     ,         - 400./81.D0 * DL1**3 - 68.069 * DL1**2 - 296.7 * DL1
     ,         - 183.8 + 33.35 * Y - 277.9 * Y**2 + 108.6 * Y*DL**2
     ,         - 49.68 * DL*DL1
      P2GQ2 = ( 64. * ( - 1./Y + 1. + 2.* Y)
     ,         + 320.* DL1 * ( 1./Y - 1. + 0.8 * Y)
     ,         + 96.* DL1**2 * ( 1./Y - 1. + 0.5 * Y) ) / 27.D0
*
      P2GQA = ( P2GQ0 + NF * (P2GQ1 + NF * P2GQ2) )
*
      RETURN
      END
*
* =====================================================================
*     P_{gg}
! ---------------------------------------------------------------------
*
*
* ..The regular piece of the gluon-gluon splitting function P_gg^(2).
*
      FUNCTION P2GGA (Y, NF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

*
      INTEGER NF
*
      DL  = LOG (Y)
      DL1 = LOG (1.-Y)
*
      P2GGA0 = 2675.8 * DL/Y + 14214./ Y - 144. * DL**4 + 72. * DL**3
     1          - 7471. * DL**2 + 274.4 * DL + 3589. * DL1 - 20852. 
     2          + 3968.* Y - 3363. * Y**2 + 4848. * Y**3 
     3          + DL*DL1 * ( 7305. + 8757. * DL )
      P2GGA1 = 157.27 * DL/Y + 182.96 / Y + 512./27.D0 * DL**4
     1          + 832./9.D0 * DL**3 + 491.3 * DL**2 + 1541. * DL
     2          - 320.0 * DL1 - 350.2 + 755.7 * Y - 713.8 * Y**2 
     3          + 559.3 * Y**3 + DL*DL1 * ( 26.15 - 808.7 * DL )
      P2GGA2 = - 680./(243.D0 * Y) - 32./27.D0 * DL**3 + 9.680 * DL**2
     1          - 3.422 * DL - 13.878 + 153.4 * Y - 187.7 * Y**2 
     2          + 52.75 * Y**3 - DL*DL1 * (115.6 - 85.25* Y + 63.23* DL)
*
      P2GGA = P2GGA0 + NF * ( P2GGA1 + NF * P2GGA2 )
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
* ..The singular piece of the gluon-gluon splitting function P_gg^(2).
*
      FUNCTION P2GGB (Y, NF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

*
      INTEGER NF
*
      P2GGB = ( 2643.521 - NF * 412.172 - NF**2 * 16./9.D0 ) / ( 1.-Y) 
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
* ..The 'local' piece of the gluon-gluon splitting function P_gg^(2).  
*
      FUNCTION P2GGC (Y, NF)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

*
      INTEGER NF
*
      DL1 = LOG (1.-Y)
*
      P2GGC =       2643.521 * DL1 + 4425.448 + 0.446
     ,       - NF * ( 412.172 * DL1 +  528.720 + 0.003 )
     ,       + NF**2 * ( - 16./9.D0 * DL1 + 6.4630)
*
      RETURN
      END
*
* =====================================================================
* =====================================================================
*     N3LO
* =====================================================================
*     
*     This code is provided as supplementary material to the MSHT20 
*     approximate N3LO PDFs.
*     
*     These approximations are not exact and allow for a level of 
*     variation by the input parameter "a" in all cases. These 
*     subroutines here use a simple method of simultaneous equations in
*     order to produce the approximations and therefore may well 
*     be superseded in the future by more robust and equally/more 
*     efficient approximations or indeed, the exact functions when 
*     available.
*
*     Please refer to arXiv:2207.XXXX for the full details of how these
*     approximations are determined.
*
*     Note that these can and will be continually updated with any  
*     extra exact calculations that are made available, which will 
*     increase the accuracy of these approximations.
*
*     When using these approximations, please reference the following.
*
*     References:
*     J. McGowan, et. al., arXiv:2207.XXXXXX
*     S. Moch, et. al., arXiv:1707.08315
*     A. Vogt, et. al., arXiv:1808.08981
*     S. Catani, et. al., hep-ph/9405388
*     L. N. Lipatov, Sov. J. Nucl. Phys. (1976)
*     E. A. Kuraev, et. al., JETP (1977)
*     I. I. Balitsky, et. al., Sov. J. Nucl. Phys. (1978)
*     V. S. Fadin, et. al., hep-ph/9802290
*     M. Ciafaloni, et. al., hep-ph/9803389
*     
*     
* ---------------------------------------------------------------------
*     P_{qq}^{NS}
* ---------------------------------------------------------------------
*
*     Contribution from approximate four-loop splitting function with
*     variational parameter (0 < aNS < 0.014). This is a subroutine and
*     outputs a regular and plus distribution part of the approximation.
*
      SUBROUTINE Pqqns3PAB(x, a, FLAVOR, PnsN3LOresult, PnsN3LOplus)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER n_dimension
      DOUBLE PRECISION N(8), C(8), N_exact(8), N_a(8), inv_A(8, 8)
      DATA PI,PI2/3.14159,9.8696/
*
      n_dimension=8
* 
*     Moments for the known exact small-x and large-x contributions (Vogt)
      N_exact(1) = -0.00021211D0 - FLAVOR*(0.1654801313006386D0 + 0.20740623526104135D0
     .-0.1502823524D0 + 0.01076711056D0 - 0.2333500675D0) 
     .- FLAVOR**2.*(-0.004654314166031352D0 -0.007842937752425193D0
     .+ 4.7723213922D-03 + 7.774039377D-03) 
     .- FLAVOR**3.*(-6.62104090137594D-05 - 0.0001312258806063389D0 
     .+ 7.653104787D-05 + 1.209052801D-04)

      N_exact(2) = -0.00021141D0 - FLAVOR*(0.11120681844515805D0 + 0.3802447646452425D0 
     .-0.2688473229D0 + 0.01076711056D0 - 0.2333500675D0) 
     .- FLAVOR**2.*(9.571697043e-3 + 7.774039377D-03 -0.0029670140315604917D0
     .-0.014378719212779521D0) 
     .- FLAVOR**3.*(-3.609651149299798D-05 -0.00024058078111162135D0 
     .+ 1.55772043D-04 + 1.209052801D-04)

      N_exact(3) = -0.00021069D0 - FLAVOR*(0.08533032993183774D0 + 0.4735775705127111D0 
     .-0.3363034306D0 + 0.01076711056D0 - 0.2333500675D0) 
     .- FLAVOR**2.*(0.01237066021 + 7.774039377D-03 
     .- 0.002236654071582081D0 -0.01790804120137086D0) 
     .- FLAVOR**3.*(-2.4727034394311423D-05 -0.0002996324273844739D0 
     .+ 2.034542201D-04 + 1.209052801D-04)

      N_exact(4) = -0.00020986D0 - FLAVOR*(0.06978655862237221D0 + 0.5377747385697003D0
     .-0.3849568093D0 + 0.01076711056D0 - 0.2333500675D0) 
     .- FLAVOR**2.*(-0.001806806533621782D0 -0.020335617172359613D0 
     .+ 0.01436839622 + 7.774039377D-03) 
     .- FLAVOR**3.*(-1.875753442493042D-05 -0.00034024996185786445D0 
     .+ 2.381022604D-04 + 1.209052801D-04)

      N_exact(5) = -0.00020809D0 - FLAVOR*(0.059311724466073266D0 + 0.5867456552285572D0 
     .-0.4234530704D0 + 0.01076711056D0 - 0.2333500675D0) 
     .- FLAVOR**2.*(-0.0015200868207649136D0 -0.02218742191946D0 
     .+ 0.01593349271 + 7.774039377D-03) 
     .- FLAVOR**3.*(-1.5092063174265756D-05 - 0.0003712338503343612 
     .+ 2.654206822D-04 + 1.209052801D-04)

      N_exact(6) = -0.0002073D0 - FLAVOR*(0.0517320184770643D0 + 0.6263413910511195D0 
     .-0.4554694245D0 + 0.01076711056D0 - 0.2333500675D0) 
     .- FLAVOR**2.*(-0.001314205205204153D0 -0.023684710035832085D0 
     .+ 0.0172249122 + 7.774039377D-03) 
     .- FLAVOR**3.*(-1.2617371403920047D-05 -0.0003962860639046623D0 
     .+ 2.879982076D-04 + 1.209052801D-04)

      N_exact(7) = -0.00020562D0 - FLAVOR*(0.045972268757607376D0 + 0.6595795697788507D0 
     .-0.4829482701D0 + 0.01076711056D0 - 0.2333500675D0) 
     .- FLAVOR**2.*(-0.0011588022146498728D0 -0.024941591085900227D0 
     .+ 0.0183264033 + 7.774039377D-03) 
     .- FLAVOR**3.*(-1.0836262424080888D-05 -0.0004173158524633704D0 + 3.072468902D-04 + 1.209052801D-04)

      N_exact(8) = -0.00020381D0 - FLAVOR*(0.04143575048971359D0 + 0.6882213832196611D0 
     .-0.5070540331D0 + 0.01076711056D0 - 0.2333500675D0) 
     .- FLAVOR**2.*(-0.0010371493846056854D0 -0.026024663442187512D0 
     .+ 0.01928783511D0 + 7.774039377D-03) 
     .- FLAVOR**3.*(-9.493881167514427D-06 -0.00043543752168996016D0 
     .+ 3.240261808D-04 + 1.209052801D-04)
* 
*     Moments for the varied function multiplied by a (ln^3(x))
      N_a(1) = -0.375D0
      N_a(2) = -0.0234375D0
      N_a(3) = -0.004629629629528962D0
      N_a(4) = -0.0014648437499999065D0
      N_a(5) = -0.0005999999999999996D0
      N_a(6) = -0.0002893518518518518D0
      N_a(7) = -0.00015618492294877136D0
      N_a(8) = -9.155273437499997D-05
*      
      DO i=1, n_dimension
        N(i) = N_exact(i) - a * N_a(i)
      ENDDO
*
      inv_A = reshape((/9.67452285D-02, 7.12205880D03, 2.55568773D03, 
     &                  4.41612871D02, -3.63251998D02, 8.54599881D03, 
     &                  -8.03082218D03, 5.66780493D01, -1.10377154D01,
     &                  -7.11199538D05, -2.58657431D05, -4.59602884D04,  
     &                  2.63108547D04, -8.24798082D05, 7.81924565D05,  
     &                  -3.25280541D03, 1.99489100D02, 1.15489637D07,  
     &                  4.23773320D06, 7.69318756D05, -3.54975874D05,  
     &                  1.31604354D07, -1.25178270D07, 4.09448382D04,  
     &                  -1.28368247D03, -6.79475656D07, -2.50856870D07,  
     &                  -4.63268153D06, 1.85137884D06, -7.66032630D07,  
     &                  7.29666144D07, -2.06388099D05, 3.82370653D03,  
     &                  1.87306734D08, 6.94570801D07, 1.30090317D07,  
     &                  -4.68326682D06, 2.09632375D08, -1.99805999D08,  
     &                  5.11458122D05, -5.75484319D03, -2.63212244D08,  
     &                  -9.79232841D07, -1.85597818D07, 6.16625773D06,  
     &                  -2.93018936D08, 2.79351504D08, -6.64139088D05,  
     &                  4.25513991D03, 1.82949771D08, 6.82322632D07,  
     &                  1.30646558D07, -4.07104678D06, 2.02835925D08,  
     &                  -1.93381943D08, 4.34087442D05, -1.22905233D03,  
     &                  -4.99408525D07, -1.86617424D07, -3.60495568D06,  
     &                  1.06569809D06, -5.51895151D07, 5.26130529D07,  
     &                  -1.12766384D05/), shape=(/8,8/)) 
* 
*     Matrix multiplication of inv_A and N
*
      DO i=1,n_dimension
        C(i) = 0.D0
        DO j=1,n_dimension
          C(i) = C(i) + inv_A(i, j) * N(j)
        ENDDO
      ENDDO
*
      PnsN3LOplus = (C(1) * (1./(1.-X))) * (4*PI)**4.
      PnsN3LOresult = (C(2) * (1.-X) * LOG(1.-X) + C(3) * (1.-X) 
     &          * (LOG(1.-X)**2.) + C(4) * (1.-X) * (LOG(1.-X)**3.) 
     &          + C(5) + C(6) * X + C(7) * (X**2.) + C(8) * (LOG(X)**2.)
     &          + (a) * (LOG(X)**3.)) * (4*PI)**4.

*
      RETURN
      END
*
* -------------------------------------------------------------------------
*
*     This is the regular exact part of P3NS.
*
      FUNCTION Pqqns3PexactA(x, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA PI,PI2/3.14159,9.8696/
*
      CA = 3.0D0
      CF = 4.0D0/3.0D0
      AL1 = LOG(1.-X)
*
      Pqqns3PexactA = (19069.80D0 * AL1 - 29733.85D0 + 25000D0*((1 - x) *
     &          (3.5242 + 8.3679D0 * x - 1.2395D0 * (x**2.) + 1.7423D0 * 
     &          (x**3.)) + 11.916D0 * x * log(x) - 0.2237D0 * x * 
     &          (log(x)**2.) - 0.0129D0 * x * (log(x)**3.) +11.937D0 * 
     &          (1. - x) * AL1 + 13.955D0 * log(x) * AL1) + 51671.33D0 * 
     &          log(x)+17120.95D0*(log(x)**2.)+2863.226D0 * (log(x)**3.) 
     &          + 297.8255D0*(log(x)**4.)+16.D0*(log(x)**5.)+(1.D0/2.D0) 
     &          * (log(x)**6.) + FLAVOR * (- 3079.761D0 * AL1 + 8115.605D0 
     &          + 2500.D0 * ((1. - x) * (-7.4077D0+4.5141D0*x-1.0069D0 * 
     &          (x**2.) + 0.7641D0 * (x**3.)) + 8.4211D0 * (1 - x) * AL1 
     &          + 7.5633D0*log(x)*AL1 + 7.5236D0 * x * log(x)+0.2208D0 
     &          * x * (log(x)**2.) + 0.05712D0 * x * (log(x)**3.)) 
     &          - 9239.374D0*log(x)-2917.312D0*(log(x)**2.)-430.5308D0 
     &          *(log(x)**3.)- 36.D0*(log(x)**4.) -(4.D0 /3.D0) * 
     &          (log(x)**5.)) + (FLAVOR**2.) * (26.68861D0 *AL1 - 376.0092D0 
     &          + 250.D0*((1. - x)*(3.0008D0 + 0.8619D0*x - 0.12411D0 * 
     &          (x**2.) +0.31595D0 * (x**3.)) - 0.37529D0 * x * log(x) 
     &          - 0.21684D0*x*(log(x)**2.) -0.02295D0 * x * (log(x)**3.)
     &          + (0.03394D0 * (1. - x) + 0.40431D0 * log(x)) * AL1) 
     &          + 393.0056D0 * log(x) + 112.5705D0 * (log(x)**2.) 
     &          + 16.52675D0 * (log(x)**3.)+0.7901235D0 * (log(x)**4.)) 
     &          + (FLAVOR**3.) * (-2.426296D0-0.8460488D0* x+(0.5267490D0* 
     &          (1./(1.-x)) - 3.687243D0+3.160494D0 * x) * log(x) -
     &          (0.1316872D0 * (10.*(1./(1. - x)) + 1.)-1.448560D0 * x) 
     &          * (log(x)**2.) - (0.2633745D0 *(1./(1.-x))-0.1316870D0 
     &          * (1.+x)) * (log(x)**3.)) - 507.152D0 * (- 1.)-1777.27D0 
     &          * (AL1**2) * (1 - x) - 204.183D0 * (AL1**3.)*(1.-x) + 
     &          3948.16D0 * (1.-x) - 2464.61D0* (2.*x - (x**2.)) 
     &          * (1.-x) - 1839.44D0* (log(x)**2.) - 402.156D0 * 
     &          (log(x)**3.) - 55.87553D0 * (log(x)**4.) - 2.831276D0
     &          * (log(x)**5.) - 0.1488340D0 * (log(x)**6.) - 2601.749D0 
     &          - 2118.867D0 * AL1 + FLAVOR * (7.33927D0 * (-1.) - 143.813D0 
     &          * AL1 * (1. - x) - 18.8803D0 *(AL1**3.) * (1. - x) 
     &          - (1116.34D0 - 1071.24D0 * x) * x * (1. - x)- 59.3041D0 
     &          * (log(x)**2.) - 8.4620D0* (log(x)**3.) + 4.658436D0 * (
     &          log(x)**4.) + 0.2798354D0 * (log(x)**5.) + 312.1643D0 
     &          + 337.9310D0 * AL1))
*
      RETURN
      END
*
* ---------------------------------------------------------------------------
*
*     This is the singular exact part of P3NS.
*
      FUNCTION Pqqns3PexactB(x, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA PI,PI2/3.14159,9.8696/
*
      CA = 3.0D0
      CF = 4.0D0/3.0D0
      AL1 = LOG(1.-X)
*
      Pqqns3PexactB = (21209.02D0 * (1. / (1. - x)) + FLAVOR * (
     &                   -5179.372D0 * (1. / (1. - x))) + (FLAVOR**2.) * 
     &                   (195.5772D0 * (1. / (1. - x))) + (FLAVOR**3.) * 
     &                   (3.272344D0 * (1. / (1. - x))) + (- 507.152D0 
     &                   + FLAVOR * 7.33927D0) * (1. / (1. - x)))
*
      RETURN
      END
*
* ----------------------------------------------------------------------------
*
*     This is the 'local' exact part of P3NS.
*
      FUNCTION Pqqns3PexactC(X, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA PI, PI2/3.14159,9.8696/
*
      CA = 3.0D0
      CF = 4.0D0/3.0D0
      AL1 = LOG(1.-X)
*
      PnsN3LOEXACTCONST = (-2405.03D0+FLAVOR*267.965D0+FLAVOR**2.*(193.8554D0
     &      +0.0037D0)+FLAVOR**3.*(3.014982D0)+(25796.09D0-1.0D0))
*
      Pqqns3PexactC = PnsN3LOEXACTCONST - (21209.02D0*(AL1) + FLAVOR * (
     &                   -5179.372D0 * AL1) + (FLAVOR**2.) * 
     &                   (195.5772D0 * AL1) + (FLAVOR**3.) * 
     &                   (3.272344D0 * AL1) + (- 507.152D0 
     &                   + FLAVOR * 7.33927D0) * AL1)
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
* ..This is the regular piece of P3_NS-. Approximation taken from 
*   arXiv:1707.08315.
*
      FUNCTION P3NSMA (Y, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      DL  = LOG (Y)
      DL1 = LOG (1.-Y)
      D81 = 1./81.D0

      PLNFminusA = FLAVOR**2. * (26.68861*DL1 - 376.0092 
     ,            + 250.0*((1.-y)*(3.2206 
     ,            + 1.7507*y + 0.13281*y**2. + 0.45969*y**3.)
     ,            + 1.5641*y*DL - 0.37902*y*DL**2. - 0.03248*y*DL**3.
     ,            + (2.7511*(1.-y) + 3.2709*DL)*DL1) + 437.8810*DL 
     ,            + 128.2948*DL**2. + 19.59945*DL**3. + 0.9876543*DL**4.)
     ,            + FLAVOR**3.*(- 2.426296
     ,            - 0.8460488*y + (0.5267490*(1./(1.-y)) - 3.687243 
     ,            + 3.160494*y)*DL - (0.1316872*(10.*(1./(1.-y)) + 1.0)
     ,            - 1.448560*y)*DL**2. - (0.2633745*(1./(1.-y)) 
     ,            - 0.1316870*(1.+y))*DL**3.)

      PL0_PL1A = 19069.80D0 * DL1 - 29733.85D0 + 25000D0*((1. - y) *
     ,          (3.5242 + 8.3679D0 * y - 1.2395D0 * (y**2.) + 1.7423D0 * 
     ,          (y**3.)) + 11.916D0 * y * log(y) - 0.2237D0 * y * 
     ,          (log(y)**2.) - 0.0129D0 * y * (log(y)**3.) +11.937D0 * 
     ,          (1. - y) * DL1 + 13.955D0 * log(y) * DL1) + 51671.33D0 * 
     ,          log(y)+17120.95D0*(log(y)**2.)+2863.226D0 * (log(y)**3.) 
     ,          + 297.8255D0*(log(y)**4.)+16.D0*(log(y)**5.)+(1.D0/2.D0) 
     ,          * (log(y)**6.) + FLAVOR * (- 3079.761D0 * DL1 + 8115.605D0 
     ,          + 2500.D0 * ((1. - y) * (-7.4077D0+4.5141D0*y-1.0069D0 * 
     ,          (y**2.) + 0.7641D0 * (y**3.)) + 8.4211D0 * (1. - y) * DL1 
     ,          + 7.5633D0*log(y)*DL1 + 7.5236D0 * y * log(y)+0.2208D0 
     ,          * y * (log(y)**2.) + 0.05712D0 * y * (log(y)**3.)) 
     ,          - 9239.374D0*log(y)-2917.312D0*(log(y)**2.)-430.5308D0 
     ,          *(log(y)**3.)- 36.D0*(log(y)**4.) -(4.D0 /3.D0) * 
     ,          (log(y)**5.))
*
      P3NSMA =   -502.481*(- 1.)
     ,            - 3997.39*DL1*(1.-y) + 511.567*DL1**3.*(1.-y)
     ,            + (4043.59 - 15386.6*y)*y*(1.-y) + 1532.96*DL**2.
     ,            + 31.6023*DL**3. + 0.4964335*(DL**6. + 18*DL**5.)
     ,            - 2601.749 - 2118.867*DL1 + FLAVOR*(7.82077*(
     ,            -1.0) - 1360.04*DL1*(1.-y) 
     ,            + 38.7337*DL1**3.*(1.-y) - (335.995*(2.+y) 
     ,            + 1605.91*y**2)*(1.-y) - 9.76627*DL**2. + 0.14218*DL**5.
     ,            + 312.1643 + 337.9310*DL1) + PL0_PL1A + PLNFminusA
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
* ..This is the singular piece of P3_NS-. Approximation taken from 
*   arXiv:1707.08315.
*
      FUNCTION P3NSMB (Y, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      DL  = LOG (Y)
      DL1 = LOG (1.-Y)
      D81 = 1./81.D0

      PLNFminusB = FLAVOR**2. * (195.5772*(1./(1.-y)))
     ,            + FLAVOR**3.*(3.272344*(1./(1.-y)))

      PL0_PL1B = 21209.02*(1./(1.-y)) + FLAVOR*(-5179.372*(1./(1.-y)))
*
      P3NSMB =   -502.481*((1./(1.-y))) + FLAVOR*(7.82077*(1./(1.-y)))
     ,            + PLNFminusB + PL0_PL1B
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
* ..This is the 'local' piece of P3_NS-. Approximation taken from 
*   arXiv:1707.08315.
*
      FUNCTION P3NSMC (Y, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      DL  = LOG (Y)
      DL1 = LOG (1.-Y)
      D81 = 1./81.D0

      PLNFminusC = FLAVOR**2. * ((193.8554 + 0.0037)) 
     ,            + FLAVOR**3.*(3.014982)

      PL0_PL1C = (25796.09 - 1.0) - FLAVOR*(5818.637 + 0.35)

      P3NSMresidual = - (21209.02D0*(DL1) + FLAVOR * (
     &                   -5179.372D0 * DL1) + (FLAVOR**2.) * 
     &                   (195.5772D0 * DL1) + (FLAVOR**3.) * 
     &                   (3.272344D0 * DL1) + (- 507.152D0 
     &                   + FLAVOR * 7.33927D0) * DL1)
*
      P3NSMC =   - 2380.255 + FLAVOR*(270.468) + PL0_PL1C + PLNFminusC
     ,            + P3NSMresidual
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
* ..This is P3_NSS, the difference of P3_NSV and P3_NS-.  
*   Approximation taken from arXiv:1707.08315.
*
      FUNCTION P3NSSA (Y, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      DL  = LOG (Y)
      DL1 = LOG (1.-Y)
      D81 = 1./81.D0

      P3NSSA = FLAVOR * (60.40*(1.-y)*DL1**2. + 4.685*(1.-y)*DL1**3. 
     ,            + (1.-y)*y*(4989.2 - 1607.73*y) + 3687.6*DL 
     ,            + 3296.6*DL**2. + 1271.11*DL**3. + 533.44*DL**4.
     ,            + 97.27*DL**5. + 4.*DL**6.)
     ,            + FLAVOR**2. * (19.70002*(1.-y)*DL1 
     ,            - 3.435474*(1.-y)*DL1**2. + 250.0*((1.-y)*(
     ,            -4.7656 + 1.6908*y + 0.1703*y**2.) - 0.41652*y*DL
     ,            + 0.90777*y*DL**2. + 0.12478*y*DL**3. 
     ,            + 0.17155*(1.-y)*DL1 + 0.17191*DL*DL1) - 647.3971*DL
     ,            - 66.41219*DL**2. - 5.353347*DL**3. - 5.925926*DL**4.
     ,            - 0.3950617*DL**5.)
*
      RETURN
      END
*
* ---------------------------------------------------------------------
*
*
*
*
      FUNCTION P3NSA(y, aNS, ipm, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      ipm2=ipm*ipm
      p3nsa = 0.d0

      CALL Pqqns3PAB(y, aNS, FLAVOR, p3nspaa, p3nspba)

      if(ipm2.ne.1.and.ipm.ne.-2) go to 1 ! G.W. 04/07/2007

      if(ipm.eq.1) p3nsa = p3nspaa + Pqqns3PexactA(y, FLAVOR)
      if(ipm.eq.-1) p3nsa = P3NSMA(y,FLAVOR)
      if(ipm.eq.-2) p3nsa = P3NSMA(y,FLAVOR) + P3NSSA (Y, 2.0D0)
      return
*
    1 print 99,ipm
      stop
   99 format(1x,'ipm=',i3)
*
      END
*
* ---------------------------------------------------------------------
*
*
*
*
      FUNCTION P3NSB(y, aNS, ipm, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      ipm2=ipm*ipm
      p3nsb = 0.d0

      CALL Pqqns3PAB(y, aNS, FLAVOR, p3nspaa, p3nspba)

      if(ipm2.ne.1.and.ipm.ne.-2) go to 1 ! G.W. 04/07/2007

      if(ipm.eq.1) p3nsb = p3nspba + Pqqns3PexactB(y, FLAVOR)
      if(ipm.eq.-1) p3nsb = P3NSMB(y,FLAVOR)
      if(ipm.eq.-2) p3nsb = P3NSMB(y,FLAVOR) ! p2nsma(y,nf) + p2nssa(y,2) ! G.W. 05/07/2007
      return
*
    1 print 99,ipm
      stop
   99 format(1x,'ipm=',i3)
*
      END
*
* ---------------------------------------------------------------------
*
*
*
*
      FUNCTION P3NSC(y, aNS, ipm, FLAVOR)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      ipm2=ipm*ipm
      p3nsc = 0.d0              ! G.W. 07/05/2007
      AL1 = LOG(1.-y)

      CALL Pqqns3PAB(y, aNS, FLAVOR, p3nspaa, p3nspba)
*
      if(ipm2.ne.1.and.ipm.ne.-2) go to 1 ! G.W. 04/07/2007
*
      if(ipm.eq.1) p3nsc = Pqqns3PexactC(y, FLAVOR) + p3nspba*(1.-y)*AL1
      if(ipm.eq.-1) p3nsc = P3NSMC(y,FLAVOR) ! p2nsmc(y,nf)
      if(ipm.eq.-2) p3nsc = P3NSMC(y,FLAVOR) ! p2nsmc(y,nf) ! G.W. 04/07/2007
      return
*
    1 print 99,ipm
      stop
   99 format(1x,'ipm=',i3)
*
      END
*
* =====================================================================
*     P_{qq}^{PS}
! ---------------------------------------------------------------------
*     
*     Contribution from approximate four-loop splitting function with
*     variational parameter (-0.7 < aPS < 0.0).
*     This is the regular part of P3PS.
*
      FUNCTION Pqqps3A(x, a)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER n_dimension
      DOUBLE PRECISION N(4), C(4),  
     &        N_exact(4), N_a(4)
      DOUBLE PRECISION inv_A(4, 4)
      DATA PI,PI2/3.14159,9.8696/
*
      n_dimension=4
* 
*     Moments for the known exact small-x contribution (Vogt)
      N_exact(1) = -1.28827235D-01
      N_exact(2) = -2.42108771D-03
      N_exact(3) = -6.81270332D-05
      N_exact(4) = 4.79640375e-05
* 
*     Moments for the varied function multiplied by a (ln(1/x)/x)*(1-x)
      N_a(1) = 0.75D0
      N_a(2) = 0.04861111111078488D0
      N_a(3) = 0.012222222222363594D0
      N_a(4) = 0.004783163265306247D0
*
      DO i=1, n_dimension
        N(i) = N_exact(i) - a * N_a(i)
      ENDDO
*
      inv_A = reshape((/-4.33423515D+01, 1.04394869D+03, 
     &                  -3.26481918D+03, 2.49215628D+03, 1.01800794D+02, 
     &                  -2.29059693D+03, 7.08126548D+03, 
     &                  -5.37933607D+03, 1.28895321D+02, 
     &                  -3.28605015D+03, 1.09724714D+04, 
     &                  -8.75612151D+03, -9.56650450D0, 2.57515782D+02, 
     &                  -9.20800163D+02, 7.84949427D+02/), 
     &                 shape=(/4,4/))

* 
*     Matrix multiplication of inv_A and N
*
      DO i=1,n_dimension
        C(i) = 0.D0
        DO j=1,n_dimension
          C(i) = C(i) + inv_A(j, i) * N(j)
        ENDDO
      ENDDO
*
      Pqqps3A = ((1.-X) * (C(1) * (1./X) + C(2) * (LOG(X)**2.) 
     &         + C(3) * (X**2.) + C(4) * (LOG(1.-X)**2.) 
     &         + a * LOG(1./X)/X) + Pqqps3exactA(X)) * (4*PI)**4.
*
      RETURN
      END
*
* ------------------------------------------------------------------------
*
*     This is the regular exact part of P3PS.
*
      FUNCTION Pqqps3exactA(x)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA PI,PI2/3.14159,9.8696/
*
      CA = 3.0D0
      CF = 4.0D0/3.0D0
      AL1 = LOG(1.-X)
*
      Pqqps3exactA = (((CA ** 2.0)*CF / (6.0D0 * (PI**4.))) * (
     &                (82.0D0/81.0D0) + 2.0 * 1.20206D0) 
     &                * (log(1.0 / x)**2.) / x)
*
      RETURN
      END
*
* =====================================================================
*     P_{qg}
! ---------------------------------------------------------------------
*
*     Contribution from approximate four-loop splitting function with
*     variational parameter (-2.5 < aQG < -0.8).
*     This is the regular part of P3QG.
*
      FUNCTION Pqg3A(x, a)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER n_dimension
      DOUBLE PRECISION N(4), C(4), N_exact(4), 
     &        N_a(4)
      DOUBLE PRECISION inv_A(4, 4)
      DATA PI,PI2/3.14159,9.8696/
*
      n_dimension=4
* 
*     Moments for the known exact small-x contribution (Vogt)
      N_exact(1) = -0.32821791D0
      N_exact(2) = -0.0079068D0
      N_exact(3) = 0.00668085D0
      N_exact(4) = 0.0112976D0
* 
*     Moments for the varied function multiplied by a (ln(1/x)/x)
      N_a(1) = 0.9999999999999999D0
      N_a(2) = 0.11111111111076043D0
      N_a(3) = 0.04000000000013782D0
      N_a(4) = 0.02040816D0
*
      DO i=1, n_dimension
        N(i) = N_exact(i) - a * N_a(i)
      ENDDO
*
      inv_A = reshape((/1.50859345D1, -1.73502577D2, 3.81181701D2, -2.25021001D2,
     &                  1.39350535D1, -2.33166358D2, 5.53133581D2, -3.37644014D2,
     &                  3.92339377D0, -8.52278788D1, 2.57493287D2, -1.79628628D2,
     &                  -1.15513419D-2, 2.74882781D-1, -9.19800075D-1, 7.20490480D-1/),
     &                  shape=(/4,4/))
*
*     Matrix multiplication of inv_A and N
*
      DO i=1,n_dimension
        C(i) = 0.D0
        DO j=1,n_dimension
          C(i) = C(i) + inv_A(j, i) * N(j)
        ENDDO
      ENDDO
*
      Pqg3A = (C(1) * LOG(X)**2. + C(2) * LOG(X) + C(3) * (X**2.) 
     &          + C(4) * LOG(1.-X)**4.
     &          + a * LOG(1./X)/X + Pqg3exactA(X)) * (4*PI)**4.
*
      RETURN
      END
*
* ---------------------------------------------------------------------------
*
*     This is the regular exact part of P3QG.
*
      FUNCTION Pqg3exactA(x)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      DATA PI,PI2/3.14159,9.8696/
*
      CA = 3.0D0
      CF = 4.0D0/3.0D0
      AL1 = LOG(1.-X)
*
      Pqg3exactA = ((CA ** 3.0)/ (6.0 * (PI**4.))) * 
     &               ((82.0D0 / 81.0D0) + 2.0 * 1.20206D0) * 
     &               (log(1.0 / x)**2.) / x
*
      RETURN
      END
*
* =====================================================================
*     P_{gq}
! ---------------------------------------------------------------------
*
*     Contribution from approximate four-loop splitting function with
*     variational parameter (-1.8 < aGQ < -1.5).
*     This is the regular part of P3GQ.
*
      FUNCTION Pgq3A(x, a)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER n_dimension
      DOUBLE PRECISION N(4), C(4), N_exact(4), 
     &        N_a(4)
      DOUBLE PRECISION inv_A(4, 4)
      DATA PI,PI2/3.14159,9.8696/
*
      n_dimension=4
*
* 
*     Moments for the known exact small-x contribution (Vogt)
      N_exact(1) = -0.78602265D0
      N_exact(2) = 0.05776039D0
      N_exact(3) = 0.04868481D0
      N_exact(4) = 0.04071129D0
* 
*     Moments for the varied function multiplied by a (ln^2(1/x)/x)
      N_a(1) = 1.9999999999999991D0
      N_a(2) = 0.07407407407348075D0
      N_a(3) = 0.016000000001040217D0
      N_a(4) = 0.0058309037900882415D0
*
      DO i=1, n_dimension
        N(i) = N_exact(i) - a * N_a(i)
      ENDDO
*
      inv_A = reshape((/-15.48054909D0, -39.22467343D0, 13.63767859D0, 
     &                 3.69946174D0, 442.68554891D0, 1038.70260338D0, 
     &                 -415.55730182D0,-113.79603728D0,-1263.91174116D0, 
     &                 -2938.1594667D0, 1302.32965951D0,362.67726047D0, 
     &                 874.63704072D0, 2025.64293423D0,-958.33374045D0, 
     &                 -272.56485302D0/), shape=(/4,4/))
*
*     Matrix multiplication of inv_A and N
      DO i=1,n_dimension
        C(i) = 0.D0
        DO j=1,n_dimension
          C(i) = C(i) + inv_A(i, j) * N(j)
        ENDDO
      ENDDO
*
      Pgq3A = (C(1) * LOG(1./X)/X + C(2) * (LOG(X)**3.) + C(3) * X 
     &          + C(4) * LOG(1.-X) + a * (LOG(1./X)**2.)/X 
     &          + Pgq3exactA(X)) * (4*PI)**4.
*
      RETURN
      END
*
* ---------------------------------------------------------------------------
*
*     This is the regular exact part of P3GQ.
*
      FUNCTION Pgq3exactA(x)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA PI,PI2/3.14159,9.8696/
*
      CA = 3.0D0
      CF = 4.0D0/3.0D0
      AL1 = LOG(1.-X)
*
      Pgq3exactA = ((CA ** 3.0) * CF / (3.0 * (PI**4.))) * 1.20206D0 
     &               * (log(1.0 / x)**3.) / x
*
      RETURN
      END
*
* =====================================================================
*     P_{gg}
! ---------------------------------------------------------------------
*
*     Contribution from approximate four-loop splitting function with
*     variational parameter (-5 < aGG < -15).
*     This is the regular part of P3GG.
*
      FUNCTION Pgg3A(x, a)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER n_dimension
      DOUBLE PRECISION N(4), C(4), N_exact(4), 
     &        N_a(4), inv_A(4, 4)
      DATA PI,PI2/3.14159,9.8696/
*
      n_dimension=4
* 
*     Moments for the known exact small-x contribution (Vogt)
      N_exact(1) = 6.94542399D0
      N_exact(2) = 0.01255781D0
      N_exact(3) = -0.24160401D0
      N_exact(4) = -0.27416992D0
* 
*     Moments for the varied function multiplied by a (ln(1/x)/x)
      N_a(1) = 0.9999999999999999D0
      N_a(2) = 0.11111111111076043D0
      N_a(3) = 0.04000000000013782D0
      N_a(4) = 0.02040816D0
*
      DO i=1, n_dimension
        N(i) = N_exact(i) - a * N_a(i)
      ENDDO
*
      inv_A = reshape((/1.57592805D1, -1.89525930D2, 4.34798303D2, -2.67019531D2,
     &                  1.51046749D1, -2.60999383D2, 6.46267166D2, -4.10596673D2,
     &                  6.05043429D0, -1.35844233D2, 4.26863395D2, -3.12298274D2,
     &                  -3.86434270D-1, 9.19582568D0, -3.07706474D1, 2.41030189D1/),
     &                  shape=(/4,4/))
*
*     Matrix multiplication of inv_A and N     
      DO i=1,n_dimension
        C(i) = 0.D0
        DO j=1,n_dimension
          C(i) = C(i) + inv_A(j, i) * N(j)
        ENDDO
      ENDDO
*
      Pgg3A = (C(1) * LOG(X)**2. + C(2) * LOG(X) + C(3) * (X**2.) 
     &          + C(4) * (LOG(1.-X)**2.) 
     &          + a * LOG(1./X)/X + Pgg3exactA(X)) * (4*PI)**4.
*
      RETURN
      END
*
* -----------------------------------------------------------------------------
*
*     This is the regular exact part of P3GG.
*
      FUNCTION Pgg3exactA(x)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA PI,PI2/3.14159,9.8696/
      DOUBLE PRECISION NF4
*
      CA = 3.0D0
      CF = 4.0D0/3.0D0
      AL1 = LOG(1.-X)
      NF4 = 4.0D0 ! Set Nf = 4 since the moments in the approximation 
                  ! are for 4 flavours.
*
      Pgg3exactA = ((CA ** 4.0D0) / (3.0D0 * (PI**4.0D0))) * 1.20206D0 
     &               * (log(1.0D0/x)**3.0D0)/ x + (1.0D0 / ((PI**4.0D0))) 
     &               * ((CA**4.0D0)*(-(1205.0D0/162.0D0)+(67.0D0/36.0D0) 
     &               * ((PI**2.0)/6.0D0)+(1.0D0/4.0D0) * (((PI**2.0) 
     &               / 6.0D0)**2.0)-(11.0D0/2.0D0)*1.20206D0) + NF4 
     &               * (CA ** 3.0) * (-(233.0D0/162.0D0)+(13.0D0/36.0D0) 
     &               * ((PI**2.0)/6.0D0)-(1.0D0/3.0D0) * 1.20206D0) + NF4 
     &               * (CA ** 2.0) * CF * ((617.0D0/243.0D0) 
     &               - (13.0D0/18.0D0)*((PI**2.0)/6.0D0) + (2.0D0/3.0D0) 
     &               * 1.20206D0))*(1.0D0/2.0D0)*((log(1.0/x)**2.0)/x)
*
      RETURN
      END
*