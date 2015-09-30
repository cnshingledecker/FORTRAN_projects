MODULE pjg76
  IMPLICIT NONE

  ! Data dictionary: Parameters for calculations
  REAL(KIND=8), PARAMETER, DIMENSION(7) :: FJVALS=(/ 0.08,0.19,0.19,0.17, &
                                                     0.11,0.16,0.1 /)
  REAL(KIND=8), PARAMETER, DIMENSION(7) :: IJVALS=(/ 12.1,16.1,16.9,18.2, &
                                                     20.3,23.0,37.0 /)
  REAL(KIND=8), PARAMETER, DIMENSION(7) :: CJVALS=(/ 1.93,2.56,2.69,2.90, &
                                                     3.23,3.66,5.89 /)
  REAL(KIND=8), PAREMETER               :: B0=0.030
  REAL(KIND=8), PARAMETER               :: B1=1.035
  REAL(KIND=8), PARAMETER               :: E0=8239 !eV
  REAL(KIND=8), PARAMETER               :: T1=68.3
  REAL(KIND=8), PARAMETER               :: G1=189.1
  REAL(KIND=8), PARAMETER               :: K = 6.55E-16
  REAL(KIND=8), PARAMETER               :: GS=13.1
  REAL(KIND=8), PARAMETER               :: TS=6.34
  REAL(KIND=8), PARAMETER               :: GAM1=5.00E5 !eV
  REAL(KIND=8), PARAMETER               :: GAM2=7.60E4 !eV
  REAL(KIND=8), PARAMETER               :: TA=2.52E3   !eV
  REAL(KIND=8), PARAMETER               :: TB=1.28E2 !eV
  REAL(KIND=8), PARAMETER               :: J=40.3 !eV
  REAL(KIND=8), PARAMETER               :: NU=3.14E-1
  REAL(KIND=8), PARAMETER               :: DLTA=132.1 !eV
  REAL(KIND=8), PARAMETER               :: CVAC=2.998E10 !cm/s
  REAL(KIND=8), PARAMETER               :: MP=1.672622E-24 !g
  REAL(KIND=8), PARAMETER               :: EV2ERG=1.602177E-12 !erg/eV

CONTAINS

FUNCTION fj(energy)
  IMPLICIT NONE

  ! Data dictionary: Calling parameters
  REAL(KIND=8), INTENT(IN)               :: energy
  REAL(KIND=8)            , DIMENSION(7) :: fj

  ! Data dictionary: Local variables
  INTEGER                                :: n
  REAL(KIND=8)                           :: e_erg
  REAL(KIND=8)                           :: num,den,den1,den2
  REAL(KIND=8)                           :: bfac

  bfac = beta(energy)
  
  DO n=1,7 
    num = ((0.5*MP*(bfac*bfac)*(CVAC*CVAC)**(NU+1))*FJVALS(n)
    den1 = (j*EV2ERG)**(NU+1)
    den2 = (0.5*MP*(bfac*bfac)*(CVAC*CVAC)**(NU+1)
    den = den1+den2
    fj(n) = num/den
  END DO
  RETURN
END FUNCTION fj

FUNCTION t0(energy)
  IMPLICIT NONE

  ! Data dictionary: Calling parameters
  REAL(KIND=8), INTENT(IN)               :: energy
  REAL(KIND=8)                           :: t0 

  ! Data dictionary: Local variables     
  REAL(KIND=8)                           :: bfac
  REAL(KIND=8)                           :: den1,den,num

  bfac = beta(energy)

  den1 = 0.5*MP*(bfac*bfac)*(CVAC*CVAC)
  den = den1 + (TB*EV2ERG)
  num = TA*EV2ERG

  t0 = TS - (num/den)
  RETURN
END FUNCTION t0

FUNCTION beta(energy_ev)
  IMPLICIT NONE

  ! Data dictionary: Calling parameters
  REAL(KIND=8), INTENT(IN)               :: energy_ev
  REAL(KIND=8)                           :: beta

  ! Data dictionary: Local variables
  REAL(KIND=8)                           :: energy_erg

  energy_erg = energy_ev*EV2ERG
  beta = SQRT((2.*energy_erg)/MP)
  RETURN
END FUNCTION beta

FUNCTION gamfac(energy)
  IMPLICIT NONE

  ! Data dictionary: Calling parameters
  REAL(KIND=8), INTENT(IN)               :: energy
  REAL(KIND=8)                           :: gamfac

  ! Data dictionary: Local variables  
  REAL(KIND=8)                           :: bfac
  REAL(KIND=8)                           :: num,den1,den2,den

  num = GAM1*EV2ERG
  den1 = 0.5*MP*(bfac*bfac)*(CVAC*CVAC)
  den2 = GAM2
  den = den1 + den2

  gamfac = GS + (num/den)
  RETURN
END FUNCTION gamfac

FUNCTION spr(pr_en,se_en,fj,gamfac)
  IMPLICIT NONE

  ! Data dictionary: Calling parameters
  REAL(KIND=8), INTENT(IN) :: pr_en !proton energy in eV
  REAL(KIND=8), INTENT(IN) :: se_en !secondary electron energy in eV
  REAL(KIND=8), INTENT(IN), DIMENSION(7) :: fj
  REAL(KIND=8), INTENT(IN)               :: gamfac
  REAL(KIND=8)             :: spr

  ! Data dictionary: Local variables
  REAL(KIND=8)             :: prefacnum,prefacden
  REAL(KIND=8)             :: prefac
  REAL(KIND=8)             :: fac1prefac,fac1insidesnum,fac1insidesden
  REAL(KIND=8)             :: fac1insides,part1_1
  REAL(KIND=8)             :: bfac
  REAL(KIND=8)             :: prenerg,seenerg
  REAL(KIND=8)             :: fac2p1,fac2p2
  INTEGER                  :: n

  bfac = beta(pr_en)

  DO n=1,7
    prefacnum = fj(n)
    prefacden = 0.5*MP*(bfac*bfac)*(CVAC*CVAC)
    prefac = prefacnum/prefacden

    fac1prefac = K*gamfac*gamfac

    fac1insidesnum = 4*MP*(bfac*bfac)*(CVAC*CVAC)*cj(n)
    fac1insidesden = 2*IJ(n)*(1-(bfac*bfac))
    fac1insides = fac1insidesnum/fac1insidesden
  END DO

END FUNCTION spr

END MODULE pjg76
