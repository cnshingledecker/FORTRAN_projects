MODULE bh80
  IMPLICIT NONE

  REAL, PARAMETER :: Z1=1.
  REAL, PARAMETER :: Z2=8.
  REAL, PARAMETER :: RHO =0.0286 !in Angstrom^-3
  DOUBLE PRECISION, PARAMETER :: ENERG =100*1E3 !1.602E-14 !in eV 
  REAL, PARAMETER :: PI=3.14159
  REAL, PARAMETER :: M1=1 !Ion mass in amu 
  REAL, PARAMETER :: M2=16 !Target mass in amu
  REAL, PARAMETER :: A0=0.529177 !Bohr radius in Angstroms 
  REAL, PARAMETER :: ECHARG2 = 14.39 !Square of the electron charge in eV*Angstroms 

  CONTAINS

    SUBROUTINE magic(eps,b,c2,s2,theta)
      IMPLICIT NONE

      !Data dictionary: Calling Parameters and Output
      DOUBLE PRECISION, INTENT(IN)  :: eps
      REAL, INTENT(IN)  :: b
      DOUBLE PRECISION, INTENT(OUT) :: c2,s2
      DOUBLE PRECISION, INTENT(OUT) :: theta

      !Data dictionary: Local Variables
      DOUBLE PRECISION              :: r, rr
      DOUBLE PRECISION              :: ex1,ex2,ex3,ex4
      DOUBLE PRECISION              :: v,v1
      DOUBLE PRECISION              :: fr,fr1
      DOUBLE PRECISION              :: q
      DOUBLE PRECISION              :: roc
      DOUBLE PRECISION              :: sqe
      DOUBLE PRECISION              :: cc,aa,ff
      DOUBLE PRECISION              :: delta
      DOUBLE PRECISION              :: co

      r     = b
      rr    = -2.7*DLOG(eps*b)
      IF ( rr .LT. b ) GOTO 1980
      rr    = -2.7*DLOG(eps*b)
      IF ( rr .LT. b ) GOTO 1980
      r     = rr
       PRINT *, 'The value of r is',r
1980  ex1   =  0.18175*EXP(-3.1998*r)
      ex2   =  0.50986*EXP(-0.94229*r)
      ex3   =  0.28022*EXP(-0.4029*r)
      ex4   = 0.028171*EXP(-0.20162*r)
      v     = (ex1 + ex2 + ex3 + ex4)/r
      v1    = -(v+3.1998*ex1+0.94229*ex2+0.4029*ex3+0.20162*ex4)/r
      fr    = b*b/r+v*r/eps-r
      fr1   = -b*b/(r*r)+(v+v1*r)/eps-1.0
      q     = fr/fr1
      r     = r-q
      IF ( ABS(q/r) .GT. 0.001) GOTO 1980
      roc   = -2.0*(eps-v)/v1
!      PRINT *, 'roc=',roc
      sqe   = SQRT(eps)
!      PRINT *, 'sqe=',sqe
      cc    = (0.011615+sqe)/(0.0071222+sqe)
      aa    = 2.0*eps*(1.0+(0.99229/sqe))*b**cc
      ff    = (SQRT(aa**2+1.0)-aa)*((9.3066+eps)/(14.813+eps))
      delta = (r-b)*aa*ff/(ff+1.0)
!      PRINT *, 'cc=',cc,'aa=',aa,'ff=',ff,'delta=',delta
      co    = (b+delta+roc)/(r+roc)
      c2    = co*co
      s2    = 1.0-c2
      PRINT *, 'co=',co,'c2=',c2,'s2=',s2
      theta = 2.0*ACOS(co)
      RETURN
    END SUBROUTINE magic

    FUNCTION b(rn,a,rho)
      IMPLICIT NONE

      !Data dictionary: Calling parameters
      REAL, INTENT(IN) :: rn 
      REAL, INTENT(IN) :: rho
      REAL, INTENT(IN) :: a
      REAL             :: b

      !Data dictionary: Local variables
      REAL             :: p

      p = SQRT(rn/(3.14159*(rho**(2./3.))))
      b = p/a
!      PRINT *, 'The value of b is',b
      RETURN
    END FUNCTION b

    FUNCTION au(z1,z2)
      IMPLICIT NONE

      !Data dictionary: Calling parameters
      REAL, INTENT(IN) :: z1,z2
      REAL             :: au

      au = (0.8853*A0)/(z1**0.23 + z2**0.23)
      PRINT *, 'The screening length is',au,'Angstroms'
      RETURN
    END FUNCTION au

    FUNCTION eps(en,z1,z2,m1,m2,au)
      !Note: All units must be Gaussian-CGS
      IMPLICIT NONE

      !Data dictionary: Calling parameters
      DOUBLE PRECISION, INTENT(IN) :: en 
      REAL, INTENT(IN) :: z1,z2
      REAL, INTENT(IN) :: m1,m2
      REAL, INTENT(IN) :: au
      DOUBLE PRECISION             :: eps

      !Data dictionary: Local variables
      DOUBLE PRECISION             :: fac1,fac2,fac3

      fac1 = au/ECHARG2
      fac2 = m2/(m1+m2)
      fac3 = 1./(z1+z2)
      eps = en*fac1*fac2*fac3 
      PRINT *, 'The center of mass energy is',en*fac2
      PRINT *, 'The reduced energy is',eps
    END FUNCTION eps

    FUNCTION mass_fac(m1,m2)
      IMPLICIT NONE

      !Data dictionary: Calling parameters
      REAL, INTENT(IN) :: m1,m2
      REAL             :: mass_fac

!      PRINT *, 'Mass 1 is',M1,'and mass 2 is',M2
      mass_fac = (4.*m1*m2)/((m1+m2)**2)
!      PRINT *, 'The mass factor is',mass_fac
      RETURN
    END FUNCTION mass_fac

    FUNCTION t(e,mf,s2)
      IMPLICIT NONE

      !Data dictionary: Calling parameters
      DOUBLE PRECISION, INTENT(IN) :: e
      REAL, INTENT(IN) :: mf
      DOUBLE PRECISION, INTENT(IN) :: s2
      DOUBLE PRECISION             :: t

      t = mf*e*s2
      RETURN
    END FUNCTION t

    FUNCTION lab_theta(cmtheta,m1,m2)
      IMPLICIT NONE

      !Data dictionary: Calling parameters
      DOUBLE PRECISION, INTENT(IN) :: cmtheta
      REAL, INTENT(IN) :: m1,m2
      REAL             :: lab_theta

      !Data dictionary: Local variables
      REAL             :: insides

      insides = SIN(cmtheta)/(COS(cmtheta)+(m1/m2))
      lab_theta = ATAN(insides)
      RETURN
    END FUNCTION lab_theta
    
  FUNCTION sneps(eps)
    IMPLICIT NONE

    !Data dictionart: Calling Parameters
    DOUBLE PRECISION, INTENT(IN) :: eps
    DOUBLE PRECISION             :: sneps

    !Data dictionary: Local variables
    DOUBLE PRECISION             :: num
    DOUBLE PRECISION             :: den1,den2,den3,den

    IF ( eps .LE. 30. ) THEN
      num = DLOG(1.+1.1383*eps)
      den1 = eps
      den2 = 0.01321*(eps**0.21226)
      den3 = 0.19593*(eps**0.5)
      den  = 2*(den1+den2+den3)
      sneps = num/den
    ELSE
      sneps = DLOG(eps)/(2*eps)
    END IF
    RETURN
  END FUNCTION sneps

  FUNCTION sne(sneps,z1,z2,m1,m2)
    IMPLICIT NONE

    !Data dicionary: Calling parameters
    DOUBLE PRECISION, INTENT(IN) :: sneps
    REAL            , INTENT(IN) :: z1,z2
    REAL            , INTENT(IN) :: m1,m2
    DOUBLE PRECISION             :: sne

    !Data dictionary: Local variables
    DOUBLE PRECISION             :: num,den
    DOUBLE PRECISION             :: den1,den2

    num = (8.462E-15)*z1*z2*m1*sneps
    den1 = m1 + m2
    den2 = z1**0.23 + z2**0.23
    den = den1*den2

    sne = num/den
    RETURN
  END FUNCTION

  FUNCTION selast(energy,sne,mf)
    IMPLICIT NONE

    !Data dictionary: Calling parameters
    DOUBLE PRECISION, INTENT(IN) :: energy
    DOUBLE PRECISION, INTENT(IN) :: sne
    REAL            , INTENT(IN) :: mf
    DOUBLE PRECISION             :: selast

    selast = (2*sne)/(mf*energy)
    RETURN
  END FUNCTION

END MODULE bh80
