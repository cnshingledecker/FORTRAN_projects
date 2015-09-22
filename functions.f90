module functions

! ZBL screening parameters
REAL, PARAMETER :: A1=0.181
REAL, PARAMETER :: A2=0.5099
REAL, PARAMETER :: A3=0.2802
REAL, PARAMETER :: A4=0.02817
REAL, PARAMETER :: E1=-3.2
REAL, PARAMETER :: E2=-0.9423
REAL, PARAMETER :: E3=-0.4029
REAL, PARAMETER :: E4=-0.2016

! Bohr radius in cm
REAL, PARAMETER :: A0=0.5

! Atomic Parameters
REAL(KIND=8), PARAMETER :: z1 = 1
REAL(KIND=8), PARAMETER :: z2 = 8

contains

real(kind=8) function f_sqrt(x)
    implicit none
    real(kind=8), intent(in) :: x

    f_sqrt = x**2 - 4.d0

end function f_sqrt


real(kind=8) function fprime_sqrt(x)
    implicit none
    real(kind=8), intent(in) :: x
    
    fprime_sqrt = 2.d0 * x

end function fprime_sqrt

FUNCTION zbl(z1,z2,r)
  IMPLICIT NONE

  !Function parameters and result
  REAL(KIND=8), INTENT(IN) :: z1,z2
  REAL(KIND=8), INTENT(IN) :: r
  REAL(KIND=8)             :: zbl
  
  !Internal variables
  REAL(KIND=8)             :: prefac
  REAL(KIND=8)             :: sfac1
  REAL(KIND=8)             :: sfac2
  REAL(KIND=8)             :: sfac3
  REAL(KIND=8)             :: sfac4
  REAL(KIND=8)             :: phi 
  REAL(KIND=8)             :: a_num,a_den,a
  REAL(KIND=8)             :: x

  !Calculate Screening length
  a_num = (0.8853*A0)
  a_den = (z1**(0.23)) + (z2**(0.23))
  a     = a_num/a_den

  !Define Reduced Length
  x = r/a 

  !Calculate Prefactor
  prefac = (z1*z2)/r

  !Calculate Universal Screening F'n
  phi = A1*EXP(E1*x) + A2*EXP(E2*x) + A3*EXP(E3*x) + A4*EXP(E4*x)

  !Calculate Interatomic Potential
  zbl = prefac*phi
  RETURN
END FUNCTION zbl

FUNCTION zbl_prime(z1,z2,r)
  IMPLICIT NONE

  !Function Parameters
  REAL(KIND=8), INTENT(IN) :: z1,z2
  REAL(KIND=8), INTENT(IN) :: r
  REAL(KIND=8)             :: zbl_prime

  !Local Variables
  REAL(KIND=8)             :: prefac1
  REAL(KIND=8)             :: prefac2
  REAL(KIND=8)             :: qfac1,qfac2,qfac3,qfac4
  REAL(KIND=8)             :: qfac
  REAL(KIND=8)             :: phi
  REAL(KIND=8)             :: phi1,phi2,phi3,phi4
  REAL(KIND=8)             :: x
  REAL(KIND=8)             :: a_num,a_den,a
  REAL(KIND=8)             :: part1,part2

  !Calculate Screening Length
  a_num   = (0.8853*A0)
  a_den   = (z1**(0.23)) + (z2**(0.23))
  a       = a_num/a_den

  !Define Reduced Length
  x       = r/a 

  !Calculate First Prefactor
  prefac1 = (z1*z2)/r

  !Calculate Exponential Component of 1
  qfac1   = (E1*A1*EXP(E1*x))/a
  qfac2   = (E2*A2*EXP(E2*x))/a
  qfac3   = (E3*A3*EXP(E3*x))/a
  qfac4   = (E4*A4*EXP(E4*x))/a
  qfac    = qfac1 + qfac2 + qfac3 + qfac4

  !Calculate Component 1 of Derivative
  part1   = prefac1*qfac

  !Calculate Prefactor of Component 2
  prefac2 = (z1*z2)/(r**2)

  !Calculate Exponential Part of Component 2
  phi1    = A1*EXP(E1*x)
  phi2    = A2*EXP(E2*x)
  phi3    = A3*EXP(E3*x)
  phi4    = A4*EXP(E4*x)
  phi     = phi1 + phi2 + phi3 + phi4

  !Calculate Component 2 of Derivative
  part2   = prefac2*phi

  !Calculate Derivative of ZBL Potential
  zbl_prime = part1 - part2
  RETURN
END FUNCTION zbl_prime

FUNCTION bh80(r)
  IMPLICIT NONE

  !Function Parameters
  REAL(KIND=8), INTENT(IN) :: r
  REAL(KIND=8)             :: bh80

  !Local Variables
  REAL(KIND=8)             :: ec
  REAL(KIND=8)             :: p
  REAL(KIND=8)             :: vr

  p = 1.5
  ec = 100000
  vr = zbl(z1,z2,r)

  bh80 = 1 - (vr/ec) - (p/r)**2
  RETURN
END FUNCTION bh80

FUNCTION bh80_prime(r)
  IMPLICIT NONE

  !Function Parameters
  REAL(KIND=8), INTENT(IN) :: r
  REAL(KIND=8)             :: bh80_prime

  !Local Variables
  REAL(KIND=8)             :: ec
  REAL(KIND=8)             :: p
  REAL(KIND=8)             :: vrp

  p = 1.5
  ec = 100000
  vrp = zbl_prime(z1,z2,r)

  bh80_prime = (2*(p**2))/(r**3) - (vrp/ec)
  RETURN
END FUNCTION bh80_prime

  

end module functions

