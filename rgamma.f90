MODULE mc_toolbox

  CONTAINS
    FUNCTION rgamma(a)
      ! This function returns a random number with a Gamma PDF using
      ! the method described by Marsaglia and Tsang 2000.
      IMPLICIT NONE

      ! Data dictionary: Calling parameters
      REAL(8), INTENT(IN) :: a
      REAL(8)             :: rgamma

      ! Data dictionary: Local variables
      REAL(8)             :: d,c,v
      REAL(8)             :: x,u

      d = a-1./3.
      c = 1./SQRT(9.*d)
      x = R8_NORMAL_01() 
      v = 1. + c*x
      DO WHILE ( v .LE. 0. )


