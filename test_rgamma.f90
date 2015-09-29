PROGRAM test_rgamma
  USE mc_toolbox
  IMPLICIT NONE

  ! Data dictionary: Declare local variables
  REAL(KIND=8) :: a
  REAL(KIND=8) :: gam

  CALL SRAND(GSEED)

  PRINT *, '*******/Starting the Program/*******'
  a = 1.394
  gam = rgamma(a)

  PRINT *, 'The rgamma function gives a value of:',gam

END PROGRAM test_rgamma
