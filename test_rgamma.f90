PROGRAM test_rgamma
  USE mc_toolbox
  IMPLICIT NONE

  ! Data dictionary: Declare local variables
  INTEGER :: iters,n
  REAL(KIND=8) :: a
  REAL(KIND=8) :: gam

  iters=1E5
  CALL SRAND(GSEED)

  PRINT *, '*******/Starting the Program/*******'
  a = 13.0 

  OPEN(UNIT=1001,FILE="granvals.csv",POSITION='APPEND', STATUS='REPLACE')

  DO n=1,iters
    gam = rgamma(a)
!    PRINT *, 'The rgamma function gives a value of:',gam
    WRITE(1001,*) gam
  END DO

  CLOSE(1001)


END PROGRAM test_rgamma
