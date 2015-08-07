PROGRAM rand_test
IMPLICIT NONE

REAL :: u,rand
INTEGER :: i

CALL RANDOM_SEED()

  DO WHILE ( u .EQ. 0.0 .AND. rand .EQ. 0.0 ) 
    CALL RANDOM_NUMBER(u)
    CALL RANDOM_NUMBER(rand)
  END DO 
  PRINT *, "u is",u,"and rand is",rand


END PROGRAM rand_test
