PROGRAM true_test
  IMPLICIT NONE

  INTEGER :: n
  LOGICAL :: q

  q = .FALSE.
  DO WHILE (q .EQV. .FALSE.)
    n = n + 1
    IF ( n .EQ. 100 ) q = .TRUE.
  END DO 

  PRINT *, n

END PROGRAM true_test

