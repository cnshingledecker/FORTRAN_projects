PROGRAM typeparam
  USE type
  IMPLICIT NONE

  TYPE(helloworld) :: greet
  TYPE(composite) :: comp
  
!  PRINT *, t(1)
!  PRINT *, t(2)
!  PRINT *, SIZE(t)

  greet%hwtext = "composite!"
  greet%n = 42

  comp%phrase = greet
  comp%val = 1947

  PRINT *, comp%phrase%hwtext
END PROGRAM
