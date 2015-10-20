MODULE type

  TYPE helloworld
    CHARACTER(len=30) :: hwtext
    INTEGER           :: n
  END TYPE helloworld

  TYPE composite
    TYPE(helloworld) :: phrase
    INTEGER :: val
  END TYPE composite

  TYPE(helloworld), PARAMETER, DIMENSION(2) :: t=(/helloworld("Hello, world.",1947),&
                                                   helloworld("Here I am",1985)/)

END MODULE type
