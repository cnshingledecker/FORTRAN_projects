MODULE type

  TYPE helloworld
    CHARACTER(len=30) :: hwtext
    INTEGER           :: n
  END TYPE helloworld

  TYPE(helloworld), PARAMETER :: t=helloworld("Hello, world.",1947) 
END MODULE type
