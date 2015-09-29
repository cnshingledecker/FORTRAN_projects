PROGRAM test_bh80
  !The goal of this program is to test the subroutines/functions in the bh80.f90
  !file. This includes the famous MAGIC subroutine described in Biersack &
  !Haggmark 1980. 
  USE bh80
  IMPLICIT NONE

  !Data dictionary
  REAL               :: rn
  DOUBLE PRECISION               :: tth,els
  REAL               :: a,p,br
  DOUBLE PRECISION               :: c2,s2
  DOUBLE PRECISION               :: cmTheta
  REAL               :: labTheta
  REAL               :: mf
  DOUBLE PRECISION               :: energy
  INTEGER            :: iters,n
  INTEGER            :: ssize
  INTEGER, DIMENSION(:),ALLOCATABLE :: seed
  DOUBLE PRECISION :: sn,snels,sigmaEl

  CALL RANDOM_SEED(SIZE=ssize)
  ALLOCATE(seed(ssize))
  seed = 718381 
  CALL RANDOM_SEED(PUT=seed)
  iters = 100
  energy = ENERG

  DO WHILE ( energy .GT. 1.6E-17 )  
  CALL RANDOM_NUMBER(rn)
  a = AU(Z1,Z2)
  br = b(rn,a,RHO)
  els = eps(energy,Z1,Z2,M1,M2,a)
  CALL magic(els,br,c2,s2,cmTheta)
  mf = mass_fac(M1,M2)
  tth = t(energy,mf,s2)
  labTheta = lab_theta(cmTheta,M1,M2)
  energy = energy - tth
!  PRINT *, 'T=',tth,'eV and E=',energy,'eV'
!  PRINT *, "T=",tth,'eV and theta=',labTheta*57.2957795,"deg"
!  PRINT *, 'e=',els,' and b=',br,' and s2=',s2
  
  snels = sneps(els)
  sn    = sne(snels,z1,z2,m1,m2)
  sigmaEl = selast(energy,sn,mf)
  PRINT *, 'Sn(E)=',sn,'sn(eps)=',snels,'sigma_el=',sigmaEl
  END DO

END PROGRAM test_bh80




