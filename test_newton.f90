program test_newton

    use newton, only: solve
    use functions

    implicit none
    real(kind=8) :: x, x0, fx
    real(kind=8) :: x0vals(7)
    integer      :: iters, itest
	  logical      :: debug         ! set to .true. or .false.


    print *, "Test routine for computing zero of f"
    debug = .false.

    ! values to test as x0:
    x0vals = (/5.d-1, 1.d-1, 1.d0, 2.d0, 1.d1, 5.d1, 100.d0 /)

    do itest=1,SIZE(x0vals)
        x0 = x0vals(itest)
        print *, ' '  ! blank line
        call solve(bh80,bh80_prime, x0, x, iters, debug)

        print 11, x, iters
11      format('solver returns x = ', e22.15, ' after', i3, ' iterations')

        fx = bh80(x) 
        print 12, fx
12      format('the value of f(x) is ', e22.15)

        if (abs(x-2.d0) > 1d-14) then
            print 13, x
13          format('*** Unexpected result: x = ', e22.15)
            endif
        enddo

end program test_newton
