      program npendulum
        use types, only: rp
        use modules, only: theta_dd_solver, initalize, zeros
        implicit none

        integer :: n, file_id, i, j
        real(rp) :: l, g
        real(rp), allocatable :: theta(:), theta_d(:), theta_dd(:)
        call initalize("initalization.dat", n, l, g, theta, theta_d, theta_dd)
!        call zeros("initalization.dat", 2)
        
        call theta_dd_solver(n, l, g, theta, theta_d, theta_dd)
!        print *, theta_dd(:)

!        print *, n, l, g
        print *, theta_dd
      end program npendulum
