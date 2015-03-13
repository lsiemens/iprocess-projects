      program npendulum
        use types, only: rp
        use modules, only: theta_dd_solver, initalize, zeros
        implicit none

        integer :: n, sets, set_size, integrator, format, file_id, i, j
        real(rp) :: l, g, dt
        real(rp), allocatable :: theta(:), theta_d(:), theta_dd(:)
        call initalize("initalization.dat", n, sets, set_size, integrator, format, l, g, dt, theta, theta_d, theta_dd)
!        call zeros("initalization.dat", 10)
        
        call theta_dd_solver(n, l, g, theta, theta_d, theta_dd)
        print *, theta_dd
      end program npendulum
