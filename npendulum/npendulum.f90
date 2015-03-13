      program npendulum
        use types, only: rp
        use modules, only: solver, save_set, initalize, zeros
        implicit none

        logical :: inital_set
        integer :: n, sets, set_size, integrator, format, file_id, i
        real(rp) :: l, g, dt
        real(rp), allocatable :: theta(:), theta_d(:), theta_dd(:), data(:, :)

        call initalize("initalization.dat", n, sets, set_size, integrator, format, l, g, dt, theta, theta_d, theta_dd)
!        call zeros("initalization.dat", 1)
        inital_set = .true.
        
        do i = 1, sets
          call solver(n, set_size, integrator, format, l, g, dt, theta, theta_d, theta_dd, data)
          call save_set("output.dat", inital_set, n, sets, set_size, integrator, format, l, g, dt, data)

          if (inital_set) then
            inital_set = .false.
          end if
        end do
      end program npendulum
