
      module modules
        use types, only: rp, sp, dp, qp
        implicit none
        private
        public :: solver, save_set, initalize
        contains

        subroutine solver(n, set_size, integrator, format, find_energy, l, g, dt, theta, theta_d, theta_dd, data, energy, hasnan)
          integer :: i
          logical, intent(inout) :: hasnan
          logical, intent(in) :: find_energy
          integer, intent(in) :: n, set_size, integrator, format
          real(rp), intent(in) :: l, g, dt
          real(rp), intent(inout), allocatable :: theta(:), theta_d(:)
          real(rp), intent(inout), allocatable :: theta_dd(:)
          real(rp), intent(out), allocatable :: data(:, :), energy(:)
          
          if ((format.lt.1).or.(format.gt.3)) then
            print *, "format must be 1, 2, or 3"
            stop "error"
          end if

          allocate(data(n, set_size * format))
          if (find_energy) then
            allocate(energy(set_size))
          end if
          
          do i=1, set_size
            if (.not.(hasnan)) then
              if (any(.not.(theta.eq.theta))) then
                hasnan = .true.
                print *, "Error: Nan detected in theta"
              end if
            end if

            if (format.eq.1) then
              data(:,i) = theta
            else if (format.eq.2) then
              data(:,format*i-1) = theta
              data(:,format*i) = theta_d
            else
              data(:,format*i-2) = theta
              data(:,format*i-1) = theta_d
              data(:,format*i) = theta_dd
            end if
            
            if (find_energy) then
              energy(i) = calculate_energy(n, l, g, theta, theta_d)
            end if

            if (.not.(hasnan)) then
              if (integrator.eq.1) then
                call velocity_verlet_integrator(n, l, g, dt, theta, theta_d, theta_dd)
              else
                print *, "integrator must be 1"
                stop "error"
              end if
            end if
          end do
        end subroutine solver

        subroutine velocity_verlet_integrator(n, l, g, dt, theta, theta_d, theta_dd)
          integer, intent(in) :: n
          real(rp), intent(in) :: l, g, dt
          real(rp), intent(inout), allocatable :: theta(:), theta_d(:)
          real(rp), intent(inout), allocatable :: theta_dd(:)
          real(rp), allocatable :: old_theta_dd(:)
          
          allocate(old_theta_dd(n))
          old_theta_dd = theta_dd(:)
          
          theta = theta + theta_d*dt+0.5_rp*theta_dd*dt*dt
          call theta_dd_solver(n, l, g, theta, theta_d, theta_dd)
          theta_d = theta_d + 0.5_rp*(old_theta_dd+theta_dd)*dt
        end subroutine velocity_verlet_integrator

        function calculate_energy(n, l, g, theta, theta_d) result(energy)
          integer :: i, j
          integer, intent(in) :: n
          real(rp), intent(in) :: l, g
          real(rp), intent(in), allocatable :: theta(:), theta_d(:)
          real(rp) :: kenetic, potential, energy
          
          kenetic = 0.0_rp
          do i=1, n
            do j=1, n
              kenetic = kenetic + 0.5*(n+1-i)*(l**2)*theta_d(i)*theta_d(j)*cos(theta(i)-theta(j))
            end do
          end do
          
          potential = 0.0_rp
          do i=1, n
            potential = potential -g*(n+1-i)*l*cos(theta(i))
          end do
          
          energy = kenetic + potential
        end function calculate_energy

        subroutine theta_dd_solver(n, l, g, theta, theta_d, theta_dd)
          integer :: i, j, solve_ok
          integer, allocatable :: pivot(:)
          real(rp) :: crossterms
          integer, intent(in) :: n
          real(rp), intent(in) :: l, g
          real(rp), intent(in), allocatable :: theta(:), theta_d(:)
          real(rp), intent(out), allocatable :: theta_dd(:)
          real(rp), allocatable :: A(:, :), b(:)

          allocate(pivot(n))
          allocate(A(n, n))
          allocate(b(n))
          allocate(theta_dd(n))

          do i = 1, n
            crossterms = 0.0_rp
            do j = 1, n
              crossterms = crossterms + bnij(n, i, j)*l*theta_d(j)*theta_d(j)*sin(theta(i)-theta(j))
            end do
            b(i) = -g*(n+1-i)*sin(theta(i)) - crossterms
          end do

          do i = 1, n
            do j = 1, n
              A(i, j) = bnij(n,i,j)*l*cos(theta(i)-theta(j))
            end do
          end do
          
          if (rp.EQ.sp) then
            call sgesv(n, 1, A, n, pivot, b, n, solve_ok)
          else if (rp.EQ.dp) then
            call dgesv(n, 1, A, n, pivot, b, n, solve_ok)
          else if (rp.EQ.qp) then
            print *, "quad precision solver is not supported."
            stop "error"
          else
            print *, "selected precision not recognized"
            stop "error"
          end if

          if (solve_ok.GT.0) then
            print *, "The matrix is singular, solution could not be found", solve_ok
            stop "error"
          else if (solve_ok.LT.0) then
            print *, "The matrix contains illegal values", solve_ok
            stop "error"
          end if
          theta_dd = b
        end subroutine theta_dd_solver

        function bnij(n, i, j) result(b)
          integer, intent(in) :: n, i, j
          integer :: b
          b = (n+1-max(i, j))
        end function bnij

        subroutine save_set(file_name, inital_set, n, sets, set_size, integrator, format, find_energy, l, g, dt, data, energy)
          use pyio, only: openpy, appendpy, writepy_logic, writepy_int, writepy_value, writepy_array, closepy
          
          logical, intent(in) :: inital_set, find_energy
          character (len=*), intent(in) :: file_name
          integer :: file_id, i, j
          integer, intent(in) :: n, sets, set_size, integrator, format
          real(rp), intent(in) :: l, g, dt
          real(rp), intent(in), allocatable :: data(:, :), energy(:)
          
          if ((format.lt.1).or.(format.gt.3)) then
            print *, "format must be 1, 2, or 3"
            stop "error"
          end if

          if (inital_set) then
            call openpy(file_id, file_name)
            call writepy_int(file_id, format)
            call writepy_logic(file_id, find_energy)
            call writepy_int(file_id, sets)
            call writepy_int(file_id, set_size)
            call writepy_int(file_id, integrator)
            call writepy_int(file_id, n)
            call writepy_value(file_id, l)
            call writepy_value(file_id, g)
            call writepy_value(file_id, dt)
          else
            call appendpy(file_id, file_name)
          end if

          do i=1, set_size
            do j=1-format, 0
              call writepy_array(file_id, n, data(:, format*i+j))
            end do
            if (find_energy) then
              call writepy_value(file_id, energy(i))
            end if
          end do

          call closepy(file_id)
        end subroutine save_set

        subroutine initalize(file_name, n, sets, set_size, integrator, format, find_energy, l, g, dt, theta, theta_d, theta_dd)
          use pyio, only: openpy, readpy_logic, readpy_int, readpy_value, readpy_array, closepy
        
          character (len=*), intent(in) :: file_name
          integer :: file_id
          logical, intent(out) :: find_energy
          integer, intent(out) :: n, sets, set_size, integrator, format
          real(rp), intent(out) :: l, g, dt
          real(rp), intent(out), allocatable :: theta(:), theta_d(:), theta_dd(:)
        
          call openpy(file_id, file_name)
          call readpy_int(file_id, format)
          call readpy_logic(file_id, find_energy)
          call readpy_int(file_id, sets)
          call readpy_int(file_id, set_size)
          call readpy_int(file_id, integrator)
          call readpy_int(file_id, n)
          call readpy_value(file_id, l)
          call readpy_value(file_id, g)
          call readpy_value(file_id, dt)
          call readpy_array(file_id, n, theta)
          call readpy_array(file_id, n, theta_d)
          call readpy_array(file_id, n, theta_dd)
          call closepy(file_id)
          
        end subroutine initalize
      end module modules
