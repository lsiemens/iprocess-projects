      module modules
        use types, only: rp
        implicit none
        public :: theta_dd_solver, initalize, zeros
        contains

        subroutine theta_dd_solver(n, l, g, theta, theta_d, theta_dd)
          implicit none
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
              crossterms = crossterms + bnij(n, i, j)*l*theta_d(j)*sin(theta(i)-theta(j))
            end do
            b(i) = -g*(n+1-i)*sin(theta(i)) - crossterms
          end do

          do i = 1, n
            do j = 1, n
              A(i, j) = bnij(n,i,j)*l*cos(theta(i)-theta(j))
            end do
          end do

          print *, A
          print *, b
          call sgesv(n, 1, A, n, pivot, b, n, solve_ok)

          if (solve_ok.GT.0) then
            print *, "The matrix is singular, solution could not be found", solve_ok
          else if (solve_ok.LT.0) then
            print *, "The matrix contains illegal values", solve_ok
          end if
          theta_dd = b
        end subroutine theta_dd_solver

        function bnij(n, i, j) result(b)
          integer, intent(in) :: n, i, j
          integer :: b
          b = (n+1-max(i, j))
        end function bnij

        subroutine initalize(file_name, n, l, g, theta, theta_d, theta_dd)
          use pyio, only: openpy, readpy_value, readpy_array, closepy
          implicit none
        
          character (len=*), intent(in) :: file_name
          integer :: file_id
          integer, intent(out) :: n
          real(rp), intent(out) :: l, g
          real(rp), intent(out), allocatable :: theta(:), theta_d(:), theta_dd(:)
        
          call openpy(file_id, file_name)
          call readpy_value(file_id, l)
          call readpy_value(file_id, g)
          call readpy_array(file_id, n, theta)
          call readpy_array(file_id, n, theta_d)
          call readpy_array(file_id, n, theta_dd)
          call closepy(file_id)
        end subroutine initalize

        subroutine zeros(file_name, n)
          use pyio, only: openpy, writepy_value, writepy_array, closepy
          implicit none
        
          character (len=*), intent(in) :: file_name
          integer :: file_id, i
          integer, intent(in) :: n
          real(rp) :: l, g
          real(rp), allocatable :: theta(:), theta_d(:), theta_dd(:)
          
          l = 0.0
          g = 0.0
          allocate(theta(n))
          allocate(theta_d(n))
          allocate(theta_dd(n))

          theta = (/ (0.0, i=1, n) /)
          theta_d = (/ (0.0, i=1, n) /)
          theta_dd = (/ (0.0, i=1, n) /)

          call openpy(file_id, file_name)
          call writepy_value(file_id, l)
          call writepy_value(file_id, g)
          call writepy_array(file_id, n, theta)
          call writepy_array(file_id, n, theta_d)
          call writepy_array(file_id, n, theta_dd)
          call closepy(file_id)
        end subroutine zeros
      end module modules
