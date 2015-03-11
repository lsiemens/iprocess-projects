      module modules
        use types, only: rp
        implicit none
        public :: initalize, zeros
        contains

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
