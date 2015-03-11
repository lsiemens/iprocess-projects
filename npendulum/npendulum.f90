      program npendulum
        use types, only: rp
        use modules, only: initalize, zeros
        implicit none

        integer :: n, file_id, i, j
        real(rp) :: l, g
        real(rp), allocatable :: theta(:), theta_d(:), theta_dd(:)
        
        call initalize("initalization.dat", n, l, g, theta, theta_d, theta_dd)
!        call zeros("initalization.dat", 10)
        
        print *, n, l, g
        print *, theta
      end program npendulum
