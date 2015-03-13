      module pyio
        use types, only: rp
        implicit none
        private :: newunit
        public :: openpy, readpy_int, readpy_value, readpy_array
        public :: writepy_int, writepy_value, writepy_array, closepy
        contains
      
        subroutine openpy(id, file_name)
          integer, intent(out) :: id
          character (len=*), intent(in) :: file_name
        
          call newunit(id)
          open(unit=id, file=file_name)
          return
        end
      
        subroutine readpy_int(id, int)
          integer, intent(in) :: id
          integer, intent(out) :: int

          read(id, *) int
          return
        end

        subroutine readpy_value(id, value)
          integer, intent(in) :: id
          real(rp), intent(out) :: value

          read(id, *) value
          return
        end

        subroutine readpy_array(id, n, array)
          integer, intent(in) :: id, n
          real(rp), intent(out), allocatable :: array(:)

          allocate(array(n))
          read(id, *) array
          return
        end
      
        subroutine writepy_int(id, int)
          integer, intent(in) :: id, int
        
          write(id, *) int
          return
        end
      
        subroutine writepy_value(id, value)
          integer, intent(in) :: id
          real(rp), intent(in) :: value
        
          write(id, *) value
          return
        end
      
        subroutine writepy_array(id, n, array)
          integer, intent(in) :: id, n
          real(rp), intent(in) :: array(n)

          write(id, *) array
          return
        end
      
        subroutine closepy(id)
          integer, intent(in) :: id
          close(id)
          return
        end
      
        subroutine newunit(unit)
          integer, intent(out) :: unit
          logical :: inuse
          integer, parameter :: nmin=10, nmax = 999
          integer :: n
        
          do n = nmin, nmax
            inquire(unit=n, opened=inuse)
            if (.not. inuse) then
              unit=n
              return
            end if
          end do
        end
      end module pyio
