      module types
        implicit none

!----    Select real presision    -------!
!        integer, parameter :: rp = selected_real_kind(6, 37)
        integer, parameter :: rp = selected_real_kind(15, 307)
!        integer, parameter :: rp = selected_real_kind(33, 4931)

!----    Known types   ------!
        integer, parameter :: sp = selected_real_kind(6, 37)
        integer, parameter :: dp = selected_real_kind(15, 307)
        integer, parameter :: qp = selected_real_kind(33, 4931)
      end module
