!!========================================================================
!!>>> service subroutines: set 5, fermi-dirac function                 <<<
!!========================================================================

!!
!! @fun fermi_dirac
!!
!! try to calculate the fermi-dirac function
!!
  function fermi_dirac(omega) result(value)
     use constants, only : dp
     use constants, only : zero, one

     use control, only : beta

     implicit none

! external arguments
! frequency point, \omega
     real(dp), intent(in) :: omega

! result value, return this
     real(dp) :: value

! check the range of omega to avoid numerical instability
     if      ( beta * omega >=  600.0_dp ) then
         value = zero
     else if ( beta * omega <= -600.0_dp ) then
         value = one
     else
         value = one / ( one + exp( beta * omega ) )
     endif ! back if ( beta * omega >=  600.0_dp ) block

     return
  end function fermi_dirac
