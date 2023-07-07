!--------------------------------------------------------------------------------!
! See the LICENSE file for license information. Please report all bugs and       !
! problems to abhishekd18 at gmail.com                                           !
!--------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------!
! Adapted from CIAO code at Institute for Combustion Technology, RWTH Aachen     !
! University. Due credit goes to the writers of the CIAO code                    !
!--------------------------------------------------------------------------------!
module precision
  implicit none
  integer, parameter  :: SP = kind(1.0)
  integer, parameter, private :: DP = kind(1.0d0)
  integer, parameter  :: WP = DP
  real(WP), private   :: sample_real_at_WP
  real(WP), parameter :: MAX_REAL_WP = HUGE(sample_real_at_WP)
  real(WP), parameter :: MIN_REAL_WP = TINY(sample_real_at_WP)
  integer, private    :: sample_int
  integer, parameter  :: MAX_INTEGER = HUGE(sample_int)
  ! for future extension
  integer, parameter  :: SPI = kind(1)
  integer, parameter  :: WPI = SPI
  integer, parameter  :: WPC = DP
  logical :: verbose=.true.
end module precision
