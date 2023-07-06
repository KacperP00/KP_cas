! =================================================================================================================================== !
! writes memory information in Mb
! to activate,
!   - function get_memusage_f() needs to be compiled in src/library using the appropriate platform-dependent preprocessor flag, e.g.
!     CFLAGS   = -DINTEL
!   - input file: Monitor preset : memory
! =================================================================================================================================== !

module monitor_memory 
  use precision
  use string
  use math
   implicit none
   
   private  ! default property
   
   ! public members
   ! =================================================
   
   public :: monitor_memory_update
 
   ! private members
   ! =================================================

   real(WP) :: mem_old_

contains

   subroutine monitor_memory_update()

      implicit none

      ! ---------------------------------
      real(WP) :: rp, vp
      integer, pointer :: ip
      character(len=str_medium) :: name
      integer :: ranks_(3),ranks(3)
      real(WP) :: mem_(5),mem(5),tot_mem(2),drmem
      ! ---------------------------------

      mem_ = 0.0_WP
#ifdef INTEL
      call get_memusage_f(mem_(1),mem_(2))

      ! real memory
      rp = mem_(1)
      ! virtual memory
      vp = mem_(2)
      
      drmem = 0.0_WP
      if (mem_old_.gt.0.0_WP) drmem = mem_(1) - mem_old_
      mem_(5) = drmem
      mem_old_ = mem_(1)

      write(*,*) 'Real (MB): ',rp,'\tVirtual (MB):',vp
#endif
      return
   end subroutine monitor_memory_update

end module monitor_memory
