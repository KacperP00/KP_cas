!--------------------------------------------------------------------------------!
! See the LICENSE file for license information. Please report all bugs and       !
! problems to abhishekd18 at gmail.com                                           !
!--------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------!
! Adapted from CIAO code at Institute for Combustion Technology                  !
! Due credit goes to writers of the CIAO code                                    !
!--------------------------------------------------------------------------------!
module string
  use, intrinsic :: iso_c_binding, only: c_char, c_null_char, c_ptr, c_loc

  implicit none

  integer, parameter :: str_short  = 8
  integer, parameter :: str_medium = 64
  integer, parameter :: str_long   = 4096
  integer, parameter :: str_control_command = 1024

  interface string_f2c
     module procedure string_f2c_std
     module procedure string_f2c_len
  end interface string_f2c

  interface string_c2f
     module procedure string_c2f_std
     module procedure string_c2f_len
  end interface string_c2f

  interface string_array_1d_f2c
     module procedure string_array_1d_f2c_std
     module procedure string_array_1d_f2c_ptr
  end interface string_array_1d_f2c

  interface string_array_1d_c2f
     module procedure string_array_1d_c2f_len
  end interface string_array_1d_c2f
  
  private :: string_f2c_std, string_f2c_len, &
             string_c2f_std, string_c2f_len, &
             string_array_1d_f2c_std, string_array_1d_f2c_ptr, &
             string_array_1d_c2f_len
             
  private :: c_char, c_null_char, c_ptr, c_loc
 
contains

  function up_case(str) result(str_out)

     implicit none

     character(len=*), intent(in) :: str
     character(len=len(str)) :: str_out

     ! ---------------------------------

     integer :: i,iachar_a,iachar_z,delta,this_iachar

     ! ---------------------------------

     iachar_a = IACHAR('a')
     iachar_z = IACHAR('z')
     delta = IACHAR('A') - iachar_a

     do i = 1,len(str)
        this_iachar = IACHAR(str(i:i))
        if ((iachar_a <= this_iachar) .and. (this_iachar <= iachar_z)) then
           str_out(i:i) = ACHAR(this_iachar + delta)
        else
           str_out(i:i) = str(i:i)
        end if
     end do

  end function up_case

  function lo_case(str) result(str_out)

     implicit none

     character(len=*), intent(in) :: str
     character(len=len(str)) :: str_out

     ! ---------------------------------

     integer :: i,iachar_a,iachar_z,delta,this_iachar

     ! ---------------------------------

     iachar_a = IACHAR('A')
     iachar_z = IACHAR('Z')
     delta = IACHAR('a') - iachar_a

     do i = 1,len(str)
        this_iachar = IACHAR(str(i:i))
        if ((iachar_a <= this_iachar).AND.(this_iachar <= iachar_z)) then
           str_out(i:i) = ACHAR(this_iachar + delta)
        else
           str_out(i:i) = str(i:i)
        end if
     end do

  end function lo_case

  function up_trim(str) result(str_out)

     implicit none

     character(len=*), intent(in) :: str
     character(len=len_trim(str)) :: str_out

     ! ---------------------------------
     
     str_out = up_case(trim(str))

  end function up_trim

  function lo_trim(str) result(str_out)

     implicit none

     character(len=*), intent(in) :: str
     character(len=len_trim(str)) :: str_out

     ! ---------------------------------
     
     str_out = lo_case(trim(str))

  end function lo_trim
  
  function string_f2c_std(f_string) result(c_string)
     implicit none

     character(len=*), intent(in) :: f_string
     character(len=1,kind=c_char) :: c_string(len_trim(f_string)+1)

     ! ---------------------------------

     integer :: n, i

     ! ---------------------------------

     n = len_trim(f_string)
     do i = 1,n
        c_string(i) = f_string(i:i)
     end do
     c_string(n + 1) = c_null_char

  end function string_f2c_std

  function string_f2c_len(f_string,len_c_str) result(c_string)
     implicit none

     character(len=*), intent(in) :: f_string
     integer, intent(in) :: len_c_str
     character(len=1,kind=c_char) :: c_string(len_c_str)

     ! ---------------------------------

     integer :: n, i

     ! ---------------------------------

     n = len_trim(f_string)
     do i = 1,n
        c_string(i) = f_string(i:i)
     end do
     c_string(n + 1:len_c_str) = c_null_char

  end function string_f2c_len


  function string_c2f_std(c_string) result(f_string)
     implicit none

     character(len=1,kind=c_char), intent(in) :: c_string(:)
     character(len=ubound(c_string,1)) :: f_string
     
     ! ---------------------------------

     integer :: n, i

     ! ---------------------------------
     
     f_string = ' '

     n = ubound(c_string,1)
     do i = 1,n
        if (c_string(i) == c_null_char) then
           f_string(i:i) = ' '
           cycle
        end if
        f_string(i:i) = c_string(i)
     end do

  end function string_c2f_std
  
  function string_c2f_medium(c_string) result(f_string)
     implicit none

     character(len=1,kind=c_char), intent(in) :: c_string(:)
     character(len=str_medium) :: f_string

     ! ---------------------------------

     integer :: n, i

     ! ---------------------------------
     
     f_string = ' '

     n = ubound(c_string,1)
     do i = 1,n
        if (c_string(i) == c_null_char) then
           f_string(i:i) = ' '
           cycle
        end if
        f_string(i:i) = c_string(i)
     end do

  end function string_c2f_medium

  function string_c2f_len(c_string,len_f_str) result(f_string)
     implicit none

     character(len=1,kind=c_char), intent(in) :: c_string(:)
     integer, intent(in) ::len_f_str
     character(len=len_f_str) :: f_string

     ! ---------------------------------

     integer :: n, i

     ! ---------------------------------
     
     f_string = ' '

     n = ubound(c_string,1)
     do i = 1,n
        if (c_string(i) == c_null_char) then
           f_string(i:i) = ' '
           cycle
        end if
        f_string(i:i) = c_string(i)
     end do

  end function string_c2f_len
  
  
  subroutine string_array_1d_c2ptr(c_string,c_string_ptr,CLEAR_STR)
     implicit none

     character(len=1,kind=c_char), intent(inout), target :: c_string(:,:)
     type(c_ptr), pointer :: c_string_ptr(:)
     logical, intent(in), optional :: CLEAR_STR
 
     ! ---------------------------------

     logical :: CLEAR_STR_
     integer :: n, i

     ! ---------------------------------

     CLEAR_STR_ = .FALSE.
     if (present(CLEAR_STR)) CLEAR_STR_ = CLEAR_STR 
      
     n = ubound(c_string,2)
     
     if (.not.associated(c_string_ptr)) then
        allocate(c_string_ptr(n))

     else if (ubound(c_string_ptr,1) /= n) then
        deallocate(c_string_ptr)
        allocate(c_string_ptr(n))

     end if
     
     if (CLEAR_STR_) then
        do i = 1,n
           call clear_string_c(c_string(:,i))
        end do
     end if

     do i = 1,n
        c_string_ptr(i) = c_loc(c_string(1,i))
     end do

  end subroutine string_array_1d_c2ptr

  subroutine string_array_1d_f2c_std(f_string,c_string)
     implicit none

     character(len=*), intent(in) :: f_string(:)
     character(len=1,kind=c_char), target :: c_string(:,:)

     ! ---------------------------------

     integer :: n, i

     ! ---------------------------------

     n = ubound(f_string,1)
     
     do i = 1,n
        call clear_string_c(c_string(:,i))
        c_string(:,i) = string_f2c(f_string(i))
     end do

  end subroutine string_array_1d_f2c_std
  
  subroutine string_array_1d_f2c_ptr(f_string,c_string,c_string_ptr)
     implicit none

     character(len=*), intent(in) :: f_string(:)
     character(len=1,kind=c_char), target :: c_string(:,:)
     type(c_ptr), pointer :: c_string_ptr(:)

     ! ---------------------------------

     integer :: n, i

     ! ---------------------------------

     n = ubound(f_string,1)
     
     if (.not.associated(c_string_ptr)) then
        allocate(c_string_ptr(n))

     else if (ubound(c_string_ptr,1) /= n) then
        deallocate(c_string_ptr)
        allocate(c_string_ptr(n))

     end if

     do i = 1,n
        call clear_string_c(c_string(:,i))
        c_string(:,i) = string_f2c(f_string(i))
        c_string_ptr(i) = c_loc(c_string(1,i))
     end do

  end subroutine string_array_1d_f2c_ptr


  subroutine string_array_1d_c2f_medium(c_string,f_string)
     implicit none

     character(len=1,kind=c_char), intent(in) :: c_string(:,:)
     character(len=str_medium) :: f_string(:)

     ! ---------------------------------

     integer :: n, i

     ! ---------------------------------

     n = ubound(c_string,2)
     
     do i = 1,n
        f_string(i) = string_c2f_medium(c_string(:,i))
     end do

  end subroutine string_array_1d_c2f_medium

  subroutine string_array_1d_c2f_len(c_string,len_f_str,f_string)
     implicit none

     character(len=1,kind=c_char), intent(in) :: c_string(:,:)
     integer, intent(in) ::len_f_str
     character(len=len_f_str) :: f_string(:)

     ! ---------------------------------

     integer :: n, i

     ! ---------------------------------

     n = ubound(c_string,2)
     
     do i = 1,n
        f_string(i) = string_c2f(c_string(:,i),len_f_str)
     end do

  end subroutine string_array_1d_c2f_len


  subroutine clear_string_c(c_string)
     implicit none

     character(len=1,kind=c_char), intent(inout) :: c_string(:)

     ! ---------------------------------

     integer ::  i

     ! ---------------------------------
     
     do i = 1,ubound(c_string,1)
        c_string(i) = ' '
     end do

  end subroutine clear_string_c

end module string
