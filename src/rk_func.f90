module rk_func
  use precision
  use math
  use rk_defs
  use spray_defs

contains

  subroutine rk_init(spray)

    implicit none

    type(spray_t), intent(inout) :: spray

    ! ---------------------------------

    type(rktvd_solver_t), pointer :: solver => null()
    character(len=128) :: scheme

    ! ---------------------------------

    solver => spray%solver%rktvd

    solver%neq = spray%solver%neq

    solver%nzo = spray%nzo

    ! Set type of RK scheme
    if (spray%solver%strang) then
       write(*,*) 'Warning! Strang splitting only implemented with Euler time stepping...'
       solver%scheme = '1'
    else
       solver%scheme = spray%solver%rktvd%scheme
    end if

    select case (trim(solver%scheme))
    case ('1','11','E11','Euler')
       solver%scheme = 'Euler-11'
    case ('2','22','TVDRK-22','TRK-22')
       solver%scheme = 'TVDRK-22'
    case ('3','33','TVDRK-33','TRK-33')
       solver%scheme = 'TVDRK-33'
    case ('4','44','TVDRK-44','TRK-44')
       solver%scheme = 'TVDRK-44'
    case default
       write(*,*) 'Error! Unknown time advancement scheme!'
       write(*,*) 'Following are the options: '
       write(*,*) 'Euler, TVDRK-22, TVDRK-33, TVDRK-44'
    end select

    select case (trim(solver%scheme))

    case ('Euler-11')
       !if (sd%irank == sd%iroot) then
          write(*,*) '====================================================================='
          write(*,*) 'A 1st. order, forward difference Euler scheme for time advancement...'
          write(*,*) '====================================================================='
       !end if

       solver%stage = 1
       solver%order = 1
       call allocate_rktvd(solver)

       solver%coeff(1,1) = 1.0_WP
       solver%coeff(1,2) = 1.0_WP

    case ('TVDRK-22')
       !if (sd%irank == sd%iroot) then
          write(*,*) '================================================================'
          write(*,*) 'A 2nd. order, 2-stage Runge-Kutta scheme for time advancement...'
          write(*,*) '- Reference: Gottlieb and Shu 1998'
          write(*,*) '================================================================'
       !end if

       solver%stage = 2
       solver%order = 2
       call allocate_rktvd(solver)

       solver%coeff(1,1) = 1.0_WP
       solver%coeff(1,2) = 1.0_WP
       
       solver%coeff(2,1) = 0.5_WP
       solver%coeff(2,2) = 0.5_WP
       solver%coeff(2,3) = 0.5_WP

    case ('TVDRK-33')
       !if (sd%irank == sd%iroot) then
          write(*,*) '================================================================'
          write(*,*) 'A 2nd. order, 2-stage Runge-Kutta scheme for time advancement...'
          write(*,*) '- Reference: Gottlieb and Shu 1998'
          write(*,*) '================================================================'
       !end if

       solver%stage = 3
       solver%order = 3
       call allocate_rktvd(solver)

       solver%coeff(1,1) = 1.0_WP
       solver%coeff(1,2) = 1.0_WP

       solver%coeff(2,1) = 0.75_WP
       solver%coeff(2,2) = 0.25_WP
       solver%coeff(2,3) = 0.25_WP

       solver%coeff(3,1) = 1.0_WP/3.0_WP
       solver%coeff(3,2) = 2.0_WP/3.0_WP
       solver%coeff(3,3) = 2.0_WP/3.0_WP

    case ('TVDRK-44')
       !if (sd%irank == sd%iroot) then
          write(*,*) '================================================================'
          write(*,*) 'A 2nd. order, 2-stage Runge-Kutta scheme for time advancement...'
          write(*,*) '- Reference: Gottlieb and Shu 1998'
          write(*,*) '================================================================'
       !end if
          
       solver%stage = 4
       solver%order = 4
       call allocate_rktvd(solver)
       
       solver%coeff(1,1) = 1.0_WP
       solver%coeff(1,2) = 0.5_WP
       
       solver%coeff(2,1) = 649.0_WP/1600.0_WP
       solver%coeff(2,2) = -10890423.0_WP/25193600.0_WP
       solver%coeff(2,3) = 951.0_WP/1600.0_WP
       solver%coeff(2,4) = 5000.0_WP/7873.0_WP

       solver%coeff(3,1) = 53989.0_WP/2500000.0_WP
       solver%coeff(3,2) = -102261.0_WP/5000000.0_WP
       solver%coeff(3,3) = 4806213.0_WP/20000000.0_WP
       solver%coeff(3,4) = -5121.0_WP/20000.0_WP
       solver%coeff(3,5) = 23619.0_WP/32000.0_WP
       solver%coeff(3,6) = 7873.0_WP/10000.0_WP

       solver%coeff(4,1) = 1.0_WP/5.0_WP
       solver%coeff(4,2) = 1.0_WP/10.0_WP
       solver%coeff(4,3) = 6127.0_WP/30000.0_WP
       solver%coeff(4,4) = 1.0_WP/6.0_WP
       solver%coeff(4,5) = 7873.0_WP/30000.0_WP
       solver%coeff(4,6) = 1.0_WP/3.0_WP
       solver%coeff(4,6) = 1.0_WP/6.0_WP

    case default
       write(*,*) 'Error: unsupported RK scheme...',trim(solver%scheme)
    end select

  end subroutine rk_init

end module rk_func
