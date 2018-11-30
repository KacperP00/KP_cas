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

    type(rk_solver_t), pointer :: solver => null()
    character(len=128) :: scheme

    ! ---------------------------------

    solver => spray%solver%rk

    ! Set type of RK scheme

    solver%scheme = spray%solver%rk%scheme

    select case (trim(solver%scheme))
    case ('1')
       solver%scheme = 'Euler-11'
    case ('2','3','4')
       solver%scheme = 'Williamson-33'
    case ('5')
       solver%scheme = 'Stanescu-25'
    case ('6')
       solver%scheme = 'Allampalli-46'
    case ('7')
       solver%scheme = 'Allampalli-47'
    case ('hu46','hu-46')
       solver%scheme = 'Hu-46'
    case ('will33','will-33')
       solver%scheme = 'Williamson-33'
    case ('wray33','wray-33')
       solver%scheme = 'Wray-33'
    case ('shu33','shu-33','tvd33','tvd-33')
       solver%scheme = 'Shu-33'
    end select

    select case (trim(solver%scheme))

    case ('Williamson-33')
       !if (sd%irank == sd%iroot) then
          write(*,*) '================================================================'
          write(*,*) 'A 3rd. order, 3-stage Runge-Kutta scheme for time advancement...'
          write(*,*) '- Reference: Williamson, JCP, 1980, vol. 35, pp. 48-56.'
          write(*,*) '- A low-storage (2-N) version'
          write(*,*) '- CFL limit with the 2nd. order FD : 1.73'
          write(*,*) '================================================================'
       !end if

       solver%stage = 3
       call allocate_rk(solver)

       solver%alpha(1) =  0.0_WP
       solver%alpha(2) = -5.0_WP/9.0_WP
       solver%alpha(3) = -153.0_WP/128.0_WP
       solver%beta(1) = 1.0_WP/3.0_WP
       solver%beta(2) = 15.0_WP/16.0_WP
       solver%beta(3) = 8.0_WP/15.0_WP
       solver%dt_ratio(1) = 0.0_WP
       solver%dt_ratio(2) = 1.0_WP/3.0_WP
       solver%dt_ratio(3) = 3.0_WP/4.0_WP

    case ('Wray-33')
       !if (sd%irank == sd%iroot) then
          write(*,*) '================================================================'
          write(*,*) 'A 3th. order, 3-stage Runge-Kutta scheme for time advancement...'
          write(*,*) '- A classic one used in several previous studies in CTR...'
          write(*,*) '- Reference: Wray, Tech. Rep. NASA AMES, 1986.'
          write(*,*) '- Tailored to a low-storage (2-N) version'
          write(*,*) '- CFL limit with the 2nd. order FD : 1.7'
          write(*,*) '================================================================'
       !end if

       solver%stage = 3
       call allocate_rk(solver)

       solver%alpha(1) =  0.0_WP
       solver%alpha(2) = -17.0_WP/25.0_WP
       solver%alpha(3) = -5.0_WP/9.0_WP
       solver%beta(1) = 8.0_WP/15.0_WP
       solver%beta(2) = 5.0_WP/12.0_WP
       solver%beta(3) = 3.0_WP/4.0_WP
       solver%dt_ratio(1) = 0.0_WP
       solver%dt_ratio(2) = 8.0_WP/15.0_WP
       solver%dt_ratio(3) = 2.0_WP/3.0_WP

    case ('Shu-33')
       !if (sd%irank == sd%iroot) then
          write(*,*) '================================================================'
          write(*,*) 'A 3th. order, 3-stage Runge-Kutta scheme for time advancement...'
          write(*,*) '- A RK3 with TVD (total variation-diminishing) condition'
          write(*,*) '- Reference: Shu & Osher, JCP, 1988, vol. 77, pp. 439-471.'
          write(*,*) '- Tailored to a low-storage (2-N) version'
          write(*,*) '- CFL limit with the 2nd. order FD : 1.73'
          write(*,*) '================================================================'
       !end if

       solver%stage = 3
       call allocate_rk(solver)

       solver%alpha(1) = 0.0_WP
       solver%alpha(2) = -3.0_WP
       solver%alpha(3) = -0.5_WP
       solver%beta(1) = 1.0_WP
       solver%beta(2) = 1.0_WP/4.0_WP
       solver%beta(3) = 2.0_WP/3.0_WP
       solver%dt_ratio(1) = 0.0_WP
       solver%dt_ratio(2) = 0.5_WP  ! may choose 1.0_WP
       solver%dt_ratio(3) = 0.5_WP

    case ('Stanescu-25')
       !if (sd%irank == sd%iroot) then
          write(*,*) '================================================================'
          write(*,*) 'A 2nd. order, 5-stage Runge-Kutta scheme for time advancement...'
          write(*,*) '- Low-dissipation and low-dispersion Runge-Kutta (LDDRK) scheme'
          write(*,*) '- Reference: Stanescu & Habashi, JCP, 1998, vol. 143, pp. 674-681.'
          write(*,*) '             Hu et al.         , JCP, 1996, vol. 124, pp. 177-191.'
          write(*,*) '- A low-storage (2-N) version'
          write(*,*) '- CFL limit with the 2nd. order FD : 3.52'
          write(*,*) '================================================================'
       !end if

       solver%stage = 5
       call allocate_rk(solver)

       solver%alpha(1) =  0.0_WP
       solver%alpha(2) = -0.6913065_WP
       solver%alpha(3) = -2.655155_WP
       solver%alpha(4) = -0.8147688_WP
       solver%alpha(5) = -0.6686587_WP
       solver%beta(1) = 0.1_WP
       solver%beta(2) = 0.75_WP
       solver%beta(3) = 0.7_WP
       solver%beta(4) = 0.479313_WP
       solver%beta(5) = 0.310392_WP
       solver%dt_ratio(1) = 0.0_WP
       solver%dt_ratio(2) = 0.1_WP
       solver%dt_ratio(3) = 0.3315201_WP
       solver%dt_ratio(4) = 0.4577796_WP
       solver%dt_ratio(5) = 0.8666528_WP

    case ('Allampalli-46')
       !if (sd%irank == sd%iroot) then
          write(*,*) '================================================================'
          write(*,*) 'A 4th. order, 6-stage Runge-Kutta scheme for time advancement...'
          write(*,*) '- High-accuracy large-step explicit Runge-Kutta (HALE-RK) scheme'
          write(*,*) '- Reference: Allampalli et al., JCP, 2009, vol. 228, pp. 3837-3850.'
          write(*,*) '- A low-storage (2-N) version'
          write(*,*) '- CFL limit with the 2nd. order FD : 4.90'
          write(*,*) '================================================================'
       !end if

       solver%stage = 6
       call allocate_rk(solver)

       solver%alpha(1) =  0.000000000000_WP
       solver%alpha(2) = -0.691750960670_WP
       solver%alpha(3) = -1.727127405211_WP
       solver%alpha(4) = -0.694890150986_WP
       solver%alpha(5) = -1.039942756197_WP
       solver%alpha(6) = -1.531977447611_WP
       solver%beta(1) = 0.122000000000_WP
       solver%beta(2) = 0.477263056358_WP
       solver%beta(3) = 0.381941220320_WP
       solver%beta(4) = 0.447757195744_WP
       solver%beta(5) = 0.498614246822_WP
       solver%beta(6) = 0.186648570846_WP
       solver%dt_ratio(1) = 0.000000000000_WP
       solver%dt_ratio(2) = 0.122000000000_WP
       solver%dt_ratio(3) = 0.269115878630_WP
       solver%dt_ratio(4) = 0.447717183551_WP
       solver%dt_ratio(5) = 0.749979795490_WP
       solver%dt_ratio(6) = 0.898555413085_WP

    case ('Allampalli-47')
       !if (sd%irank == sd%iroot) then
          write(*,*) '================================================================'
          write(*,*) 'A 4th. order, 7-stage Runge-Kutta scheme for time advancement...'
          write(*,*) '- High-accuracy large-step explicit Runge-Kutta (HALE-RK) scheme'
          write(*,*) '- Reference: Allampalli et al., JCP, 2009, vol. 228, pp. 3837-3850.'
          write(*,*) '- A low-storage (2-N) version'
          write(*,*) '- CFL limit with the 2nd. order FD : 5.67'
          write(*,*) '================================================================'
       !end if

       solver%stage = 7
       call allocate_rk(solver)

       solver%alpha(1) =  0.000000000000_WP
       solver%alpha(2) = -0.647900745934_WP
       solver%alpha(3) = -2.704760863204_WP
       solver%alpha(4) = -0.460080550118_WP
       solver%alpha(5) = -0.500581787785_WP
       solver%alpha(6) = -1.906532255913_WP
       solver%alpha(7) = -1.450000000000_WP
       solver%beta(1) = 0.117322146869_WP
       solver%beta(2) = 0.503270262127_WP
       solver%beta(3) = 0.233663281658_WP
       solver%beta(4) = 0.283419634625_WP
       solver%beta(5) = 0.540367414023_WP
       solver%beta(6) = 0.371499414620_WP
       solver%beta(7) = 0.136670099385_WP
       solver%dt_ratio(1) = 0.000000000000_WP
       solver%dt_ratio(2) = 0.117322146869_WP
       solver%dt_ratio(3) = 0.294523230758_WP
       solver%dt_ratio(4) = 0.305658622131_WP
       solver%dt_ratio(5) = 0.582864148403_WP
       solver%dt_ratio(6) = 0.858664273599_WP
       solver%dt_ratio(7) = 0.868664273599_WP

    case ('Hu-46')
       !if (sd%irank == sd%iroot) then
          write(*,*) '================================================================'
          write(*,*) 'A 4th. order, 6-stage Runge-Kutta scheme for time advancement...'
          write(*,*) '- Low-dissipation and low-dispersion Runge-Kutta (LDDRK) scheme'
          write(*,*) '- Reference: Stanescu & Habashi, JCP, 1998, vol. 143, pp. 674-681.'
          write(*,*) '             Hu et al.         , JCP, 1996, vol. 124, pp. 177-191.'
          write(*,*) '- A low-storage (2-N) version'
          write(*,*) '- CFL limit with the 2nd. order FD : 1.65'
          write(*,*) '- A better accuracy with reduced stability...'
          write(*,*) '================================================================'
       !end if

       solver%stage = 6
       call allocate_rk(solver)

       solver%alpha(1) =  0.0000000_WP
       solver%alpha(2) = -0.4919575_WP
       solver%alpha(3) = -0.8946264_WP
       solver%alpha(4) = -1.5526678_WP
       solver%alpha(5) = -3.4077973_WP
       solver%alpha(6) = -1.0742640_WP
       solver%beta(1) = 0.1453095_WP
       solver%beta(2) = 0.4653797_WP
       solver%beta(3) = 0.4675397_WP
       solver%beta(4) = 0.7795279_WP
       solver%beta(5) = 0.3574327_WP
       solver%beta(6) = 0.1500000_WP
       solver%dt_ratio(1) = 0.0000000_WP
       solver%dt_ratio(2) = 0.1453095_WP
       solver%dt_ratio(3) = 0.3817422_WP
       solver%dt_ratio(4) = 0.6367813_WP
       solver%dt_ratio(5) = 0.7560744_WP
       solver%dt_ratio(6) = 0.9271047_WP

    case ('Euler-11')
       !if (sd%irank == sd%iroot) then
          write(*,*) '================================================================'
          write(*,*) 'The 1st. order explicit Euler scheme for time advancement...'
          write(*,*) '- CFL limit with the 2nd. order FD : 1.0'
          write(*,*) '================================================================'
       !end if

       solver%stage = 4
       call allocate_rk(solver)

       solver%alpha(1) =  0.0_WP
       solver%beta(1) = 1.0_WP
       solver%dt_ratio(1) = 0.0_WP

    case default
       write(*,*) 'Error: unsupported RK scheme...',trim(solver%scheme)
    end select

  end subroutine rk_init

end module rk_func
