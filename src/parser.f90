module parser
   use precision
   use spray_defs
contains

   subroutine read_inputs(spray)
      implicit none

      ! ---------------------------------
      type(spray_t), intent(inout) :: spray

      ! ---------------------------------
      logical :: exist_file
      character(len=128) :: line
      integer :: ioerr, i, idx
      character(len=128) :: param

      inquire(file=spray%inp_fname, exist=exist_file)
      
      if(exist_file) then
         open(unit=100,file=spray%inp_fname,form="formatted",status="old",action="read")
         write(*,*) 'Reading input parameters...'
         do while (.true.)
            read(unit=100,fmt='(a)',iostat=ioerr) line
            
            if (ioerr .ne. 0) then
               exit
            end if

            if (index(line,'!') > 0 .or. index(line,'#') > 0) then
               cycle
            end if

            if (index(line,':') > 0) then
               idx = index(line,':')
               param = trim(line(1:idx-1))
               if(param .eq. 'Ambient pressure') then
                  read(line(idx+1:len(line)),*) spray%P_a
               else if(param .eq. 'Ambient temperature') then
                  read(line(idx+1:len(line)),*) spray%T_a
               else if(param .eq. 'Ambient molecular weight') then
                  read(line(idx+1:len(line)),*) spray%MW_a
               else if(param .eq. 'Compressibility factor') then
                  read(line(idx+1:len(line)),*) spray%Z_a
               else if(param .eq. 'Fuel') then
                  read(line(idx+1:len(line)),*) spray%Fuel
               else if(param .eq. 'Fuel temperature') then
                  read(line(idx+1:len(line)),*) spray%T_fuel
               else if(param .eq. 'Fuel properties table') then
                  read(line(idx+1:len(line)),*) spray%FPTname
               else if(param .eq. 'Fuel molecular weight') then
                  read(line(idx+1:len(line)),*) spray%MW_f
               else if(param .eq. 'Use nozzle flow model') then
                  read(line(idx+1:len(line)),*) spray%noz_flow_model
               else if(param .eq. 'Injection pressure') then
                  read(line(idx+1:len(line)),*) spray%P_inj
               else if(param .eq. 'Nozzle exit diameter') then
                  read(line(idx+1:len(line)),*) spray%noz_D
               else if(param .eq. 'Nozzle LD ratio') then
                  read(line(idx+1:len(line)),*) spray%noz_LD
               else if(param .eq. 'Nozzle rD ratio') then
                  read(line(idx+1:len(line)),*) spray%noz_rD
               else if(param .eq. 'Nozzle sac diameter') then
                  read(line(idx+1:len(line)),*) spray%noz_Dsac
               else if(param .eq. 'Nozzle K inlet') then
                  read(line(idx+1:len(line)),*) spray%K_in
               else if(param .eq. 'Nozzle K expansion') then
                  read(line(idx+1:len(line)),*) spray%K_exp
               else if(param .eq. 'Constant injection velocity') then
                  read(line(idx+1:len(line)),*) spray%const_inj_vel
               else if(param .eq. 'Injected mass') then
                  read(line(idx+1:len(line)),*) spray%inj_mass
               else if(param .eq. 'Area of nozzle') then
                  read(line(idx+1:len(line)),*) spray%Anoz
               else if(param .eq. 'Discharge coefficient of nozzle') then
                  read(line(idx+1:len(line)),*) spray%C_d
               else if(param .eq. 'Number of nozzles') then
                  read(line(idx+1:len(line)),*) spray%num_noz
               else if(param .eq. 'Number of grid points') then
                  read(line(idx+1:len(line)),*) spray%Nz
               else if(param .eq. 'Length') then
                  read(line(idx+1:len(line)),*) spray%Lz
               else if(param .eq. 'Maximum CFL') then
                  read(line(idx+1:len(line)),*) spray%MaxCFL
               else if(param .eq. 'Maximum time') then
                  read(line(idx+1:len(line)),*) spray%final_time
               else if(param .eq. 'Initial DSD') then
                  read(line(idx+1:len(line)),*) spray%init_dsd_name
               else if(param .eq. 'Number of diameter classes') then
                  read(line(idx+1:len(line)),*) spray%nd
               else if(param .eq. 'Droplet grid factor') then
                  read(line(idx+1:len(line)),*) spray%dgf
               else if(param .eq. 'Initial Dm') then
                  read(line(idx+1:len(line)),*) spray%init_dm
               else if(param .eq. 'Initial D2') then
                  read(line(idx+1:len(line)),*) spray%init_d2
               else if(param .eq. 'Initial D3') then
                  read(line(idx+1:len(line)),*) spray%init_d3
               else if(param .eq. 'Initial Dvar') then
                  read(line(idx+1:len(line)),*) spray%init_dvar
               else if(param .eq. 'Use entrainment model') then
                  read(line(idx+1:len(line)),*) spray%ent_model
               else if(param .eq. 'Use drag model') then
                  read(line(idx+1:len(line)),*) spray%drag_model
               else if(param .eq. 'Use breakup model') then
                  read(line(idx+1:len(line)),*) spray%bre_model
               else if(param .eq. 'Use evaporation model') then
                  read(line(idx+1:len(line)),*) spray%evap_model
               else if(param .eq. 'B0') then
                  read(line(idx+1:len(line)),*) spray%B0
               else if(param .eq. 'B1') then
                  read(line(idx+1:len(line)),*) spray%B1
               else if(param .eq. 'C3') then
                  read(line(idx+1:len(line)),*) spray%C3
               else if(param .eq. 'Crel') then
                  read(line(idx+1:len(line)),*) spray%Crel
               else if(param .eq. 'Cevap') then
                  read(line(idx+1:len(line)),*) spray%Cevap
               else if(param .eq. 'C_chi') then
                  read(line(idx+1:len(line)),*) spray%C_chi
               else if(param .eq. 'Tolerance for Newton solver') then
                  read(line(idx+1:len(line)),*) spray%solver%nr%tol
               else if(param .eq. 'Initial relaxation coefficient for Newton solver') then
                  read(line(idx+1:len(line)),*) spray%solver%nr%relax_coeff
               else if(param .eq. 'Alpha for Newton solver') then
                  read(line(idx+1:len(line)),*) spray%solver%nr%alpha
               else if(param .eq. 'Maximum iterations for Newton solver') then
                  read(line(idx+1:len(line)),*) spray%solver%nr%max_count
               else if(param .eq. 'Spray half-cone angle model') then
                  read(line(idx+1:len(line)),*) spray%spray_angle_model
               else if(param .eq. 'C_theta') then
                  read(line(idx+1:len(line)),*) spray%C_theta
               else if(param .eq. 'Spray half-cone angle') then
                  read(line(idx+1:len(line)),*) spray%theta
               else if(param .eq. 'Tangent of spray half-cone angle') then
                  read(line(idx+1:len(line)),*) spray%beta
               else if(param .eq. 'ROI profile') then
                  read(line(idx+1:len(line)),*) spray%roi_file
               else if(param .eq. 'Strang splitting') then
                  read(line(idx+1:len(line)),*) spray%solver%strang
               else if(param .eq. 'RK scheme') then
                  read(line(idx+1:len(line)),*) spray%solver%rktvd%scheme
               else if(param .eq. 'Conv scheme') then
                  read(line(idx+1:len(line)),*) spray%solver%scheme
               else if(param .eq. 'Data File') then
                  read(line(idx+1:len(line)),*) spray%datafilename
               else if(param .eq. 'Datafile frequency') then
                  read(line(idx+1:len(line)),*) spray%datafreq
               else if(param .eq. 'PL File') then
                  read(line(idx+1:len(line)),*) spray%plfilename
               else if(param .eq. 'Output frequency') then
                  read(line(idx+1:len(line)),*) spray%outfreq
               else if(param .eq. 'Turbulence model') then
                  read(line(idx+1:len(line)),*) spray%turb_model
               else if(param .eq. 'Combustion model') then
                  read(line(idx+1:len(line)),*) spray%combustion_model
               else if(param .eq. 'Flamelet interval') then
                  read(line(idx+1:len(line)),*) spray%flm_int
               else if(param .eq. 'Number of flamelets') then
                  read(line(idx+1:len(line)),*) spray%nflamelet
               else if(param .eq. 'Y_O2') then
                  read(line(idx+1:len(line)),*) spray%Y_O2
               else if(param .eq. 'Fixed Re') then
                  read(line(idx+1:len(line)),*) spray%fixed_Re
               else if(param .eq. 'Fixed We') then
                  read(line(idx+1:len(line)),*) spray%fixed_We
               else if(param .eq. 'Fixed DRa') then
                  read(line(idx+1:len(line)),*) spray%fixed_DRa
               else if(param .eq. 'Fixed DRv') then
                  read(line(idx+1:len(line)),*) spray%fixed_DRv
               else if(param .eq. 'Fixed VRa') then
                  read(line(idx+1:len(line)),*) spray%fixed_VRa
               else if(param .eq. 'Fixed VRv') then
                  read(line(idx+1:len(line)),*) spray%fixed_VRv
               else if(param .eq. 'Fixed De') then
                  read(line(idx+1:len(line)),*) spray%fixed_De
               else if(param .eq. 'f_rho_l') then
                  read(line(idx+1:len(line)),*) spray%f_rho_l
               else if(param .eq. 'f_mu_l') then
                  read(line(idx+1:len(line)),*) spray%f_mu_l
               else if(param .eq. 'f_lambda_l') then
                  read(line(idx+1:len(line)),*) spray%f_lambda_l
               else if(param .eq. 'f_Lv') then
                  read(line(idx+1:len(line)),*) spray%f_Lv
               else if(param .eq. 'f_C_l') then
                  read(line(idx+1:len(line)),*) spray%f_C_l
               else if(param .eq. 'f_Pv') then
                  read(line(idx+1:len(line)),*) spray%f_Pv
               else if(param .eq. 'f_sigma') then
                  read(line(idx+1:len(line)),*) spray%f_sigma
               else if(param .eq. 'f_rho_v') then
                  read(line(idx+1:len(line)),*) spray%f_rho_v
               else if(param .eq. 'f_mu_v') then
                  read(line(idx+1:len(line)),*) spray%f_mu_v
               else if(param .eq. 'f_lambda_v') then
                  read(line(idx+1:len(line)),*) spray%f_lambda_v
               else if(param .eq. 'f_Cp_v') then
                  read(line(idx+1:len(line)),*) spray%f_Cp_v
               else if(param .eq. 'f_G_v') then
                  read(line(idx+1:len(line)),*) spray%f_G_v
               else
               end if
            end if
         end do
         close(unit=100)
      end if

    end subroutine read_inputs

   subroutine read_param(fname,input,val,default)
      implicit none

      ! ---------------------------------
      character(len=*), intent(in) :: fname, input, default
      character(len=*), intent(out) :: val

      ! ---------------------------------
      logical :: exist_file
      character(len=128) :: line
      integer :: ioerr, i, idx
      character(len=128) :: param

      inquire(file=fname, exist=exist_file)
      
      if(exist_file) then
         open(unit=100,file=fname,form="formatted",status="old",action="read")
         write(*,*) 'Reading ',trim(input),' ...'
         do while (.true.)
            read(unit=100,fmt='(a)',iostat=ioerr) line
            
            if (ioerr .ne. 0) then
               exit
            end if

            if (index(line,'!') > 0 .or. index(line,'#') > 0) then
               cycle
            end if

            if (index(line,':') > 0) then
               idx = index(line,':')
               param = trim(line(1:idx-1))
               if(param .eq. trim(input)) then
                  read(line(idx+1:len(line)),*) val
               else
                  val = default
               end if
            end if
         end do
         close(unit=100)
      end if

    end subroutine read_param

   subroutine modify_inputs(spray)
      implicit none

      ! ---------------------------------
      type(spray_t), intent(inout) :: spray

      ! ---------------------------------
      logical :: exist_file
      character(len=128) :: fname, line, param, tmp_char
      integer :: ioerr, i, idx

      fname = 'modify.in'

      inquire(file=fname, exist=exist_file)
      
      if(exist_file) then
         open(unit=100,file=fname,form="formatted",status="old",action="read")
         do while (.true.)
            read(unit=100,fmt='(a)',iostat=ioerr) line
            
            if (ioerr .ne. 0) then
               exit
            end if

            if (index(line,'!') > 0 .or. index(line,'#') > 0) then
               cycle
            end if

            ! End simulation
            if(trim(line) .eq. 'END') then
               spray%end = .true.
               write(*,*) 'Spray simulation ended by user...'
               exit
            end if

            ! Save data file 
            if(trim(line) .eq. 'SAVE_DATA_FILE') then
               spray%saveDataFile = .true.
               write(*,*) 'Data file written...'
               exit
            end if

            if (index(line,':') > 0) then
               idx = index(line,':')
               param = trim(line(1:idx-1))
               if(param .eq. 'Maximum CFL') then
                  write(tmp_char,*) spray%MaxCFL
                  write(*,*) 'Modifying '//param//' from '//tmp_char//' to '//line(idx+1:len(line))
                  read(line(idx+1:len(line)),*) spray%MaxCFL
               else if(param .eq. 'Maximum time') then
                  write(tmp_char,*) spray%final_time
                  write(*,*) 'Modifying '//param//' from '//tmp_char//' to '//line(idx+1:len(line))
                  read(line(idx+1:len(line)),*) spray%final_time
               else if(param .eq. 'Number of diameter classes') then
                  write(tmp_char,*) spray%nd
                  write(*,*) 'Modifying '//param//' from '//tmp_char//' to '//line(idx+1:len(line))
                  read(line(idx+1:len(line)),*) spray%nd
               else if(param .eq. 'Tolerance for Newton solver') then
                  write(tmp_char,*) spray%solver%nr%tol
                  write(*,*) 'Modifying '//param//' from '//tmp_char//' to '//line(idx+1:len(line))
                  read(line(idx+1:len(line)),*) spray%solver%nr%tol
               else if(param .eq. 'Initial relaxation coefficient for Newton solver') then
                  write(tmp_char,*) spray%solver%nr%relax_coeff
                  write(*,*) 'Modifying '//param//' from '//tmp_char//' to '//line(idx+1:len(line))
                  read(line(idx+1:len(line)),*) spray%solver%nr%relax_coeff
               else if(param .eq. 'Alpha for Newton solver') then
                  write(tmp_char,*) spray%solver%nr%alpha
                  write(*,*) 'Modifying '//param//' from '//tmp_char//' to '//line(idx+1:len(line))
                  read(line(idx+1:len(line)),*) spray%solver%nr%alpha
               else if(param .eq. 'Maximum iterations for Newton solver') then
                  write(tmp_char,*) spray%solver%nr%max_count
                  write(*,*) 'Modifying '//param//' from '//tmp_char//' to '//line(idx+1:len(line))
                  read(line(idx+1:len(line)),*) spray%solver%nr%max_count
               else if(param .eq. 'Data File') then
                  write(tmp_char,*) spray%datafilename
                  write(*,*) 'Modifying '//param//' from '//tmp_char//' to '//line(idx+1:len(line))
                  read(line(idx+1:len(line)),*) spray%datafilename
               else if(param .eq. 'PL File') then
                  write(tmp_char,*) spray%plfilename
                  write(*,*) 'Modifying '//param//' from '//tmp_char//' to '//line(idx+1:len(line))
                  read(line(idx+1:len(line)),*) spray%plfilename
               else if(param .eq. 'Output frequency') then
                  write(tmp_char,*) spray%outfreq
                  write(*,*) 'Modifying '//param//' from '//tmp_char//' to '//line(idx+1:len(line))
                  read(line(idx+1:len(line)),*) spray%outfreq
               else
                  write(*,*) param//' is not modified. Continuing with the old value...'
               end if
            else
               write(*,*) 'Unknown keyword '//trim(line)//' !!! Nothing changed...'
            end if
         end do
         close(unit=100)

         call system('rm '//trim(fname))
      end if

    end subroutine modify_inputs

end module parser
