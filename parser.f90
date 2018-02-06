module parser
   use precision
   use spray_defs
contains

   subroutine read_inputs(spray,fname)
      implicit none

      ! ---------------------------------
      type(spray_t), intent(inout) :: spray
      character(len=128), intent(in) :: fname

      ! ---------------------------------
      logical :: exist_file
      character(len=128) :: line
      integer :: ioerr, i, idx
      character(len=128) :: param

      inquire(file=fname, exist=exist_file)
      
      if(exist_file) then
         open(unit=100,file=fname,form="formatted",status="old",action="read")
         write(*,*) 'Reading input parameters...'
         do while (.true.)
            read(unit=100,fmt='(a)',iostat=ioerr) line
            
            if (ioerr .ne. 0) then
               exit
            end if

            if (index(line,'!') .or. index(line,'#')) then
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
               else if(param .eq. 'Liquid fuel properties table') then
                  read(line(idx+1:len(line)),*) spray%LFPTname
               else if(param .eq. 'Normal boiling point') then
                  read(line(idx+1:len(line)),*) spray%NBP
               else if(param .eq. 'Melting point') then
                  read(line(idx+1:len(line)),*) spray%MP
               else if(param .eq. 'Vapor fuel properties table') then
                  read(line(idx+1:len(line)),*) spray%VFPTname
               else if(param .eq. 'Fuel molecular weight') then
                  read(line(idx+1:len(line)),*) spray%MW_f
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
               else if(param .eq. 'Initial Dm') then
                  read(line(idx+1:len(line)),*) spray%init_dm
               else if(param .eq. 'Initial D2') then
                  read(line(idx+1:len(line)),*) spray%init_d2
               else if(param .eq. 'Initial D3') then
                  read(line(idx+1:len(line)),*) spray%init_d3
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
               else if(param .eq. 'RK scheme') then
                  read(line(idx+1:len(line)),*) spray%solver%rk%scheme
               else if(param .eq. 'Conv scheme') then
                  read(line(idx+1:len(line)),*) spray%solver%scheme
               else if(param .eq. 'Data File') then
                  read(line(idx+1:len(line)),*) spray%datafilename
               else if(param .eq. 'PL File') then
                  read(line(idx+1:len(line)),*) spray%plfilename
               else if(param .eq. 'Output frequency') then
                  read(line(idx+1:len(line)),*) spray%outfreq
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
               else
               end if
            end if
         end do
         close(unit=100)
      end if

    end subroutine read_inputs

end module parser
