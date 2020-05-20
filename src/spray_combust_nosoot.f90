! Adapted from knock-indication routine in CIAO from Marco Davidovic
module spray_combust
   use spray_defs
   use combust_mduc_data
   use combust_mduc_interface
   use MPI

   implicit none

   !include 'mpif.h'


   private ! default property
   
   ! public members
   public spray_combust_init
   public spray_combust_poststep

   type(mduc_t), target                 :: mduc_  ! mduc solver memory
   real(WP), dimension(:), pointer      :: x  ! flamelet grid
   integer                              :: nspec
   integer                              :: nprog = 4
   integer, dimension(:), allocatable   :: iprog_species
   real(WP)                             :: PROG_THRESH 
   logical                              :: USE_KNOCK
   

contains

   subroutine spray_combust_init(spray)
      use parser
      implicit none
   
      type(spray_t), intent(inout) :: spray
   
      ! ----------------------------------
   
      character(len=str_medium)         :: mduc_input_file
      type(mduc_t), pointer             :: mduc
      character(len=str_medium), dimension(:), allocatable :: species,prog_species
      real(WP)                          :: Pavg,Hmin,Hmax,nu
      integer                           :: s,c
      ! ----------------------------------

      if(spray%combustion_model == 'None' .or. spray%combustion_model == 'noname') then
         write(*,*) 'No combustion model specified...'
         write(*,*) 'Computing spray without combustion...'
         return
      end if

      mduc => mduc_ 
   
      call read_param(spray%inp_fname,'RIF Input File',mduc_input_file,'mduc.inp')

      ! Initialize mduc solver
      call mduc_init(mduc,trim(mduc_input_file),0,spray%nranks,spray%comm)
   
      ! Check if times of ACAS and MDUC are consistent 
      if (abs(mduc_get_time(mduc,0)-spray%ndtime*spray%tau).gt.0.1_WP*spray%dt*spray%tau) then
         write(*,'(A,E15.6,A,E15.6)') 'Error: RIF and Spray are at different times. RIF: ', &
                                      mduc_get_time(mduc,0),'; spray: ',spray%time
      end if
   
      ! Progress variable
      allocate(prog_species(nprog))
      allocate(iprog_species(nprog))
      prog_species(1) = 'CO'
      prog_species(2) = 'CO2'
      prog_species(3) = 'H2'
      prog_species(4) = 'H2O'

      ! get number of species in flamelet
      nspec = mduc_get_nnons_species(mduc)

      ! get species from flamelet 
      allocate(species(nspec))
      call mduc_get_species_names(mduc,species,len(species))
 
      ! find the progress variable species
      do c=1,nprog
         do s=1,nspec
            if (trim(species(s)).eq.trim(prog_species(c))) then
               iprog_species(c) = s
               write(*,*) trim(prog_species(c)), ' has index' , s
            end if
         end do
      end do
 
   
      ! allocate flamelet grid
      allocate(x(mduc%nx1))

      ! get the flamelet grid
      call mduc_get_grid_1D(mduc,x)

      spray%zz => x

      Hmin = 0.0_WP
      Hmax = 1.0_WP
    
      ! initialize flow data in flamelet solver
      call mduc_init_flow(mduc,(spray%ndtime-spray%dt)*spray%tau,spray%P_a,Hmin,Hmax)
   
      if (spray%step.eq.1) then
         ! overwrite mduc pressure 
         call mduc_set_pressure(mduc,-1,spray%P_a)
      end if
   
      ! update the time step in flamelet solver
      call mduc_update_timestep(spray%step,mduc)
   
      deallocate(species)
      return
   
   end subroutine

   subroutine spray_combust_poststep(spray)
      implicit none

      type(spray_t), intent(inout) :: spray
   
      ! ----------------------------------
  
      type(mduc_t), pointer     :: mduc
      real(WP)                  :: Pavg,Hmin,Hmax,Nu
      real(WP)                  :: fT(mduc_%nx1), &
                                   fY(mduc_%nx1,nspec), &
                                   fS(mduc_%nx1,4), &
                                   chi(mduc_%nx1),&
                                   prog_fl(mduc_%nx1)
      real(WP)                  :: Zmix_st, Zmix_max
      real(WP)                  :: Temp, sootM00, SootNumDen, &
                                   sootM10, SootVolume, &
                                   SootMass, gasDen, cellVol
      integer                   :: i,c, ispec, k
      character(len=str_medium), dimension(:), pointer :: species
      real(WP), dimension(:), pointer :: Yspec, bpdf
      character(len=512) :: fname
      integer, dimension(2) :: bound
      ! ----------------------------------

      if(spray%combustion_model == 'None' .or. spray%combustion_model == 'noname') return

      mduc => mduc_

      !if(spray%chi_st(1) == 0.0_WP) return

      Zmix_st = 1.0_WP/(1.0_WP+spray%stoic_coeff/spray%Y_O2)
      Zmix_max = maxval(spray%Zmix_g) + 2.0_WP*sqrt(maxval(spray%zvar_g))
      Zmix_max = min(1.0,Zmix_max)

      chi = spray%C_chi*spray%chi_st(1)*(x/Zmix_st)**2*(log(x/Zmix_max)/log(Zmix_st/Zmix_max))

      where(chi/=chi) chi = 0.0_WP
      where(chi<0.0_WP) chi = 0.0_WP

      !write(*,*) 'chi ', chi

      Pavg = spray%P_a
      Hmin = 0.0_WP
      Hmax = 1.0_WP

      ! Advance flamelet
      !call mduc_advance_flamelet_1d(mduc,spray%dt*spray%tau,Pavg,Hmin,Hmax,chi,0,.false.,fY,fT,fS)
      call mduc_advance_flamelet_1d(mduc,spray%dt*spray%tau,Pavg,Hmin,Hmax,chi,0,.false.,fY,fT)

      if(spray%irank .eq. 0) then

      if (mod(spray%step,spray%datafreq) == 0 .or. spray%end) then

         ! get number of species in flamelet
         nspec = mduc_get_nnons_species(mduc)

         ! get species from flamelet 
         allocate(species(nspec))
         allocate(Yspec(nspec))
         allocate(bpdf(mduc_%nx1))
         call mduc_get_species_names(mduc,species,len(species))

         write(fname,"(A,A1,I0.6)") trim('species.out'), '_', spray%step

         open(unit=100,file=trim(fname),form="formatted",status="replace",action="write")

         !write(100,'(<nspec+5>(A))') '#z<1> ','Temperature<2> ','SootNumDen<3> ','SootVolume<4> ','SootMass<5> ',(species(i), i=1,nspec)
         write(100,'(<nspec+2>(A))') '#z<1> ','Temperature<2> ',(species(i), i=1,nspec)
         do k = spray%kmin,spray%kmax

            call betaPDF(spray%Zmix_g(k),spray%Zvar_g(k),mduc_%nx1,x,bpdf,bound)

            ! Convolute
            do ispec = 1,nspec

               Yspec(ispec) = sum(fY(bound(1):bound(2),ispec)*bpdf(bound(1):bound(2)))

            end do

!!$            gasDen = spray%rho(k)*spray%rho_l*spray%Y_g(k)
!!$            cellVol = spray%dz*Pi*spray%b(k)**2*spray%D_eff**3
!!$            sootM00 = sum(fS(bound(1):bound(2),1)*bpdf(bound(1):bound(2)))
!!$            SootNumDen = sootM00*gasDen
!!$            sootM10 = sum(fS(bound(1):bound(2),2)*bpdf(bound(1):bound(2)))
!!$            SootVolume = sootM10*gasDen*24.0_WP/1800.0_WP*cellVol
!!$            SootMass = 1800.0_WP*SootVolume

            Temp = sum(fT(bound(1):bound(2))*bpdf(bound(1):bound(2)))

            !write(100, '(<nspec+5>(ES15.5E3))') spray%z(k),Temp,SootNumDen,SootVolume,SootMass,(Yspec(i), i=1,nspec)
            write(100, '(<nspec+2>(ES15.5E3))') spray%z(k), Temp, (Yspec(i), i=1,nspec)

         end do
         close(unit=100)

         !deallocate(species,Yspec,bpdf)
      end if
      end if
      
      prog_fl = 0.0_WP
      do c=1,nprog
         do i = 1,mduc_%nx1
            prog_fl(i) = prog_fl(i) + fY(i,iprog_species(c))
         end do
      end do

      return

   end subroutine

end module 
