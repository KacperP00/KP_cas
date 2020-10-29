! Adapted from knock-indication routine in CIAO from Marco Davidovic
module spray_combust
   use spray_defs
   use solver_defs
   use combust_mduc_data
   use combust_mduc_interface
   use MPI

   implicit none

   !include 'mpif.h'


   private ! default property
   
   ! public members
   public spray_combust_init
   public spray_combust_poststep_mrif

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
      type(srif_t), pointer             :: mrif, next_mrif
      character(len=str_medium), dimension(:), allocatable :: species,prog_species
      real(WP)                          :: Pavg,Hmin,Hmax,nu
      integer                           :: s,c,i
      ! ----------------------------------

      if(spray%combustion_model == 'None' .or. spray%combustion_model == 'noname') then
         write(*,*) 'No combustion model specified...'
         write(*,*) 'Computing spray without combustion...'
         return
      end if
   
      call read_param(spray%inp_fname,'RIF Input File',mduc_input_file,'mduc.inp')

      if(spray%combustion_model == 'MRIF') then
         call init_flamelets(spray)
      else
         allocate(mduc)
      end if

      mrif=>spray%first_rif
      do while (associated(mrif))
         mduc=>mrif%mduc
         !mduc => mduc_

         ! Initialize mduc solver
         call mduc_init(mduc,trim(mduc_input_file),0,spray%nranks,spray%comm)
   
         ! Check if times of ACAS and MDUC are consistent 
         if (abs(mduc_get_time(mduc,0)-spray%ndtime*spray%tau).gt.0.1_WP*spray%dt*spray%tau) then
            write(*,'(A,E15.6,A,E15.6)') 'Error: RIF and Spray are at different times. RIF: ', &
                 mduc_get_time(mduc,0),'; spray: ',spray%time
         end if
   
         ! Progress variable
         allocate(mrif%prog_species(nprog))
         allocate(mrif%iprog_species(nprog))
         mrif%prog_species(1) = 'CO'
         mrif%prog_species(2) = 'CO2'
         mrif%prog_species(3) = 'H2'
         mrif%prog_species(4) = 'H2O'

         ! get number of species in flamelet
         nspec = mduc_get_nnons_species(mduc)

         ! get species from flamelet 
         allocate(mrif%species(nspec))
         call mduc_get_species_names(mduc,mrif%species,len(mrif%species))
 
         ! find the progress variable species
         do c=1,nprog
            do s=1,nspec
               if (trim(mrif%species(s)).eq.trim(mrif%prog_species(c))) then
                  mrif%iprog_species(c) = s
                  write(*,*) trim(mrif%prog_species(c)), ' has index' , s
               end if
            end do
         end do
 
         ! allocate flamelet grid
         allocate(x(mduc%nx1))

         ! get the flamelet grid
         call mduc_get_grid_1D(mduc,x)

         !spray%zz => x
         mrif%zz => x
         ! allocate flamelet solution variables
         allocate(mrif%fY_1d(mduc%nx1,nspec))
         allocate(mrif%fT_1d(mduc%nx1))
         allocate(mrif%chi(mduc%nx1))
#ifdef SOOT
         allocate(mrif%fS_1d(mduc%nx1,4))
         allocate(mrif%sootM00(spray%nzo)); mrif%sootM00= 0.0_WP
         allocate(mrif%sootM10(spray%nzo)); mrif%sootM10 = 0.0_WP
#endif
         ! get species from flamelet
         allocate(mrif%Yspec(nspec,spray%nzo)); mrif%Yspec = 0.0_WP
         allocate(mrif%bpdf(mduc%nx1))
         allocate(mrif%Temp(spray%nzo)); mrif%Temp = 0.0_WP

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

         next_mrif => mrif%next
         mrif => next_mrif

      end do

   end subroutine

   subroutine init_flamelets(spray)
     implicit none

     ! ---------------------------------
     type(spray_t), intent(inout) :: spray

     ! ---------------------------------
     type(srif_t), pointer :: mrif, next_mrif
     integer, pointer :: kmin, kmax, kmino, kmaxo
     integer :: k, ierr, i
     real(WP) :: max_mfr,doi,tot_inj_mass

     max_mfr = spray%rho_l*pi/4*spray%D_eff**2*spray%U_inj
     doi = spray%final_time

     k = 1
     tot_inj_mass = 0.0_WP

     do k=2,size(spray%roi,1)
        tot_inj_mass = tot_inj_mass + max_mfr &
                                    * 0.5_WP*(spray%roi(k+1,2)+spray%roi(k,2)) &
                                    * (spray%roi(k+1,1)-spray%roi(k,1))*spray%tau
     end do

     spray%inj_mass = tot_inj_mass

     ! Allocate scalars for flamelets
     allocate(spray%Sc(spray%nsc,spray%nzo)); spray%Sc = 0.0_WP
     allocate(spray%Wt(spray%nsc,spray%nzo)); spray%Wt = 0.0_WP

     allocate(spray%first_rif)
     mrif => spray%first_rif
     allocate(mrif%mduc)
     do i = 1,spray%nsc-1
        call list_insert(mrif)
        next_mrif => mrif%next
        mrif => next_mrif
     end do

   end subroutine init_flamelets

   subroutine list_insert( elem )
     type(srif_t), pointer  :: elem

     type(srif_t), pointer :: next

     allocate(next)
     allocate(next%mduc)

     next%next => elem%next
     elem%next => next
   end subroutine list_insert

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
      integer                   :: i,c, ispec, k, nsc
      character(len=str_medium), dimension(:), pointer :: species
      real(WP), dimension(:), pointer :: Yspec, bpdf
      character(len=512) :: fname
      integer, dimension(2) :: bound
      ! ----------------------------------

      if(spray%combustion_model == 'None' .or. spray%combustion_model == 'noname') return

      mduc => mduc_

      !if(spray%chi_st(1) == 0.0_WP) return

      Zmix_st = 1.0_WP/(1.0_WP+spray%stoic_coeff/spray%Y_O2)

      ! Warning!!! If you change anything in the following two lines,
      ! remember to change it also in getScalarDissipationRate spray_func.f90
      ! Zmix_max < 1.0 results in too high scalar dissipation rate and 
      ! it becomes hard to ignite. Therefore, Zmix_max = 1.0 is used here.
      Zmix_max = maxval(spray%Zmix_g) + 2.0_WP*sqrt(maxval(spray%zvar_g))
      Zmix_max = 1.0_WP !min(1.0,Zmix_max)

      chi = spray%C_chi*spray%chi_st(1)*(x/Zmix_st)**2*(log(x/Zmix_max)/log(Zmix_st/Zmix_max))

      where(chi/=chi) chi = 0.0_WP
      where(chi<0.0_WP) chi = 0.0_WP

      !write(*,*) 'x, chi '

      Pavg = spray%P_a
      Hmin = 0.0_WP
      Hmax = Zmix_max

      ! Advance flamelet
#ifdef SOOT
         call mduc_advance_flamelet_1d(mduc,spray%dt*spray%tau,Pavg,Hmin,Hmax,chi,0,.false.,fY,fT,fS)
#else
         call mduc_advance_flamelet_1d(mduc,spray%dt*spray%tau,Pavg,Hmin,Hmax,chi,0,.false.,fY,fT)
#endif

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

#ifdef SOOT
               write(100,'(<nspec+5>(A))') '#z<1> ','Temperature<2> ','SootNumDen<3> ','SootVolume<4> ','SootMass<5> ',(species(i), i=1,nspec)
#else
               write(100,'(<nspec+2>(A))') '#z<1> ','Temperature<2> ',(species(i), i=1,nspec)
#endif

               do k = spray%kmin,spray%kmax

                  call betaPDF(spray%Zmix_g(k),spray%Zvar_g(k),mduc_%nx1,x,bpdf,bound)

                  ! Convolute
                  do ispec = 1,nspec

                     Yspec(ispec) = sum(fY(bound(1):bound(2),ispec)*bpdf(bound(1):bound(2)))

                  end do

#ifdef SOOT
                  gasDen = spray%rho(k)*spray%rho_l*spray%Y_g(k)
                  cellVol = spray%dz*Pi*spray%b(k)**2*spray%D_eff**3
                  sootM00 = sum(fS(bound(1):bound(2),1)*bpdf(bound(1):bound(2)))
                  SootNumDen = sootM00*gasDen
                  sootM10 = sum(fS(bound(1):bound(2),2)*bpdf(bound(1):bound(2)))
                  SootVolume = sootM10*gasDen*24.0_WP/1800.0_WP*cellVol
                  SootMass = 1800.0_WP*SootVolume
#endif

                  Temp = sum(fT(bound(1):bound(2))*bpdf(bound(1):bound(2)))

#ifdef SOOT
                  write(100, '(<nspec+5>(ES15.5E3))') spray%z(k),Temp,SootNumDen,SootVolume,SootMass,(Yspec(i), i=1,nspec)
#else
                  write(100, '(<nspec+2>(ES15.5E3))') spray%z(k), Temp, (Yspec(i), i=1,nspec)
#endif

               end do
               close(unit=100)

               if(spray%combustion_model == 'MRIF') then
                  write(fname,"(A,A1,I0.6)") trim('flamelets.out'), '_', spray%step

                  open(unit=200,file=trim(fname),form="formatted",status="replace",action="write")

                  nsc = spray%nsc
                  !write(200,'(<nsc+1>(A))') '#z<1> ',(i, i=2,nsc)

                  do k = spray%kmin,spray%kmax
                     write(200, '(<nsc+1>(ES15.5E3))') spray%z(k), (spray%Sc(i,k), i=1,nsc)
                  end do
                  close(unit=200)

                  write(fname,"(A,A1,I0.6)") trim('weights.out'), '_', spray%step

                  open(unit=300,file=trim(fname),form="formatted",status="replace",action="write")

                  nsc = spray%nsc
                  !write(300,'(<nsc+1>(A))') '#z<1> ',(i, i=2,nsc)

                  do k = spray%kmin,spray%kmax
                     write(300, '(<nsc+1>(ES15.5E3))') spray%z(k), (spray%Wt(i,k), i=1,nsc)
                  end do
                  close(unit=300)
               end if

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

    end subroutine spray_combust_poststep

   subroutine spray_combust_poststep_mrif(spray)
      implicit none

      type(spray_t), intent(inout) :: spray
   
      ! ----------------------------------
      type(srif_t), pointer     :: mrif, next_mrif
      type(mduc_t), pointer     :: mduc
      real(WP)                  :: Pavg,Hmin,Hmax,Nu
      real(WP), dimension(:,:), pointer :: fY,fS
      real(WP), dimension(:), pointer :: fT, chi
      real(WP)                  :: Zmix_st, Zmix_max
      real(WP)                  :: Temp, sootM00, SootNumDen, &
                                   sootM10, SootVolume, &
                                   SootMass, gasDen, cellVol
      integer                   :: i,c, ispec, k, nsc, nflm
      character(len=str_medium), dimension(:), pointer :: species
      real(WP), dimension(:), pointer :: Yspec, bpdf
      character(len=512) :: fname
      integer, dimension(2) :: bound
      real(WP) :: curr_vap_mass, factor
      ! ----------------------------------

      if(spray%combustion_model == 'None' .or. spray%combustion_model == 'noname') return

      !if(spray%chi_st(1) == 0.0_WP) return

      Zmix_st = 1.0_WP/(1.0_WP+spray%stoic_coeff/spray%Y_O2)

      ! Warning!!! If you change anything in the following two lines,
      ! remember to change it also in getScalarDissipationRate spray_func.f90
      ! Zmix_max < 1.0 results in too high scalar dissipation rate and 
      ! it becomes hard to ignite. Therefore, Zmix_max = 1.0 is used here.
      Zmix_max = maxval(spray%Zmix_g) + 2.0_WP*sqrt(maxval(spray%zvar_g))
      Zmix_max = 1.0_WP !min(1.0,Zmix_max)

      Pavg = spray%P_a
      Hmin = 0.0_WP
      Hmax = Zmix_max

      curr_vap_mass = sum(spray%rho*spray%rho_l*spray%Y_v*spray%dz*Pi*spray%b**2*spray%D_eff**3)
      factor = int(curr_vap_mass/(spray%inj_mass/spray%nsc))

      ! Advance flamelet
      nflm = 1
      nullify(mrif)
      mrif=>spray%first_rif
      do while (associated(mrif))
         mduc=>mrif%mduc
         !mduc => mduc_
         if(factor<nflm-1) then
            nflm = nflm + 1
            next_mrif => mrif%next
            mrif => next_mrif
            cycle
         end if
         x => mrif%zz
         fY => mrif%fY_1d
         fT => mrif%fT_1d
         chi => mrif%chi

         chi = spray%C_chi*spray%chi_st_m(nflm)*(x/Zmix_st)**2*(log(x/Zmix_max)/log(Zmix_st/Zmix_max))
         
         where(chi/=chi) chi = 0.0_WP
         where(chi<0.0_WP) chi = 0.0_WP
         
         !write(*,*) 'x, chi '

#ifdef SOOT
         fS => mrif%fS_1d
         call mduc_advance_flamelet_1d(mduc,spray%dt*spray%tau,Pavg,Hmin,Hmax,chi,0,.false.,fY,fT,fS)
#else
         call mduc_advance_flamelet_1d(mduc,spray%dt*spray%tau,Pavg,Hmin,Hmax,chi,0,.false.,fY,fT)
#endif
         if(spray%irank .eq. 0) then

            if (mod(spray%step,spray%datafreq) == 0 .or. spray%end) then

               do k = spray%kmin,spray%kmax

                  call betaPDF(spray%Zmix_g(k),spray%Zvar_g(k),mduc%nx1,x,mrif%bpdf,bound)

                  ! Convolute
                  do ispec = 1,nspec

                     mrif%Yspec(ispec,k) = spray%Wt(nflm,k)*sum(fY(bound(1):bound(2),ispec)*mrif%bpdf(bound(1):bound(2)))

                  end do

#ifdef SOOT
                  mrif%sootM00(k) = spray%Wt(nflm,k)*sum(fS(bound(1):bound(2),1)*mrif%bpdf(bound(1):bound(2)))
                  mrif%sootM10(k) = spray%Wt(nflm,k)*sum(fS(bound(1):bound(2),2)*mrif%bpdf(bound(1):bound(2)))
#endif

                  mrif%Temp(k) = spray%Wt(nflm,k)*sum(fT(bound(1):bound(2))*mrif%bpdf(bound(1):bound(2)))

               end do

            end if

         end if

         nflm = nflm + 1
         next_mrif => mrif%next
         mrif => next_mrif
      end do

      if(spray%irank .eq. 0) then

         if (mod(spray%step,spray%datafreq) == 0 .or. spray%end) then

            write(fname,"(A,A1,I0.6)") trim('species.out'), '_', spray%step

            open(unit=100,file=trim(fname),form="formatted",status="replace",action="write")

#ifdef SOOT
            write(100,'(<nspec+5>(A))') '#z<1> ','Temperature<2> ','SootNumDen<3> ','SootVolume<4> ','SootMass<5> ',(species(i), i=1,nspec)
#else
            write(100,'(<nspec+2>(A))') '#z<1> ','Temperature<2> ',(spray%first_rif%species(i), i=1,nspec)
#endif

            do k = spray%kmin,spray%kmax
               nullify(mrif)
               mrif=>spray%first_rif
               ! get number of species in flamelet
               nspec = mduc_get_nnons_species(mrif%mduc)
               allocate(Yspec(nspec)); Yspec = 0.0_WP
               Temp = 0.0_WP
               do while (associated(mrif))
                  mduc=>mrif%mduc
                  !mduc => mduc_
                  ! Sum
                  do ispec = 1,nspec
                     Yspec(ispec) = Yspec(ispec) + mrif%Yspec(ispec,k)
                  end do

#ifdef SOOT
                  gasDen = spray%rho(k)*spray%rho_l*spray%Y_g(k)
                  cellVol = spray%dz*Pi*spray%b(k)**2*spray%D_eff**3
                  SootNumDen = SootNumDen + mrif%sootM00(k)*gasDen
                  SootVolume = SootVolume + mrif%sootM10(k)*gasDen*24.0_WP/1800.0_WP*cellVol
                  SootMass = 1800.0_WP*SootVolume
#endif

                  Temp = Temp + mrif%Temp(k)
                  next_mrif => mrif%next
                  mrif => next_mrif
               end do
#ifdef SOOT
               write(100, '(<nspec+5>(ES15.5E3))') spray%z(k),Temp,SootNumDen,SootVolume,SootMass,(Yspec(i), i=1,nspec)
#else
               write(100, '(<nspec+2>(ES15.5E3))') spray%z(k), Temp, (Yspec(i), i=1,nspec)
#endif

            end do
            close(unit=100)

            if(spray%combustion_model == 'MRIF') then
               write(fname,"(A,A1,I0.6)") trim('flamelets.out'), '_', spray%step

               open(unit=200,file=trim(fname),form="formatted",status="replace",action="write")

               nsc = spray%nsc
               !write(200,'(<nsc+1>(A))') '#z<1> ',(i, i=2,nsc)

               do k = spray%kmin,spray%kmax
                  write(200, '(<nsc+1>(ES15.5E3))') spray%z(k), (spray%Sc(i,k), i=1,nsc)
               end do
               close(unit=200)

               write(fname,"(A,A1,I0.6)") trim('weights.out'), '_', spray%step

               open(unit=300,file=trim(fname),form="formatted",status="replace",action="write")

               nsc = spray%nsc
               !write(300,'(<nsc+1>(A))') '#z<1> ',(i, i=2,nsc)

               do k = spray%kmin,spray%kmax
                  write(300, '(<nsc+1>(ES15.5E3))') spray%z(k), (spray%Wt(i,k), i=1,nsc)
               end do
               close(unit=300)
            end if

            !deallocate(species,Yspec,bpdf)
         end if
      end if

      return

    end subroutine spray_combust_poststep_mrif

end module 
