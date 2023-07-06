program fuel_properties
  use precision
  use pc_defs
  use pc_func

  implicit none

  ! ---------------------------------
  type(pc_t), dimension(:), pointer :: pc => null()
  type(pc_t), pointer :: pc_l=> null()
  real(WP), dimension(:), pointer :: T,Cp_l,p_vap,L_f,T_sat,sigma,rho_l,rho_v,visc_l,visc_v,lambda_l,lambda_v,Cp_v, G_v, p_vap_cc
  integer :: kmin, kmax
  integer :: k, nz
  character(len=512) :: headfmt_l, rowfmt_l, headfmt_v, rowfmt_v
  character(len=128) :: Fuel
  real(WP) :: R_gas, Pamb
  
  nz = 100
  R_gas = 8.3144598_WP
  Pamb = 50.2e5 !46.5E+05_WP
  Fuel = '1-octanol'

  allocate(T(nz))
  allocate(Cp_l(nz))
  allocate(p_vap(nz))
  allocate(p_vap_cc(nz))
  allocate(L_f(nz))
  allocate(T_sat(nz))
  allocate(sigma(nz))
  allocate(rho_l(nz))
  allocate(rho_v(nz))
  allocate(visc_l(nz))
  allocate(visc_v(nz))
  allocate(lambda_l(nz))
  allocate(lambda_v(nz))
  allocate(Cp_v(nz))
  allocate(G_v(nz))

  kmin = 1; kmax = nz

  do k=kmin,kmax
     !T(k) = 0.0_WP + (k-1)*10.0_WP
     T(k) = 20.0_WP + (k-1)*10.0_WP
  end do
  allocate(pc(nz))

  do k = kmin,kmax
     pc_l => pc(k)

     call allocate_pc(pc_l)

     pc(k) = pc_l

     deallocate(pc_l)

  end do

  do k=kmin,kmax
     pc_l => pc(k)

     call choosePC(pc_l,Fuel)

     ! Set temperature of liquid and vapor phaseson the grid
     pc_l%T = T(k)
     pc_l%p = Pamb

     ! Compute fuel properties: Liquid phase
     call computeLiqDensity(pc_l)
     call computeLiqViscosity(pc_l)
     call computeLiqHeatCapacity(pc_l)
     call computeVapPressure(pc_l)
     call computeHeatOfVap(pc_l)
     call computeSurfaceTension(pc_l)
     call computeLiqThermalConductivity(pc_l)

     call computeIG_vapDensity(pc_l)
     call computeVapViscosity(pc_l)
     call computeVapThermalConductivity(pc_l)
     call computeIG_HeatCapacity(pc_l)
     !call computeDiffusionCoefficientTheory(pc_l)
     !call computeDiffusionCoefficientFuller(pc_l)
     call computeDiffusionCoefficientWilkeLee(pc_l)

     ! Set fuel properties to spray
     Cp_l(k) = pc_l%liqHeatCapacity
     p_vap(k) = pc_l%vapPressure%val
     L_f(k) = pc_l%HeatOfVap
     if(L_f(k) == 0.0_WP) then
        T_sat(k) = pc_l%Tcrit
     else
        T_sat(k) = min(pc_l%Tcrit,1.0_WP/max(1E-16_WP,(1.0_WP/pc_l%NormalBoilingPoint - R_gas*log(pc_l%p/101325.0_WP)/pc_l%MolecularWeight*1000.0_WP/L_f(k))))
     end if

     !p_vap_cc(k) = 101325.0_WP*exp((L_f(k)*pc_l%MolecularWeight/1000.0_WP/R_gas)*(1.0_WP/pc_l%NormalBoilingPoint-1.0_WP/T(k))) 
     !p_vap_cc(k) = 101325.0_WP*exp((L_f(1)*pc_l%MolecularWeight/1000.0_WP/R_gas)*(1.0_WP/pc_l%NormalBoilingPoint-1.0_WP/T(k))) 
     !p_vap_cc(k) = pc_l%Pcrit*exp((L_f(k)*pc_l%MolecularWeight/1000.0_WP/R_gas)*(1.0_WP/pc_l%Tcrit-1.0_WP/T(k))) 
     !p_vap_cc(k) = pc_l%Pcrit*exp((L_f(1)*pc_l%MolecularWeight/1000.0_WP/R_gas)*(1.0_WP/pc_l%Tcrit-1.0_WP/T(k))) 
     !p_vap_cc(k) = pc_l%Pcrit*exp((0.5_WP*L_f(k)*pc_l%MolecularWeight/1000.0_WP/R_gas)*(1.0_WP/pc_l%Tcrit-1.0_WP/T(k))) 
     p_vap_cc(k) = 101325.0_WP*exp((0.5_WP*(2.60538E+005_WP+L_f(k))*pc_l%MolecularWeight/1000.0_WP/R_gas)*(1.0_WP/pc_l%NormalBoilingPoint-1.0_WP/T(k))) 

     ! Set fuel properties to spray (constant throughout spray simulation)
     sigma(k) = pc_l%SurfaceTension%val
     rho_l(k) = pc_l%liqDensity
     visc_l(k) = pc_l%liqViscosity%val
     lambda_l(k) = pc_l%liqThermalConductivity%val
     ! Set fuel vapor properties in spray object
     rho_v(k) = pc_l%IG_vapDensity
     visc_v(k) = pc_l%vapViscosity%val
     lambda_v(k) = pc_l%vapThermalConductivity%val
     Cp_v(k) = pc_l%IG_HeatCapacity
     G_v(k) = pc_l%DiffusionCoefficientWilkeLee

     nullify(pc_l)
  end do

  open(unit=100,file='LiquidFuelProperties.dat',form="formatted",status="unknown",action="write")
  headfmt_l = "(A20,A20,A20,A20,A20,A20,A20,A20,A20,A20)"
  rowfmt_l = "(ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3)"
  open(unit=101,file='VaporFuelProperties.dat',form="formatted",status="unknown",action="write")
  headfmt_v = "(A20,A20,A20,A20,A20,A20)"
  rowfmt_v = "(ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3, ES15.5E3)"

  write(100,headfmt_l) "Temperature<1> ","Density<2> ","Viscosity<3> ","ThermCond<4> ","SurfaceTension<5> ","Cp<6> ","HeatOfVap<7> ","Pvap<8> ","Tsat<9> ","Pvap_cc<10>"
  write(101,headfmt_v) "Temperature<1> ","Density<2> ","Viscosity<3> ","ThermCond<4> ","Cp<5> ","DiffCoeff<6>"

  do k = kmin,kmax
     write(100,FMT=rowfmt_l) T(k),rho_l(k),visc_l(k),lambda_l(k),sigma(k),Cp_l(k),L_f(k),p_vap(k),T_sat(k),p_vap_cc(k)
     write(101,FMT=rowfmt_v) T(k),rho_v(k),visc_v(k),lambda_v(k),Cp_v(k),G_v(k)
  end do

  close(unit=100)
  close(unit=101)

  deallocate(pc,T,Cp_l,p_vap,p_vap_cc,L_f,T_sat,sigma,rho_l,rho_v,visc_l,visc_v,lambda_l,lambda_v,Cp_v,G_v)
  
end program fuel_properties
