subroutine choosePC(pc,fuel)
  implicit none

  ! ---------------------------------------------
  type(pc_t), pointer, intent(inout) :: pc
  character(len=128), intent(in) :: fuel
  ! ---------------------------------------------

  select case(fuel)
   case('1butanol','1-butanol','1-Butanol','n-butanol','n-Butanol, Butan-1-ol')
         
         call pc_1butanol(pc)

   case('1decanol','1-decanol','1-Decanol','n-decanol','n-Decanol', 'decyl alcohol')
         
         call pc_1decanol(pc)

   case('1dodecanol','1-dodecanol','1-Dodecanol','n-dodecanol','Lauryl 24')
         
         call pc_1dodecanol(pc)

   case('1heptanol','1-heptanol','1-Heptanol','n-heptanol','n-Heptanol', 'heptyl alcohol')
         
         call pc_1heptanol(pc)

   case('1hexanol','1-hexanol','1-Hexanol','1Hexanol','n-Hexanol','n-hexanol',)
         
         call pc_1hexanol(pc)

   case('1nonanol','1-nonanol','1Nonanol','1-Nonanol','Alcohol C-9','n-Nonyl Alcohol')
         
         call pc_1nonanol(pc)

   case('1octanol','1-octanol','1-Octanol','n-octanol','n-Octanol, octanol')
         
         call pc_1octanol(pc)

   case('1pentanal','1-pentanal','1-Pentanal','1Pentanal','n-Pentanal','n-pentanal')
         
         call pc_1pentanal(pc)

   case('1pentanol','1-pentanol','1-Pentanol','1Pentanol','n-Pentanol','n-pentanol')
         
         call pc_1pentanol(pc)

   case('1propanol','1-propanol','1Ppropanol','n-propanol','n-Propanol', 'propanol-1')
         
         call pc_1propanol(pc)

   case('1tetradecanol','n-tetradecanol','n-Tetradecanol','1Tetradecanol','1-tetradecanol','1-Tetradecanol')
         
         call pc_1tetradecanol(pc)

   case('1tridecanol','n-tridecanol','n-Tridecanol','1-Tridecanol')
         
         call pc_1tridecanol(pc)

   case('1undecanol','1-undecanol','1-Undecanol','1Undecanol','n-Undecanol','n-undecanol')
         
         call pc_1undecanol(pc)

   case('2butanol','2-butanol','2-Butanol','2Butanol','Butane')
         
         call pc_2butanol(pc)

   case('2butanone','2-butanone','2-Butanone','2-butanone','butanone')
         
         call pc_2butanone(pc)

   case('2heptanone','2-heptanone','2Heptanone','2-Heptanone')
         
         call pc_2heptanone(pc)

   case('2hexanol','2-hexanol','2-Hexanol','2Hexanol')
         
         call pc_2hexanol(pc)

   case('2hexanone','2-hexanone','2-Hexanone','2Hexanone')
         
         call pc_2hexanone(pc)

   case('2methyl1butanol','2methyl-1butanol','amylalcohol')
         
         call pc_2methyl1butanol(pc)

   case('2methyl1propanol','2methyl-1propanol','Fermation Butyl Alcohol','Isobutyl Alcohol')
         
         call pc_2methyl1propanol(pc)

   case('2methyl2propanol','2methyl-2propanol','t-Butanol')
         
         call pc_2methyl2propanol(pc)

   case('2methylnonane','2-methylnonane')
         
         call pc_2methylnonane(pc)

   case('2methyloctane','2-methyloctane')
         
         call pc_2methyloctane(pc)

   case('2nonanone','2-nonanone','2-Nonanone')
         
         call pc_2nonanone(pc)

   case('2octanone','2-octanone','methyl hexyl ketone')
         
         call pc_2octanone(pc)

   case('2pentanol','2-pentanol','pentanol-2')
         
         call pc_2pentanol(pc)

   case('2pentanone','2-pentanone','Ethylacetone')
         
         call pc_2pentanone(pc)

   case('2propanol','2-propanol','dimethyl carbinol','isohol')
         
         call pc_2propanol(pc)

   case('3heptanone','3-heptanone','butyl ethyl ketone','heptane-3-one')
         
         call pc_3heptanone(pc)

   case('3hexanone','3-hexanone','ethyl-n-propylketone')
         
         call pc_3hexanone(pc)

   case('3pentanone','3-pentanone','metacetone','ethyl ketone')
         
         call pc_3pentanone(pc)

   case('3petanol','3-petanol','pentanol-3','diethyl carbinol')
         
         call pc_3petanol(pc)

   case('acetone','Acetone','Allylic Alcohol','Dimethylketal','3-Hydroxypropene','Propanone','beta-Ketopropane')
         
         call pc_acetone(pc)

   case('diethylether','aether','ether','ethyl ether')
         
         call pc_diethylether(pc)

   case('OME1','Methylal','ome1','methylal')
         
         call pc_dimethoxymethane(pc)

   case('dimethylether','methyl ether','ether','ethyl ether')
         
         call pc_dimethylether(pc)

   case('dinbutylether','di-n-butyl-ether','butyl oxide','1,1-oxybisbutane')
         
         call pc_dinbutylether(pc)

   case('dinpropylether','di-n-propyl-ether','n-propyl ether')
         
         call pc_dinpropylether(pc)

   case('ethane','bimethyl','ethyl hydride','methylmethane')
         
         call pc_ethane(pc)

   case('ethanol','Ethanol')
         
         call pc_ethanol(pc)

   case('gammavalerolactone','gamma-pentalactone','4-pentanolide')
         
         call pc_gammavalerolactone(pc)

   case('isobutane','iso-butane','i-butane')
         
         call pc_isobutane(pc)

   case('isoheptane','iso-heptane')
         
         call pc_isoheptane(pc)

   case('isohexane','iso-hexane','2-methylpentane')
         
         call pc_isohexane(pc)

   case('isooctane','iso-octane','heptane')
         
         call pc_isooctane(pc)

   case('isopentane','iso-pentane','butane','2-methylbutane')
         
         call pc_isopentane(pc)

   case('methane','Methane')
         
         call pc_methane(pc)

   case('methanol','Methanol','methyl alcohol')
         
         call pc_methanol(pc)

   case('methyltertbutylether','methyl tert-buty lether','propane')
         
         call pc_methyltertbutylether(pc)

   case('nbutane','n-butane','methylethylmethane')
         
         call pc_nbutane(pc)

   case('nbutylformate','n-butyl formate','butyl methanoate')
         
         call pc_nbutylformate(pc)

   case('ndecane','n-decane','decyl hydride')
         
         call pc_ndecane(pc)

   case('ndodecane','n-dodecane','duodecane','adakane 12')
         
         call pc_ndodecane(pc)

   case('nheptane','n-heptane','dipropymethane')
         
         call pc_nheptane(pc)

   case('nhexane','n-hexane','hexyl hydride')
         
         call pc_nhexane(pc)

   case('nnonane','n-nonane','nonyl hydride')
         
         call pc_nnonane(pc)

   case('noctane','n-octane','octyl hydride')
         
         call pc_noctane(pc)

   case('npentane','n-pentane','pental hydride')
         
         call pc_npentane(pc)

   case('ntetradecane','n-tetradecane','tetradecyl hydride')
         
         call pc_ntetradecane(pc)

   case('ntridecane','n-tridecane','tridecyl hydride')
         
         call pc_ntridecane(pc)

   case('propane','dimethylmethane','propyl hydride')
         
         call pc_propane(pc)

   case('tetrahydrofuran','hydrofuran','tetramethylene oxide')
         
         call pc_tetrahydrofuran(pc)

   case('tetrahydrofurfurylalcohol','tetrahydrofuryl carbinol','oxolan-2-methanol')
         
         call pc_tetrahydrofurfurylalcohol(pc)

   case default

         write(*,*) 'Error:: Fuel not available in the database'

  end select
