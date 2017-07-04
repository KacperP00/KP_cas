ls pc*|tr '\t' '\n'|cut -d'.' -f1>fuels.list
echo"">fuels.f90
for name in $(cat fuels.list)
do
echo $name
sed "/function/c\subroutine $name(pc)" $name".m" > tmp.f90
sed -i "/subroutine/a\    implicit none\n\n    ! ---------------------------------------------\n    type(pc_t), pointer, intent(inout) :: pc\n    ! ---------------------------------------------" tmp.f90
sed -i "s/end/end subroutine $name/g" tmp.f90
sed -i "/%%/c\ " tmp.f90
sed -i "s/%/!/g" tmp.f90
sed -i "s/out\./pc%/g" tmp.f90
sed -i "s/.A /%A /g" tmp.f90
sed -i "s/.B /%B /g" tmp.f90
sed -i "s/.C /%C /g" tmp.f90
sed -i "s/.D /%D /g" tmp.f90
sed -i "s/.E /%E /g" tmp.f90
sed -i "s/.H /%H /g" tmp.f90
sed -i "s/.O /%O /g" tmp.f90
sed -i "s/.range/%range/g" tmp.f90
sed -i "s/.eqn/%eqn/g" tmp.f90
sed -i "s/\[/(\/ /g" tmp.f90
sed -i "s/\]/ \/)/g" tmp.f90
sed -i "s/NaN/MAX_REAL_WP/g" tmp.f90
sed -i "s/\[/(\/ /g" tmp.f90
sed -i "s/range/range(:,1)/g" tmp.f90
sed -i "s/ChemicalFormular/ChemicalFormula/g" tmp.f90
sed -i "s/pc%Names/pc%name/g" tmp.f90
grep -a "pc%name"
cat fuels.f90 tmp.f90 > tmpp.f90
mv tmpp.f90 fuels.f90
done

rm fuelselect.f90
echo "subroutine choosePC(pc,fuel)
  implicit none

  ! ---------------------------------------------
  type(pc_t), pointer, intent(inout) :: pc
  character(len=128), intent(in) :: fuel
  ! ---------------------------------------------

  select case(fuel)" > fuelselect.f90

for name in $(cat fuels.list)
do
echo "   case($(grep -a out.Names $name.m|cut -d= -f2|cut -d'{' -f2|cut -d'}' -f1))
	 
	 call $name(pc)
" >> fuelselect.f90
done
echo "   case default

	 write(*,*) 'Error:: Fuel not available in the database'

  end select

end subroutine choosePC">> fuelselect.f90

awk '/range/{ print; print; }' pc_ndodecane.m
