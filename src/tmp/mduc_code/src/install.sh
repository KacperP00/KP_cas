#!/bin/bash


MDUCSRC=($PWD)

cd ..
mkdir -p bin lib obj dep
cd $MDUCSRC

SUNDIALS_PATH="$MDUCSRC/../../sundials/sundials-3.2.1/install"
echo "SDIAL: $SUNDIALS_PATH"
# link mechanisms
#cd mechanism
## v3
#ln -s ../../../../run/01_ConvertMech/NXC12H26_822Spec_v3/NXC12H26_v3.c . 
#ln -s ../../../../run/01_ConvertMech/NXC12H26_822Spec_v3/NXC12H26_v3.h .
#
## v4
#ln -s ../../../../run/01_ConvertMech/NXC12H26_822Spec_v4/NXC12H26C.c NXC12H26.c
#ln -s ../../../../run/01_ConvertMech/NXC12H26_822Spec_v4/NXC12H26.h NXC12H26.h
#
#cd ..
#
## first compile for verification cases
#make clean
#rm -rf ../lib.verify 
#rm -rf ../bin.verify
#
#cp Makefile.in.verify Makefile.in
#sed -i "/SDIAL/c\SDIAL=$SUNDIALS_PATH" Makefile.in
#make mduc.debug 2>&1 | tee mduc.verify.log
#make mducPlot.debug 2>&1 | tee mducPlot.verify.log
#
#cp -r ../bin ../bin.verify
#cp -r ../lib ../lib.verify

# now compile for all other cases
make clean 
cp Makefile.in.production Makefile.in
sed -i "/SDIAL/c\SDIAL=$SUNDIALS_PATH" Makefile.in
make all.opt 2>&1 | tee mduc.all.log
