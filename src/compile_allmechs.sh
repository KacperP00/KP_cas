#!/bin/bash 

sed 's/MECH/ndodecane_Cai/g' Makefile.template Makefile
make veryclean;make;make allopt;cd ../bin/;mv cas_opt carts_opt_ndo; cd ../src/

sed 's/MECH/ome0to5/g' Makefile.template > Makefile
make veryclean;make;make allopt;cd ../bin/;mv cas_opt carts_opt_ome0to5; cd ../src/

sed 's/MECH/dinbutylether/g' Makefile.template > Makefile
make veryclean;make;make allopt;cd ../bin/;mv cas_opt carts_opt_dinbutylether; cd ../src/

sed 's/MECH/1octanol/g' Makefile.template > Makefile
make veryclean;make;make allopt;cd ../bin/;mv cas_opt carts_opt_1octanol; cd ../src/

sed 's/MECH/diethoxymethane/g' Makefile.template > Makefile
make veryclean;make;make allopt;cd ../bin/;mv cas_opt carts_opt_diethoxymethane; cd ../src/

sed 's/MECH/1_3_dioxolane/g' Makefile.template > Makefile
make veryclean;make;make allopt;cd ../bin/;mv cas_opt carts_opt_1_3_dioxolane; cd ../src/
