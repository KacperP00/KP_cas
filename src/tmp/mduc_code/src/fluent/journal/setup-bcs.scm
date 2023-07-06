; ----- set a boundary condition
;(ti"\ndefine/boundary-conditions/wall cyl-head 0.008 n 0 n y temperature n 373.15 n n n n 0 n 0.5"


(ti-menu-load-string "\ndefine/user-defined/compiled load /scratch2/doe1pal/ACCESS/M271/libudf/mduc")
(ti-menu-load-string "\ndefine/user-defined/compiled load /scratch2/doe1pal/ACCESS/M271/libudf/spray")
;
; ----- defined the UDS/UDM variables
(ti-menu-load-string "define/user-defined/user-defined-scalars 6" )
(ti-menu-load-string "define/user-defined/user-defined-memory 14" )
;
; set solution discretization 
(ti-menu-load-string "\nsolve/set/discretization-scheme uds-0 1")
(ti-menu-load-string "\nsolve/set/discretization-scheme uds-1 1")
(ti-menu-load-string "\nsolve/set/discretization-scheme uds-2 1")
(ti-menu-load-string "\nsolve/set/discretization-scheme uds-3 1")
(ti-menu-load-string "\nsolve/set/discretization-scheme uds-4 1")
(ti-menu-load-string "\nsolve/set/discretization-scheme uds-5 1")
;
; ----- define the diffusivity for the scalars
(ti-menu-load-string "\ndefine/materials/change-create surrogate-1 surrogate-1
                   n n n n n n y defined-per-uds
                   0 user-defined \"Diff_Z1::mduc\"
                   1 user-defined \"Diff_Z1::mduc\"
                   2 user-defined \"Diff_Z2::mduc\"
                   3 user-defined \"Diff_Z2::mduc\"
                   4 user-defined \"Diff_Z::mduc\"
                   5 user-defined \"Diffusivity::mduc\" n ")
;

; ----- define the source terms in the fluid zones 
; setup fluid boundary conditions
(define (set-mduc-srcs)
  (for-each
   (lambda (name)
     (display (format #f " setting MDUC scalar sources for zone: ~a\n" name) )
     (ti-menu-load-string (format #f "\ndefine/boundary-conditions/fluid ~a
       mixture y surrogate-1 y ~a ~a ~a
       1 n y \"Z1_src::mduc\"
       1 n y \"Z1var_src::mduc\"
       1 n y \"Z2_src::mduc\"
       1 n y \"Z2var_src::mduc\"
       1 n y \"Zvar_src::mduc\"
       1 n y \"Htot_src::mduc\"" name (dummy-pass 6)
       (dummy-pass (- n-species 1)) (dummy-pass 1) ))
     ); end lambda
   (get-avail-zones)
   ); end for-each
  ); end defined
;
; execute the command
(set-mduc-srcs)
;
; ----- define source terms for species mass fraction from Wiebe function
(define (set-mduc-wiebe-srcs)
  (for-each
   (lambda (name)
     (display (format #f " setting Wiebe scalar sources for zone: ~a\n" name) )
     (ti-menu-load-string (format #f "\ndefine/boundary-conditions/fluid ~a
       mixture y surrogate-1 y \n \n \n \n \n \n 
       1 n y \"Wiebe_Y_fuel\"
       \n \n \n \n
       1 n y \"Wiebe_Y_O2::mduc\"
       1 n y \"Wiebe_Y_H2O::mduc\"
       1 n y \"Wiebe_Y_CO2::mduc\"
       \n \n \n \n \n \n
       1 n y \"Wiebe_Enthalpy::mduc\"
       \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n" name) )
     ); end lambda
   (get-avail-zones)
   ); end for-each
  ); end defined
;
; execute the command
;(set-mduc-wiebe-srcs)
;
;
(define (set-bc-wall-temp walls type val)
  (for-each 
   (lambda (name)
     (display (format #f " setting boundary condition for wall: ~a\n" name) )
     (ti-menu-load-string 
      (format #f "\ndefine/boundary-conditions/wall ~a ~a
        y ~a n ~d" name (dummy-pass 4) type val ))
     ); end lambda
   walls
   ); end for-each
  ); end defined
; (set-bc-wall-temp '(wall-99) "heat-flux" 10) 
;
(define (set-bc-wall-dpm walls type params)
  (for-each 
   (lambda (name)
     (display (format #f " setting wall film on wall: ~a\n" name) )
     (ti-menu-load-string 
      (format #f "\ndefine/boundary-conditions/wall ~a ~a
        y ~a ~a" name (dummy-pass 14) type params))
     ); end lambda
   walls
   ); end for-each
  ); end defined
; (set-bc-wall-dpm '(wall-99) "wall-film" "4")
;
; ----- set the boundary conditions
(define bc-wall-cyl      '(cyl-tri wall-99 wall-103))
(define bc-wall-cyl-head '(cyl-head))
(define bc-wall-piston   '(piston))
(define bc-wall-in-port  '(invalve-1-port))
(define bc-wall-ex-port  '(exvalve-1-port))
(define bc-wall-in-valve '(invalve1-ch invalve1-ib invalve1-ob invalve1-seat invalve1-stem))
(define bc-wall-ex-valve '(exvalve1-ch exvalve1-ib exvalve1-ob exvalve1-seat exvalve1-stem))
(define bc-T-cyl-wall 400)
(define bc-T-cyl-head 670)
(define bc-T-piston   450.0)
(define bc-T-in-port  300.0)
(define bc-T-ex-port  415.0)
(define bc-T-in-valve 400.0)
(define bc-T-ex-valve 450.0)

(set-bc-wall-temp bc-wall-cyl      "temperature" bc-T-cyl-wall)
(set-bc-wall-temp bc-wall-head     "temperature" bc-T-cyl-head)
(set-bc-wall-temp bc-wall-piston   "temperature" bc-T-piston)
(set-bc-wall-temp bc-wall-in-port  "temperature" bc-T-in-port)
(set-bc-wall-temp bc-wall-ex-port  "temperature" bc-T-ex-port)
(set-bc-wall-temp bc-wall-in-valve "temperature" bc-T-in-valve)
(set-bc-wall-temp bc-wall-ex-valve "temperature" bc-T-ex-valve)



