; =========================================================================== ;
; File to set up DPM model in Fluent.      ;
;                                                                             ;
; =========================================================================== ;
;
; ----- load the desired material and turn on species transport
(define udf-mat-dbase "/scratch2/doe1pal/ACCESS/M271/common/material_surrogate_1.scm")
(ti-menu-load-string
 (format #f "\ndefine/materials/data-base database-type user-defined \"~a\""
         udf-mat-dbase))
(ti-menu-load-string "\ndefine/models/species/species-transport y surrogate-1")
;
; ----- turn on the spray model
(ti-menu-load-string "\ndefine/models/dpm/interaction couple-calculations y")
(ti-menu-load-string "\ndefine/models/dpm/interaction update-dpm-sources y")
(ti-menu-load-string "\ndefine/models/dpm/interaction no.-of-cont 200")
(ti-menu-load-string "\ndefine/models/dpm/unsteady-tracking y y")
;
; =========================================================================== ;

; ----- change the initial times
(rpsetvar 'flow-time 0.0)
(rpsetvar 'dynamesh/in-cyn/crank-angle 420.0)
(rpsetvar 'time-step 0)

(ti-menu-load-string "\ndefine/dynamic-mesh/dynamic-mesh y y 1500 420 720")
;
; set solver bounds
(ti-menu-load-string "\nsolve/set/set-all-species-together")
(ti-menu-load-string "\nsolve/set/limits 800 1000e5 60 5000 1e-14 1e-22 1e5")
