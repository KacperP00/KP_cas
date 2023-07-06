; =========================================================================== ;
; File to preview the mesh motion for a given case                            ;
;   Called from: case.scm                                                     ;
; =========================================================================== ;
;
(define mp-case "M271_CR13_sym_coarse")
(define common-dir "/scratch2/doe1pal/ACCESS/M271/common")
(define mp-start-cad 0)
(define mp-end-cad   180)
(define mp-case-freq 10)
(define mp-step-size 0.25)
(define mp-time-steps (/ (- mp-end-cad mp-start-cad) mp-step-size))
;
(define mp-events  (format #f "~a/M271_events_no_valves.scm" common-dir))
(define mp-invalve (format #f "~a/intake_zero_lift.ttab" common-dir))
(define mp-exvalve (format #f "~a/exhaust_zero_lift.ttab" common-dir))
;
;(define mp-events  "M271_TM02_083_events.scm")
;(define mp-invalve "M271_TM02_083_intake_lift.ttab")
;(define mp-exvalve "M271_TM02_083_exhaust_lift.ttab")
;
; ----- read in the case file
(ti-read-case (format #f "~a_~03dCAD" mp-case mp-start-cad))
(ti-menu-load-string "views/camera/projection orthographic")
;
; ----- generate the planes/geometry for this setup
(ti-read-scheme (format #f "~a/M271_create_planes.scm" common-dir))
;
; ----- scripts for creating view
(ti-menu-load-string "view/restore-view right")
(ti-menu-load-string "view/camera/roll -90")
(ti-menu-load-string "view/camera/target 0 0 -0.04")
(ti-menu-load-string "view/camera/field 0.22 0.1")
(ti-menu-load-string "view/save-view a_mesh_motion")
(ti-menu-load-string "\nview/restore-view a_mesh_motion")
;(ti-menu-load-string "\nview/restore-view a_valve_motion")
(ti-menu-load-string "\ndisplay/surface-grid '(a_intake)")
;
;(ti-menu-load-string "view/restore-view right")
;(ti-menu-load-string "view/camera/roll -90")
;(ti-menu-load-string "view/camera/target 0 0 -0.024")
;(ti-menu-load-string "view/camera/field 0.22 0.165")
;(ti-menu-load-string "view/save-view a_mesh_motion")
;
; ----- if generating blank start files, use zero lift profiles
(ti-menu-load-string (format #f "file/read-transient-table ~a" mp-invalve))
(ti-menu-load-string (format #f "file/read-transient-table ~a" mp-exvalve))
;
; ----- read the event file without any valve motion
(ti-menu-load-string (format #f "define/dynamic-mesh/events import \"~a\""
                             mp-events))
;
; ----- set up pictures
(ti-menu-load-string "display/set/picture/invert y")
(ti-menu-load-string "display/set/picture/driver png")
(ti-menu-load-string "display/set/picture/color-mode color")
(ti-menu-load-string "display/set/picture/x-resolution 960")
(ti-menu-load-string "display/set/picture/y-resolution 720")
;
; ----- loop to write a case file at a specified crank angle interval
(define mp-stop 0)
(define mp-cad (rpgetvar 'dynamesh/in-cyn/crank-angle))
(define mp-cad-write (+ mp-cad mp-case-freq))
(do ((i 1 (+ i 1)))((= mp-stop 1))
  (do ((j 1 (+ j 1)))((>= mp-cad mp-cad-write))
    (ti-menu-load-string "solve/mesh-motion 1 y 1 y n y")
    (set! mp-cad (rpgetvar 'dynamesh/in-cyn/crank-angle))
    )
  (ti-write-case (format #f "~a_~03dCAD" mp-case mp-cad))
  (set! mp-cad-write (+ mp-cad-write mp-case-freq))
  (if (= mp-cad mp-end-cad) (set! mp-stop 1))
); end mp-stop
;(do ((i 1 (+ i 1)))((= i (/ mp-time-steps mp-case-freq)))
;  (ti-menu-load-string
;   (format #f "solve/mesh-motion ~d y 1 y n y" mp-case-freq ))
;  (ti-write-case (format #f "~a_~03.0dCAD" mp-case
;                         (rpgetvar 'dynamesh/in-cyn/crank-angle)))
;); end i
;(ti-menu-load-string "solve/mesh-motion 2880 y 1 y n y")
;
; ===== end mesh-motion.scm =====
