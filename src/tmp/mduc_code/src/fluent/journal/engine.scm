;; ========================================================================== ;;
;; File to setup engine settings for a case                                   ;;
;;   Called from: run.scm                                                     ;;
;; ========================================================================== ;;
;;
;; ***** define functions
(define get-engine-zones
  (lambda ()
    (let ( (cur-fluids '())
           (fluid-tmp "")
           )
      (set! cur-fluids (map thread-name (sort-threads-by-name
                                         (get-threads-of-type 'fluid))))
      (do ((i 0 (+ i 1))) ((= i (length cur-fluids )))
        (set! fluid-tmp (symbol->string (list-ref cur-fluids i)))
        (cond ( (string=? "fluid-ch" (substring fluid-tmp 0 8))
                (set! eng-cyl (append eng-cyl (list (list-ref cur-fluids i))))
                )
              ( (string=? "fluid-in" (substring fluid-tmp 0 8 ))
                (set! eng-int (append eng-int (list (list-ref cur-fluids i))))
                )
              ( (string=? "fluid-ex" (substring fluid-tmp 0 8 ))
                (set! eng-exh (append eng-exh (list (list-ref cur-fluids i))))
                )
              (else (display (format #f "~a not an engine zone" fluid-tmp) ))
              ) ;; end cond
        ) ;; end do
      ) ;; end let
    )
  ) ;; end get-engine-zones
;; ***** end function definition
;;
(define mduc-reset-engine-egr
  (lambda ()
    ;; make sure engine zones are up-to-date
    (ti-menu-load-string "(get-engine-zones)")
    ;; cylinder 
    (ti-menu-load-string "\nsolve/patch eng-cyl uds-2 1.0")
    (ti-menu-load-string "\nsolve/patch eng-cyl uds-3 0.0")
    ;; intake
    (ti-menu-load-string "\nsolve/patch eng-int uds-2 0.0")
    (ti-menu-load-string "\nsolve/patch eng-int uds-3 0.0")
    ;; exhaust
    (ti-menu-load-string "\nsolve/patch eng-exh uds-2 1.0")
    (ti-menu-load-string "\nsolve/patch eng-exh uds-3 0.0")
    ) ;; lambda
  ) ;; end get-engine-zones
;;
(define mduc-reset-engine-fuel
  (lambda ()
    ;; make sure engine zones are up-to-date
    (ti-menu-load-string "(get-engine-zones)")
    ;; cylinder 
    (ti-menu-load-string "\nsolve/patch eng-cyl uds-0  0.0")
    (ti-menu-load-string "\nsolve/patch eng-cyl uds-1  0.0")
    ;; intake
    (ti-menu-load-string "\nsolve/patch eng-int uds-0  0.0")
    (ti-menu-load-string "\nsolve/patch eng-int uds-1  0.0")
    ;; exhaust
    (ti-menu-load-string "\nsolve/patch eng-exh uds-0  0.0")
    (ti-menu-load-string "\nsolve/patch eng-exh uds-1  0.0")
    ) ;; lambda
  ) ;; end get-engine-zones
;;
;; ----- initialize the engine zones 
(define eng-cyl '())
(define eng-int '())
(define eng-exh '())
(ti-menu-load-string "(get-engine-zones)")
(display "Engine zones set to:\n")
(display (format #f " cylinder: ~a\n" eng-cyl))
(display (format #f " intake:   ~a\n" eng-int))
(display (format #f " exhaust   ~a\n" eng-exh))
;;
;; ----- read in the ICE specific functions necessary to run
(if (or (string=? (format #f "~a" (inquire-release)) "(14 0 0)") 
        (string=? (format #f "~a" (inquire-release)) "(14 5 0)") )
    (ti-read-scheme (format #f "~a/journal/WB-ICE-Solver-Setup.scm" case-root))
   )
;;
;; ----- read in boundary profiles
(ti-menu-load-string (format #f "/file/read-transient-table ~a/~a/init/~a"
                             case-root case-spec bc-manifold))
;; ----- read in valve profiles
(ti-menu-load-string (format #f "/file/read-transient-table ~a/~a/init/~a"
                             case-root case-spec invalve-lift))
(ti-menu-load-string (format #f "/file/read-transient-table ~a/~a/init/~a"
                             case-root case-spec exvalve-lift))
;; ----- read in events
(ti-menu-load-string
 (format #f "define/dynamic-mesh/events import \"~a/~a/init/~a\""
             case-root case-spec event-file))
;;
;; ----- set up in-cylinder settings
(ti-menu-load-string
 (format #f "\ndefine/dynamic-mesh/dynamic-mesh yes yes ~d ~d ~d" 
	op-speed ca-start ca-period) )
(ti-menu-load-string
 (format #f "\ndefine/dynamic-mesh/controls/in-cylinder-output y 1 ~a
   ~a/monitor/cylinder-swirl.~a cg 0 0 1 1 0 0 0 1 0"
         eng-cyl case-dir case-sub-id))
;;
;; ========================================================================== ;;