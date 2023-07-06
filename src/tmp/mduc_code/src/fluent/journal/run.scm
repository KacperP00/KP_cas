;; ========================================================================== ;;
;; File to run the Fluent case.                                               ;;
;;   Called from: case.scm                                                    ;;
;; ========================================================================== ;;
;;
;; ----- run time file names -----
(define case-dir    (format #f "~a/~a/rd_~a" case-root case-spec case-id))
(display (format #f "Case directory: ~a\n" case-dir))
;;
(define file-prefix "")
(if (> (string-length case-name) 0)
    (set! file-prefix case-name) )
(if (and (> (string-length case-name) 0) (> (string-length case-map) 0))
    (set! file-prefix (string-append file-prefix  "-" )) )
(if (> (string-length case-map) 0) 
    (set! file-prefix (string-append file-prefix (format #f "~a_" case-map))) )
(if (> (string-length case-spec) 0)
    (set! file-prefix (string-append file-prefix (format #f "~a" case-spec))) )
(if (> (string-length case-id) 0)
    (set! file-prefix (string-append file-prefix (format #f "-~a" case-id))) )
(display (format #f "File prefix:    ~a\n" file-prefix))
;;
;; load some matrix stuff
(ti-read-scheme (format #f "~a/journal/matrix.scm" case-root ))
;;
(define transcript-file (format #f "~a/~a.out" case-dir file-prefix))
(define udf-at-end "")
;;
;; ----- set up directory structure -----
(ti-read-scheme (format #f "~a/journal/directory.scm" case-root ))
;;
;; ----- start the transcript file
(if (file-exists? transcript-file)
    (system (format #f "mv ~a ~a.pre" transcript-file transcript-file )) )
(ti-menu-load-string (format #f "file/start-transcript ~a/~a.out"
                             case-dir file-prefix ))
;;
;; ----- read in the custom functions used by MDUC library
(ti-read-scheme (format #f "~a/journal/custom-funcs.scm" case-root ))
;;
;; ----- read initial case and data file -----
(if case-restart
    (if (> (string-length restart-file) 0)
        (define init-file (format #f "~a/~a/~a" case-dir fl-data restart-file))
        (define init-file (format #f "~a/~a/~a_~a" case-dir fl-data file-prefix
                                  (mduc-file-suffix fl-freq-type nt-start)))
        )
    (define init-file (format #f "~a/~a/init/~a" case-root case-spec
                              init-name))
)
(ti-read-case-data (format #f "~a" init-file) )
;;
;; check manifold conditions
(if manifold-case
    (manifold-check-restart) )
;; 
;; ----- setup data export -----
(if (> fl-freq 0) 
    (begin ;; data will be exported
      (if (string=? fl-freq-type "fluent")
          (begin ;; fluent data export)
            (ti-menu-load-string
             "file/auto-save/case-frequency if-case-is-modified")
            (ti-menu-load-string (format #f "file/auto-save/data-frequency ~d"
                                         fl-freq))
            (ti-menu-load-string
             (format #f " file/auto-save/root-name \"~a/~a/~a\""
                     case-dir fl-data file-prefix))
            )
          (begin ;; custom export
            (ti-menu-load-string "file/auto-save/data-frequency 0")
            (ti-read-scheme (format #f "~a/journal/data-export.scm" case-root ))
            ) ;; end custom export
          )
      )
    (display "Run: no data files will be saved")
    ) ;; end fl-freq
;;
;; ----- force an update on the solver thread names for UDF libraries -----
(ti-menu-load-string "(update-solver-thread-names)")
;;
;; ----- set up engine parameters -----
(if engine-case
    (ti-read-scheme (format #f "~a/journal/engine.scm" case-root )) )
;;
;; ----- set up spray parameters -----
(if inj-on
    (ti-read-scheme (format #f "~a/journal/spray.scm" case-root )) )
;;
;; ----- set up combustion module -----
(if mduc-on
    (ti-read-scheme (format #f "~a/journal/mduc.scm" case-root )) )
;;
;; ----- set up monitors -----
(ti-read-scheme (format #f "~a/journal/monitor.scm" case-root ))
;;
;; ----- set up data export (if specified) -----
(if exp-sol
    (ti-read-scheme (format #f "~a/journal/export.scm" case-root )) )
(if exp-udf
    (ti-read-scheme (format #f "~a/journal/export-udf.scm" case-root )) )
;;
;; ----- hook the necessary functions -----
(ti-menu-load-string
 (format #f "define/user-defined/function-hooks/execute-at-end/ ~a" udf-at-end))
;;
(if (> plot-freq 0)
    (ti-read-scheme (format #f "~a/journal/plots.scm" case-root ))
)
;; ----- read in any custom functions for the specific case
(if (> (string-length case-custom) 0)
    (ti-read-scheme (format #f "~a/~a/init/~a" case-root case-spec case-custom))
    )
;;
;; ----- have Fluent check the run directory for run-time files
(ckpt/time-step? #t)
(set! checkpoint/check-filename (format #f "~a/check-fluent" case-dir ))
(set! checkpoint/exit-filename  (format #f "~a/exit-fluent"  case-dir ))
(rpsetvar 'checkpoint/filename (format #f  "~a/fluent/restart-file_" case-dir ))
;;
;; ----- write the initilization case before starting the solver
;(if (not case-restart)
;    (ti-menu-load-string (format #f "wcd ~a/~a/~a-init ok" case-dir fl-data
;;                                 file-prefix)) )
;;
;; ----- turn off warnings
(ti-menu-load-string "/solve/set/flow-warnings n")
;;
;; ----- run the solver -----
(if (> t-step-size 0)
    (if (string=? fl-freq-type "time-cad")
	(begin ; TRUE	
          (ti-menu-load-string
           (format #f "solve/set/time-step ~d" (/ t-step-size (* 6 op-speed))) )
    	)
	(begin ; FALSE
          (ti-menu-load-string (format #f "solve/set/time-step ~d" t-step-size) )
          )
        )
    )
;;
(if (string=? fl-freq-type "time-step")
    (begin ;TRUE
	(ti-menu-load-string (format #f "solve/dual-time-iterate ~d ~d"
                                 (- nt-stop nt-start) max-its) )
    )
    (begin ;FALSE
    (ti-menu-load-string (format #f "solve/dual-time-iterate 100000 ~d"
                                 max-its))
    )
)
;;
(ti-menu-load-string "file/stop-transcript")
;;
;; ========================================================================== ;;
