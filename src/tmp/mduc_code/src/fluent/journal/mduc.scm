;; ========================================================================== ;;
;; File to setup the MDUC combustion module.                                  ;;
;;   Called from: run.scm                                                     ;;
;; ========================================================================== ;;
;;
;; ***** define all rpvars that the mduc library could use *****
(rp-var-define 'mduc/directory               "" 'string #f)
(rp-var-define 'mduc/file-prefix             "" 'string #f)
(rp-var-define 'mduc/case-id                 "" 'string #f)
(rp-var-define 'mduc/inputfile               "" 'string #f)
(rp-var-define 'mduc/combustion-model        "" 'string #f)
(rp-var-define 'mduc/sdm-model               #f 'boolean #f)
(rp-var-define 'mduc/restart                 #f 'boolean #f)
(rp-var-define 'mduc/scalars                 "" 'string #f)
(rp-var-define 'mduc/egr-species             "" 'string #f)
(rp-var-define 'mduc/egr-mass-fractions      "" 'string #f)
(rp-var-define 'mduc/drop-ref-temp         -1.0 'real #f)
(rp-var-define 'mduc/sdm/nist-sat-data       "" 'string #f)
(rp-var-define 'mduc/sdm/xi-max            -1.0 'real #f)
(rp-var-define 'mduc/initial-pressure-rate -1.0 'real #f)
(rp-var-define 'mduc/initial-temperature   -1.0 'real #f) 
(rp-var-define 'mduc/cylinder-mass         -1.0 'real #f) 
(rp-var-define 'mduc/avg-z	           -1.0 'real #f) 
(rp-var-define 'mduc/avg-egr                -1.0 'real #f) 
(rp-var-define 'mduc/scale-z                -1.0 'real #f) 
(rp-var-define 'mduc/avg-mixture              #f 'boolean #f) 
(rp-var-define 'mduc/transformation          -1 'integer #f)
;; injection properties
(rp-var-define 'mduc/init-pilot-time       -1.0 'real #f)
(rp-var-define 'mduc/init-main-time        -1.0 'real #f)
(rp-var-define 'mduc/inj2/n-grid             -1 'integer #f)
(rp-var-define 'mduc/inj2/n-fuel             -1 'integer #f)
(rp-var-define 'mduc/inj2/fuel               "" 'string #f)
(rp-var-define 'mduc/inj2/mass-fraction      "" 'string #f)
(rp-var-define 'mduc/inj2/temp             -1.0 'real #f)
;; scalar dissipation rate
(rp-var-define 'mduc/chi-1d-mode             "" 'string #f)
(rp-var-define 'mduc/chi-1d-form             "" 'string #f)
(rp-var-define 'mduc/chi-1d-scaling          "" 'string #f)
(rp-var-define 'mduc/chi-2d-mode             "" 'string #f)
(rp-var-define 'mduc/chi-2d-form             "" 'string #f)
(rp-var-define 'mduc/chi-2d-scaling          "" 'string #f)
(rp-var-define 'mduc/variance-scale         1.0 'real #f)
(rp-var-define 'mduc/dissipation-scale      1.0 'real #f)
;; distribution functions
(rp-var-define 'mduc/pdf-1d                  "" 'string #f)
(rp-var-define 'mduc/pdf-2d                  "" 'string #f)
;; rif model
(rp-var-define 'mduc/rif/dimensions          -1 'integer #f)
(rp-var-define 'mduc/rif/x-coord      "default" 'string #f)
(rp-var-define 'mduc/rif/y-coord      "default" 'string #f)
(rp-var-define 'mduc/enthalpy-buffer       -1.0 'real #f)
(rp-var-define 'mduc/solve-zones             "" 'string #f)
(rp-var-define 'mduc/intake-zones            "" 'string #f)
(rp-var-define 'mduc/exhaust-zones           "" 'string #f)
(rp-var-define 'mduc/solve-thread-names      "" 'string #f)
(rp-var-define 'mduc/solve-thread-ids        "" 'string #f)
(rp-var-define 'mduc/thread-names            "" 'string #f)
(rp-var-define 'mduc/thread-ids              "" 'string #f)
(rp-var-define 'mduc/init-h-from-table       -1 'integer #f)
;; monitor
(rp-var-define 'mduc/monitor/zone-names      "" 'string #f)
(rp-var-define 'mduc/monitor/zone-threads    "" 'string #f)
(rp-var-define 'mduc/monitor/wall-names      "" 'string #f)
(rp-var-define 'mduc/monitor/wall-threads    "" 'string #f)
(rp-var-define 'mduc/monitor/surf-names      "" 'string #f)
(rp-var-define 'mduc/monitor/surf-threads    "" 'string #f)
(rp-var-define 'mduc/monitor/surf-ref-names  "" 'string #f)
(rp-var-define 'mduc/monitor/surf-ref-ids    "" 'string #f)
;; wiebe function
(rp-var-define 'mduc/mass-fuel             -1.0 'real #f)
(rp-var-define 'mduc/mass-fuel-2           -1.0 'real #f)
(rp-var-define 'wiebe/coeff-a              -1.0 'real #f)
(rp-var-define 'wiebe/coeff-m              -1.0 'real #f)
(rp-var-define 'wiebe/soc                   0.0 'real #f)
(rp-var-define 'wiebe/soc-alpha            -1.0 'real #f)
(rp-var-define 'wiebe/soc-beta             -1.0 'real #f)
(rp-var-define 'wiebe/rxn-c                -1.0 'real #f)
(rp-var-define 'wiebe/rxn-h                -1.0 'real #f)
(rp-var-define 'mduc/lower-heating-value   -1.0 'real #f)
;; premixed variables
(rp-var-define 'mduc/prog-var-species           "" 'string #f)
(rp-var-define 'mduc/premixed/table             "" 'string #f)
(rp-var-define 'mduc/premixed/table-dim          3 'integer #f)
(rp-var-define 'mduc/premixed/pressure-table    "" 'string #f)
(rp-var-define 'mduc/premixed/ignition-time   -1.0 'real #f)
(rp-var-define 'mduc/premixed/ignition-geom     "" 'string #f)
(rp-var-define 'mduc/premixed/ignition-loc      "" 'string #f)
(rp-var-define 'mduc/premixed/ignition-radius -1.0 'real #f)
(rp-var-define 'mduc/premixed/ignition-size     "" 'string #f)
(rp-var-define 'mduc/premixed/turb-flame-speed  "zimont" 'string #f)
(rp-var-define 'mduc/premixed/prog-threshold    0.1 'real #f)
(rp-var-define 'mduc/premixed/prog-src-scale    1.0 'real #f)
(rp-var-define 'mduc/premixed/ignition-c-scale  1.0 'real #f)
;;
;; ***** end variable definition *****
;;
;; ----- set user defined scalars and memory -----
;;(ti-menu-load-string "define/user-defined/user-defined-scalars 6" )
;;(ti-menu-load-string "define/user-defined/user-defined-memory 14" )
;;
;; ----- load and hook library-----
;;(ti-menu-load-string "\ndefine/user-defined/compiled load libudf/mduc")
(ti-menu-load-string
 "define/user-defined/function-hooks/adjust/\"adjust::mduc\"" )
(set! udf-at-end (string-append udf-at-end " \"solve::mduc\""))
;;
;; ----- set the additional quantities to monitor
;;(set! mon-m-int (append mon-m-int '(z1-mean z2-mean)))
;;(set! mon-m-avg (append mon-m-avg '(z1-variance z2-variance z-variance h-mean
;;                                                chi-ref-0 chi-ref-1 chi-ref-2)))
;;
;; ----- set up mduc variables that are implicitly passed on
(rpsetvar 'mduc/directory   case-dir)
(rpsetvar 'mduc/file-prefix file-prefix)
(rpsetvar 'mduc/case-id     case-sub-id)
(rpsetvar 'mduc/inputfile
          (format #f "~a/~a/init/~a" case-root case-spec mduc-input))
(rpsetvar 'mduc/combustion-model  mduc-model)
(rpsetvar 'mduc/sdm-model         mduc-sdm)
(rpsetvar 'mduc/restart           case-restart)
;;
;; ----- read/set case specific mduc variables -----
(ti-read-scheme (format #f "~a/~a/init/~a" case-root case-spec mduc-fl) )
;;
;; ----- need to get the ids for the cx surfaces here and pass them to library
(define trim-chars  (apply char-set '(#\( #\) #\space)))
(define mduc-surf (str-split (rpgetvar 'mduc/monitor/surf-threads) trim-chars) )
(define mduc-surf-nm "")
(define mduc-surf-id "")
;;
(for-each
 (lambda (surf)
   (define tmp-id (surface-name/id? (string->symbol surf)))
   (set! mduc-surf-id (string-append mduc-surf-id
                                     (format #f "~d " (list-ref tmp-id 0))))
   (set! mduc-surf-nm (string-append mduc-surf-nm (format #f "~a " surf)))
   ); end lambda
   mduc-surf
 ); end for-each
(display (format #f "\nMDUC surface summary:\n"))
(display (format #f "surface-ref-names: ~a\n" mduc-surf-nm))
(display (format #f "surface-ref-ids:   ~a\n" mduc-surf-id))
(rpsetvar 'mduc/monitor/surf-ref-names mduc-surf-nm)
(rpsetvar 'mduc/monitor/surf-ref-ids   mduc-surf-id)
;;
;; ----- make sure the facets on the surfaces are initialized
;(ti-menu-load-string (format #f "display/surface-mesh ~a" mduc-mon-surfaces) )
;(ti-menu-load-string "surface/list-surfaces")
;;
;; ----- initialize library -----
(ti-menu-load-string "define/user-defined/execute-on-demand/\"init::mduc\"")
;;
;; ========================================================================== ;;
