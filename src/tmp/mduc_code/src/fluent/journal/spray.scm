;; ========================================================================== ;;
;; File to setup the spray libraries                                          ;;
;;   Called from: run.scm                                                     ;;
;; ========================================================================== ;;
;;
;; ***** define all rpvars that the spray library needs
(rp-var-define 'spray/injection-data     "" 'string #f)
(rp-var-define 'spray/injection-massflow "" 'string #f)
(rp-var-define 'spray/injection-velocity "" 'string #f)
(rp-var-define 'spray/injection-dropsize "" 'string #f)
(rp-var-define 'spray/monitor-file       "" 'string #f)
(rp-var-define 'spray/nozzles            -1 'integer #f)
(rp-var-define 'spray/injection-ids      "" 'string #f)
(rp-var-define 'spray/dumouchel-alpha  -1.0 'real #f)
;; ***** end variable definition
;;
;; ----- read nozzle geometry and injection parameters
(ti-read-scheme (format #f "~a/setup/~a" case-root inj-geom) )
(ti-read-scheme
 (format #f "~a/~a/injection/~a" case-root case-spec inj-settings) )
;;
;; ----- set variables in solver 
;; specify mass flow/velocity/droplet distribution files
(rpsetvar 'spray/injection-data
          (format #f "~a/~a/~a/" case-root case-spec inj-data))
(rpsetvar 'spray/injection-massflow inj-mdot-file)
(rpsetvar 'spray/injection-velocity inj-vel-file)
(rpsetvar 'spray/injection-dropsize inj-drop-file)
;; nozzle info
(rpsetvar 'spray/nozzles  (length noz-x0))
(rpsetvar 'spray/injection-ids  (format #f "~a" inj-name))
;; if computing a Dumouchel distribution, need additional parameter
(rpsetvar 'spray/dumouchel-alpha  inj-alpha)
;;
;; ----- define injections 
(if (not case-restart)
    (ti-read-scheme (format #f "~a/journal/injector.scm" case-root) )
    ;; delete the dummy injection
    (ti-menu-load-string "define/injections/delete-injection injection-0")
)
;;
;; ----- load the library -----
;;(ti-menu-load-string "\ndefine/user-defined/compiled load libudf/spray")
(set! udf-at-end (string-append udf-at-end " \"writeStats::spray\""))
;; set the monitor file
(rpsetvar 'spray/monitor-file
          (format #f "~a/monitor/spray.~a" case-dir case-sub-id))
;;
;; ----- load the mass flow/velocity profiles
(ti-menu-load-string "(ti-udf-on-demand) \"loadProfile::spray\"")
;;
;; ========================================================================== ;;
