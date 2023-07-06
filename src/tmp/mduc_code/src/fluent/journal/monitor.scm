;; ========================================================================== ;;
;;                                                                            ;;
;; File to set up monitors for Fluent variables.                              ;;
;; Note that some variables require that mduc-cff.scm (MDUC library)          ;;
;; be loaded                                                                  ;;
;; Requires: case.scm, mon_zones                                              ;;
;;                                                                            ;;
;; ========================================================================== ;;
;;
;
;; by-default, monitor all volume zones.  If a specific set is necessary, define
;; a list according to (define xxx-vols '(zone1 zone2 ...))
(define mon-zones (get-avail-zones))
;;
;; Volume average monitors
(do ((n 1 (+ n 1) )) ((> n (length mon-v-avg)))
  (ti-menu-load-string
   (format #f "\nsolve/monitors/volume/set-monitor ~a \"Volume-Average\" ~a '~a n n y ~a/monitor/~a.~a 1 y flow-time"
           (list-ref mon-v-avg (- n 1)) (list-ref mon-v-avg (- n 1)) mon-zones
           case-dir (list-ref mon-v-avg (- n 1)) case-sub-id) )
  )
;;
;; Mass average monitors
(do ((n 1 (+ n 1) )) ((> n (length mon-m-avg)))
  (ti-menu-load-string
   (format #f "\nsolve/monitors/volume/set-monitor ~a \"Mass-Average\" ~a '~a n n y ~a/monitor/~a.~a 1 y flow-time"
           (list-ref mon-m-avg (- n 1)) (list-ref mon-m-avg (- n 1)) mon-zones
           case-dir (list-ref mon-m-avg (- n 1)) case-sub-id) )
  )
;;
;; add species to mass integral monitors
(set! mon-m-int (append mon-m-int mon-spec))
;; Mass integral monitors
(do ((n 1 (+ n 1) )) ((> n (length mon-m-int)))
  (display (- n 1))
  (ti-menu-load-string
   (format #f "\nsolve/monitors/volume/set-monitor ~a \"Mass Integral\" ~a '~a n n y ~a/monitor/~a.~a 1 y flow-time"
           (list-ref mon-m-int (- n 1)) (list-ref mon-m-int (- n 1)) mon-zones
           case-dir (list-ref mon-m-int (- n 1)) case-sub-id) )
  )
;;
;; ----- read in scheme residuals
(if (not mduc-on)
    (ti-menu-load-string "\nsolve/execute-commands disable my_res_out") )
;;
;; ----- read in custom monitors if defined
(if (> (string-length mon-custom-file) 0)
    (ti-read-scheme (format #f "~a/setup/~a" case-root mon-custom-file)) )
;;
;; ========================================================================== ;;
