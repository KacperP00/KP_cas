;; ========================================================================== ;;
;; File to set up data to be exported during Fluent run using custom UDF.     ;;
;;  Called from: run.scm                                                      ;;
;;  Requires:    case.scm                                                     ;;
;; ========================================================================== ;;
;;
;; ***** define all export variables *****
(rp-var-define 'export/directory             "" 'string  #f)
(rp-var-define 'export/file-prefix           "" 'string  #f)
(rp-var-define 'export/case-name             "" 'string  #f)
(rp-var-define 'export/dump-sdm              #f 'boolean #f)
(rp-var-define 'export/dump-particle         #f 'boolean #f)
(rp-var-define 'export/init                  #f 'boolean #f)
(rp-var-define 'export/binary                #t 'boolean #f)
(rp-var-define 'export/freq                  -1 'integer #f)
(rp-var-define 'export/species               "" 'string  #f)
(rp-var-define 'export/zones                 "" 'string  #f)
(rp-var-define 'export/zone-threads          "" 'string  #f)
(rp-var-define 'export/wall-surfaces         "" 'string  #f)
(rp-var-define 'export/wall-surface-threads  "" 'string  #f)
(rp-var-define 'export/surfaces              "" 'string  #f)
(rp-var-define 'export/surface-threads       "" 'string  #f)
(rp-var-define 'export/surface-thread-ids    "" 'string  #f)
(rp-var-define 'export/surface-ref-names     "" 'string  #f)
(rp-var-define 'export/surface-ref-ids       "" 'string  #f)
; ***** end variable definition *****
;
;; ----- define export functions to check frequency
;; function to be called at every time step to check whether to write
(define (exp-udf-write-data)
  ;; if exporting surfaces, need to reactivate surfaces to run
  (if (> (string-length (rpgetvar 'export/surface-ref-ids)) 0)
      (ti-menu-load-string (format #f "report/surface-integrals area '(~a)"
                                   (rpgetvar 'export/surface-ref-ids) ) ) )
  (ti-menu-load-string "(client-fill-node-values 'density)")
  (ti-menu-load-string "(ti-udf-on-demand) \"writeData::export\"")
) ; end exp-udf-write-data
;;
(define (exp-udf-set-surf-cx-ids)
  (define trim-chars  (apply char-set '(#\( #\) #\space)))
  (define tmp-surf    (str-split (rpgetvar 'export/surface-threads) trim-chars) )
  (define tmp-surf-nm "")
  (define tmp-surf-id "")
  ;;
  (for-each
   (lambda (surf)
     (define tmp-id (surface-name/id? (string->symbol surf)))
     (display tmp-id)
     (set! tmp-surf-id
           (string-append tmp-surf-id (format #f "~d " (list-ref tmp-id 0))))
     (set! tmp-surf-nm (string-append tmp-surf-nm (format #f "~a " surf)))
     ); end lambda
   tmp-surf
   ); end for-each
  (display (format #f "\nSurface summary:\n"))
  (display (format #f "surface-ref-names: ~a\n" tmp-surf-nm))
  (display (format #f "surface-ref-ids:   ~a\n" tmp-surf-id))
  (rpsetvar 'export/surface-ref-names tmp-surf-nm)
  (rpsetvar 'export/surface-ref-ids   tmp-surf-id)
  )
;;
;; function to check whether data should be exported at this time step
(define (exp-udf-check-freq)
  (display "\n ----- Export Check -----\n")
  (set! exp-cur-t (mduc-get-time exp-freq-type))
  (if (or (string=? fl-freq-type "cad") (string=? fl-freq-type "time-cad"))
      (begin
        (display (format #f "~8a ~8a ~8a ~8a ~8a ~8a\n"
                         "CAD" "current" "freq" "last" "next" "mod") )
        (display (format #f "~8a ~8.3F ~8.3F ~8.3F ~8.3F ~8.3F\n"
                         "" exp-cur-t exp-freq exp-pre-t (+ exp-pre-t exp-freq)
                         (mduc-mod exp-cur-t exp-freq) ))
        )
      (begin
        (display (format #f "~8a ~8a ~8a ~8a ~8a ~8a\n"
                         "time(s)" "current" "freq" "last" "next" "mod") )
        (display (format #f "~8a ~8.3E ~8.3E ~8.3E ~8.3E ~8.3E\n"
                         "" exp-cur-t exp-freq exp-pre-t (+ exp-pre-t exp-freq)
                         (mduc-mod exp-cur-t exp-freq) ))
        )
      )
  (if (or (< (mduc-mod exp-cur-t exp-freq) (/ (rpgetvar 'physical-time-step) 2))
          (> (+ exp-cur-t 1.0e-9) (+ exp-pre-t exp-freq)) )
      (begin
        (set! exp-pre-t exp-cur-t)
        (ti-menu-load-string "(exp-udf-write-data)")
        ) )
  (display "\n ---------------------------------\n")
  ); end exp-udf-check-freq
;;
;; initialize the export variables
(define exp-cur-t (mduc-get-time exp-freq-type))
(define exp-pre-t (mduc-get-time exp-freq-type))
;;
;; ----- load the library
(if (not case-restart)
    (ti-menu-load-string
     (format #f "\ndefine/user-defined/compiled load libudf/export") )
    );; if restart
;;
;; ----- add the command to check export at every time-step
(ti-menu-load-string "\nsolve/execute-commands/add-edit export-udf 1 
                        \"time-step\" \"(exp-udf-check-freq)\"")
;;
;; set file name information
(rpsetvar 'export/directory   (format #f "~a/~a" case-dir exp-data))
(rpsetvar 'export/file-prefix (format #f "~a" file-prefix))
;;
;; pass flags to UDF library
(rpsetvar 'export/case-name     case-name)
(rpsetvar 'export/dump-sdm      exp-sdm)
(rpsetvar 'export/dump-particle exp-particle)
(rpsetvar 'export/init          case-restart)
(rpsetvar 'export/binary        #t)
(rpsetvar 'export/freq          exp-freq)
;;
;; ----- call any custom setup for the species, zones, surfaces, etc...
(if (> (string-length exp-setup) 0)
    (ti-read-scheme (format #f "~a/~a/init/~a" case-root case-spec exp-setup) ) )
;;
;; ----- need to get the ids for the cx surfaces here and pass them to library
(if (> (string-length (rpgetvar 'export/surfaces)) 0)
    (ti-menu-load-string "(exp-udf-set-surf-cx-ids)") )
;;
;; ----- initialize the library
(ti-menu-load-string "(ti-udf-on-demand) \"init::export\"")
;;
;; ----- check whether this time-step should be exported
(ti-menu-load-string "(exp-udf-check-freq)")
;;
;; ========================================================================== ;;
