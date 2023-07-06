;; ========================================================================== ;;
;; File to run specify the data frequency                                     ;;
;;   Called from: case.scm                                                    ;;
;; ========================================================================== ;;
;;
;; ***** define functions
(define (mduc-write-case-data suffix)
  (define tmp-fn (format #f "~a/~a/~a_~a" case-dir fl-data file-prefix suffix ))
  (display (format #f "\nMDUC writing case file: ~a\n" tmp-fn))
  (if (file-exists? (format #f "~a.cas" tmp-fn))
      (ti-menu-load-string (format #f "wcd ~a y\n" tmp-fn))
      (ti-menu-load-string (format #f "wcd ~a\n"   tmp-fn))
      )
  ); end mduc-file-name
;;
(define (mduc-resubmit)
  (display "Restarting the simulation\n")
  (mduc-write-case-data (mduc-file-suffix fl-freq-type -1))
  (system (format #f "~a/journal/resubmit.bash ~a ~d ~d ~a" case-root
                  case-script (%rpgetvar 'time-step) 
                  (mduc-get-time fl-freq-type)
                  (format #f "resubmit-~.0f" 
                          (mduc-mod (mduc-get-time fl-freq-type) 720.0) )))
  (exit)
  )
;;
(define (mduc-check-simulation)
  (define mduc-wall-t-now  (/ (- (time) mduc-wall-t-start) 3600))
  (display "\n------ checking simulation -----\n")
  (display (format #f "\n Wall time: ~8.3F --> Remaining ~8.3F\n\n"
                   mduc-wall-t-now (- max-wall-time mduc-wall-t-now) ) )
  (set! mduc-time (mduc-get-time fl-freq-type))
  (if (or (string=? fl-freq-type "cad") (string=? fl-freq-type "time-cad"))
      (begin
        (display (format #f "~8a ~8a ~8a ~8a ~8a ~8a\n"
                         "CAD" "current" "freq" "last" "next" "mod") )
        (display (format #f "~8a ~8.3F ~8.3F ~8.3F ~8.3F ~8.3F\n"
                         "" mduc-time fl-freq mduc-pre-t (+ mduc-pre-t fl-freq)
                         (mduc-mod mduc-time fl-freq) ))
        )
      (begin
        (display (format #f "~8a ~8a ~8a ~8a ~8a ~8a\n"
                         "time(s)" "current" "freq" "last" "next" "mod") )
        (display (format #f "~8a ~8.3E ~8.3E ~8.3E ~8.3E ~8.3E\n"
                         "" mduc-time fl-freq mduc-pre-t (+ mduc-pre-t fl-freq)
                         (mduc-mod mduc-time fl-freq) ))
        )
      )
  ;; check to see if the simulation has reached the maximum wall time (- 20min)
  (if (> (+ mduc-wall-t-now 0.3) max-wall-time)
      (begin ;;
        (display "\nMaximum wall time exceeded --> stopping simulation\n")
        (mduc-write-case-data (mduc-file-suffix fl-freq-type -1))
        ;; resubmit the job if not at the stop time
        (if (< mduc-time nt-stop)
            (begin
              (display "\nResubmitting simulation with command:\n")
              (display (format #f "~a/journal/resubmit.bash ~a ~d ~a" case-root
                               case-script (%rpgetvar 'time-step)
                               (mduc-get-time fl-freq-type) ))
              (system (format #f "~a/journal/resubmit.bash ~a ~d ~a" case-root
                              case-script (%rpgetvar 'time-step)
                              (mduc-get-time fl-freq-type) ))
              )
            )
        (exit)
        )
      ); end wall time check
  ;; check to see if a data file needs to be written
  (if (or (< (mduc-mod mduc-time fl-freq) (/ (rpgetvar 'physical-time-step) 2))
          (> (+ mduc-time 1.0e-9) (+ mduc-pre-t fl-freq)) )
      (begin ;;
        (display "\nWriting Fluent data file\n")
        (set! mduc-pre-t mduc-time)
        (mduc-write-case-data (mduc-file-suffix fl-freq-type -1))
        )
      ); end data write
  ;; check to see if final time has been reached 
  (if (>= mduc-time nt-stop)
      (begin ;; if
        (mduc-write-case-data 
         (format #f "final_~a" (mduc-file-suffix fl-freq-type -1)))
        (exit)
        )
      ); end data write
  (display "\n ---------------------------------\n")
  ); end mduc-check-simulation
;; ***** end function definition
;;
;; ----- initialize variables
(define mduc-wall-t-start (time)) ; initial wall time of simulation
(define mduc-wall-t-stop 0) ; flag for noting when simulation should stop
(display (format #f "Starting wall clock: ~a" mduc-wall-t-start))
(define mduc-time (mduc-get-time fl-freq-type))
(display (format #f "Simulation start time: ~a" mduc-time))
(define mduc-pre-t (mduc-get-time fl-freq-type))
;;
;; ----- add the command to execute
(ti-menu-load-string "\nsolve/execute-commands/add-edit mduc-check 1 
                        \"time-step\" \"(mduc-check-simulation)\"")
;;
;; ========================================================================== ;;
