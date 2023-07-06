;;=========================================================================== ;;
;; File to set up data directory structure for Fluent cases                   ;;
;;   Called from: run.scm                                                     ;;
;;   Requires:    case.scm                                                    ;;
;; ========================================================================== ;;
;; 
;; define function to check directory path
(display "setting up directory structure...")
(define (check-dir path-name)
  (if (file-exists? path-name )
      (display (format #f "~a exists\n" path-name))
      (begin
        (display (format #f "creating ~a\n" path-name))
        (system (format #f "mkdir ~a" path-name )) )
      )
  )
;; make sure that the specified directory is there
(check-dir (format #f "~a/~a" case-root case-spec) )
(check-dir (format #f "~a" case-dir) )
;;
;; ----- fluent data directory -----
(check-dir (format #f "~a/~a" case-dir fl-data) )
;; link libudf folder there
(if (file-exists? (format #f "~a/~a/libudf" case-dir fl-data))
    (begin ;; if
      (pp "link to libudf directory exists\n")
      ) ;; end if
    (begin ;; else
      (pp "linking libudf directory to run data\n")
      (system (format #f "ln -s ../../../libudf ~a/~a/libudf" case-dir fl-data))
      ) ;; end else
    )
;;
;; ----- monitor directory -----
(check-dir (format #f "~a/monitor" case-dir) )
;;
;; ----- plots directory -----
(check-dir (format #f "~a/plots" case-dir) )
;;;
;; ----- pdf directory -----
(check-dir (format #f "~a/pdf" case-dir) )
;;
;; ----- mduc directory -----
(if mduc-on
    (begin
      (check-dir (format #f "~a/mduc" case-dir) )
      ;; copy the run files to the case directory
      (system (format #f "cp -f ~a/~a/init/~a ~a/~a" case-root case-spec
                      mduc-input case-dir mduc-input))
      (system (format #f "cp -f ~a/~a/init/~a ~a/~a" case-root case-spec mduc-fl
                      case-dir mduc-fl))
      )
    )
;; ----- create export directory if necessary -----
(if exp-sol
    (begin ;; if
      (check-dir (format #f "~a/~a" case-dir exp-data) )
      (check-dir (format #f "~a/~a/export" case-dir exp-data) )
      (if exp-particle
          (check-dir (format #f "~a/particle" case-dir) )
          (display "No particle data will be exported\n" ) )
      )
    (display "Fluent will not export any data\n" )
    )
(if exp-udf
    (begin ;; if
      (check-dir (format #f "~a/~a" case-dir exp-data) )
      (check-dir (format #f "~a/~a/data" case-dir exp-data) )
      (if exp-particle
          (check-dir (format #f "~a/~a/particle" case-dir exp-data) )
          (display "No particle data will be exported\n" ) )
      )
    (display "No ensight will be exported\n" )
    )
;;
;; ========================================================================== ;;
