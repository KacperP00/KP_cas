;; ========================================================================== ;;
;; File to automatically plot images during the simulation                    ;;
;;  Called from: run.scm                                                      ;;
;;  Requires:    case.scm, custom-funcs.scm                                   ;;
;; ========================================================================== ;;
;;
;; setup the plot view (intake for now)
;; ----- define plot functions
(define plot-surf '(plane-122))
(define (mduc-plot-image)
  (define plot-file-name (format #f "~a/plots/plot_scalar_~a.png" case-dir
                                 (mduc-file-suffix plot-freq-type -1)))
  (ti-menu-load-string "/display/set/contours/surfaces a_intake\n")
  (ti-menu-load-string "/views/restore-view a_plot_view\n")
  (ti-menu-load-string "/display/contour egr 0 1\n")
  (if (file-exists? plot-file-name)
      (ti-menu-load-string
       (format #f "/display/save-picture ~a y\n" plot-file-name))
      (ti-menu-load-string
       (format #f "/display/save-picture ~a\n" plot-file-name))
      )
  ) ;; end mduc-plot-image
;;
(define (mduc-plot-freq)
  (if (< (mduc-mod (+ (mduc-get-time plot-freq-type) 1.0e-12) plot-freq) 1.0e-9)
      (begin ;; if
        (display "\n ----- plotting image ----- \n")
        (ti-menu-load-string "(mduc-plot-image)")
        (display "\n ---------------------------------\n")
        ) ;; end if
      ) ;; end if
  ); end exp-udf-check-freq
;;
(define plot-cur-t (mduc-get-time plot-freq-type))
;;
;; ----- set up the parameters for 
(ti-menu-load-string "/display/set/picture/invert y") 
(ti-menu-load-string "/display/set/picture/driver png")
(ti-menu-load-string "/display/set/picture/color-mode color")
(ti-menu-load-string "/display/set/picture/x-resolution 960")
(ti-menu-load-string "/display/set/picture/y-resolution 720") 
(ti-menu-load-string "/display/set/contours/clip-to-range no")
(ti-menu-load-string "/display/set/contours/filled-contours yes")
(ti-menu-load-string "/display/set/contours/node-values yes")
;;
;; ----- set up the view for plotting (intake for now)
(ti-menu-load-string
 (format #f "/views/read-views/ \"~a/setup/plot_views.vw\"" case-root))
;;
;; ----- hook plot function to be called at the right frequency 
(ti-menu-load-string "\nsolve/execute-commands/add-edit mduc-plot 1 
                        \"time-step\" \"(mduc-plot-freq)\"")
;;
;; ========================================================================== ;;
