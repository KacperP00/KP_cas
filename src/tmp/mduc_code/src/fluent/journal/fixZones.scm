;;;;     SCHEME COMMAND TO CREATE VOLUME AND SURFACE MONITORS      ;;;;
;;;;Created By Sourabh Shrivastava : sourabh.shrivastava@ansys.com ;;;;
;;;User has to provide following data
;;;1). List of old surface boundaries and cell zones that needs to be included in the new monitors
;;;2). Number of monitors to be defined
;;;3). Type of Monitors (for volume and surface monitors separately)
;;;4). Quantity to monitor
;;;5). Name for monitor files (please make sure to assign different names for monitors compare to the existing one)
;;;6). Frequency of data entry
;;;
;;;Steps to use this file
;;;1)-Read case file
;;;2)-Read the Scheme file
;;;3)-Define 2 new events and define 'Execute Command' as Type for these events
;;;4)-For first event define any crank angle such the this occures before the 'Insert Interior Layer' event and write command as '(mesh-replace-initial) under Execute Command Definition fo the event
;;;5)-For second event define crank angle same as angle at which 'Insert Interior Layer' event is executed and Define command as '(mesh-replace-final) under Execute Command Definition fo the event
;;;

(define initial-walls '())
(define initial-fluids '())
(define final-walls '())
(define final-fluids '())
(define old-cells '(fluid-ch)) ; write cell volume that already exist and that needs to be included in the current monitors

(define old-boundaries '(cyl-head cyl-tri)) ; write boundaries that already exist and that needs to be included in the current monitors

(define cyl-zone "")
(define wall-zone "")
(define dumy1 "")
(define dumy2 "")
(define dumy3 "")
(define dumy4 "")
(define mesh-replace-initial
	(lambda ()
    (set! initial-walls
          (map thread-name (sort-threads-by-name (get-threads-of-type 'wall))))
    ; (pp walls)
    (set! initial-fluids
          (map thread-name (sort-threads-by-name (get-threads-of-type 'fluid))))
    (set! dumy1 "")
    (set! dumy2 "")
    (set! dumy3 "")
    (set! dumy4 "")
    (ti-menu-load-string "solve/monitors/volume/clear-monitors")
    (ti-menu-load-string "solve/monitors/surface/clear-monitors")
    ; (ti-menu-load-string (format #f "\nfile/transient-export/delete cfd-post-1 "))
;    (ti-menu-load-string
;     (format #f "\nsolve/execute-commands/disable export-udf-r14"))
;    (ti-menu-load-string
;     (format #f "\nsolve/execute-commands/disable export-udf-write"))
;    (ti-menu-load-string
;     (format #f "\nsolve/execute-commands/disable plot-pdf "))
    (ti-menu-load-string "\ndefine/dynamic-mesh/controls/ in-cylinder-out n")
    (set-bc-wall-hflux initial-walls #f)
    ; (pp fluids)
    ) ; end lambda
) ; end mesh-replace-initial

; ----- function to call when 'Insert Interior Layer is called'
(define mesh-replace-final
  (lambda ()
    (let ((wall-new '()) (fluid-new '()))
      (define count 0)
      (set! final-walls (map thread-name (sort-threads-by-name (get-threads-of-type 'wall))))
      (set! final-fluids (map thread-name (sort-threads-by-name (get-threads-of-type 'fluid))))
      ; ----- Boundary operations
      (do ((i 0 (+ i 1))) ((= i (length final-walls)))
        (set! count 0)
        ; (pp i)
        (do ((j 0 (+ j 1))) ((= j (length initial-walls)))
          
          (if (string-ci=? (symbol->string (list-ref final-walls i))
                           (symbol->string (list-ref initial-walls j)))
              (set! count 1)
              )
          ;(pp "loop2")
          )
        (if (= count 0)
            (set! wall-new (append wall-new (list (list-ref final-walls i))))
            )
        )
      (set! wall-new (append wall-new old-boundaries)) 
      (pp wall-new)
      (do ((i 0 (+ i 1))) ((= i (length wall-new)))
        (set! dumy3 (symbol->string (list-ref wall-new i)))
        (set! dumy4 (string-append dumy4 " "dumy3))
        ) 
      ; ----- Cell Zone operations
      (do ((i 0 (+ i 1))) ((= i (length final-fluids)))
                                        ;(pp "fluid loop")
        (set! count 0)
        (do ((j 0 (+ j 1))) ((= j (length initial-fluids)))
          (if (string-ci=? (symbol->string (list-ref final-fluids i))
                           (symbol->string (list-ref initial-fluids j)))
              (set! count 1)
              )
          )
        (if (= count 0)
            (set! fluid-new (append fluid-new (list (list-ref final-fluids i))))
            )
        )
      (pp fluid-new)
      (set! fluid-new (append fluid-new old-cells)) 
      (do ((i 0 (+ i 1))) ((= i (length fluid-new)))
        (set! dumy1 (symbol->string (list-ref fluid-new i)))
        (set! dumy2 (string-append dumy2 " "dumy1))
        )
      
	    ; delete all monitor plots, cylinder-output and solution export and
      ; re-define them
      (set! wall-zone (string-append "(" dumy4 ")"))
      (set! cyl-zone (string-append "(" dumy2 ")"))
      (ti-read-scheme (format #f "~a/journal/monitor.scm" case-root ))
      (ti-menu-load-string "\ndefine/dynamic-mesh/controls/in-cylinder-out y")
      (set-bc-wall-hflux final-walls #t)
;      (rpsetvar 'export/init #t)
;      (ti-menu-load-string "(ti-udf-on-demand) \"init::export\"")
;      (ti-menu-load-string
;       (format #f "\nsolve/execute-commands/enable export-udf-r14 ~d" exp-freq))
;      (ti-menu-load-string
;       (format #f "\nsolve/execute-commands/enable export-udf-write ~d" exp-freq))
;      (ti-menu-load-string
;       (format #f "\nsolve/execute-commands/enable plot-pdf ~d"exp-freq))
       ;(ti-read-scheme (format #f "~a/journal/export_hot.scm" case-root ))
      )
    )
  ) ; end mesh-replace-final

;(mesh-replace-initial)
;(mesh-replace-final)
