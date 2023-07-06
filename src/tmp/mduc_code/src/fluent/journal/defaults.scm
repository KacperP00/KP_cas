;; ========================================================================== ;;
;;                             Fluent Run Script                              ;;
;;                                                                            ;;
;;                       (ti-read-scheme "case.scm")                          ;;
;; ========================================================================== ;;
;; ----- Fluent  -----
(define case-root  (get-pwd))
(define case-script "")  ;; name of current file
;;
;; ----- set up case naming -----
(define case-name    "") ;; name of top level directory
(define case-map     "") ;; test map directory
(define case-spec    "") ;; operating point directory
(define case-id      "") ;; specify arbitrary label here
(define case-sub-id  "") ;; for monitor files
(define init-name    "") ;; name of initial file (non-restart)
(define restart-file "") ;; set to empty if restarting from non-timestep file
(define case-custom  "") ;; a file to read (in /init) for custom settings
;;
;; ----- set up run options -----
(define case-restart    #f) ;; is this case restarting?
(define nt-init          0) ;; initial time step number
(define nt-start         0) ;; start from
(define nt-stop      10e-5) ;; run til
(define t-step-size    -1)  ;; < 0 uses what is specified in case file
(define max-its        100)
(define max-wall-time   24) ;; in hours
;;
;; ----- set up data/export options -----
(define fl-data       "fluent")
(define fl-freq       0)
(define fl-freq-type  "flow-time")
;;
(define exp-freq      0)
(define exp-freq-type "flow-time")
(define exp-data      "vtk")
(define exp-setup     "")
(define exp-sol       #f)    ; export the solution according to exp-data format
(define exp-udf       #f)    ; export basic data to ensight
(define exp-particle  #f)    ; export particle history data
(define exp-sdm       #f)    ; export single droplet parameters
;;
;; plotting
(define plot-freq      0)
(define plot-freq-type "time-step")
;;
;; ----- monitor info -----
(define mon-v-avg '())
(define mon-m-avg '())
(define mon-m-int '())
(define mon-spec  '())
(define mon-custom-file "")
;;
;; ----- engine ------
(define engine-case #f)
(define manifold-case #f)
(define op-speed     2000)
(define ca-period 	 720)
(define ca-start	   0)
(define event-file   "events-sym.scm")
(define bc-manifold  "manifold_pressure.ttab")
(define invalve-lift "intake_lift.ttab")
(define exvalve-lift "exhaust_lift.ttab")
;;
;; ----- define spray option files -----
(define inj-on        #f)
(define wall-film-on  #f)
(define inj-data      "injection") ; directory below case-spec
(define inj-geom      "inj-geom.scm")
(define inj-settings  "inj-setup.scm")
(define khrt	      #f) ; wave model is default
(define p-dep-boiling #f)
;;
;; ----- set up the combustion case (located in case directory) -----
(define mduc-on      #f)
(define mduc-model   "none")
(define mduc-sdm     #f)
(define mduc-restart #f)
(define mduc-input   "mduc_rif.inp")
(define mduc-fl      "mduc-setup.scm")
;;
;; ----- start the run -----
;;(ti-read-scheme (format #f "~a/journal/run.scm" case-root ) )
;;
;; ========================================================================== ;;
