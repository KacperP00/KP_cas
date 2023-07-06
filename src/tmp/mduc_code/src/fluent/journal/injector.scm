;; ========================================================================== ;;
;; File creates injection parameters based on                                 ;;
;; file inj-settings.scm                                                      ;;
;; ========================================================================== ;;
;;
(pp (format #f "setting up ~d injections\n" (length inj-name)))
;; ***** define functions
(define (inj-calc-xv n-noz z-ref r-ref cone-angle phi-offset sym)
  ;; note, will operate on x0 and v0, so they need to be defined
  ;; requires matrix.scm, custom-funcs.scm
  (define x '())
  (define v '())
  (define noz-x)
  (define noz-v)
  (define phi)
  ;; loop over nozzles
  ;;  (do ((i 0 (+ i 1)) ) ((= i n-noz))
  (do ((i 0 (+ i 1)) ) ((= i (* n-noz (- 1 (* sym 0.5)))))
    (set! phi (+ (deg->rad phi-offset) (* (* 2 3.141592653589793) (/ i n-noz))))
    (pp (format #f "current angle: ~d\n" (rad->deg phi)))
    ;; calculate the hole location and vector of this hole
    (set! noz-x (list (* r-ref (cos phi)) (* r-ref (sin phi)) z-ref) )
    (set! noz-v (list (cos phi) (sin phi)
                      (/ r-ref (tan (deg->rad (/ cone-angle 2)))) ) )
    ;; add to the global arrays
    (set! x (append x (list noz-x)))
    (set! v (append v (list noz-v)))
    ) ;; end loop over nozzles
  (list x v) 
  ) ;; end function to calculate inj
;;
;; ----- function to draw a nozzle vector
;; requires matrix.scm
(define (draw-noz-vec k x0 v)
  (define name "")
  (if (= k 0)
      (set! name "injector")
      (set! name "nozzle")
      )
  ;; create a point at the hole origin
  (ti-menu-load-string
   (format #f "/surface/point-surface b-~a-pt-~d ~d ~d ~d" name k
           (list-ref x0 0) (list-ref x0 1) (list-ref x0 2) ) )
  ;; create a vector ~10mm long (assumes unit vector)
  (define x1 (vec-add x0 (vec-scale 0.01 v)))
  (ti-menu-load-string
   (format #f "/surface/line-surface b-~a-vec-~d ~d ~d ~d ~d ~d ~d" name k
           (list-ref x0 0) (list-ref x0 1) (list-ref x0 2)
           (list-ref x1 0) (list-ref x1 1) (list-ref x1 2) ))
) ;; end draw-noz-vec
;;
(define (mduc-set-injection inj noz x v)
  (cond 
   ( (string=? (format #f "~a" (inquire-release)) "(15 0 0)") 
     (begin
       (pp "setting injection for 15.0\n")
       (ti-menu-load-string (format #f "define/injections/create-injection
        ~a-~d y ~a y ~a ~a ~d n n y \"initDPM::spray\" n n ~a y n ~d ~d ~d ~d
        ~d ~d 0 360.0 ~d ~d ~d constant ~d constant ~d ~d constant ~d
        ~d ~d ~d ~d ~d"
        (list-ref inj-name inj) noz inj-particle inj-type inj-material
        inj-n-part (list-ref inj-evap-spec inj) 
        (list-ref x 0) (list-ref x 1) (list-ref x 2)
        inj-part-temp (list-ref inj-t-start inj) (list-ref inj-t-stop inj)
        (list-ref v 0) (list-ref v 1) (list-ref v 2)
        inj-vel-mag inj-hoca inj-rad inj-mflow inj-dmin inj-dmax inj-dmean
        inj-spread inj-n-diam))
       )
     )
   ( (string=? (format #f "~a" (inquire-release)) "(14 5 0)")
     (begin
       (pp "setting injection for 14.5\n")
       (ti-menu-load-string (format #f "define/injections/create-injection
        ~a-~d y ~a y ~a ~a ~d n n y \"initDPM::spray\" n n ~a y n ~d ~d ~d ~d
        ~d ~d 0 360.0 ~d ~d ~d constant ~d constant ~d ~d constant ~d
        ~d ~d ~d ~d ~d"
        (list-ref inj-name inj) noz inj-particle inj-type inj-material
        inj-n-part (list-ref inj-evap-spec inj)
        (list-ref x 0) (list-ref x 1) (list-ref x 2)
        inj-part-temp (list-ref inj-t-start inj) (list-ref inj-t-stop inj)
        (list-ref v 0) (list-ref v 1) (list-ref v 2)
        inj-vel-mag inj-hoca inj-rad inj-mflow inj-dmin inj-dmax inj-dmean
        inj-spread inj-n-diam))
       )
     )
   (else
    (begin ;; default < 14.5
      (pp "setting injection for 14.0 or earlier\n")
      (ti-menu-load-string (format #f "define/injections/create-injection
       ~a-~d y ~a y ~a ~a ~d n n y \"initDPM::spray\" n ~a y n ~d ~d ~d ~d
       ~d ~d ~d ~d ~d ~d ~d ~d ~d ~d ~d ~d ~d ~d"
       (list-ref inj-name inj) noz inj-particle inj-type inj-material
       inj-n-part (list-ref inj-evap-spec inj)
       (list-ref x 0) (list-ref x 1) (list-ref x 2)
       inj-part-temp (list-ref inj-t-start inj) (list-ref inj-t-stop inj)
       (list-ref v 0) (list-ref v 1) (list-ref v 2)
       inj-vel-mag inj-hoca inj-rad inj-mflow inj-dmin inj-dmax inj-dmean
       inj-spread inj-n-diam))
      )
    )
   ) ;; end cond
  ) ;; end mduc-set-injection
;;
;; ==================== end functions ====================
;;
;; transform the nozzle holes defined in local co-ordinate system to global
;; first
;; build the transformation matrix
(define inj-trans (build-euler-trans inj-phi inj-theta inj-psi))
(display (format #f "transformation matrix: ~a\n" inj-trans))
;;
(define noz-x '())
(define noz-v '())
(define noz-vp '())
(define noz-offset '(0.0 0.0 0.0))
(define noz-x-g '())
(define noz-v-g '())
(do ((i 0 (+ i 1)) ) ((= i (length noz-x0)))
  ;; make sure the vector is a unit vector
  (vec-norm (list-ref noz-v0 i))
  ;; transform the nozzle hole vector to the global co-ordinate system
  (set! noz-v   (mat-vec-prod inj-trans (list-ref noz-v0 i)))
  (set! noz-v-g (append noz-v-g (list noz-v)))
  ;; calculate the nozzle hole location in the global co-ordinate system
  (set! noz-vp (mat-vec-prod inj-trans (list-ref noz-x0 i)))
  (set! noz-x (vec-add inj-zero noz-vp))
  ;; move the point slightly along the unit vector to ensure that it is located
  ;; in the domain
  (set! noz-offset (vec-scale 0.00001 noz-v))
  ;; set the hole location/vector in the global arrays
  (set! noz-x-g (append noz-x-g (list (vec-add noz-x noz-offset))))
  );; end loop over nozzles
;;
;; loop over number of injection events
(do ((k 0 (+ k 1)) ) ((= k (length inj-name)))
  ;; loop over nozzles
  (do ((i 0 (+ i 1)) ) ((= i (length noz-x-g)))
    ;; compute the current angle
    (set! noz-x (list-ref noz-x-g i))
    (set! noz-v (list-ref noz-v-g i))
    (pp (format "nozzle position: ~a\n" noz-x))
    (pp (format "nozzle vector:   ~a\n" noz-v))
    ;; create the injection
    (mduc-set-injection k (+ i 1) noz-x noz-v)
    ;; create the visualization aids
    (if (and (not case-restart) (= k 0))
        (begin
          ;; injection plane
          (ti-menu-load-string (format #f "\nsurface/plane-bounded aa_nozzle_~02d
                       ~d ~d 0.02 ~d ~d 0.02 ~d ~d -2.0 n" 
                       (+ i 1) (list-ref noz-x 0) (list-ref noz-x 1)
                       (* (list-ref noz-v 0) 10) (* (list-ref noz-v 1) 10)
                       (* (list-ref noz-v 0) 10) (* (list-ref noz-v 1) 10)))
          ;; draw a point/line for the nozzle origin/vector
          (draw-noz-vec (+ i 1) noz-x noz-v)
          )
        ) ;; end if
    ) ;; end loop over nozzles
  ) ;; end loop over injections
;;
;; build some planes downstream from the injector
(define inj-z-plane-dist '(5 10 20 30)) ; in mm
(define inj-z-plane-x '())
(define inj-z-plane-v (mat-vec-prod inj-trans '(0 0 1)))
(for-each
 (lambda (z)
   (set! inj-z-plane-x (vec-add inj-zero (vec-scale (* z 0.001) inj-z-plane-v)))
   (display inj-z-plane-x)
   (ti-menu-load-string
    (format #f "\nsurface/plane-point-n-normal aa_inj_z_~02dmm
                ~d ~d ~d ~d ~d ~d" z (list-ref inj-z-plane-x 0)
                (list-ref inj-z-plane-x 1) (list-ref inj-z-plane-x 2) 
                (list-ref inj-z-plane-v 0) (list-ref inj-z-plane-v 1) 
                (list-ref inj-z-plane-v 2) ))
   ); end lambda
   inj-z-plane-dist
 ); end for-each
;;
;; Pressure dependent boiling
(if p-dep-boiling
	(pp "Enable pressure dependent boiling\n")
	(ti-menu-load-string
		(format #f "/define/models/dpm/options> vaporization-options yes")))
;; ----- setting up droplet breakup parameters 
;; Break up always enabled! 
(pp "Enable break up\n")
(ti-menu-load-string
	(format #f "/define/models/dpm/spray-model enable-breakup y"))
;; Choose break up model -> wave break up is default
(if khrt
 	(begin ;TRUE -> Kelvin Helmholtz Rayleigh Taylor breakup
        (pp "KHRT breakup model\n")
		(ti-menu-load-string	
 			(format #f "define/models/dpm/spray-model khrt-model ~d ~d ~d ~d ~d"
        		inj-wave-b0 inj-wave-b1 inj-khrt-c-l inj-khrt-c-tau inj-khrt-c-rt))
	)
	(begin ;FALSE -> Wave breakup
        (pp "Wave breakup model\n")
	(ti-menu-load-string	
 		(format #f "define/models/dpm/spray-model wave-model ~d ~d"
         	inj-wave-b0 inj-wave-b1))
	)
)
;; Collision
(pp "Disable collision model\n")
(ti-menu-load-string
	(format #f "define/models/dpm/spray-model droplet-collision ~a" inj-collision))
;;
;; =========================================================================== ;
