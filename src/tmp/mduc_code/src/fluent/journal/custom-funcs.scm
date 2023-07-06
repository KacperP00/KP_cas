;; ========================================================================== ;;
;;                                                                            ;;
;; Define some custom functions used without Scheme codes.                    ;;
;;                                                                            ;;
;; ========================================================================== ;;
;;
;; remainder for float values less than 1
(define (mduc-mod a b) (- a (* (floor (/ a b)) b) )); end mduc-mod
;;
;; angle conversion 
(define (deg->rad a) (/ (* a 3.141592653589793) 180.0) )
(define (rad->deg r) (/ (* r 180) 3.141592653589793) )
;;
(define (mduc-get-time type)
  (cond ((string=? type "time") (rpgetvar  'flow-time) )
        ((string=? type "cad")  (%rpgetvar 'dynamesh/in-cyn/crank-angle) )
        ((string=? type "time-cad")
         (+ (* (rpgetvar 'flow-time) (* 6 op-speed)) ca-start)
         )
        (else                   (rpgetvar  'time-step) )
        ) ;; end cond
  ) ; end mduc-get-time
;;
(define (mduc-file-suffix type cur-time)
  (if (> 0 cur-time) (set! cur-time (mduc-get-time type)) )
  (cond ( (string=? type "time") (format #f "~.3Es"   cur-time)   ); time
        ( (or (string=? type "cad") (string=? type "time-cad"))
          (begin
            (set! cur-time (+ cur-time 0.001))
            (cond ((< cur-time 10)   (format #f "000~1.0fCAD" cur-time))
                  ((< cur-time 100)  (format #f "00~2.0fCAD"  cur-time))
                  ((< cur-time 1000) (format #f "0~3.0fCAD"  cur-time))
                  (else              (format #f "~4.0fCAD"   cur-time))
                  ) ;; end cond
            )
          ); cad
        ( else                   (format #f "~05d"    cur-time)   ); time-step
        ) ;; end cond
  ); end mduc-file-name
;;
(define (dummy-pass n)
  (define tmp-str "")
  (do ((i 1 (+ i 1))) ((> i n))
    (set! tmp-str (string-append tmp-str "\n "))
    ); end do
  tmp-str
  ); end dummy-pass
;;
(define i-uds-enthalpy 4)
;;
(define (turn-uds-on uds-on?)
  (do ((i 0 (+ i 1))) ((= i (rpgetvar 'uds/n-uds)))
    (if uds-on?
        (ti-menu-load-string (format #f "solve/set/equations uds-~d y" i))
        (ti-menu-load-string (format #f "solve/set/equations uds-~d n" i))
        )
    )
  )
;;
(define (set-bc-wall-hflux walls hflux-on?)
  ;;(define dummy-dpm 7); for reflect boundary conditions
  (define dummy-dpm 2); for wall-film boundary conditions
  (for-each 
   (lambda (name)
     (if hflux-on?
         (begin
           (display (format #f " turning on heat flux for wall: ~a\n" name) )
           (ti-menu-load-string 
            (format #f "\ndefine/boundary-conditions/wall ~a ~a ~a ~a ~a ~a
             y y \"udf\" \"H_FLUX::mduc\"" name (dummy-pass 14)
             (dummy-pass dummy-dpm)
             (dummy-pass (- (length (inquire-species-names)) 1))
             (dummy-pass (rpgetvar 'uds/n-uds))
             (dummy-pass (* 2 i-uds-enthalpy)) ))
           )
         (begin ;else
           (display (format #f " turning off heat flux for wall: ~a\n" name) )
           (ti-menu-load-string 
            (format #f "\ndefine/boundary-conditions/wall ~a ~a ~a ~a ~a ~a
             n 0" name (dummy-pass 14) (dummy-pass dummy-dpm)
             (dummy-pass (- (length (inquire-species-names)) 1))
             (dummy-pass (rpgetvar 'uds/n-uds))
             (dummy-pass (* 2 i-uds-enthalpy)) ))
           )
         )
     ); end lambda
   walls
   ); end for-each
  ); end setbc-wall-hflux
;; (set-bc-wall-hflux '(cyl-head) #t) 
;;
;; ===== some string functions
(define (str-split str ch)
  (let ((len (string-length str)))
    (letrec
        ((split
          (lambda (a b)
            (cond
             ((>= b len) (if (= a b) '() (cons (substring str a b) '())))
             ((char-set-member? (string-ref str b) ch) (if (= a b)
                                                 (split (+ 1 a) (+ 1 b))
                                                 (cons (substring str a b) (split b b))))
             (else (split a (+ 1 b)))))))
      (split 0 0))))
;;
(define (char-set-member? c-check c-set)
  (define found-in-set #f)
  (for-each
   (lambda (cc)
     (if (char=? c-check cc) (set! found-in-set #t) )
   ); end lambda
   c-set
   ); end for-each
  found-in-set
  )
;;
(define (set-bc-wall-pressure outlets wall?)
  ;; (define dummy-dpm 2); for wall-film boundary conditions
  (define back-flow-temp 800)
  (define back-turb-intensity 2)
  (define back-hyd-diam 0.04)
  (for-each 
   (lambda (name)
     (if wall?
         (begin
           (display (format #f " changing to a wall: ~a\n" name) )
           (ti-menu-load-string
            (format #f "\ndefine/boundary-conditions zone-type ~a wall" name) )
           (ti-menu-load-string 
            (format #f "\ndefine/boundary-conditions/wall ~a" name) )
           )
         (begin ;else
           (display (format #f " changing to a pressure outlet: ~a\n" name) )
           (ti-menu-load-string
            (format #f "\ndefine/boundary-conditions
              zone-type ~a pressure-outlet" name) )
           (ti-menu-load-string 
            (format #f "\ndefine/boundary-conditions/pressure-outlet
              ~a y n \"pressure\" \"~a\" n ~d n y n n n y ~d ~d n y n 1.0"
              name name back-flow-temp back-turb-intensity back-hyd-diam) )
           )
         )         
     ); end lambda
   outlets
   ); end for-each
  ); end setbc-wall-hflux
;; ========================================================================== ;;
;;
(define (manifold-check-restart)
  (display "\n----- checking manifold restart -----\n")
  (define manifold-time)
  (define manifold-ca)
  (define tmp-ivo)
  (define tmp-ivc)
  (define man-outlet-name '(dummy))
  ;; get the simulation time
  (set! manifold-time (mduc-get-time fl-freq-type))
  ;; convert to crank angle if necessary
  (if (or (string=? fl-freq-type "cad") (string=? fl-freq-type "time-cad"))
      (set! manifold-ca manifold-time)
      (set! manifold-ca (* manifold-time (* 6 op-speed)))
      )
  ;; output what's going on to the screen
  (display (format #f " Current CAD:   ~d\n" manifold-ca))
  ;; make sure the value is between 0-720
  (set! manifold-ca (mduc-mod manifold-ca ca-period))
  ;;
  (display (format #f " Reference CAD: ~d\n" manifold-ca))
  ;; check whether any valves should be open
  (do ((i 0 (+ i 1) )) ((= i (length manifold-ca-ivo)))
    (list-set! manifold-valve-open i #f) ;; assume closed
    (display (format #f " Checking cylinder ~d\n" (+ i 1)))
    ;; check to see if closing is on next cycle
    (if (< (list-ref manifold-ca-ivc i) (list-ref manifold-ca-ivo i))
        (begin
          (if (> manifold-ca (list-ref manifold-ca-ivo i))
              (begin
                (display "   IVO is on next cycle\n")
                (set! tmp-ivo 0)
                (set! tmp-ivc 720)
                )
              (begin
                (set! tmp-ivo -720)
                (set! tmp-ivc 0)
                )
              )
          )
        (begin
          (set! tmp-ivo 0)
          (set! tmp-ivc 0)
          )
        )
    (display (format #f "~d ~d\n" manifold-ca (+ (list-ref manifold-ca-ivo i)
                                                 tmp-ivo)))
    (display (format #f "~d ~d\n" manifold-ca (+ (list-ref manifold-ca-ivc i)
                                                 tmp-ivc)))
    (display (format #f "~a\n" (list-ref manifold-valve-open i)))
    (if (and (>= manifold-ca (+ (list-ref manifold-ca-ivo i) tmp-ivo))
             (<  manifold-ca (+ (list-ref manifold-ca-ivc i) tmp-ivc)) )
        (begin
          (display (format #f " manifold cylinder ~d: open\n" (+ i 1)))
          (list-set! manifold-valve-open i #t)
          )
        )
    ) ;; end loop over ivo
  ;;
  (display "manifold valve states\n")
  (display manifold-valve-open)
  (display "\n------------------------------\n")
  ); end manifold-check-restart
;; ========================================================================== ;;
;;
(define (manifold-check-valves)
  (display "\n----- checking manifolds -----\n")
  (define manifold-time)
  (define manifold-ca)
  (define tmp-ivc)
  (define man-outlet-name '(dummy))
  ;; get the simulation time
  (set! manifold-time (mduc-get-time fl-freq-type))
  ;; convert to crank angle if necessary
  (if (or (string=? fl-freq-type "cad") (string=? fl-freq-type "time-cad"))
      (set! manifold-ca manifold-time)
      (set! manifold-ca (* manifold-time (* 6 op-speed)))
      )
  ;; output what's going on to the screen
  (display (format #f " Current CAD:   ~d\n" manifold-ca))
  ;; make sure the value is between 0-720
  (set! manifold-ca (mduc-mod manifold-ca ca-period))
  ;;
  (display (format #f " Reference CAD: ~d\n" manifold-ca))
  ;; check whether any valve should be opened
  (do ((i 0 (+ i 1) )) ((= i (length manifold-ca-ivo)))
    (display (format #f " Checking cylinder ~d\n" (+ i 1)))
    ;; check to see if closing is on next cycle
    (if (< (list-ref manifold-ca-ivc i) (list-ref manifold-ca-ivo i))
        (begin
          (display "   IVO is on next cycle\n")
          (set! tmp-ivc 720)
          )
        (begin
          (set! tmp-ivc 0)
          )
        )
    (if (and (and (>= manifold-ca (list-ref manifold-ca-ivo i))
                  (< manifold-ca (+ (list-ref manifold-ca-ivc i) tmp-ivc)) )
             (not (list-ref manifold-valve-open i) ) )
        (begin
          (display (format #f " Opening manifold cylinder ~d\n" (+ i 1)))
          (list-set! man-outlet-name 0 (format #f "cyl0~d_outlet" (+ i 1)))
          (set-bc-wall-pressure man-outlet-name #f)
          (list-set! manifold-valve-open i #t)
          )
        )
    (if (and (and (>= manifold-ca (list-ref manifold-ca-ivc i))
                  (< manifold-ca
                     (- (+ (list-ref manifold-ca-ivo i) 720) tmp-ivc)) )
             (list-ref manifold-valve-open i) )
        (begin
          (display (format #f " Closing manifold cylinder ~d\n" (+ i 1)))
          (list-set! man-outlet-name 0 (format #f "cyl0~d_outlet" (+ i 1)))
          (set-bc-wall-pressure man-outlet-name #t)
          (list-set! manifold-valve-open i #f)
          )
        )
    ) ;; end loop over ivo
  ;;
  (display "\n------------------------------\n")
  ); end manifold-check-valves
;;
;; ***** end function definition
