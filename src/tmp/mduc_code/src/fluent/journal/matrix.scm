;; ========================================================================== ;;
;; Define some matrix operations, mostly for injector configuration           ;;
;; file matrix.scm                                                            ;;
;; ========================================================================== ;;
;;
;; add components of two vectors
(define (vec-add a b)
  (define result '())
  (define tmp 0.0)
  (do ((i 0 (+ i 1))) ((= i (length a)))
    (set! tmp (+ (list-ref a i) (list-ref b i)))
    (set! result (append result (list tmp)))
    ) ;; end do
  result
  ) ;; end vec-add
;;
;; scale magnitude of vector
(define (vec-scale sc v)
  (define result '())
  (define sc-tmp 0.0)
  (do ((i 0 (+ i 1))) ((= i (length v)))
    (set! sc-tmp (* sc (list-ref v i)))
    (set! result (append result (list sc-tmp)))
    ) ;; end do
  result
  ) ;; end vec-scale
;;
;; dot product of two vectors
(define (dot-prod a b)
  (define result 0.0)
  (do ((i 0 (+ i 1))) ((= i (length a)))
    (set! result (+ result (* (list-ref a i) (list-ref b i))))
    ) ;; end do
  result
  ) ;; end dot-prod
;;
;; product of matrix (m) and vector (v)
(define (mat-vec-prod m v)
  (define result '())
  (if (not (= (length (list-ref m 0)) (length v)))
      (begin
        (display "matrix vector sizes not consistent")
        )
      (begin
        (do ((i 0 (+ i 1))) ((= i (length v)))
;          (set! result (append result '(0.0)))
;          (list-set! result i (dot-prod (list-ref m i) v))
          (set! result (append result (list (dot-prod (list-ref m i) v))))
          ) ;; end row
        )
      )
  result
  ) ;; end dot-prod
;;
;; transpose matrix m
(define (mat-transpose m)
  (define result '())
  (define row '())
  ;; loop over columns
  (do ((j 0 (+ j 1))) ((= j (length (list-ref m 0))))
    (set! row '())
    ;; loop over rows
    (do ((i 0 (+ i 1))) ((= i (length m)))
      (set! row (append row (list (list-ref (list-ref m i) j))))
      )
    (set! result (append result (list row)))
    )
  result
  ) ;; end mat-transpose
;;
;; product of two matrices
(define (mat-prod m1 m2)
  (define result '())
  (define row '())
  (define m2t (mat-transpose m2))
  (do ((j 0 (+ j 1))) ((= j (length m1)))
    (set! row '())
    (do ((i 0 (+ i 1))) ((= i (length m2t)))
      (set! row (append row (list (dot-prod (list-ref m1 j) (list-ref m2t i)))))
      )
    (set! result (append result (list row)))
    )
  result
  ) ;; end dot-prod
;;
;; normalize the vector a by it's magnitude, i.e. create unit vector
;; NOTE: this operates directly on the a matrix
(define (vec-norm a)
  ;; calculate the magnitude
  (define mag (sqrt (dot-prod a a)))
  (do ((i 0 (+ i 1))) ((= i (length a)))
    (list-set! a i (/ (list-ref a i) mag))
    ) ;; end do
  ) ;; end vec-norm
;;
;; build a transformation matrix by using Euler angles
(define (build-euler-trans phi-a theta-a psi-a)
  (define phi   (deg->rad phi-a))
  (define theta (deg->rad theta-a))
  (define psi   (deg->rad psi-a))
  (define t1  (list (list         (cos phi)   (sin phi)  0.0)
                    (list (* -1.0 (sin phi))  (cos phi)  0.0)
                    (list          0.0         0.0       1.0) ) )
  (define t2  (list (list   1.0            0.0           0.0)
                    (list   0.0           (cos theta)   (sin theta))
                    (list   0.0   (* -1.0 (sin theta))  (cos theta)) ) )
  (define t3  (list (list         (cos psi)   (sin psi)  0.0)
                    (list (* -1.0 (sin psi))  (cos psi)  0.0)
                    (list          0.0         0.0       1.0) ) )
  (define t12 (mat-prod t2 t1))
  (mat-prod t3 t12)
  )
;;
;; ========================================================================== ;;
