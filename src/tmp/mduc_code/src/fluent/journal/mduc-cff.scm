; ============================================================================ ;
; File defining custom field function definitions based on variables in the    ;
; MDUC library
; ============================================================================ ;
;
(define cff-cmd "define/custom-field-functions define")
(display "---------- Setting up custom field functions ----------\n\n")
; ----- set scalar names -----
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "z1-mean" "uds_0_scalar"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "z1-variance" "uds_1_scalar"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "z2-mean" "uds_2_scalar"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "z2-variance" "uds_3_scalar"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "z-variance" "uds_4_scalar"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "h-mean" "uds_5_scalar"))
;(ti-menu-load-string
; (format #f "~a \"~a\" ~a\n" cff-cmd "fresh-gas" "uds_6_scalar"))
; ----- set memory names -----
; general
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "chi-ref-0" "udm_0"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "chi-ref-1" "udm_1"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "chi-ref-2" "udm_2"))
; sdm components
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "fuel-volume-1" "udm_3"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "fuel-volume-2" "udm_4"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "z1-var-src" "udm_5"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "z2-var-src" "udm_6"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "z-var-src" "udm_7"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "z1-var-diss" "udm_8"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "z2-var-diss" "udm_9"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "z-var-diss" "udm_10"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "fdt-1" "udm_11"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "fdt-2" "udm_12"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "fdt" "udm_13"))
; ----- set computed values -----
; general
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "chi-mean-1" 
         "2*turb_diss_rate/(turb_kinetic_energy+0.0000001)*uds_1_scalar"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "chi-mean-2" 
         "2*turb_diss_rate/(turb_kinetic_energy+0.000001)*uds_3_scalar"))
; sdm components
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "drop-diameter-1"
         "(6/PI*fuel_volume_1)^(1/3)"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "drop-diameter-2"
         "(6/PI*fuel_volume_2)^(1/3)"))
(ti-menu-load-string
 (format #f "~a \"~a\" ~a\n" cff-cmd "drop-diameter"
         "(6/PI*(fuel_volume_1+fuel_volume_2))^(1/3)"))
; swirl
;;(ti-menu-load-string
;; (format #f "~a \"~a\" ~a\n" cff-cmd "x-swirl" "-244*x_coordinate"))
;;(ti-menu-load-string
;; (format #f "~a \"~a\" ~a\n" cff-cmd "y-swirl" "244*y_coordinate"))
