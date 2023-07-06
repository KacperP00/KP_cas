;; ----- command to comment/uncomment text
(defun mduc-comment-dwim (arg)
  "Comment or uncomment current line or region in a smart way"
  (interactive "*P")
  (require 'newcomment)
  (let (
        (comment-start "#") (comment-end "")
        )
    (comment-dwim arg))
  )

;; ----- define several class of keywords
;; RED
(setq mduc-keywords '("begin" "end") )
;; GREEN
(setq mduc-types
      '(
        "heat-release" "rho" "density" "viscosity" "dry-nox"
        "scalar-dissipation" "specific-heat" "cp" "heat-loss" "qdot"
        "gamma" "convection-coefficient"
        )
      )
;; PURPLE
(setq mduc-constants
      '("none" "from-file" "beta" "beta-pdf"
        "value" "derivative"
        ;; time modes
        "steady" "input" "file-list" "scalar-file" "engine"
        ;; pressure modes
        "isobaric" "isochoric" "constant-pressure" "constant-volume"
        "compression" "implicit" "explicit"
        ;; scalar dissipation rate modes
        "erfc" "erfc-form" "log" "log-form" "equation"
        "linear" "least-squares" "laplace"
        ;; chemistry/soot
        "HMOM" "compiled" 
        ;; solver settings
        "cvode" "ida"
        "gmres" "bicgstab" "tfqmr" "band" "banded"
        ;; engine types
        "woschni"
        "HCCI" "SI" "SACI" "hcci" "si" "saci"
        )
      )
;; PINK
(setq mduc-events
      '("FLAMELET" "ENDFLAMELET" "OxidizerSide" "FuelSide1" "FuelSide2"))
;; BLUE
(setq mduc-functions
      '("title" "author" "date" "mechanism" "number-of-flamelets"
        "solver-mode" "pressure-mode" "mixing-mode" "time-mode"
        "chi-mode" "chi-interp" "chi-scale" "pdf-mode" 
        "diffusion-mode" "convection-mode" "bounds-filter"
        "enthalpy-dimension" "heat-loss" "dns-mode"
        "sundials-solver" "sundials-method" "spatial-order" "energy-equation"
        "absolute-tolerance" "relative-tolerance" "viscosity-update"
        "start-time" "time-step-number" "time-step-size" "final-time" "CFL"
        "time-steps" "output-variables" "output-species" "monitor-species"
        "data-directory" "read-data-file" "read-data" "read-flow-file"
        "write-data-file" "flow-data-file"
        "frequency-type" "write-frequency" "print-frequency"
        "output-var" "output-total-enthalpy" "output-heat-release"
        "output-chi" "output-rho" "output-pdf" "output-cfl" "output-dry-nox"
        "convection-coefficient" "post-process-file" "flow-prefix"
        "file-number-style" "file-number-precision" "file-number-start"
        "file-number-stop" "monitor-file" "error-file" "ensight-directory"
        "MaxGridPoints" "MaxGridPoints2" "Z2-transformation"
        "pressure"
        "Z-ref" "Z0-ref" "Z1-ref" "Z2-ref"
        "Z-mean" "Z0-mean" "Z1-mean" "Z2-mean"
        "Z-variance" "Z0-variance" "Z1-variance" "Z2-variance"
        "chi-ref" "chi0-ref" "chi1-ref" "chi2-ref"
        "strain-rate" "strain-rate-0" "strain-rate-1" "strain-rate-2"
        "Min-enthalpy" "Max-enthalpy"
        "reverse-reaction-rates" "reverse-reaction-debug"
        "read-thermo-data" "read-transport-data" "transport-data-file"
        "soot-model" "soot-moment-order"
        "fuel-species" "fuel-mass-fractions" "fuel-mole-fractions"
        "fuel-lower-heating-value"
        "engine-setup-file" "engine-type" "engine-mass" "engine-bore"
        "engine-stroke" "engine-rod-length" "engine-comp-ratio"
        "engine-volume-tdc" "engine-speed"
        "engine-start-cad" "engine-stop-cad" "engine-step-size"
        "engine-manifold-pressure"
        "intake-valve-open" "intake-valve-close" "intake-port-diameter"
        "intake-stem-diameter" "intake-head-diameter" "intake-seat-diameter"
        "intake-seat-width" "intake-seat-angle" "intake-valve-lift"
        "intake-manifold-pressure" "intake-mass-flow"
        "exhaust-valve-open" "exhaust-valve-close" "exhaust-port-diameter"
        "exhaust-stem-diameter" "exhaust-head-diameter" "exhaust-seat-diameter"
        "exhaust-seat-width" "exhaust-seat-angle" "exhaust-valve-lift"
        "exhaust-manifold-pressure" "exhaust-mass-flow"
        "injection-profile" "injection-timing" "injection-duration"
        "engine-spark-model" "engine-spark-advance"
        "engine-spark-product-species" "engine-spark-product-fraction"
        "engine-spark-diameter"
        "wiebe-a" "wiebe-m" "wiebe-soc" "wiebe-alpha" "wiebe-beta"
        "wiebe-spark-a" "wiebe-spark-m" "wiebe-spark-soc" "wiebe-spark-alpha"
        "wiebe-spark-beta" 
        "engine-burn-rate" "burn-rate-scale"
        "engine-wall-temp"
        "woschni-a" "woschni-b" "woschni-c" "woschni-d" "woschni-c2"
        "woschni-alpha" "woschni-ref-temp" "woschni-ref-pressure"
        "dimensions"
        "gridPoints1" "gridPivot1" "gridPercent1"
        "gridPoints2" "gridPivot2" "gridPercent2"
        )
      )


;; create the regex string for each class of keywords
(setq mduc-keywords-regexp (regexp-opt mduc-keywords 'words))
(setq mduc-type-regexp (regexp-opt mduc-types 'words))
(setq mduc-constant-regexp (regexp-opt mduc-constants 'words))
(setq mduc-event-regexp (regexp-opt mduc-events 'words))
(setq mduc-functions-regexp (regexp-opt mduc-functions 'words))

;; create the list for font-lock.
;; each class of keyword is given a particular face
(setq mduc-font-lock-keywords
  `(
    (,mduc-type-regexp . font-lock-type-face)
    (,mduc-constant-regexp . font-lock-constant-face)
    (,mduc-event-regexp . font-lock-builtin-face)
    (,mduc-functions-regexp . font-lock-function-name-face)
    (,mduc-keywords-regexp . font-lock-keyword-face)
    ;; note: order above matters. ''mduc-keywords-regexp'' goes last because
    ;; otherwise the keyword ''state'' in the function ''state_entry''
    ;; would be highlighted.
))

;; syntax table
(defvar mduc-syntax-table
  (let ((synTable (make-syntax-table)))
    
    ;; bash style comment
    (modify-syntax-entry ?# "< b" synTable)
    (modify-syntax-entry ?\n "> b" synTable)

    ;; bash style comment
    (modify-syntax-entry ?! "< b" synTable)
    (modify-syntax-entry ?\n "> b" synTable)
    
    synTable)
  "Syntax table for `mduc-mode'")

(setq mduc-indent-start '("FLAMELET" "begin"))
(setq mduc-indent-stop '("ENDFLAMELET" "end"))
(setq mduc-indent-start-regexp (regexp-opt mduc-indent-start 'words))
(setq mduc-indent-stop-regexp (regexp-opt mduc-indent-stop 'words))

(defun mduc-indent-line ()
  "Indent current line as MDUC code"
  (interactive)
  (beginning-of-line)

  (if (bobp)
      (indent-line-to 0)
    (let ((not-indented t) cur-indent)
      (if (looking-at "^[ \t]*\\(\}\\|END\\|end\\)") ; Check for rule 2
          (progn
            (save-excursion
              (forward-line -1)
              (setq cur-indent (- (current-indentation) default-tab-width)))
            (if (< cur-indent 0)
                (setq cur-indent 0)))
        (save-excursion 
          (while not-indented
            (forward-line -1)
            (if (looking-at "^[ \t]*\\(\}\\|END\\|end\\)") ; Check for rule 3
                (progn
                  (setq cur-indent (current-indentation))
                  (setq not-indented nil))
              ; Check for rule 4
              (if (looking-at "^[ \t]*\\(\{\\|FLAMELET\\|begin\\)")
                  (progn
                    (setq cur-indent (+ (current-indentation) default-tab-width)
                          )
                    (setq not-indented nil))
                (if (bobp) ; Check for rule 5
                    (setq not-indented nil)))))))
      (if cur-indent
          (indent-line-to cur-indent)
        ;; If we didn't see an indentation hint, then allow no indentation
        (indent-line-to 0))))
  )

;; define the major mode.
(define-derived-mode mduc-mode fundamental-mode
  "mduc-mode is a major mode for editing MDUC input files"
  (set-syntax-table mduc-syntax-table)
  
  (setq font-lock-defaults '(mduc-font-lock-keywords))
  (setq mode-name "mduc")

  ;; modify the keymap
  (define-key mduc-mode-map [remap comment-dwim] 'mduc-comment-dwim)

  (set (make-local-variable 'indent-line-function) 'mduc-indent-line)

  ;; clear memory
  (setq mduc-keywords-regexp nil)
  (setq mduc-types-regexp nil)
  (setq mduc-constants-regexp nil)
  (setq mduc-events-regexp nil)
  (setq mduc-functions-regexp nil)
)

(provide 'mduc-mode)