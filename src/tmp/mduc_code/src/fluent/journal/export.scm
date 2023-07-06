; ============================================================================ ;
; File to set up data to be exported during Fluent run.                        ;
; Requires: case.scm, mduc-cff.scm int-zones, inj-names                        ;
; ============================================================================ ;
;
; solution data (to ensight-gold format)
(if dump-ensight
    (begin ;; if
      (define int-zones '(all_walls))
      (define ens-basic "absolute-pressure temp")
      (define ens-mduc
        "z1-mean z1-variance z2-mean z2-variance z-variance h-mean")
      (define ens-sdm "z1-var-src z2-var-src z-var-src z1-var-diss z2-var-diss z-var-diss fuel-volume-1 fuel-volume-2 fdt-1 fdt-2 fdt xi-1 xi-2 xi alpha-b1 alpha-b2 alpha-b spalding-no")
      (define ens-species "n-c7h16 o2 n2 co2 h2o co")

      (if dump-sdm
          (define ens-vars (format #f "~a ~a ~a ~a" ens-basic ens-mduc
            ens-species ens-sdm))
          (define ens-vars (format #f "~a ~a ~a" ens-basic ens-mduc
            ens-species))
          )

      (ti-menu-load-string (format #f "file/transient-export/ensight-gold-transient ~a/~a/~a- '~a ~a q n y export-standard ~d y" case-dir ens-data file-prefix
                                   int-zones ens-vars dump-freq) )

    )
    (begin ;; else
      (display "\nsolution data not exported\n")
      )
)
;
; particle data (temperature)
;
; NOTE: Right now the particle data export has a bug related to the maximum file name size.
;       Typically, if an absolute path is given, the filename is too long and the particle
;       data filenames will be garbage, as well as eventually causing Fluent to crash.
;       Therefore, dump-particle currently assumes that the Fluent working directory
;       is the case directory (i.e. /case-root/case-name), rather than the case-root.
(if dump-particle
    (begin ;; if
      (define inj-names '(pilot-1 pilot-2 pilot-3 pilot-4 pilot-5 pilot-6
                                  pilot-7 main-1 main-2 main-3 main-4 main-5
                                  main-6 main-7))
      (define inj-pilot '(nozzle-1))
      (define inj-main  '(main-1 main-2 main-3 main-4 main-5 main-6 main-7))

      ; dump the ensight case file (just geometry and velocity)
      (ti-menu-load-string (format #f "file/transient-export/ensight-gold-transient particle/~a- '~a q n y export-standard ~d y" file-prefix int-zones dump-freq))

      ; particle temperature of main/pilot
      (ti-menu-load-string (format #f "file/transient-export/particle-history-data particle-temperature '~a ensight particle/~a- y ~d particle/~a-pilot-temp" inj-pilot file-prefix dump-freq file-prefix))
;      (ti-menu-load-string (format #f "file/transient-export/particle-history-data particle-temperature '~a ensight ~a/particle/~a- y ~d ~a/particle/~a-main-temp" inj-main case-id file-prefix dump-freq case-id file-prefix))
      ; particle diameter of main/pilot
      (ti-menu-load-string (format #f "file/transient-export/particle-history-data particle-diameter '~a ensight particle/~a- y ~d particle/~a-pilot-diam" inj-pilot file-prefix dump-freq file-prefix))
;      (ti-menu-load-string (format #f "file/transient-export/particle-history-data particle-diameter '~a ensight ~a/particle/~a- y ~d ~a/particle/~a-main-diam" inj-main case-id file-prefix dump-freq case-id file-prefix))
      )
    (begin ;; else
      (display "\nparticle data not exported\n")
      )
)