(define read-var)
(define wb-setup-ic-case)
(define initial-fluids '())
(define initial-mesh-walls '()); EMD
(define final-fluids '())
(define fluid-ch-original 0)
(define final-mesh-walls '()); EMD
(define wb-ice-report-dir-fluent)
(define wb-ice-report-dir)
(define wb-setup-ic-solver)
(define wb-ice-inlet-outlet-bc)
(define wb-ice-wall-bc)
(define wb-ice-viscous-model-setup)
(define wb-ice-initialization)
(define wb-ice-create-planes)
(define wb-ice-create-default-images-on-planes)
(define wb-ice-create-view-for-planes)
(define wb-ice-report-mesh)
(define wb-ice-report-urf)
(define wb-ice-plot-profile)
(define wb-ice-swirl-data)
(define wb-ice-postproc-monitors-default)
(define wb-ice-postproc-monitors-user)
(define patching-values)
(define ice-panel-read-data)
(define ice-read-data)
(define ice-apply-call-back)
(define wb-ice-driver)
(define wb-ice-solver-setup)
(define wb-ice-urf-setting)
(define set-smoothing-events)
(define define-deforming-zone-event)
(define delete-dynamic-zone-event)
(define define-stationary-zone-event)
(define define-rigidbody-zone-event)
(define read-valve-profile)
(define new-interiors '())
(define initial-interiors '())
;; Engine Parameters
(define valves-check)
(define V-adk-old 1)
(define smoothing-on #f)
(define setup-valve-list '((1 "iname" "iprof" 340 580 0.015 (0 0 1) (0 0 0) #f #f 5)) )
(define inlet-outlet-bc-list)
(define wall-bc-list)
(define rpm 2160)
(define dCA 0.25)
(define crank_radius 105e-3)
(define piston_pin_offset 0)
(define Connecting_rod_length 145e-3)
(define Piston_stroke_cutoff 13e-3)
(define Min_valve_lift 0.35e-3)
(define increment 5)
(define auto_activate_deactivate #f)
(define auto_activate_deactivate_s "#f")
;
(define non_const_URFs #f)
(define non_const_dt #f)
(define non_const_URFs-s "#f")
(define non_const_dt-s "#f")
(define valveStore 1)

(define valve-value)
(define symmetryengine #f)
(define point_on_symmetry (list 0 0 0))
(define symm_normal (list 0 1 0))



;; Valve Parameters
(define no_of_valves 1)
(define V-adk 1)
(define valve_root_name)
(define valve_profile_name)
(define valve_axis)
(define open_valve)
(define close_valve)
(define valve_mar_rad 0.015)
(define valve_origin (list 0e-3 0e-3 0e-3))
(define valve_axis (list 0 0 1))

;;Piston Parameters
(define cylinder_axis_origin (list 0e-3 0e-3 0e-3))
(define cylinder_axis (list 0 0 1))
(define cylinder_radius 46.0e-3)
(define piston_type 1)
(define start-crank-angle 0)
(define insert_angle_power_stroke 30)
(define insert_angle_intake_stroke 40)
(define withbowl #f)
(define hybrid/layering "h")
(define hybrid #t)
(define layering #f)
(define hybrid-s "#t")
(define layering-s "#f")

(define conf/nonconf "c")
(define conformal #t)
(define nonconformal #f)

(define ICDyMeshSetFile "IC-parameters")
(define solver-set-file "ic-solver-set-file")
(define valve-profile-path "WBICValveProfFile.prof")
(define piston-full-prof  "**piston-full**")
(define piston-limit-prof "**piston-limit**")
(define wb-ic-register-post-journal)
(define wb-ice-post-iteration-method)
(define wb-ice-post-iteration-journal #f)
(define wb-ice-pre-iteration-journal #f)
(define wb-ice-user-bc-file #f)
(define wb-ic-write-custom-files)
(define wb-ic-update-time-steps)



(define ic_setup_panel
	(let ((panel)
	     (rpm-w)
 	     (starting-ca)
       (d-ca)
       (crank-radius)
       (piston-pin-offset)
       (connecting-rod-length)
       (piston-stroke-cutoff)
       (min-valve-lift)
	     (valves)
	     (no-of-valves)
	     (valve-root-name)
       (valve-profile-name)
       (valve-axis)
	     (open-valve)
       (close-valve)
	     (valve-mar-rad)
	     (cylinder-radius)
	     (piston-type)
	     (insert-angle)
	     (insert-angle2)
	     (valve-hybrid)
	     (valve-layering)
	     (piston-hybrid)
	     (piston-layering)
	     (row_space_ep)
	     (row_space_piston)
	     (row_space_piston_1)
	     (row_space_piston_2)
	     (row_space_valve)
	     (meshing_strategy)
	     (toggle_hybrid)
	     (toggle_layering)
	     (symmetry-engine)
	     (point-on-symmetry)
	     (symm-normal)
	     (frame_main)
	     (frame_engine_parameters)
	     (frame_meshing_approach)
	     (frame_valve_motion) 
	     (v_frame_valve_axis) (a1) (a2) (a3)
	     (v_frame_valve_origin) (b1) (b2) (b3)
	     (frame_piston_motion)
	     (v_frame_piston_axis) (r1) (r2) (r3)
	     (v_frame_piston_origin) (q1) (q2) (q3)
	     (v_frame_piston_motion) (vframe_piston_pt) (vframe_ia)
	     (vframe_valve_pt1) (vframe_valve_pt2) (vframe_valve_pt3) (comment-1) (comment-2) (vframe7) 
	     (vframe8) (vframe9) (vframe10) (with_bowl)
	     (vframe8-1) (pos1) (pos2) (pos3)
	     (vframe8-2) (sn1) (sn2) (sn3)
	     (increment_)
	     (toggle_read_variable_time)
	     (toggle_read_urf)
	     (toggle_read_activate)
			
	     (layering_setup)
	     (toggle_con)
	     (toggle_non)
	      
	     (apply-button)	
	     )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  START OF UPDATE CALL-BACK FUNCTION
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

       (define close-callback
              (lambda()
									(set! rpm (cx-show-real-entry rpm-w))
									(set! dCA (cx-show-real-entry d-ca))
									(set! crank_radius (cx-show-real-entry crank-radius))
									(set! Connecting_rod_length (cx-show-real-entry connecting-rod-length))
									(set! Piston_stroke_cutoff (cx-show-real-entry  piston-stroke-cutoff))
									(set! Min_valve_lift (cx-show-real-entry min-valve-lift))
			
									(set! auto_activate_deactivate_s (cx-show-toggle-button toggle_read_activate))
											(if auto_activate_deactivate_s
												(set! auto_activate_deactivate #t)
												(set! auto_activate_deactivate #f)
											) 
									(if (cx-show-toggle-button symmetry-engine)
										(set! symmetryengine  #t)
										(set! symmetryengine #f)
									)

								(set! symmetryengine (cx-show-toggle-button symmetry-engine))
								(set! point_on_symmetry (list  (cx-show-real-entry pos1)
												  (cx-show-real-entry pos2)
													(cx-show-real-entry pos3)
     		                      )
													)	
								(set! symm_normal (list  (cx-show-real-entry sn1)
                                      (cx-show-real-entry sn2)
                                      (cx-show-real-entry sn3)
                          	     )
																)
								(if (string-ci=? hybrid/layering "h")
										(set! hybrid #t)
										(set! layering #t)
								)
								(if (string-ci=? hybrid/layering "l")
										(set! hybrid #f)
										(set! layering #f)
								)
								(if (string-ci=? conf/nonconf "c")
										(set! conformal #t)
										(set! nonconformal #t)
								)
								
								(if (string-ci=? conf/nonconf "n")
										(set! conformal #f)
										(set! nonconformal #f)
								)

								
	
							(set! piston_type (cx-show-integer-entry piston-type))
                        
								(if (cx-show-toggle-button toggle_hybrid)
									(begin 
					  		
										(if (= (cx-show-integer-entry piston-type) 1)
				 
											(begin
											;;PISTON CALLBACK
											(set! Piston_stroke_cutoff (cx-show-real-entry piston-stroke-cutoff))
											(set! cylinder_radius (cx-show-real-entry cylinder-radius))
											(set! cylinder_axis (list (cx-show-real-entry r1)
													(cx-show-real-entry r2)
													(cx-show-real-entry r3)))
											(set! cylinder_axis_origin (list (cx-show-real-entry q1)
											  (cx-show-real-entry q2)
											  (cx-show-real-entry q3)))
													;; VALVE CALLBACK

									)
								)
		
								(if (= (cx-show-integer-entry piston-type) 2)
				
										(begin
										;; PISTON CALLBACK
										(set! Piston_stroke_cutoff (cx-show-real-entry piston-stroke-cutoff))
										(set! cylinder_radius (cx-show-real-entry cylinder-radius))
										(set! insert_angle_power_stroke (cx-show-real-entry insert-angle))
										(set! insert_angle_intake_stroke (cx-show-real-entry insert-angle2))
										(set! withbowl (cx-show-toggle-button with_bowl)) 
										(set! cylinder_axis (list (cx-show-real-entry r1)
													(cx-show-real-entry r2)
													(cx-show-real-entry r3)))
										(set! cylinder_axis_origin (list (cx-show-real-entry q1)
												  (cx-show-real-entry q2)
												  (cx-show-real-entry q3)))
														;; VALVE CALLBACK
										)
								)

								(if (= (cx-show-integer-entry piston-type) 3)
				
										(begin
										;; PISTON CALLBACK
										(set! Piston_stroke_cutoff (cx-show-real-entry piston-stroke-cutoff))
										(set! cylinder_radius (cx-show-real-entry cylinder-radius))
										(set! insert_angle_power_stroke (cx-show-real-entry insert-angle))

										(set! cylinder_axis (list (cx-show-real-entry r1)
													(cx-show-real-entry r2)
													(cx-show-real-entry r3)))
										(set! cylinder_axis_origin (list (cx-show-real-entry q1)
												  (cx-show-real-entry q2)
												  (cx-show-real-entry q3)))
										(set! cylinder_radius (cx-show-real-entry cylinder-radius))
										;; VALVE CALLBACK
							
									)
								)
							)
						)
				
						(if (cx-show-toggle-button toggle_layering )

								(begin
								;;PISTON CALLBACK
								;(set! Piston_stroke_cutoff (cx-show-real-entry piston-stroke-cutoff))

								(set! cylinder_axis (list (cx-show-real-entry r1)
											(cx-show-real-entry r2)
											(cx-show-real-entry r3)))
								;; VALVE CALLBACK
								)

						)
							   (set! valveStore (- V-adk 1))
                 (save-data  (- V-adk 1)  (length setup-valve-list)  #f )
					 )
       )


	(define update-callback
		(lambda ()
		(cx-set-real-entry rpm-w rpm)
  	(cx-set-real-entry d-ca dCA)
   	(cx-set-real-entry crank-radius crank_radius)
   	(cx-set-real-entry connecting-rod-length Connecting_rod_length)
   	(cx-set-real-entry min-valve-lift Min_valve_lift)

		(cx-set-toggle-button toggle_read_activate auto_activate_deactivate)
		(cx-set-toggle-button symmetry-engine  symmetryengine) 
		
    (cx-set-real-entry pos1 (car  point_on_symmetry))
	  (cx-set-real-entry pos2 (cadr  point_on_symmetry))
	  (cx-set-real-entry pos3 (caddr  point_on_symmetry))

   	(cx-set-real-entry sn1 (car  symm_normal))
   	(cx-set-real-entry sn2 (cadr  symm_normal))
   	(cx-set-real-entry sn3 (caddr  symm_normal))
		
		(symmetry-on)
		
			
		(if (string-ci=? hybrid/layering "h")
		(begin
		(cx-set-toggle-button toggle_hybrid hybrid)

		(cx-set-integer-entry piston-type piston_type)
		(cx-set-real-entry cylinder-radius cylinder_radius)
		(cx-set-real-entry piston-stroke-cutoff Piston_stroke_cutoff)
		(cx-set-real-entry insert-angle insert_angle_power_stroke)
		(cx-set-real-entry insert-angle2 insert_angle_intake_stroke)
		(cx-set-toggle-button  with_bowl     withbowl)
		(cx-set-real-entry r1 (car cylinder_axis))
   	(cx-set-real-entry r2 (cadr cylinder_axis))
   	(cx-set-real-entry r3 (caddr cylinder_axis))
		(cx-set-real-entry q1 (car cylinder_axis_origin))
   	(cx-set-real-entry q2 (cadr cylinder_axis_origin))
   	(cx-set-real-entry q3 (caddr cylinder_axis_origin))

		)

		(begin
		(cx-set-toggle-button toggle_layering layering)
		(if (string-ci=? conf/nonconf "c")
		(cx-set-toggle-button toggle_con conformal)
		(cx-set-toggle-button toggle_non nonconformal)
		)
		;(cx-set-integer-entry piston-type piston_type)
		
		(cx-set-real-entry r1 (car cylinder_axis))
   	(cx-set-real-entry r2 (cadr cylinder_axis))
   	(cx-set-real-entry r3 (caddr cylinder_axis))
		
		)
		)
		
		
		(check1)

		(cx-set-integer-entry no-of-valves no_of_valves)
                (cx-set-integer-entry valves V-adk)

		
		   (do ((i 0 (+ i 1))) ((= i (length setup-valve-list)))

                        (if (= (list-ref (list-ref setup-valve-list i) 0)   V-adk  )
                             (begin
                               ;(set! V-adk (list-ref(list-ref setup-valve-list i) 0) )
                               (set! valve_root_name (list-ref (list-ref setup-valve-list i) 1) )
                               (set! valve_profile_name (list-ref (list-ref setup-valve-list i) 2) )
                               (set! open_valve (list-ref (list-ref setup-valve-list i) 3) )
                               (set! close_valve (list-ref (list-ref setup-valve-list i) 4) )
                               (set! valve_mar_rad (list-ref (list-ref setup-valve-list i) 5) )
        			                       
	
			       (set! valve_axis (list-ref (list-ref setup-valve-list i) 6) )
                               
			       (set! valve_origin (list-ref (list-ref setup-valve-list i) 7) )
			       (set! non_const_URFs (list-ref (list-ref setup-valve-list i) 8) )
			       (set! non_const_dt (list-ref (list-ref setup-valve-list i) 9))
			       (set! increment (list-ref (list-ref setup-valve-list i ) 10))
                               (update-valve-panel)
			       (check2)
			       (read-valve-profile (list-ref setup-valve-list i))

                             )                
                        )  
                    )
		
		   (if  (and valve-profile-path (file-exists? valve-profile-path))
		     (begin
		     (ti-menu-load-string (format #f "/file/read-profile  ~s " valve-profile-path))
		     ))
		;(check2)
	        (check3)


		)
	)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  End of UPDATE CALL-BACK FUNCTION
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  START OF UPDATE VALVE-PANEL
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	(define update-valve-panel
		(lambda ()
		(cx-set-integer-entry no-of-valves no_of_valves)
		(cx-set-integer-entry valves V-adk)
		
				
		(cx-set-real-entry a1 (car valve_axis))
    (cx-set-real-entry a2 (cadr valve_axis))
    (cx-set-real-entry a3 (caddr valve_axis))

		(cx-set-real-entry b1 (car valve_origin))
    (cx-set-real-entry b2 (cadr valve_origin))
    (cx-set-real-entry b3 (caddr valve_origin))
	        	
	  (cx-set-text-entry valve-root-name valve_root_name)			

	  (cx-set-text-entry valve-profile-name valve_profile_name)
    (cx-set-real-entry open-valve open_valve)
    (cx-set-real-entry close-valve close_valve)
    (cx-set-real-entry valve-mar-rad valve_mar_rad )
			
		(cx-set-toggle-button toggle_read_urf non_const_URFs)
		(cx-set-toggle-button toggle_read_variable_time non_const_dt)
		(cx-set-real-entry increment_ increment)		

		)
	)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  END OF UPDATE VALVE-PANEL
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; CHECK1 Meshing Strategies (Hybrid Vs Layering) (ACTIVATE-CALLBACK)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		(define check1
			(lambda	()

				(if (cx-show-toggle-button toggle_hybrid)
					(begin
					(set! hybrid/layering "h")
					(set! layering #f)
					(set! hybrid #t)
					(cx-hide-item layering_setup)
					(cx-show-item vframe_ia)
					(cx-show-item piston-type)
					(cx-show-item vframe10)
					(cx-show-item v_frame_piston_axis)
					(cx-show-item v_frame_piston_origin)
					(cx-show-item cylinder-radius)
					(cx-show-item piston-stroke-cutoff)

					(if (= (cx-show-integer-entry piston-type) 1)
				 
							(begin
							;;PISTON CALLBACK
							(cx-hide-item with_bowl)
							(cx-hide-item insert-angle)
							(cx-hide-item insert-angle2)
							(cx-hide-item comment-1)
							;; VALVE CALLBACK
							(cx-show-item vframe_valve_pt1)
							(cx-hide-item comment-2)
							(cx-show-item v_frame_valve_axis)
							(cx-show-item v_frame_valve_origin)
							(cx-show-item valve-root-name)
							(cx-show-item valve-profile-name)
							(cx-show-item open-valve)
							(cx-show-item close-valve)
							(cx-show-item valve-mar-rad)
							(cx-show-item toggle_read_urf)
							(cx-show-item toggle_read_variable_time)
							(check2)
							)
						)
						(if (= (cx-show-integer-entry piston-type) 2)
				
							(begin
							;; PISTON CALLBACK
							(cx-show-item with_bowl)
							(cx-show-item insert-angle)
							(cx-show-item insert-angle2)
							(cx-hide-item comment-1)
							;; VALVE CALLBACK
							(cx-show-item vframe_valve_pt1)
							(cx-hide-item comment-2)
							(cx-show-item v_frame_valve_axis)
							(cx-show-item v_frame_valve_origin)
							(cx-show-item valve-root-name)
							(cx-show-item valve-profile-name)
							(cx-show-item open-valve)
							(cx-show-item close-valve)
							(cx-show-item valve-mar-rad)
							(cx-show-item toggle_read_urf)
							(cx-show-item toggle_read_variable_time)
							(check2)
							)
						)

						(if (= (cx-show-integer-entry piston-type) 3)
				
							(begin
							;; PISTON CALLBACK
							(cx-hide-item with_bowl)
							(cx-show-item insert-angle)
							(cx-hide-item insert-angle2)
							(cx-hide-item comment-1)
							;; VALVE CALLBACK
							(cx-show-item vframe_valve_pt1)
							(cx-hide-item comment-2)
							(cx-show-item v_frame_valve_axis)
							(cx-show-item v_frame_valve_origin)
							(cx-show-item valve-root-name)
							(cx-show-item valve-profile-name)
							(cx-show-item open-valve)
							(cx-show-item close-valve)
							(cx-show-item valve-mar-rad)
							(cx-show-item toggle_read_urf)
							(cx-show-item toggle_read_variable_time)
							(check2)
							)
						)

					)
				)
				(if (cx-show-toggle-button toggle_layering)
					(begin
					(set! hybrid/layering "l")
					(set! layering #t)
					(set! hybrid #f)
					
					
					(if (cx-show-toggle-button toggle_con)
					    (begin
						(set! conf/nonconf "c")
						(set! conformal #t)
						(set! nonconformal #f)
					     )
					)
					(if (cx-show-toggle-button toggle_non)
					    (begin
						(set! conf/nonconf "n")
						(set! conformal #f)
						(set! nonconformal #t)
					     )
					)

							(begin
							;;PISTON CALLBACK
						  (cx-hide-item vframe_ia)
							(cx-hide-item piston-type)
							(cx-show-item layering_setup)
							(cx-show-item vframe10)
							(cx-show-item v_frame_piston_axis)
							(cx-hide-item v_frame_piston_origin)
							(cx-hide-item cylinder-radius)
							(cx-hide-item piston-stroke-cutoff)
							(cx-hide-item with_bowl)
							(cx-hide-item insert-angle)
							(cx-hide-item insert-angle2)
							(cx-hide-item comment-1)
							;; VALVE CALLBACK
							(cx-show-item vframe_valve_pt1)
							(cx-hide-item comment-2)
							(cx-hide-item v_frame_valve_axis)
							(cx-hide-item v_frame_valve_origin)
							(cx-show-item valve-root-name)
							(cx-show-item valve-profile-name)
							(cx-show-item open-valve)
							(cx-show-item close-valve)
							(cx-hide-item valve-mar-rad)
							(cx-show-item toggle_read_urf)
							(cx-show-item toggle_read_variable_time)
				
							(check2)
							)
		
					)
				)
			)
		)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; END OF CHECK1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; START OF CHECK2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

			
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; END OF CHECK2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define check2
        (lambda ()
        
        (if (or
            (cx-show-toggle-button toggle_read_urf) 
            (cx-show-toggle-button toggle_read_variable_time)
            )
            (cx-show-item increment_)
            (cx-hide-item increment_)  
        )
    )
    )



(define symmetry-on
	(lambda ()
		(if (cx-show-toggle-button symmetry-engine)
			(begin
			(cx-show-item vframe8-1)
			(cx-show-item vframe8-2)
			)

			(begin
                        (cx-hide-item vframe8-1)
			(cx-hide-item vframe8-2)
			)



		)
	)
   )	

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Start of Check3
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define check3
(lambda ()

(if (cx-show-toggle-button with_bowl)
    (begin
      (set!    withbowl #t) 
    )
    (begin    
      (set!    withbowl #f) 
    )
)



)
)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;End of Check3
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define valve-axis-replace
  (lambda (master pos new-mem)
    (map 
      (lambda (lis)
        (let ((len)(head)(tail))
          (set! len (length lis))
          (set! head (reverse (memq (list-ref (reverse lis) (- len pos)) (reverse lis))))
          (set! tail (memq (list-ref lis pos) lis))
          (set-car! tail new-mem)
          (append head tail)
        )
      ) master
    )
  )
)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; START OF APPLY-CALLBACK FUNCTION
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	(set! ice-apply-call-back
		(lambda()
                        
      (save    V-adk)
			(cx-hide-panel ic_setup_panel)

			(set! rpm (cx-show-real-entry rpm-w))
			(set! dCA (cx-show-real-entry d-ca))
      (set! crank_radius (cx-show-real-entry crank-radius))
      (set! Connecting_rod_length (cx-show-real-entry connecting-rod-length))
      (set! Piston_stroke_cutoff (cx-show-real-entry  piston-stroke-cutoff))
      (set! Min_valve_lift (cx-show-real-entry min-valve-lift))
			
			(set! auto_activate_deactivate_s (cx-show-toggle-button toggle_read_activate))
      				  (if auto_activate_deactivate_s
       				     (set! auto_activate_deactivate #t)
       				     (set! auto_activate_deactivate #f)
     			) 

			(if (cx-show-toggle-button symmetry-engine)

				(set! symmetryengine  #t)
				(set! symmetryengine #f)
      			)

			(set! symmetryengine (cx-show-toggle-button symmetry-engine))
			(set! point_on_symmetry (list  (cx-show-real-entry pos1)
				      (cx-show-real-entry pos2)
			              (cx-show-real-entry pos3)
     		                      )
     			)	

      			(set! symm_normal (list  (cx-show-real-entry sn1)
                                      (cx-show-real-entry sn2)
                                      (cx-show-real-entry sn3)
                          	     )
       			)
			(if (string-ci=? hybrid/layering "h")
				(set! hybrid #t)
				(set! layering #t)
			)
			(if (string-ci=? hybrid/layering "l")
				(set! hybrid #f)
				(set! layering #f)
			)
			(if (string-ci=? conf/nonconf "c")
				(set! conformal #t)
				(set! nonconformal #t)
			)
			(if (string-ci=? conf/nonconf "n")
				(set! conformal #f)
				(set! nonconformal #f)
			)
	

			(set! piston_type (cx-show-integer-entry piston-type))
                        
				(if (cx-show-toggle-button toggle_hybrid)
					(begin 
					  		
						(if (= (cx-show-integer-entry piston-type) 1)
				 
							(begin
							;;PISTON CALLBACK
							(set! Piston_stroke_cutoff (cx-show-real-entry piston-stroke-cutoff))
							(set! cylinder_radius (cx-show-real-entry cylinder-radius))
							(set! cylinder_axis (list (cx-show-real-entry r1)
											(cx-show-real-entry r2)
											(cx-show-real-entry r3)))
							(set! cylinder_axis_origin (list (cx-show-real-entry q1)
											  (cx-show-real-entry q2)
											  (cx-show-real-entry q3)))
													;; VALVE CALLBACK
							(set! valve_axis (list (cx-show-real-entry a1)
										       (cx-show-real-entry a2)
										       (cx-show-real-entry a3)))

							(set! valve_origin (list (cx-show-real-entry b1)
											 (cx-show-real-entry b2)
											 (cx-show-real-entry b3)))
							(set! valve_root_name (cx-show-text-entry valve-root-name))
							(set! valve_profile_name (cx-show-text-entry valve-profile-name))
							(set! open_valve (cx-show-real-entry open-valve))
							(set! close_valve (cx-show-real-entry close-valve))
							(set! valve_mar_rad (cx-show-real-entry valve-mar-rad))
							(set! non_const_URFs (cx-show-toggle-button toggle_read_urf))
							(set! non_const_dt (cx-show-toggle-button toggle_read_variable_time))
							(set! increment (cx-show-real-entry increment_) )

							)
						)
		
						(if (= (cx-show-integer-entry piston-type) 2)
				
							(begin
							;; PISTON CALLBACK
							(set! Piston_stroke_cutoff (cx-show-real-entry piston-stroke-cutoff))
							(set! cylinder_radius (cx-show-real-entry cylinder-radius))
							(set! insert_angle_power_stroke (cx-show-real-entry insert-angle))
							(set! insert_angle_intake_stroke (cx-show-real-entry insert-angle2))
							(Set! withbowl (cx-show-toggle-button with_bowl)) 


							(set! cylinder_axis (list (cx-show-real-entry r1)
											(cx-show-real-entry r2)
											(cx-show-real-entry r3)))
							(set! cylinder_axis_origin (list (cx-show-real-entry q1)
											  (cx-show-real-entry q2)
											  (cx-show-real-entry q3)))
														;; VALVE CALLBACK
							(set! valve_axis (list (cx-show-real-entry a1)
										       (cx-show-real-entry a2)
										       (cx-show-real-entry a3)))

							(set! valve_origin (list (cx-show-real-entry b1)
											 (cx-show-real-entry b2)
											 (cx-show-real-entry b3)))
							(set! valve_root_name (cx-show-text-entry valve-root-name))
							(set! valve_profile_name (cx-show-text-entry valve-profile-name))
							(set! open_valve (cx-show-real-entry open-valve))
							(set! close_valve (cx-show-real-entry close-valve))
							(set! valve_mar_rad (cx-show-real-entry valve-mar-rad))
							(set! non_const_URFs (cx-show-toggle-button toggle_read_urf))
							(set! non_const_dt (cx-show-toggle-button toggle_read_variable_time))
							(set! increment (cx-show-real-entry increment_) )
							)
						)

						(if (= (cx-show-integer-entry piston-type) 3)
				
							(begin
							;; PISTON CALLBACK
							(set! Piston_stroke_cutoff (cx-show-real-entry piston-stroke-cutoff))
							(set! cylinder_radius (cx-show-real-entry cylinder-radius))

							(set! insert_angle_power_stroke (cx-show-real-entry insert-angle))


							(set! cylinder_axis (list (cx-show-real-entry r1)
											(cx-show-real-entry r2)
											(cx-show-real-entry r3)))
							(set! cylinder_axis_origin (list (cx-show-real-entry q1)
											  (cx-show-real-entry q2)
											  (cx-show-real-entry q3)))
							(set! cylinder_radius (cx-show-real-entry cylinder-radius))
							;; VALVE CALLBACK
							
							(set! valve_axis (list (cx-show-real-entry a1)
										       (cx-show-real-entry a2)
										       (cx-show-real-entry a3)))

							(set! valve_origin (list (cx-show-real-entry b1)
											 (cx-show-real-entry b2)
											 (cx-show-real-entry b3)))
							(set! valve_root_name (cx-show-text-entry valve-root-name))
							(set! valve_profile_name (cx-show-text-entry valve-profile-name))
							(set! open_valve (cx-show-real-entry open-valve))
							(set! close_valve (cx-show-real-entry close-valve))
							(set! valve_mar_rad (cx-show-real-entry valve-mar-rad))
							(set! non_const_URFs (cx-show-toggle-button toggle_read_urf))
							(set! non_const_dt (cx-show-toggle-button toggle_read_variable_time))
							(set! increment (cx-show-real-entry increment_) )
							)
						)
					)
				)
				
				(if (cx-show-toggle-button toggle_layering )

							(begin
							;;PISTON CALLBACK
							
							;(set! Piston_stroke_cutoff (cx-show-real-entry piston-stroke-cutoff))

							(set! cylinder_axis (list (cx-show-real-entry r1)
											(cx-show-real-entry r2)
											(cx-show-real-entry r3)))
							;; VALVE CALLBACK
							;(set! valve_axis (list (cx-show-real-entry a1)
							;				(cx-show-real-entry a2)
							;				(cx-show-real-entry a3)))
	
							(set! valve_root_name (cx-show-text-entry valve-root-name))
							(set! valve_profile_name(cx-show-text-entry valve-profile-name))
							(set! open_valve (cx-show-real-entry open-valve))
							(set! close_valve (cx-show-real-entry close-valve))
							(set! non_const_URFs (cx-show-toggle-button toggle_read_urf))
							(set! non_const_dt (cx-show-toggle-button toggle_read_variable_time))
							(set! increment (cx-show-real-entry increment_) )
							)

				)

		       (set! no_of_valves (cx-show-integer-entry no-of-valves))
		       (wb-ice-driver)
                        
			
		

		)
	)

	(set! wb-ice-driver
	  (lambda ()

	    (format  "\n\n\t\t ___________Setting Up Dynamic Mesh _____________\t\n\n")  	
	    	(if (string-ci=? hybrid/layering "l")
				(begin
				(valve-axis-replace setup-valve-list 6 cylinder_axis)
				(valve-axis-replace setup-valve-list 7 '(0.0 0.0 0.0))
				(journal-l)
				)  	
				(journal-h)
			)
	   (format  "\n\n\t\t ___________Dynamic Mesh Setup Done_____________\t\n\n")  


	    )
	  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;END OF APPLY CALL-BACK FUNCTION
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Valves Call Function (ACTIVATE-CALLBACK)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

       (set!  valves-check
               (lambda ()
                    (set! no_of_valves (cx-show-integer-entry no-of-valves))
                    (set! V-adk (cx-show-integer-entry valves))
                (if (> V-adk no_of_valves)
                    (begin
                           (cx-set-integer-entry valves no_of_valves)
                           (set! V-adk no_of_valves)
                           (retrive (- V-adk 1))
                           (set! V-adk-old V-adk)
                    )
                )

                   (if (< (length setup-valve-list) no_of_valves)  
                      (begin
                               (do (( counter (length setup-valve-list) (+ counter 1))) ((=  counter no_of_valves)) 
                                 (set! valveStore counter)
                                 (save-data  1 1  #t )
                               ) 
                      )
                    )
               )
       )


	(define valves-call
		(lambda ()
              
                (set! no_of_valves (cx-show-integer-entry no-of-valves))


                 (if (< (length setup-valve-list) no_of_valves)
                      (begin
                               (do (( counter (length setup-valve-list) (+ counter 1))) ((=  counter no_of_valves))
                                 (set! valveStore counter)
                                 (save-data  1 1  #t )
                               )
                      )
                 )




		(set! V-adk (cx-show-integer-entry valves))	
                (save    V-adk-old)

                (if (> V-adk no_of_valves)
                    (begin
                            (pp "You can not go ahead ")
			    (cx-set-integer-entry valves no_of_valves)
                            (set! V-adk no_of_valves)
                            (set! V-adk-old V-adk)
                            
                    )
                )
                (retrive (- V-adk 1))
        )
	)


(define save
  (lambda (par1)
     (set! valveStore  (- par1 1))
     (save-data  (- par1 1)  (length setup-valve-list)  #f )
  )
)

(define retrive
   (lambda (par1 )

                   (do ((i 0 (+ i 1))) ((= i (length setup-valve-list)))

                        (if (= i  par1 )
                             (begin

                               (set! valve_root_name (list-ref (list-ref setup-valve-list i) 1) )
                               (set! valve_profile_name (list-ref (list-ref setup-valve-list i) 2) )

                               (set! open_valve (list-ref (list-ref setup-valve-list i) 3) )

                               (set! close_valve (list-ref (list-ref setup-valve-list i) 4) )

                               (set! valve_mar_rad (list-ref (list-ref setup-valve-list i) 5) )

                               (set! valve_axis (list-ref (list-ref setup-valve-list i) 6) )
                               (set! valve_origin (list-ref (list-ref setup-valve-list i) 7) )
                               (set! non_const_URFs (list-ref (list-ref setup-valve-list i) 8))
                               (set! non_const_dt (list-ref (list-ref setup-valve-list i ) 9 ))
                               (set! increment (list-ref (list-ref setup-valve-list i) 10))

                               (update-valve-panel)
                               (check2)
                             )
                        )
                    )



   )
)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; START OF SAVE-DATA FUNCTION to list setup-valve-list
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define save-data
(lambda (par1 till bool )
(let
    (  (tmp (list) )   )
(do ((i 0 (+ i 1))) ((= i till ))
   (if  (or
          (= i  par1)   
          bool
        )
   (begin
    (set! tmp (append tmp 
             (list (list
                (+ valveStore 1)
                ;V-adk-old 
                ;(cx-show-integer-entry valves )
	             	(cx-show-text-entry valve-root-name )			
	            	(cx-show-text-entry valve-profile-name )
                (cx-show-real-entry open-valve )
                (cx-show-real-entry close-valve )
                (cx-show-real-entry valve-mar-rad      )    
	        (list			
	          	(cx-show-real-entry a1 (car valve_axis))
          		(cx-show-real-entry a2 (cadr valve_axis))
   	        	(cx-show-real-entry a3 (caddr valve_axis))
                ) 
          (list
	          	(cx-show-real-entry b1 (car valve_origin))
       	    	(cx-show-real-entry b2 (cadr valve_origin))
       	    	(cx-show-real-entry b3 (caddr valve_origin))
	        )
		
	      	(cx-show-toggle-button toggle_read_urf)
      		(cx-show-toggle-button toggle_read_variable_time)
		
		;(if (or 
		;	(cx-show-toggle-button toggle_read_urf)
		;	(cx-show-toggle-button toggle_read_variable_time))
		
	      	(cx-show-real-entry increment_)
		;)
                ))
           ))
   )
   (begin
        (set! tmp (append tmp  (list(list-ref setup-valve-list i))))
   )
   )
)

      (if bool (set! tmp (append setup-valve-list tmp)) )
;commented by nag
;(set! setup-valve-list tmp)
(set! V-adk-old V-adk)
)
)
)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; END OF SAVE DATA FUNCTION 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; START OF READ_DATA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(set! ice-panel-read-data
  (lambda ()

    (ice-read-data)
    (update-callback)
    ))
(set! ice-read-data

(lambda ()

(pp "in-ice-read-data")
(set! rpm                        (read-var ICDyMeshSetFile 'rpm))
(set! dCA                        (read-var ICDyMeshSetFile 'crank_angle_step_size))
(set! crank_radius               (read-var ICDyMeshSetFile 'crank_radius))
(set! piston_pin_offset          (read-var ICDyMeshSetFile 'piston_pin_offset))
(set! Connecting_rod_length      (read-var ICDyMeshSetFile 'connecting_rod_length))
(set! Min_valve_lift             (read-var ICDyMeshSetFile 'min_valve_lift))
(set! auto_activate_deactivate   (read-var ICDyMeshSetFile 'auto_activate_deactivate))
(set! symmetryengine             (read-var ICDyMeshSetFile 'symmetry_engine))
(set! point_on_symmetry          (read-var ICDyMeshSetFile 'point_on_symmetry))
(set! symm_normal                (read-var ICDyMeshSetFile 'symm_normal))
(set! hybrid/layering            (read-var ICDyMeshSetFile 'hybrid/layering))
(set! hybrid                     (read-var ICDyMeshSetFile 'hybrid))
(set! layering                   (read-var ICDyMeshSetFile 'layering))

(set! conf/nonconf               (read-var ICDyMeshSetFile 'conf/nonconf))
(set! conformal                  (read-var ICDyMeshSetFile 'conformal))
(set! nonconformal               (read-var ICDyMeshSetFile 'nonconformal))

(set! piston_type                (read-var ICDyMeshSetFile 'piston_type))
(set! Piston_stroke_cutoff       (read-var ICDyMeshSetFile 'piston_stroke_cutoff))
(set! cylinder_radius            (read-var ICDyMeshSetFile 'cylinder_radius))

(set! withbowl                   (read-var ICDyMeshSetFile 'withbowl))
(if  (read-var ICDyMeshSetFile 'start_crank_angle)
(set! start-crank-angle          (read-var ICDyMeshSetFile 'start_crank_angle))
)
(set! insert_angle_power_stroke  (read-var ICDyMeshSetFile 'insert_angle_power_stroke))
(set! insert_angle_intake_stroke (read-var ICDyMeshSetFile 'insert_angle_intake_stroke))

(set! cylinder_axis              (read-var ICDyMeshSetFile 'cylinder_axis))
(set! cylinder_axis_origin       (read-var ICDyMeshSetFile 'cylinder_axis_origin))
(set! no_of_valves               (read-var ICDyMeshSetFile 'no_of_valves))
(set! V-adk no_of_valves)
(set! setup-valve-list           (read-var ICDyMeshSetFile 'setup_valve_list))
(set! inlet-outlet-bc-list           (read-var ICDyMeshSetFile 'inlet_outlet_bc_list))

(set! patching-values           (read-var ICDyMeshSetFile 'patching_values))
(set! valve-profile-path         (read-var ICDyMeshSetFile 'profile-path))

(if (read-var ICDyMeshSetFile 'piston-full)
  (set! piston-full-prof (read-var ICDyMeshSetFile 'piston-full))
  )

(if (read-var ICDyMeshSetFile 'piston-limit)
  (set! piston-limit-prof (read-var ICDyMeshSetFile 'piston-limit))
  )



;;(set! monitor-def  (read-var solver-set-file 'monitor_def))
(set! V-adk-old no_of_valves)
(if (and valve-profile-path (file-exists? valve-profile-path))
		     (begin
		     (ti-menu-load-string (format #f "/file/read-profile  ~s " valve-profile-path))
		     ))
)
)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; END OF READ_DATA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; RELATED TO READ & WRITE FILE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(set! read-var
  (lambda (par-file var-name)
    (let ((port (open-input-file par-file)))
      (let loop ((x (read port)) (y (read port)))
             (if (eof-object? x) (close-input-port port)
               (if (equal? x var-name)
                   (begin (close-input-port port) y)
                   (loop (read port) (read port))))))
  )
)

(define writedata
(lambda ()
(set! ICDyMeshSetFile (string-append "IC-motion-parameters" ".par"))
(cx-file-dialog "data_files" "*.par" ICDyMeshSetFile applycallback_filedialog_wd #f)
)
)

(define readdata
(lambda ()
(set! ICDyMeshSetFile (string-append "IC-motion-parameters" ".par"))
(cx-file-dialog "data_files" "*.par" ICDyMeshSetFile applycallback_filedialog_rd #f)
)
)

(define applycallback_filedialog_wd
(lambda (widget name)
(set! ICDyMeshSetFile name)

    (if (file-exists? ICDyMeshSetFile)

	(begin

             (let
                (
                (panel)
                (t1)
                (t2) 
                )
                (set! panel (cx-create-panel " Warning " 'modal #t    'ice-apply-call-back  write_data    ))


            (let ((help (cx-get-panel-button-id panel "Help")))
            (cx-hide-item help))


                (set! t1 (cx-create-text panel (string-append "The file \"" ICDyMeshSetFile "\" already exists.")))
                (set! t2 (cx-create-text panel "OK to overwrite? " 'below t1  ))

                (cx-show-panel panel)
             )
         )
      (write_data)
   )
)
)

(define applycallback_filedialog_rd
 (lambda (widget  name )
 (set! ICDyMeshSetFile name)
 (ice-panel-read-data)
 )
)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; START OF WRITE DATA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define write_data
(lambda ()



     (set! valveStore (- V-adk 1))
     (save-data  (- V-adk 1)  (length setup-valve-list)  #f )

     (if (string-ci=? hybrid/layering "l")
       (begin
       (valve-axis-replace setup-valve-list 6 (list (cx-show-real-entry r1) (cx-show-real-entry r2) (cx-show-real-entry r3)))
       (valve-axis-replace setup-valve-list 7 '(0.0 0.0 0.0))
       )
     )

(define outputport (open-output-file  ICDyMeshSetFile))
(let ((write-par))
  (set! write-par (lambda (key val)
                    (display (format #f "~30@a ~a~%" key val) outputport)))

  ; Engine Parameters
  (write-par 'rpm (cx-show-real-entry rpm-w))
  (write-par 'crank_angle_step_size (cx-show-real-entry d-ca))
  (write-par 'crank_radius (cx-show-real-entry crank-radius))
  (write-par 'connecting_rod_length (cx-show-real-entry connecting-rod-length))
  (write-par 'min_valve_lift (cx-show-real-entry min-valve-lift))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  (write-par 'auto_activate_deactivate (cx-show-toggle-button toggle_read_activate))
  (write-par 'symmetry_engine (cx-show-toggle-button symmetry-engine))
  (write-par 'point_on_symmetry (list (cx-show-real-entry pos1)
                                      (cx-show-real-entry pos2)
                                      (cx-show-real-entry pos3)))
  (write-par 'symm_normal (list (cx-show-real-entry sn1)
                                (cx-show-real-entry sn2)
                                (cx-show-real-entry sn3)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; Piston Parameters
  (write-par 'hybrid/layering (format #f "~s" hybrid/layering))
  (write-par 'hybrid (cx-show-toggle-button toggle_hybrid))
  (write-par 'layering (cx-show-toggle-button toggle_layering))
  (write-par 'conf/nonconf (format #f "~s" conf/nonconf))
  (write-par 'conformal conformal)
  (write-par 'nonconformal nonconformal)
  (write-par 'piston_type (cx-show-integer-entry piston-type))
  (write-par 'piston_stroke_cutoff (cx-show-real-entry  piston-stroke-cutoff))
  (write-par 'cylinder_radius (cx-show-real-entry cylinder-radius))

  (write-par 'withbowl (cx-show-toggle-button with_bowl))
  (write-par 'insert_angle_power_stroke (cx-show-real-entry insert-angle))
  (write-par 'insert_angle_intake_stroke (cx-show-real-entry insert-angle2))

  (write-par 'cylinder_axis (list (cx-show-real-entry r1)
                                  (cx-show-real-entry r2)
                                  (cx-show-real-entry r3)))

  (write-par 'cylinder_axis_origin (list (cx-show-real-entry q1)
                                         (cx-show-real-entry q2)
                                         (cx-show-real-entry q3)))

  ; Valves Parameters
  (write-par 'no_of_valves (cx-show-integer-entry no-of-valves))

  (display (format #f "~30@a (~%" 'setup_valve_list) outputport)
  (map (lambda (x)
         (display (format #f "~@33a~a " "(" (list-ref x 0)) outputport)
         (display (format #f "~s " (list-ref x 1)) outputport)
         (display (format #f "~s " (list-ref x 2)) outputport)
         (display (format #f "~a " (list-ref x 3)) outputport)
         (display (format #f "~a " (list-ref x 4)) outputport)
         (display (format #f "~a " (list-ref x 5)) outputport)
         (display (format #f "~a " (list-ref x 6)) outputport)
         (display (format #f "~a " (list-ref x 7)) outputport)

         (display (format #f "~a " (list-ref x 8)) outputport)
         (display (format #f "~a " (list-ref x 9)) outputport)
         (display (format #f "~a)~%" (list-ref x 10)) outputport)
         ;(display (format #f "~a " (cx-show-toggle-button toggle_read_urf)) outputport)
         ;(display (format #f "~a " (cx-show-toggle-button toggle_read_variable_time)) outputport)
         ;(display (format #f "~a)~%" (cx-show-real-entry increment_)) outputport)
       ) setup-valve-list)
  (display (format #f "~@32a" ")") outputport)
)
(close-output-port outputport)
)
)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; END OF WRITE_DATA FUNCTION
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	(lambda()
		(if (not panel)
		(begin
		   (set! ic_setup_panel  (cx-create-panel "IC Setup" 'update-callback update-callback 'close-callback close-callback  ))

			  (let ((help (cx-get-panel-button-id ic_setup_panel "Help")))
          		  (cx-hide-item help))
  			  (cx-add-panel-button ic_setup_panel "OK"  ice-apply-call-back)
   			  (cx-add-panel-button ic_setup_panel "Read" readdata)
    			  (cx-add-panel-button ic_setup_panel "Write" writedata)
		

		  (set! frame_main (cx-create-frame ic_setup_panel ""   'tabbed #t))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	                      
                                     ;tab1                       Engine Parameters
                   
		    (set! frame_engine_parameters (cx-create-frame frame_main "Engine Parameters " 'right-of 0 'left-of 0))
		   

		  	(set! row_space_ep (cx-create-text frame_engine_parameters "" 'below 0 ))  
		   	(set! rpm-w (cx-create-real-entry frame_engine_parameters
								"Crank Shaft Speed (rpm)"
								'below row_space_ep 'right-of 0
								'value 2000 'width 20
								))
		   	(set! d-ca (cx-create-real-entry frame_engine_parameters
		                         			"Crank Angle Step Size (deg)" 
                                       				'below rpm-w  'right-of 0
                                          			'value 0.25 'width 20 
                                        			))
      
              	 	(set! crank-radius (cx-create-real-entry frame_engine_parameters 
                                                   		"Piston Stroke (m)"
                                                 		'below d-ca 'right-of 0
                                                		'value 90e-3 'width 20
                                             		        ))

            	  	(set! connecting-rod-length (cx-create-real-entry frame_engine_parameters
                                                             	"Connecting Rod Length (m)"
                                                             	'below crank-radius 'right-of 0
                                                             	'value 145e-3 'width 20
                                                             	))
             
            	  	(set! min-valve-lift (cx-create-real-entry frame_engine_parameters
                                                      		"Minimum Valve Lift (m)"
                                                      		'below connecting-rod-length 'right-of 0
                                                       		'value 0.1e-3 'width 20
                                                       		))
			(set! vframe7 (cx-create-frame frame_engine_parameters
                                                "   "
                                                'below 0 'right-of d-ca
                                                ))
			(set! toggle_read_activate (cx-create-toggle-button vframe7 "Auto Fluid Zone Activate/Deactivate"
                                                           'below 0
                                                           'right-of 0
                                                           ))
			(cx-enable-item toggle_read_activate #f)

			(set! vframe8 (cx-create-frame frame_engine_parameters
                                                "      "
                                                'below vframe7 'right-of d-ca
                                                ))

			(set! symmetry-engine (cx-create-toggle-button vframe8  "Symmetry Engine " 'activate-callback symmetry-on))


			(set! vframe8-1 (cx-create-frame vframe8 "Point on Symmetry" 'below symmetry-engine  'right-of 0))
			(set! vframe8-2 (cx-create-frame vframe8 "Symmetry Normal" 'below symmetry-engine  'right-of vframe8-1  ))

            (set! pos1 (cx-create-real-entry vframe8-1
                                           "X (m)"
                                           'below 0  'right-of 0
                                           'value 0.0
                                           ))

            (set! pos2 (cx-create-real-entry vframe8-1
                                           "Y (m)"
                                           'below pos1  'right-of 0
                                           'value 0
                                           ))

            (set! pos3 (cx-create-real-entry vframe8-1
                                           "Z (m)"
                                           'below pos2 'right-of 0
                                           'value  1
                                           ))

            (set! sn1 (cx-create-real-entry vframe8-2
                                           "X"
                                           'below  0  'right-of 0
                                           'value 0.0
                                           ))

            (set! sn2 (cx-create-real-entry vframe8-2
                                           "Y"
                                           'below sn1  'right-of 0
                                           'value 0
                                           ))

            (set! sn3 (cx-create-real-entry vframe8-2
                                           "Z"
                                           'below sn2 'right-of 0
                                           'value  1
                                           ))


	(cx-hide-item vframe8-1)
	(cx-hide-item vframe8-2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
		                        ;tab2                       Piston Motion Setup
                                               
		        (set! frame_piston_motion (cx-create-frame frame_main "Piston Motion Setup "  'right-of 0 'left-of 0))
			
	      		(set! meshing_strategy (cx-create-button-box frame_piston_motion "Meshing Strategy" 'below 0 'right-of 0 'radio-mode #t ))
						
		      	(set! toggle_hybrid (cx-create-toggle-button meshing_strategy "Hybrid Approach" 'activate-callback check1))
		      	(set! toggle_layering (cx-create-toggle-button meshing_strategy "Layering Approach" 'activate-callback check1))
			
		      	(set! vframe_piston_pt (cx-create-frame frame_piston_motion "" 'below meshing_strategy 'right-of 0))

		      	(set! layering_setup (cx-create-button-box vframe_piston_pt "Layering Setup Type" 'below 0 'right-of 0 'radio-mode #t ))
		
		      	(set! toggle_con (cx-create-toggle-button layering_setup "Conformal Setup    " 'activate-callback check1))
		      	(set! toggle_non (cx-create-toggle-button layering_setup "Non-conformal Setup" 'activate-callback check1))
	; initial state
		      	(cx-set-toggle-button toggle_con #t)
			
		      	(set! vframe_ia (cx-create-frame vframe_piston_pt "" 'below 0 'right-of 0 'border #f ))
			
		      	(set! piston-type (cx-create-integer-entry vframe_ia "Piston Type     " 'label-position 'left 'below 0 'right-of 0 'activate-callback check1 
		          'maximum 3 'minimum 1 'width 5 ))
		      	(set! comment-1 (cx-create-text vframe_ia "For Layering, This Piston Type Not Available" 'below piston-type 'right-of 0))
			
		      	(set! with_bowl (cx-create-toggle-button  vframe_ia  "With Bowl" 'below 0  'right-of piston-type  'activate-callback check3))

		      	(set! piston-stroke-cutoff (cx-create-real-entry vframe_ia "Piston Stroke Cutoff (m)" 'below piston-type 'right-of 0 'value 20e-3 'width 20))
		      	(set! cylinder-radius (cx-create-real-entry vframe_ia "Cylinder Radius" 'below piston-stroke-cutoff 'right-of 0 'value 0 'width 20))
		      	(set! insert-angle (cx-create-real-entry vframe_ia "Power Stroke"   'below cylinder-radius 'right-of 0 'value 0 'width 20))
		      	(set! insert-angle2 (cx-create-real-entry vframe_ia "Intake Stroke" 'below insert-angle 'right-of 0 'value 0 'width 20))
			
		      	(set! vframe10 (cx-create-frame frame_piston_motion "" 'below meshing_strategy 'right-of vframe_piston_pt))			

		      	(set! v_frame_piston_axis (cx-create-frame vframe10
					                        		"Cylinder Axis Direction"
						                          'below 0
						                         	'right-of 0
						))
				        (set! r1 (cx-create-real-entry v_frame_piston_axis
                                         		"X"
                                        		'below 0 'right-of 0
                                         		'value 0
                                           		))
            
         	  	  (set! r2 (cx-create-real-entry v_frame_piston_axis
                                           		"Y"
                                          		'below r1 'right-of 0
                                           		'value 0
                                           		))
            
        		  	(set! r3 (cx-create-real-entry v_frame_piston_axis
                                           		"Z"
                                           		'below r2 'right-of 0
                                           		'value 1
                                           		))

			       (set! v_frame_piston_origin (cx-create-frame vframe10 
                                        		"Cylinder Axis Origin "
                                        		'below 0 
                        	               		'right-of v_frame_piston_axis
                                	     		))
            
           			(set! q1 (cx-create-real-entry v_frame_piston_origin 
                                         		"X (m)"
                                        		'below 0 'right-of 0
                                       		'value 0.0
                                      		))
           
        		   	(set! q2 (cx-create-real-entry v_frame_piston_origin
                                          		"Y (m)"
                                          		'below q1 'right-of 0
                                          		'value 0 
                                          		))
            
          			(set! q3 (cx-create-real-entry v_frame_piston_origin 
                                          		"Z (m)"
                                          		'below q2 'right-of 0
                                          		'value  1
                                          		))
			
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

                                      ;tab3                      Valve Motion Setup                                                  
			
		        (set! frame_valve_motion (cx-create-frame frame_main "Valve Motion Setup " 'right-of 0 'left-of 0))
			
		      	(set! comment-2 (cx-create-text frame_valve_motion " For Layering, This Piston Type Not Available" 'below 0 'right-of 0))

		      	(set! row_space_valve (cx-create-text frame_valve_motion "" 'below 0 'right-of 0))

		      	(set! vframe_valve_pt1 (cx-create-frame frame_valve_motion
					                      		""
						                      	'below row_space_valve
						))
		       	(set! no-of-valves (cx-create-integer-entry vframe_valve_pt1 
					          		 "Number of Valves " 'label-position
							           'left 'below 0 'right-of 0 
                         'activate-callback valves-check
                         'value 1 'maximum 10 'minimum 1   'width 4 
							))

			       (set! valves (cx-create-integer-entry vframe_valve_pt1
						          	 "Valve Number       " 'label-position
							           'left 'below no-of-valves 'right-of 0 
							           'activate-callback valves-call 
							           'value 1  'maximum 10 'minimum 1  'width 4
							))
			
		        	(set! vframe_valve_pt2 (cx-create-frame frame_valve_motion
						                    	""
				                    			'below vframe_valve_pt1
							))

	         		(set! valve-root-name (cx-create-text-entry vframe_valve_pt2 "Valve Name" 'value "iname"
				                     			'below 0 'right-of 0 'width 22
							))

			        (set! valve-profile-name (cx-create-text-entry vframe_valve_pt2
							                   "Valve Profile Name" 'value "iprof" 
							                   'below valve-root-name
						                    'right-of 0 'width 22
							))
			
		         	(set! open-valve (cx-create-real-entry vframe_valve_pt2
						                  	" Open Valve (deg)"
						                  	'below valve-profile-name 'right-of 0
				                   			'value 340 'width 22
							))

		         	(set! close-valve (cx-create-real-entry vframe_valve_pt2 
						                  	" Close Valve (deg)"
							                  'below open-valve 'right-of 0
							                  'value 580 'width 22
							))
			
		         	(set! valve-mar-rad (cx-create-real-entry vframe_valve_pt2 
						                  	"Valve Margin Radius (m)"
						                  	'below close-valve 'right-of 0
						                    'value 0.015 'width 22
							))
		        	(set! vframe_valve_pt3 (cx-create-frame frame_valve_motion
						                  	""
						                  	'below vframe_valve_pt1
						                  	'right-of vframe_valve_pt2
							))
		        	(set! toggle_read_urf (cx-create-toggle-button vframe_valve_pt3  "Variable URF" 
                                                           'below 0            
                                                           'right-of 0
                                                           'activate-callback check2))

          		(set! toggle_read_variable_time (cx-create-toggle-button vframe_valve_pt3  "Variable Crank Angle Step Size" 
                                                           'below toggle_read_urf
                                                           'right-of 0
                                                           'activate-callback check2))

		         	(set!  increment_ (cx-create-real-entry vframe_valve_pt3
                                                            "Duration (deg)"
                                                            'below toggle_read_variable_time
                                                            'right-of 0
                                                            'value 5
                                                            'activate-callback check2 ))

			        (set! v_frame_valve_axis (cx-create-frame frame_valve_motion
						                       	"Valve Axis Direction  "
						                      	'below vframe_valve_pt3
						                      	'right-of vframe_valve_pt2
							))
				
				      (set! a1 (cx-create-real-entry v_frame_valve_axis
                                          		 "X"
                                          		 'below 0 'right-of 0
                                          		 'value 0
                                         		 ))
            
         		  (set! a2 (cx-create-real-entry v_frame_valve_axis
                                           		"Y"
                                           		'below a1 'right-of 0
                                           		'value 0
                                           		))
            
        			 (set! a3 (cx-create-real-entry v_frame_valve_axis
                                           		"Z"
                                           		'below a2 'right-of 0
                                           		'value 1
                                           		))
			
		          (set! v_frame_valve_origin (cx-create-frame frame_valve_motion 
						                  	"Valve Axis Origin  "
						                  	'below vframe_valve_pt3
						                    'right-of v_frame_valve_axis
							))
			       	(set! b1 (cx-create-real-entry v_frame_valve_origin
					                		"X (m)"
					    	            	'below 0 'right-of 0
					    	            	'value 0.0
					    ))
				
			       	(set! b2 (cx-create-real-entry v_frame_valve_origin
					                		"Y (m)"
					    	            	'below b1 'right-of 0
					                		'value 0.0
					    ))
				
		      		(set! b3 (cx-create-real-entry v_frame_valve_origin
					                 		"Z (m)"
					    	            	'below b2 'right-of 0
					    	             	'value 0.0
					   	))
								
		)
		)
(cx-show-panel ic_setup_panel)
(check1)

	)
       
     

)
)
;(cx-add-item "User-Defined"
;	     "In-Cylinder Mesh Motion Setup..."
;	     #f
;	     #f
;	     (lambda () (and (cx-client?) ))
;	     ic_setup_panel
;)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; JOURNAL FOR VALVE AND PISTON SETUP
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define journal-l
 (lambda()
 (define valve-value 0)
 (define ib_cell_height)
 (define ob_cell_height)
 (define no_invalve 2)
 (define no_exvalve 2)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Change length units, otherwise the limit setup will not be correct because SI is assumed for limits
;
(ti-menu-load-string "/def/units/length m")
(ti-menu-load-string "/def/units/temperature k")
(ti-menu-load-string "/def/units/pressure pascal")
;
; Solver setup
;
(ti-menu-load-string "/def/mod/unsteady-1st-order yes")
(ti-menu-load-string "/def/material/cc air , y ideal-gas y piecewise-polynomial ,,,,,,,,,,,,,,,,,,,,,,, ,,,,,")
(ti-menu-load-string "/dis/set/hardcopy/driver tiff /dis/set/hardcopy/color/color")
(ti-menu-load-string "/def/mod/viscous/ke-standard y")
(ti-menu-load-string "/sol/mon/residual/con-c 0.1 ,,,,, 2e-3")
(ti-menu-load-string "/sol/set/dis-s/press 14")
(ti-menu-load-string "/sol/set/p-v-coupling 22 /sol/set/p-v-controls 0 1 yes")
(ti-menu-load-string "/sol/set/un-rel/press 0.5 /sol/set/un-rel/k 0.4 /sol/set/un-rel/epsilon 0.4")
(ti-menu-load-string "/sol/set/limits 800 , 60 ,,,,")
(rpsetvar 'partition/encapsulate-si? #f)

(define pre-event-list)
(define cur-event-list)
(define event01)
(define event02)
(define event03)
(define event04)
(define event05)
(define event06)
(define event07)
(define event08)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Zone NAME CHECK (Comman zones for all valve)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
(define mylist '())
(define mylist-conformal '())
(define mylist-nonconformal '())



	    (format  "\n\n\t\t ___________Validating Zone Names _____________\t\n\n")  
(set! mylist (list 
			(string-append "fluid-" "piston" ) ;fluid-piston from fluid-bowl
			(string-append "fluid-ch-" "lower")	
			(string-append "piston") ;piston from bowl
			(string-append "intf-piston" "-ch")
			(string-append "intf-piston" "-bowl")
			)
)

(set! mylist-nonconformal (list 
				(string-append "intf-deck-fluid-ch-lower")
				(string-append "intf-deck-fluid-ch-upper")
			  )
)



;(if (string-ci=? conf/nonconf "n") 
;	(set! mylist (append mylist 
;			(list (string-append "intf-deck-fluid-ch-lower")
;			      (string-append "intf-deck-fluid-ch-upper")
;				)
;			)
;
;	)	
;)




(define mylist-prev)
(do ((valve-value 0 (+ valve-value 1))) ((= valve-value (length setup-valve-list)))


(set! mylist-prev (list
	(string-append "fluid-ch-" (list-ref (list-ref setup-valve-list valve-value) 1) )

	(string-append "fluid-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib")
	(string-append "fluid-" (list-ref (list-ref setup-valve-list valve-value) 1) "-vlayer")
	
	(string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ib")
	(string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ob-quad")
	(string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ob-tri")
	(string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch-lower")
	(string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch-upper")
	(string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch-valve")
	(string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-vlayer")
	(string-append "intf-piston-" (list-ref (list-ref setup-valve-list valve-value) 1) )
	(string-append (list-ref (list-ref setup-valve-list valve-value) 1) "-ch" )
	(string-append (list-ref (list-ref setup-valve-list valve-value) 1) "-ib" )
	(string-append (list-ref (list-ref setup-valve-list valve-value) 1) "-ob" )
	(string-append (list-ref (list-ref setup-valve-list valve-value) 1) "-seat")
	
	
	
	)
)
;(if (string-ci=? conf/nonconf "c")
;	(set! mylist-prev (append mylist-prev (list 
;						(string-append "int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib" )
;						(string-append "int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob" )
;					      )
;			  )
;	)
;)

	(set! mylist-conformal (append mylist-conformal (list 
						(string-append "int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib" )
						(string-append "int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob" )
					      )
			  )
	)
		
;(if (string-ci=? conf/nonconf "n")
;	(set! mylist-prev (append mylist-prev (list 
;						(string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ib")
;						(string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-port")
;						(string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-port")
;						(string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-vlayer")
;					      )
;			   )
;	)
;)

	(set! mylist-nonconformal (append mylist-nonconformal (list 
						(string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ib")
						(string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-port")
						(string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-port")
						(string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-vlayer")
					      )
			   )
	)


(set! mylist (append mylist mylist-prev))

)

( if (string-ci=? conf/nonconf "c")
	(set! mylist (append mylist mylist-conformal))
	
	(set! mylist (append mylist mylist-nonconformal))
	
)	





(define all_surfaces (inquire-surface-names))
  (set!  all_surfaces (append all_surfaces (inquire-zone-names)))
  (define surface_counter 0)
  (define final_list (list))
  (define matched #f)
(for-each
          (LAMBDA (surfacename)
					 (if surfacename
						(begin      
	         (set! matched #f)
		  (do
			    (
			    (i 0 (+ i 1))
			    )
			    (( = i (length  all_surfaces)))

                            (if (string-ci=?  surfacename  (symbol->string  (list-ref all_surfaces i)))
                            (begin 
                               (set! surface_counter (+ surface_counter 1))
                               (set! matched #t)
                            )
                            ) 
                          
                           
                   )


                    (if (not matched)
                    (set! final_list (append final_list (list    surfacename     )))
                    )
						)
           )    
				 )
mylist)

(if (= 0 (length final_list))

(begin
	    (format  "\n\n\t\t ___________Zone Names Validation PASS _____________\t\n\n")  

; Change bc
;
	    (format  "\n\n\t\t ___________Setting the Zone Types _____________\t\n\n")  
;
(set! valve-value 0)

(do ((valve-value 0 (+ valve-value 1))) ((= valve-value (length setup-valve-list)))

(if (string-ci=? conf/nonconf "c")
	(begin
		(ti-menu-load-string (string-append  "/def/bc/zt int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib interior"))
		(ti-menu-load-string (string-append  "/def/bc/zt int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob interior"))
	)
)
;
(ti-menu-load-string (string-append  "/def/bc/zt intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ib interface"))
(ti-menu-load-string (string-append  "/def/bc/zt intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ob-quad interface"))
(ti-menu-load-string (string-append  "/def/bc/zt intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ob-tri interface"))
(ti-menu-load-string (string-append  "/def/bc/zt intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch-upper interface"))
(ti-menu-load-string (string-append  "/def/bc/zt intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch-lower interface"))
(ti-menu-load-string (string-append  "/def/bc/zt intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-vlayer wall"))
(ti-menu-load-string (string-append  "/def/bc/zt intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch-valve interface"))
;
(ti-menu-load-string (string-append  "/def/bc/zt intf-piston-" (list-ref (list-ref setup-valve-list valve-value) 1) " interface"))

(if (string-ci=? conf/nonconf "n")
	(begin

(ti-menu-load-string (string-append  "/def/bc/zt intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ib interface"))
(ti-menu-load-string (string-append  "/def/bc/zt intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-port interface"))
(ti-menu-load-string (string-append  "/def/bc/zt intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-port interface"))
(ti-menu-load-string (string-append  "/def/bc/zt intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-vlayer interface"))

	)
)


 (format  "\n\n\t\t ___________ Creating Grid Interfaces _____________\t\n\n")  
 
;
; Create non-conformal interface
;
(rpsetvar 'nonconformal/virtual-polygons? #t)
(ti-menu-load-string (string-append  "/def/gi/create intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ib () intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ob-tri intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ob-quad () no"))
(ti-menu-load-string (string-append  "/def/gi/create " (list-ref (list-ref setup-valve-list valve-value) 1) " intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch-valve () intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch-upper intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch-lower () no"))


(if (string-ci=? conf/nonconf "n")
	
	(begin

		(ti-menu-load-string (string-append  "/def/gi/create intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ib () intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-port () no"))
		(ti-menu-load-string (string-append  "/def/gi/create intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-port () intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-vlayer () no"))

	)

)

 (format  "\n\n\t\t ___________ Setting Dynamic Mesh Controls _____________\t\n\n")  
;
; Dynamic mesh parameters setup
;
(ti-menu-load-string (string-append "/def/dynamic-mesh/dynamic-mesh y y " (number->string rpm) " 0 720"))
(ti-menu-load-string (string-append "/def/dynamic-mesh/controls/in-cylinder-parameter/crank-angle-step " (number->string dCA)))

(if start-crank-angle
(rpsetvar 'dynamesh/in-cyn/crank-start-angle start-crank-angle)
(rpsetvar 'dynamesh/in-cyn/crank-start-angle 0)
)

(ic-update-flow-time-from-angle ( rpgetvar 'dynamesh/in-cyn/crank-start-angle))
(ti-menu-load-string (string-append "/def/dynamic-mesh/controls/in-cylinder-parameter/piston-data " (number->string crank_radius) " " (number->string Connecting_rod_length)" " (number->string piston_pin_offset)))
(ti-menu-load-string (string-append "/def/dynamic-mesh/controls/in-cylinder-parameter/piston-stroke-cutoff " (number->string Piston_stroke_cutoff)))
(ti-menu-load-string (string-append "/def/dynamic-mesh/controls/in-cylinder-parameter/min-lift " (number->string Min_valve_lift)))
;
(ti-menu-load-string "/def/dynamic-mesh/controls/layering y")
(ti-menu-load-string "/def/dynamic-mesh/controls/smoothing y")
(ti-menu-load-string "/def/dynamic-mesh/controls/layering-parameter/constant-height no")
(ti-menu-load-string "/def/dynamic-mesh/controls/layering-parameter/collapse-factor 0.2")
(ti-menu-load-string "/def/dynamic-mesh/controls/remeshing n")
;
; Rigid body dynamic zone setup
;

 (format  "\n\n\t\t ___________ Setting Dynamic Mesh Zones _____________\t\n\n")  
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create fluid-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib rigid-body " (list-ref (list-ref setup-valve-list valve-value) 2) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 0)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 1)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 2))))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create fluid-" (list-ref (list-ref setup-valve-list valve-value) 1) "-vlayer rigid-body " (list-ref (list-ref setup-valve-list valve-value) 2) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 0)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 1)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 2))))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create " (list-ref (list-ref setup-valve-list valve-value) 1) "-ib rigid-body " (list-ref (list-ref setup-valve-list valve-value) 2) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 0)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 1)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 2)) " ,, 0"))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create " (list-ref (list-ref setup-valve-list valve-value) 1) "-ob rigid-body " (list-ref (list-ref setup-valve-list valve-value) 2) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 0)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 1)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 2)) " ,, 0"))
;(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create " (list-ref (list-ref setup-valve-list valve-value) 1) "-ch rigid-body " (list-ref (list-ref setup-valve-list valve-value) 2) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 0)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 1)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 2)) " ,, " (number->string (+ (* 0.6 (car (%get-min-max-length-scale-on-zone (thread-name->id  (string-append (list-ref (list-ref setup-valve-list valve-value) 1) "-ch")) #f))) (* 0.6 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  (string-append (list-ref (list-ref setup-valve-list valve-value) 1) "-ch")) #f)))))) ))
;setting invalve-ch new layer height as 2*vlayer height for smoothing-on as false otherwise set it to 1mm
(set! smoothing-on (cadr (assoc 'smoothing (list-ref setup-valve-list valve-value))))
(if smoothing-on (ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create " (list-ref (list-ref setup-valve-list valve-value) 1) "-ch rigid-body " (list-ref (list-ref setup-valve-list valve-value) 2) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 0)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 1)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 2)) " ,, 0.001" ))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create " (list-ref (list-ref setup-valve-list valve-value) 1) "-ch rigid-body " (list-ref (list-ref setup-valve-list valve-value) 2) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 0)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 1)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 2)) " ,, " (number->string (+ (* 0.6 (car (%get-min-max-length-scale-on-zone (thread-name->id  (string-append (list-ref (list-ref setup-valve-list valve-value) 1) "-ch")) #f))) (* 0.6 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  (string-append (list-ref (list-ref setup-valve-list valve-value) 1) "-ch")) #f)))))) ))
)
;
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create fluid-ch-" (list-ref (list-ref setup-valve-list valve-value) 1) " rigid-body " piston-full-prof "  " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 0)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 1)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 2))))
;
; Stationary dynamic zone setup
;

(if (string-ci=? conf/nonconf "c")
	(begin

		(set! ib_cell_height  (+ (* 0.5 (car (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib")) #f))) (* 0.5 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib")) #f))))))
		(set! ob_cell_height  (+ (* 0.25 (car (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob")) #f))) (* 0.25 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob")) #f))))))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib stationary " (string-append "fluid-" (list-ref (list-ref setup-valve-list valve-value) 1)  "-ib") " ," (number->string ib_cell_height) " ,, 0"))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob stationary " (string-append "fluid-" (list-ref (list-ref setup-valve-list valve-value) 1)  "-vlayer") " ," (number->string ob_cell_height) " ,, 0"))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create " (list-ref (list-ref setup-valve-list valve-value) 1) "-seat stationary ,, " (number->string ob_cell_height)))
	
	)

)

(if (string-ci=? conf/nonconf "n")
	(begin
	
		(set! ib_cell_height  (+ (* 0.5 (car (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ib")) #f))) (* 0.5 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ib")) #f))))))
		(set! ob_cell_height  (+ (* 0.25 (car (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-vlayer")) #f))) (* 0.25 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-vlayer")) #f))))))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ib stationary " (string-append "fluid-" (list-ref (list-ref setup-valve-list valve-value) 1)  "-ib") " ," (number->string ib_cell_height)))
(set! smoothing-on #f)

(if (assoc 'smoothing (list-ref setup-valve-list valve-value)) (set! smoothing-on (cadr (assoc 'smoothing (list-ref setup-valve-list valve-value)))))

(if smoothing-on 
  (begin
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-vlayer stationary " (string-append "fluid-" (list-ref (list-ref setup-valve-list valve-value) 1)  "-vlayer") " , 0 " ))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create " (list-ref (list-ref setup-valve-list valve-value) 1) "-seat stationary ,, 0 " ))
)
(begin
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-vlayer stationary " (string-append "fluid-" (list-ref (list-ref setup-valve-list valve-value) 1)  "-vlayer") " , " (number->string ob_cell_height) ))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create " (list-ref (list-ref setup-valve-list valve-value) 1) "-seat stationary ,, " (number->string ob_cell_height) ))

  )

)

)
)

 (format  "\n\n\t\t ___________ Setting Dynamic Mesh Events _____________\t\n\n")  

;
; Dynamic mesh events setup
;
;
; Open and close valve
;
(set! event01 (list '(active . #t) (cons 'name (string-append "open-" (list-ref (list-ref setup-valve-list valve-value) 1))) (cons 'start  (list-ref (list-ref setup-valve-list valve-value) 3) ) (list 'action 2 (list-ref (list-ref setup-valve-list valve-value) 1) (list (thread-name->id (string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch-valve")) (thread-name->id (string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-vlayer"))) (list (thread-name->id (string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch-lower")) (thread-name->id (string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch-upper"))) 0 0)))

(set! event02 (list '(active . #t) (cons 'name (string-append "close-" (list-ref (list-ref setup-valve-list valve-value) 1))) (cons 'start  (list-ref (list-ref setup-valve-list valve-value) 4) ) (list 'action 2 (list-ref (list-ref setup-valve-list valve-value) 1) (list (thread-name->id (string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch-valve"))) (list (thread-name->id (string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch-lower")) (thread-name->id (string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch-upper"))) 0 0)))
;
; Events setup
;
(set! cur-event-list (list event01 event02)) 
;

; Append current events to previously set up events if any
;
(set! pre-event-list (rpgetvar 'dynamesh/in-cyn/event-list))

(rpsetvar 'dynamesh/in-cyn/event-list (append pre-event-list cur-event-list))
;
;



;
; Events for non-constant URFs and dt 
;
;

(if (list-ref (list-ref setup-valve-list valve-value) 8)
(begin
; Change URFs
;
(set! event01 (list '(active . #t) (cons 'name (string-append "reduce-urf-due-to-open-" (list-ref(list-ref setup-valve-list valve-value) 1))) (cons 'start  (- (list-ref (list-ref setup-valve-list valve-value) 3) 1) ) (list 'action 12 (list '(1 epsilon 0.2) '(1 k 0.2) '(1 mom 0.4) '(1 pressure 0.2)))))
(set! event02 (list '(active . #t) (cons 'name (string-append "increase-urf-due-to-open-" (list-ref(list-ref setup-valve-list valve-value) 1))) (cons 'start  (+ (list-ref (list-ref setup-valve-list valve-value) 3)  increment)) (list 'action 12 (list '(1 epsilon 0.4) '(1 k 0.4) '(1 mom 0.7) '(1 pressure 0.5)))))
(set! event03 (list '(active . #t) (cons 'name (string-append "reduce-urf-due-to-close-" (list-ref(list-ref setup-valve-list valve-value) 1))) (cons 'start (- (list-ref (list-ref setup-valve-list valve-value) 4) 1) ) (list 'action 12 (list '(1 epsilon 0.2) '(1 k 0.2) '(1 mom 0.4) '(1 pressure 0.2)))))
(set! event04 (list '(active . #t) (cons 'name (string-append "increase-urf-due-to-close-" (list-ref(list-ref setup-valve-list valve-value) 1))) (cons 'start  (+ (list-ref (list-ref setup-valve-list valve-value) 4)  increment)) (list 'action 12 (list '(1 epsilon 0.4) '(1 k 0.4) '(1 mom 0.7) '(1 pressure 0.5)))))    

;(set! pre-event-list (rpgetvar 'dynamesh/in-cyn/event-list))

(set! cur-event-list (append cur-event-list  (list event01 event02 event03 event04)))

;(rpsetvar 'dynamesh/in-cyn/event-list (append pre-event-list cur-event-list))

)
)
;



(if (list-ref (list-ref setup-valve-list valve-value) 9)

; Change dt
;
(begin
(set! event05 (list '(active . #t) (cons 'name (string-append "reduce-dt-due-to-open-" (list-ref(list-ref setup-valve-list valve-value) 1))) (cons 'start (list-ref (list-ref setup-valve-list valve-value) 3) ) (list 'action 5 (* 0.5 dCA))))
(set! event06 (list '(active . #t) (cons 'name (string-append "increase-dt-due-to-open-" (list-ref(list-ref setup-valve-list valve-value) 1))) (cons 'start (+ (list-ref (list-ref setup-valve-list valve-value) 3)  increment)) (list 'action 5 dCA)))    
(set! event07 (list '(active . #t) (cons 'name (string-append "reduce-dt-due-to-close-" (list-ref(list-ref setup-valve-list valve-value) 1))) (cons 'start (list-ref (list-ref setup-valve-list valve-value) 4) ) (list 'action 5 (* 0.5 dCA))))
(set! event08 (list '(active . #t) (cons 'name (string-append "increase-dt-due-to-close-" (list-ref(list-ref setup-valve-list valve-value) 1))) (cons 'start (+ (list-ref (list-ref setup-valve-list valve-value) 4)  increment)) (list 'action 5 dCA)))
;
;(set! pre-event-list (rpgetvar 'dynamesh/in-cyn/event-list))

(set! cur-event-list (append cur-event-list  (list event05 event06 event07 event08)))

;(rpsetvar 'dynamesh/in-cyn/event-list (append pre-event-list cur-event-list))

)
)

; Find the existing event list
;
;(set! pre-event-list (rpgetvar 'dynamesh/in-cyn/event-list))
;
; Set the events
;
(rpsetvar 'dynamesh/in-cyn/event-list (append pre-event-list cur-event-list))


;;set the smoothing events
 
(format "\n valve value =~d setup valve list =~s" valve-value setup-valve-list) 
(set-smoothing-events (list-ref setup-valve-list valve-value))

)

;activate spring smoothing on all shapes

(ti-menu-load-string (format #f "/define/dynamic-mesh/controls/smoothing-parameters/spring-on-all-shapes? y"))




;===============================================================

;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Piston setup
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;

 (format  "\n\n\t\t ___________  Creating Grid Interfaces for Piston _____________\t\n\n")  
(define min_scale)
(define max_scale)
(define ave_scale)
(define piston_layer_height)
;
; Change length unit to meter
;
(ti-menu-load-string "/def/units/length m")
;
; Change bc
;
(ti-menu-load-string (string-append  "/def/bc/zt intf-piston-ch interface"))
(ti-menu-load-string (string-append  "/def/bc/zt intf-piston-bowl interface"))

(if (string-ci=? conf/nonconf "n")

	(begin

	(ti-menu-load-string (string-append  "/def/bc/zt intf-deck-fluid-ch-lower interface"))
	(ti-menu-load-string (string-append  "/def/bc/zt intf-deck-fluid-ch-upper interface"))
	
	)

)

;
; Define interface for the piston 
;

(define for-piston-intf "")
(do
    ((i 0 (+ i 1)))
    ((= i (length setup-valve-list)))
  (begin
    (set! for-piston-intf (string-append for-piston-intf (string-append "intf-piston-" (list-ref (list-ref setup-valve-list i) 1))))
    (set! for-piston-intf (string-append for-piston-intf " "))
    )
)
(set! for-piston-intf (string-append "/def/gi/create intf-piston intf-piston-bowl () intf-piston-ch " for-piston-intf "() no"))
(ti-menu-load-string for-piston-intf)

(if (string-ci=? conf/nonconf "n")

(ti-menu-load-string (string-append  "/def/gi/create intf-deck intf-deck-fluid-ch-upper () intf-deck-fluid-ch-lower () no"))

)
;
; Rigid body dynamic zone setup
;



(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create fluid-piston rigid-body  " piston-full-prof "  " (number->string (list-ref cylinder_axis 0)) " " (number->string (list-ref cylinder_axis 1)) " " (number->string (list-ref cylinder_axis 2))))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create piston rigid-body " piston-full-prof "  " (number->string (list-ref cylinder_axis 0)) " " (number->string (list-ref cylinder_axis 1)) " " (number->string (list-ref cylinder_axis 2)) " ,, 0")) ;piston from bowl

(set! ave_scale (+ (* 0.5 (car (%get-min-max-length-scale-on-zone (thread-name->id  "intf-piston-ch") #f))) (* 0.5 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  "intf-piston-ch") #f))))))
;(set! piston_layer_height (* 1.2 ave_scale))
;setting fluid-ch-lower new layer height to 2*vlayer layer height
(if smoothing-on (set! piston_layer_height 0.001) (set! piston_layer_height (* 2 ob_cell_height)))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create intf-piston-ch rigid-body  " piston-full-prof "  " (number->string (list-ref cylinder_axis 0)) " " (number->string (list-ref cylinder_axis 1)) " " (number->string (list-ref cylinder_axis 2)) " ,, " (number->string piston_layer_height)))
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;;
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create fluid-ch-lower rigid-body  " piston-full-prof "  "  (number->string (list-ref cylinder_axis 0)) " " (number->string (list-ref cylinder_axis 1)) " " (number->string (list-ref cylinder_axis 2))))

(set! ave_scale  (caddr (%get-min-max-length-scale-on-zone (thread-name->id  "fluid-ch-lower") #f))) 
;(set! piston_layer_height (* 1.2 ave_scale))
;setting fluid-ch-lower new layer height to 2*vlayer layer height
(if smoothing-on (set! piston_layer_height 0.001) (set! piston_layer_height (* 2 ob_cell_height)))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create intf-deck-fluid-ch-lower stationary , , "  (number->string piston_layer_height)))



)

(begin
(pp (format #f "~%"))
(pp "The following zones did not match")
(pp final_list)
(if (equal? final_list mylist-conformal)
	(format "~% Error in Layering Setup: \"Conformal Setup\" selected for Non-conformal Setup mesh ~%")
)

(if (equal? final_list mylist-nonconformal)
	(format "~% Error in Layering Setup: \"Non-conformal Setup\" selected for Conformal Setup mesh ~%")
)


)

)

)


)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Hybrid Setup Journal
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; JOURNAL FOR VALVE AND PISTON SETUP HYBRID APPROACH
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 (define journal-h
 (lambda()
 (define valve-value 0)
 (define ib_cell_height)
 (define ob_cell_height)
 (define pre-event-list)
 (define cur-event-list)
 (define event01)
 (define event02)
 (define event03)
 (define event04)
 (define event05)
 (define event06)
 (define event07)
 (define event08)

 ; Change length units, otherwise the limit setup will not be correct because SI is assumed for limits
;
(ti-menu-load-string "/def/units/length m")
(ti-menu-load-string "/def/units/temperature k")
(ti-menu-load-string "/def/units/pressure pascal")
;
; Solver setup
;
(ti-menu-load-string "/def/mod/unsteady-1st-order yes")
(ti-menu-load-string "/def/material/cc air , y ideal-gas y piecewise-polynomial ,,,,,,,,,,,,,,,,,,,,,,, ,,,,,")
(ti-menu-load-string "/dis/set/hardcopy/driver tiff /dis/set/hardcopy/color/color")
(ti-menu-load-string "/def/mod/viscous/ke-standard y")
(ti-menu-load-string "/sol/mon/residual/con-c 0.1 ,,,,, 2e-3")
(ti-menu-load-string "/sol/set/dis-s/press 14")
(ti-menu-load-string "/sol/set/p-v-coupling 22 /sol/set/p-v-controls 0 1 yes")
(ti-menu-load-string "/sol/set/un-rel/press 0.5 /sol/set/un-rel/k 0.4 /sol/set/un-rel/epsilon 0.4")
(ti-menu-load-string "/sol/set/limits 800 , 60 ,,,,")
(rpsetvar 'partition/encapsulate-si? #f)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Zone NAME CHECK ( Comman zones for all valves)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define mylist '())

	    (format  "\n\n\t\t ___________Validating Zone Names _____________\t\n\n")  
(if symmetryengine
	(begin
		(set! mylist (list 
				(string-append "fluid-" "ch" )
				(string-append "cyl-" "tri")
				(string-append "symm-" "cyl-" "tri")
	 		      )
		)	

	)

	(begin 
		(set! mylist (list 
				(string-append "fluid-" "ch" )
				(string-append "cyl-" "tri")
	       		     )

		)

	)

)

	(if (= piston_type 1)
	
		(set! mylist (append mylist (list (string-append "fluid" "-piston" "-layer")
						  (string-append "int" "-piston")
						  (string-append "piston")
					    )
			     )
		)
	)
	
	(if (= piston_type 2)
		(if withbowl	
		(set! mylist (append mylist (list (string-append "piston")
						  (string-append "bowl")
						  (string-append "fluid-" "piston") ;piston from bowl
						  (string-append "intf-" "bowl")
					    )
			     )
		)
		(set! mylist (append mylist (list (string-append "piston")
					    )
			     )
		)

		)
	)

	(if (= piston_type 3)
	
		(set! mylist (append mylist (list (string-append "piston")
					    )
			     )
		)
	)
	
		
(define mylist-prev)

(do ((valve-value 0 (+ valve-value 1))) ((= valve-value (length setup-valve-list)))
	
	(set! mylist-prev(list
		
					(string-append "fluid-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib")
	    	;(string-append "fluid-" (list-ref (list-ref setup-valve-list valve-value) 1) "-port")
      		(string-append "fluid-" (list-ref (list-ref setup-valve-list valve-value) 1) "-vlayer")
      		;(string-append "int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib")
      		;(string-append "int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob")
		;-------------------------
		(string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib" "-fluid" "-ib")
		(string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib" "-fluid" "-port")
		(string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob" "-fluid" "-port")
		(string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob" "-fluid" "-vlayer")
      		;-------------------------
		(string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib" "-fluid" "-ib")
      		(string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib" "-fluid" "-ob" "-quad")
      		(string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib" "-fluid" "-ob" "-tri")
      		(string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob" "-fluid" "-ch")
      		(string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob" "-fluid" "-vlayer") ;vlayer
      		(string-append (list-ref (list-ref setup-valve-list valve-value) 1) "-seat")
      		(string-append (list-ref (list-ref setup-valve-list valve-value) 1) "-ch")
      		(string-append (list-ref (list-ref setup-valve-list valve-value) 1) "-ib")
      		(string-append (list-ref (list-ref setup-valve-list valve-value) 1) "-ob")
		)
	)
(set! mylist (append mylist mylist-prev))
)

(define all_surfaces (inquire-surface-names))  ;this gives the list of surfaces in a given case)
(set!  all_surfaces (append all_surfaces (inquire-zone-names))) ;this gives the list of zones in a given case)
(define surface_counter 0)
(define final_list (list))
(define matched #f)
(for-each
             (LAMBDA (surfacename)
                  (set! matched #f)
		  (do
			    (
			    (i 0 (+ i 1))
			    )
			    (( = i (length  all_surfaces)))

                            (if (string-ci=?  surfacename  (symbol->string  (list-ref all_surfaces i)))
                            (begin 
                               (set! surface_counter (+ surface_counter 1))
                               (set! matched #t)
                            )
                            ) 
                          
                           
                   )


                    (if (not matched)
                    (set! final_list (append final_list (list    surfacename     )))
                    )
              )
mylist)


(if (= 0 (length final_list))

(begin

 (format  "\n\n\t\t ___________Zone Names Validation PASS _____________\t\n\n")  

; Change bc
;
	    (format  "\n\n\t\t ___________Setting the Zone Types _____________\t\n\n")  
	    
(set! valve-value 0)
(do ((valve-value 0 (+ valve-value 1))) ((= valve-value (length setup-valve-list)))

; Change bc
;
;(ti-menu-load-string (string-append  "/def/bc/zt int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib interior"))
;(ti-menu-load-string (string-append  "/def/bc/zt int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob interior"))
;-------------------------
(ti-menu-load-string (string-append  "/def/bc/zt intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ib interface"))
(ti-menu-load-string (string-append  "/def/bc/zt intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-port interface"))
(ti-menu-load-string (string-append  "/def/bc/zt intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-port interface"))
(ti-menu-load-string (string-append  "/def/bc/zt intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-vlayer interface"))
;-------------------------
(ti-menu-load-string (string-append  "/def/bc/zt intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ib interface"))
(ti-menu-load-string (string-append  "/def/bc/zt intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ob-quad interface"))
(ti-menu-load-string (string-append  "/def/bc/zt intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ob-tri interface"))
(ti-menu-load-string (string-append  "/def/bc/zt intf-" (list-ref (list-ref setup-valve-list valve-value) 1)"-ob-fluid-ch wall"))
(ti-menu-load-string (string-append  "/def/bc/zt intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-vlayer wall")) ;vlayer

 (format  "\n\n\t\t ___________ Creating Grid Interfaces _____________\t\n\n")  
 
; Create non-conformal interface
;
(rpsetvar 'nonconformal/virtual-polygons? #t)
(ti-menu-load-string (string-append  "/def/gi/create intf-" (list-ref (list-ref setup-valve-list valve-value) 1) " intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ib () intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ob-quad intf-" (list-ref(list-ref setup-valve-list valve-value) 1) "-ib-fluid-ob-tri () no "))
;-------------------------
(ti-menu-load-string (string-append  "/def/gi/create intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ib () intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-port () no"))

(ti-menu-load-string (string-append  "/def/gi/create intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-port () intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-vlayer () no"))
;-------------------------

(format  "\n\n\t\t ___________ Setting Dynamic Mesh Controls _____________\t\n\n")  
; Dynamic mesh parameters setup
;
(ti-menu-load-string (string-append "/def/dynamic-mesh/dynamic-mesh y y " (number->string rpm) " 0 720"))
(ti-menu-load-string (string-append "/def/dynamic-mesh/controls/in-cylinder-parameter/crank-angle-step " (number->string dCA)))

(if start-crank-angle
(rpsetvar 'dynamesh/in-cyn/crank-start-angle start-crank-angle)
(rpsetvar 'dynamesh/in-cyn/crank-start-angle 0)
)

(ic-update-flow-time-from-angle ( rpgetvar 'dynamesh/in-cyn/crank-start-angle))
(ti-menu-load-string (string-append "/def/dynamic-mesh/controls/in-cylinder-parameter/piston-data " (number->string crank_radius) " " (number->string Connecting_rod_length)" " (number->string piston_pin_offset)))
(ti-menu-load-string (string-append "/def/dynamic-mesh/controls/in-cylinder-parameter/piston-stroke-cutoff " (number->string Piston_stroke_cutoff)))
(ti-menu-load-string (string-append "/def/dynamic-mesh/controls/in-cylinder-parameter/min-lift " (number->string Min_valve_lift)))
;
(ti-menu-load-string "/def/dynamic-mesh/controls/layering y")
(ti-menu-load-string "/def/dynamic-mesh/controls/layering-parameter/constant-height no")
(ti-menu-load-string "/def/dynamic-mesh/controls/layering-parameter/collapse-factor 0.2")
(ti-menu-load-string "/def/dynamic-mesh/controls/remeshing y")
(ti-menu-load-string (string-append "/def/dynamic-mesh/controls/remeshing-parameter/length-min " (number->string (* 0.4 (caddr (%get-min-max-length-scale-on-zone (thread-name->id  "fluid-ch") #f))))))
;(ti-menu-load-string (string-append "/def/dynamic-mesh/controls/remeshing-parameter/length-min " (number->string (* 0.95 (car (%get-min-max-length-scale-on-zone (thread-name->id  "fluid-ch") #f))))))
(ti-menu-load-string (string-append "/def/dynamic-mesh/controls/remeshing-parameter/length-max "  (number->string (* 1.4 (caddr (%get-min-max-length-scale-on-zone (thread-name->id  "fluid-ch") #f))))))
;(ti-menu-load-string (string-append "/def/dynamic-mesh/controls/remeshing-parameter/length-max "  (number->string (* 1.1 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  "fluid-ch") #f)))))))
;(ti-menu-load-string (string-append "/def/dynamic-mesh/controls/remeshing-parameter/size-remesh-interval 5"))
(ti-menu-load-string (string-append "/def/dynamic-mesh/controls/remeshing-parameter/size-remesh-interval 1"))
(rpsetvar 'dynamesh/remesh/methods '(#t #f #f #t #f #f))
;
; Rigid body dynamic zone setup
;
 (format  "\n\n\t\t ___________ Setting Dynamic Mesh Zones _____________\t\n\n") 

(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create fluid-" (list-ref (list-ref setup-valve-list valve-value) 1) "-vlayer rigid-body " (list-ref (list-ref setup-valve-list valve-value) 2) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 0)) " " (number->string (list-ref(list-ref (list-ref setup-valve-list valve-value) 6) 1)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 2))))

(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create fluid-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib rigid-body " (list-ref (list-ref setup-valve-list valve-value) 2) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 0)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 1)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 2))))

(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create " (list-ref (list-ref setup-valve-list valve-value) 1) "-ib rigid-body " (list-ref (list-ref setup-valve-list valve-value) 2) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 0)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 1)) " " (number->string (list-ref(list-ref (list-ref setup-valve-list valve-value) 6) 2)) " ,, 0"))

(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create " (list-ref (list-ref setup-valve-list valve-value) 1) "-ob rigid-body " (list-ref (list-ref setup-valve-list valve-value) 2) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 0)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 1)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 2)) " ,, 0"))

(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create " (list-ref (list-ref setup-valve-list valve-value) 1) "-ch rigid-body " (list-ref (list-ref setup-valve-list valve-value) 2) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 0)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 1)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 2)) " ,, " (number->string (+ (* 0.5 (car (%get-min-max-length-scale-on-zone (thread-name->id  (string-append (list-ref (list-ref setup-valve-list valve-value) 1) "-ch")) #f))) (* 0.5 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  (string-append (list-ref (list-ref setup-valve-list valve-value) 1) "-ch")) #f)))))) ))


; Stationary dynamic zone setup
;
(define ib_cell_height)
(define ob_cell_height 0.5e-3)


;(set! ib_cell_height  (+ (* 0.5 (car (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "int-" (list-ref (list-ref setup-valve-list valve-value) 1)"-ib")) #f))) (* 0.5 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "int-" (list-ref (list-ref setup-valve-list valve-value) 1)"-ib")) #f))))))
;(set! ob_cell_height  (+ (* 0.25 (car (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "int-" (list-ref (list-ref setup-valve-list valve-value) 1)"-ob")) #f))) (* 0.25 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "int-" (list-ref (list-ref setup-valve-list valve-value) 1)"-ob")) #f))))))
;(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create int-" (list-ref (list-ref setup-valve-list valve-value) 1)"-ib stationary " (string-append "fluid-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib") " ," (number->string ib_cell_height) " ,, 0"))
;(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob stationary " (string-append "fluid-" (list-ref (list-ref setup-valve-list valve-value) 1)  "-vlayer") " ," (number->string ob_cell_height) " ,, 0"))
;-------------------------

(define smoothing-on #f)

(set! ib_cell_height  (+ (* 0.5 (car (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ib")) #f))) (* 0.5 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ib")) #f))))))
(set! ob_cell_height  (+ (* 0.25 (car (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-vlayer")) #f))) (* 0.25 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-vlayer")) #f))))))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ib-fluid-ib stationary " (string-append "fluid-" (list-ref (list-ref setup-valve-list valve-value) 1)  "-ib") " ," (number->string ib_cell_height)))


(if (assoc 'smoothing (list-ref setup-valve-list valve-value)) (set! smoothing-on (cadr (assoc 'smoothing (list-ref setup-valve-list valve-value)))))

(if smoothing-on
(begin

(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-vlayer stationary " (string-append "fluid-" (list-ref (list-ref setup-valve-list valve-value) 1)  "-vlayer") " , 0 " ))
;-------------------------

(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create " (list-ref (list-ref setup-valve-list valve-value) 1)"-seat stationary ,, 0 " ))

)
(begin


(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create intf-int-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-vlayer stationary " (string-append "fluid-" (list-ref (list-ref setup-valve-list valve-value) 1)  "-vlayer") " ," (number->string ob_cell_height)))
;-------------------------

(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create " (list-ref (list-ref setup-valve-list valve-value) 1)"-seat stationary ,, " (number->string ob_cell_height)))


  )

)
;
; Deforming dynamic zone setup
;
(define min_scale)
(define max_scale)
(define ave_scale)
;(set! ave_scale (+ (* 0.5 (car (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch")) #f))) (* 0.5 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch")) #f))))))
 
(set! ave_scale (caddr  (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch")) #f)))

(set! min_scale (* 0.4 ave_scale))
(set! max_scale (* 1.4 ave_scale))

(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch deforming cylinder " (number->string (list-ref (list-ref setup-valve-list valve-value) 5)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 7) 0)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 7) 1)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 7) 2)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 0)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 1)) " " (number->string (list-ref (list-ref (list-ref setup-valve-list valve-value) 6) 2)) " no no yes yes no " (number->string min_scale) " " (number->string max_scale) " 0.6"))


;
;

; Symmetry plane setup
;

(if symmetryengine
(begin

;(set! ave_scale (+ (* 0.5 (car (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "symm-cyl-tri")) #f))) (* 0.5 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "symm-cyl-tri")) #f))))))
;(set! min_scale (* 0.5 ave_scale))
;(set! max_scale (* 1.5 ave_scale))

(set! ave_scale (caddr  (%get-min-max-length-scale-on-zone (thread-name->id  (string-append "symm-cyl-tri")) #f)))

(set! min_scale (* 0.4 ave_scale))
(set! max_scale (* 1.4 ave_scale))

;
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create symm-cyl-tri deforming plane " (number->string (list-ref point_on_symmetry 0)) " " (number->string (list-ref point_on_symmetry 1)) " " (number->string (list-ref point_on_symmetry 2)) " " (number->string (list-ref symm_normal 0)) " "  (number->string (list-ref symm_normal 1)) " " (number->string (list-ref symm_normal 2)) " no no yes yes no " (number->string min_scale) " " (number->string max_scale) " 0.65"))
(ti-menu-load-string "/def/bc/zt symm-cyl-tri symmetry")

)
)
;
; Dynamic mesh events setup

;Open and Close Valve

 (format  "\n\n\t\t ___________ Setting Dynamic Mesh Events _____________\t\n\n") 
 
(set! event01 (list '(active . #t) (cons 'name (string-append "open-" (list-ref (list-ref setup-valve-list valve-value) 1))) (cons 'start  (list-ref (list-ref setup-valve-list valve-value) 3) ) (list 'action 2 (list-ref (list-ref setup-valve-list valve-value) 1) (list (thread-name->id (string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-ch"))) (list (thread-name->id (string-append "intf-" (list-ref (list-ref setup-valve-list valve-value) 1) "-ob-fluid-vlayer")) ) 0 0))) ;vlayer


(set! event02 (list '(active . #t) (cons 'name (string-append "close-" (list-ref (list-ref setup-valve-list valve-value) 1))) (cons 'start (list-ref (list-ref setup-valve-list valve-value) 4)) (list 'action 3 (list-ref (list-ref setup-valve-list valve-value) 1))))


;Events Setup

(set! cur-event-list (list event01 event02)) 

;Append current events to previsoutly setup events if any

(set! pre-event-list (rpgetvar 'dynamesh/in-cyn/event-list))

(rpsetvar 'dynamesh/in-cyn/event-list (append pre-event-list cur-event-list))

;Events for non-constant URFs and dt

(if (list-ref (list-ref setup-valve-list valve-value) 8)

(begin

;Change URFs

(set! event01 (list '(active . #t) (cons 'name (string-append "reduce-urf-due-to-open-" (list-ref (list-ref setup-valve-list valve-value) 1))) (cons 'start  (- (list-ref (list-ref setup-valve-list valve-value) 3) 1) ) (list 'action 12 (list '(1 epsilon 0.2) '(1 k 0.2) '(1 mom 0.4) '(1 pressure 0.2)))))
(set! event02 (list '(active . #t) (cons 'name (string-append "increase-urf-due-to-open-" (list-ref (list-ref setup-valve-list valve-value) 1))) (cons 'start  (+ (list-ref (list-ref setup-valve-list valve-value) 3) increment)) (list 'action 12 (list '(1 epsilon 0.4) '(1 k 0.4) '(1 mom 0.7) '(1 pressure 0.5)))))
(set! event03 (list '(active . #t) (cons 'name (string-append "reduce-urf-due-to-close-" (list-ref (list-ref setup-valve-list valve-value) 1))) (cons 'start  (- (list-ref (list-ref setup-valve-list valve-value) 4) 1) ) (list 'action 12 (list '(1 epsilon 0.2) '(1 k 0.2) '(1 mom 0.4) '(1 pressure 0.2)))))
(set! event04 (list '(active . #t) (cons 'name (string-append "increase-urf-due-to-close-" (list-ref (list-ref setup-valve-list valve-value) 1))) (cons 'start  (+ (list-ref (list-ref setup-valve-list valve-value) 4) increment)) (list 'action 12 (list '(1 epsilon 0.4) '(1 k 0.4) '(1 mom 0.7) '(1 pressure 0.5)))))   

(set! cur-event-list (append cur-event-list  (list event01 event02 event03 event04)))

)

)

; Change dt
;
(if (list-ref (list-ref setup-valve-list valve-value) 9)
 
(begin
(set! event05 (list '(active . #t) (cons 'name (string-append "reduce-dt-due-to-open-" (list-ref (list-ref setup-valve-list valve-value) 1))) (cons 'start (list-ref (list-ref setup-valve-list valve-value) 3) ) (list 'action 5 (* 0.5 dCA))))
(set! event06 (list '(active . #t) (cons 'name (string-append "increase-dt-due-to-open-" (list-ref (list-ref setup-valve-list valve-value) 1))) (cons 'start (+ (list-ref (list-ref setup-valve-list valve-value) 3) increment)) (list 'action 5 dCA)))    
(set! event07 (list '(active . #t) (cons 'name (string-append "reduce-dt-due-to-close-" (list-ref (list-ref setup-valve-list valve-value) 1))) (cons 'start (list-ref (list-ref setup-valve-list valve-value) 4) ) (list 'action 5 (* 0.5 dCA))))
(set! event08 (list '(active . #t) (cons 'name (string-append "increase-dt-due-to-close-" (list-ref (list-ref setup-valve-list valve-value) 1))) (cons 'start (+ (list-ref (list-ref setup-valve-list valve-value) 4) increment)) (list 'action 5 dCA)))
;

(set! cur-event-list (append cur-event-list  (list event05 event06 event07 event08)))

)

)
;Find the existing event list

;(set! pre-event-list (rpgetvar 'dynamesh/in-cyn/event-list))

;Set the events

(rpsetvar 'dynamesh/in-cyn/event-list (append pre-event-list cur-event-list))

;;set the smoothing events
 
(format "\n valve value =~d setup valve list =~s" valve-value setup-valve-list) 
(set-smoothing-events (list-ref setup-valve-list valve-value))

)

;activate spring smoothing on all shapes

(ti-menu-load-string (format #f "/define/dynamic-mesh/controls/smoothing-parameters/spring-on-all-shapes? y"))




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;Piston Setup
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  PISTON-1
 (format  "\n\n\t\t ___________  Creating Grid Interfaces for Piston _____________\t\n\n")  
 
(if (= piston_type 1)
(begin

(define min_scale)
(define max_scale)
(define ave_scale)
(define piston_layer_height)

; Change length unit to meter
;
(ti-menu-load-string "/def/units/length m")
;
; Change bc
;
(ti-menu-load-string (string-append  "/def/bc/zt int-piston interior"))
;
; Rigid body dynamic zone setup
;
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create fluid-piston-layer rigid-body  " piston-full-prof "  " (number->string (list-ref cylinder_axis 0)) " " (number->string (list-ref cylinder_axis 1)) " " (number->string (list-ref cylinder_axis 2))))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create piston rigid-body  " piston-full-prof "  " (number->string (list-ref cylinder_axis 0)) " " (number->string (list-ref cylinder_axis 1)) " " (number->string (list-ref cylinder_axis 2)) " ,, 0"))
(set! ave_scale (+ (* 0.5 (car (%get-min-max-length-scale-on-zone (thread-name->id  "int-piston") #f))) (* 0.5 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  "int-piston") #f))))))
(set! piston_layer_height (* 1.2 ave_scale))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create int-piston rigid-body " piston-limit-prof "  "(number->string (list-ref cylinder_axis 0)) " " (number->string (list-ref cylinder_axis 1)) " " (number->string (list-ref cylinder_axis 2)) " fluid-ch " "," (number->string ave_scale) " fluid-piston-layer " "," (number->string piston_layer_height)))
;
; Deforming dyanmic zone setup
;
(set! ave_scale (+ (* 0.5 (car (%get-min-max-length-scale-on-zone (thread-name->id  "cyl-tri") #f))) (* 0.5 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  "cyl-tri") #f))))))
(set! min_scale (* 0.5 ave_scale))
(set! max_scale (* 1.5 ave_scale))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create cyl-tri deforming cylinder " (number->string cylinder_radius) " " (number->string (list-ref cylinder_axis_origin 0)) " " (number->string (list-ref cylinder_axis_origin 1)) " " (number->string (list-ref cylinder_axis_origin 2)) " " (number->string (list-ref cylinder_axis 0)) " " (number->string (list-ref cylinder_axis 1)) " " (number->string (list-ref cylinder_axis 2)) " no no yes yes no " (number->string min_scale) " " (number->string max_scale) " 0.6"))
;
) 

)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  PISTON-2
(if  (= piston_type 2)

(begin

(define min_scale)
(define max_scale)
(define ave_scale)
;
; Change length unit to meter
;
(ti-menu-load-string "/def/units/length m")
;
; Rigid body dynamic zone setup
;
(set! ave_scale (+ (* 0.5 (car (%get-min-max-length-scale-on-zone (thread-name->id  "piston") #t))) (* 0.5 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  "piston") #t))))))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create piston rigid-body  " piston-full-prof "  " (number->string (list-ref cylinder_axis 0)) " " (number->string (list-ref cylinder_axis 1)) " " (number->string (list-ref cylinder_axis 2)) " ,, " (number->string ave_scale)))
;
; Deforming dynamic zone setup
;
(set! ave_scale (+ (* 0.5 (car (%get-min-max-length-scale-on-zone (thread-name->id  "cyl-tri") #t))) (* 0.5 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  "cyl-tri") #t))))))
(set! min_scale (* 0.5 ave_scale))
(set! max_scale (* 1.5 ave_scale))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create cyl-tri deforming cylinder " (number->string cylinder_radius) " " (number->string (list-ref cylinder_axis_origin 0)) " " (number->string (list-ref cylinder_axis_origin 1)) " " (number->string (list-ref cylinder_axis_origin 2)) " " (number->string (list-ref cylinder_axis 0)) " " (number->string (list-ref cylinder_axis 1)) " " (number->string (list-ref cylinder_axis 2)) " no no yes yes no " (number->string min_scale) " " (number->string max_scale) " 0.6"))
;
; If with a bowl
;
(if withbowl
(begin
;
; Change bc
;
(ti-menu-load-string "/def/bc/zt intf-bowl interface")

(ti-menu-load-string "/def/bc/zt piston interface")
;
; Create non-conformal interface
;
(rpsetvar 'nonconformal/virtual-polygons? #t)
(ti-menu-load-string "/def/gi/create intf-bowl piston () intf-bowl () no no")
;
; Set up dynamic zones
;
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create fluid-piston rigid-body  " piston-full-prof "  " (number->string (list-ref cylinder_axis 0)) " " (number->string (list-ref cylinder_axis 1)) " " (number->string (list-ref cylinder_axis 2)))) ;fluid-piston from fluid-bowl
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create bowl rigid-body " piston-full-prof "  "  (number->string (list-ref cylinder_axis 0)) " " (number->string (list-ref cylinder_axis 1)) " " (number->string (list-ref cylinder_axis 2)) " ,, 0"))

)
)
; 
; Events setup
;
(define pre-event-list)
(define cur-event-list)
(define event01)
(define event02)
(define event03)
(define event04)
;
(set! event01 (list '(active . #t) (cons 'name "insert-piston-boundary-layer-power-stroke") (cons 'start  (+ 0   insert_angle_power_stroke)) (list 'action 6 (thread-name->id "piston") (thread-name->id "cyl-tri"))))
(set! event02 (list '(active . #t) (cons 'name "delete-piston-boundary-layer-power-stroke") (cons 'start  (- 360 insert_angle_power_stroke)) (list 'action 7 (thread-name->id "piston"))))
(set! event03 (list '(active . #t) (cons 'name "insert-piston-boundary-layer-intake-stroke") (cons 'start  (+ 360 insert_angle_intake_stroke)) (list 'action 6 (thread-name->id "piston") (thread-name->id "cyl-tri"))))
(set! event04 (list '(active . #t) (cons 'name "delete-piston-boundary-layer-intake-stroke") (cons 'start  (- 720 insert_angle_intake_stroke)) (list 'action 7 (thread-name->id "piston"))))
;
(set! cur-event-list (list event01 event02 event03 event04))
;
(set! pre-event-list (rpgetvar 'dynamesh/in-cyn/event-list))
(rpsetvar 'dynamesh/in-cyn/event-list (append pre-event-list cur-event-list))
;
;
)
)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  PISTON-3

(if  (= piston_type 3)
(begin

(define min_scale)
(define max_scale)
(define ave_scale)
;
; Change length unit to meter
;
(ti-menu-load-string "/def/units/length m")
;
; Rigid body dynamic zone setup
;
(set! ave_scale (+ (* 0.5 (car (%get-min-max-length-scale-on-zone (thread-name->id  "piston") #t))) (* 0.5 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  "piston") #t))))))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create piston rigid-body  " piston-full-prof "  "  (number->string (list-ref cylinder_axis 0)) " " (number->string (list-ref cylinder_axis 1)) " " (number->string (list-ref cylinder_axis
 2)) " ,, " (number->string ave_scale)))
;


; Deforming dynamic zone setup
;
(set! ave_scale (+ (* 0.5 (car (%get-min-max-length-scale-on-zone (thread-name->id  "cyl-tri") #t))) (* 0.5 (car (cdr (%get-min-max-length-scale-on-zone (thread-name->id  "cyl-tri") #t))))))
(set! min_scale (* 0.5 ave_scale))
(set! max_scale (* 1.5 ave_scale))
(ti-menu-load-string (string-append  "/define/dynamic-mesh/zones/create cyl-tri deforming cylinder " (number->string cylinder_radius) " " (number->string (list-ref cylinder_axis_origin 0)) " " (number->string (list-ref cylinder_axis_origin 1))
 " " (number->string (list-ref cylinder_axis_origin 2)) " " (number->string (list-ref cylinder_axis 0)) " " (number->string (list-ref cylinder_axis 1)) " " (number->string (list-ref cylinder_axis 2)) " no no yes yes no " (number->string min_scale) " " (number->string max_scale) " 0.6"))
;
; Events setup
;
(define pre-event-list)
(define cur-event-list)
(define event01)
(define event02)
(define event03)
(define event04)
(define meshreplace1)
(define meshreplace2)
(define meshreplace3)
(define meshreplace4)
;
(set! meshreplace1 (list '(active . #t) (cons 'name "call-mesh-replace-initial") (cons 'start (- insert_angle_power_stroke 1)) (list 'action 15 (format #f "(mesh-replace-initial)")))) 
(set! event01 (list '(active . #t) (cons 'name "insert-piston-interior-layer1") (cons 'start  insert_angle_power_stroke) (list 'action 8 (thread-name->id "piston") (thread-name->id "cyl-tri") "zone1" "zone2")))
(set! meshreplace2 (list '(active . #t) (cons 'name "call-mesh-replace-final") (cons 'start insert_angle_power_stroke) (list 'action 15 (format #f "(mesh-replace-final)")))) 
(set! event02 (list '(active . #t) (cons 'name "delete-piston-interior-layer1") (cons 'start  (- 360 insert_angle_power_stroke)) (list 'action 9 "zone1" "zone2")))

(set! meshreplace3 (list '(active . #t) (cons 'name "call-mesh-replace-initial") (cons 'start (+ 360 (- insert_angle_power_stroke 1))) (list 'action 15 (format #f "(mesh-replace-initial)"))))
(set! event03 (list '(active . #t) (cons 'name "insert-piston-interior-layer2") (cons 'start  (+ 360 insert_angle_power_stroke)) (list 'action 8 (thread-name->id "piston") (thread-name->id "cyl-tri") "zone1" "zone2")))
(set! meshreplace4 (list '(active . #t) (cons 'name "call-mesh-replace-final") (cons 'start (+ 360 insert_angle_power_stroke)) (list 'action 15 (format #f "(mesh-replace-final)"))))
 
(set! event04 (list '(active . #t) (cons 'name "delete-piston-interior-layer2") (cons 'start  (- 720 insert_angle_power_stroke)) (list 'action 9 "zone1" "zone2")))
;
(set! cur-event-list (list meshreplace1 event01 meshreplace2 event02 meshreplace3 event03 meshreplace4 event04))
;
(set! pre-event-list (rpgetvar 'dynamesh/in-cyn/event-list))
(rpsetvar 'dynamesh/in-cyn/event-list (append pre-event-list cur-event-list))

)
)


)

(begin 
(pp (format #f "~%"))
(pp "The following zones did not match")
(pp final_list)
)


)

)

)


;;Retrieve the zones from the valveId and set the smoothing events from the information in valveData
(set!  set-smoothing-events
	  (lambda (valveData)
	    (let 
	      ((valve-id (list-ref valveData 1))
	      (smoothing-on #f)
	      (valve-open (list-ref valveData 3))
	      (valve-close (list-ref valveData 4))
	      (smoothing-end 0)
	      (smoothing-start 0)
	      (zone-name "")
	      (event-name "")
	      (adj-cell-zone "")
	      (event-list '())
	      (cell-height 0)
	      (pre-event-list (rpgetvar 'dynamesh/in-cyn/event-list))
	      )

	      (pp "\n in smoothing events")
	      (pp valve-id)
	      

	      (if (assoc 'smoothing valveData) (set! smoothing-on (cadr (assoc 'smoothing valveData))))

	      (format "valve data =  ~s smoothing=~s " valveData smoothing-on)

	      (if smoothing-on
		(begin
		   (if (assoc 'end-valve-open-smoothing valveData) (set! smoothing-end (cadr (assoc 'end-valve-open-smoothing valveData)))
		     (set! smoothing-on #f))
		    (if (assoc 'begin-valve-close-smoothing valveData) (set! smoothing-start (cadr (assoc 'begin-valve-close-smoothing valveData)))
		     (set! smoothing-on #f))
		    )
		)

	      (format "valve data =  ~s smoothing=~s  smoothing-end ~d smoothin-start=~d" valveData smoothing-on smoothing-end smoothing-start)
	      ( if smoothing-on
		   (begin 

		     
		     (set! zone-name (string-append "intf-" valve-id "-ob-fluid-vlayer"))
		     (set! event-name (string-append "make-" zone-name "-opening-deforming"))
		     (set! event-list (list-add event-list (define-deforming-zone-event (thread-name->id  zone-name) valve-open event-name ))) 
		      
		     (set! zone-name (string-append "intf-" valve-id "-ib-fluid-ob-quad"))
		     (set! event-name (string-append "make-" zone-name "-opening-deforming"))
		     (set! event-list (list-add event-list (define-deforming-zone-event (thread-name->id  zone-name) valve-open event-name )) )

		      (set! zone-name (string-append "intf-int-" valve-id "-ob-fluid-vlayer"))
		      (set! adj-cell-zone (string-append "fluid-" valve-id "-vlayer"))
		     (set! event-name (string-append "reduce-stationary-" zone-name "-cell-height"))
		     (set! event-list(list-add event-list (define-stationary-zone-event (thread-name->id  zone-name) valve-open event-name (thread-name->id  adj-cell-zone) 0) ) )
		    
		     (set! zone-name (string-append valve-id "-seat"))
		      (set! adj-cell-zone (string-append "fluid-" valve-id "-vlayer"))
		     (set! event-name (string-append "reduce-stationary-" zone-name "-cell-height"))
		     (set! event-list(list-add event-list (define-stationary-zone-event (thread-name->id  zone-name) valve-open event-name (thread-name->id  adj-cell-zone) 0) ) )

		       (set! zone-name (string-append "fluid-" valve-id "-vlayer"))
		       (set! event-name (string-append "delete-" zone-name ))
		     (set! event-list(list-add event-list (delete-dynamic-zone-event (thread-name->id  zone-name) valve-open event-name )) )



		     (set! zone-name (string-append "intf-" valve-id "-ob-fluid-vlayer"))
		     (set! event-name (string-append "delete-" zone-name "-deforming"))
		     (set! event-list(list-add event-list (delete-dynamic-zone-event (thread-name->id  zone-name) smoothing-end event-name )))

		     (set! zone-name (string-append "intf-" valve-id "-ib-fluid-ob-quad"))
		     (set! event-name (string-append "delete-" zone-name "deforming"))
		     (set! event-list(list-add event-list (delete-dynamic-zone-event (thread-name->id  zone-name) smoothing-end event-name )))

		     (set! zone-name (string-append "fluid-" valve-id "-vlayer"))
		     (set! event-name (string-append "make-" zone-name "-rigid-body" ))
		     (set! event-list(list-add event-list (define-rigidbody-zone-event (thread-name->id  zone-name) smoothing-end event-name (list-ref valveData 2) (list-ref valveData 6))) )

		     (set! zone-name (string-append "intf-int-" valve-id "-ob-fluid-vlayer"))
		     (set! adj-cell-zone (string-append "fluid-" valve-id "-vlayer"))
		     (set! event-name (string-append "increase-stationary-" zone-name "-cell-height"))
		     ;(set! cell-height (* 0.5 (caddr (%get-min-max-length-scale-on-zone (thread-name->id  zone-name) #f))))
			 ; setting the vlayer height equal to 2mm /4; considering we are doing smoothing till 2mm
			 (set! cell-height 0.0005)
		     (set! event-list(list-add event-list (define-stationary-zone-event (thread-name->id  zone-name) smoothing-end event-name (thread-name->id  adj-cell-zone) cell-height ) ) )

		     (set! zone-name (string-append valve-id "-seat"))
		     (set! adj-cell-zone (string-append "fluid-" valve-id "-vlayer"))
		     (set! event-name (string-append "increase-stationary-" zone-name "-cell-height"))
		     (set! event-list(list-add event-list (define-stationary-zone-event (thread-name->id  zone-name) smoothing-end event-name (thread-name->id  adj-cell-zone) cell-height) )) 



		     (set! zone-name (string-append "intf-" valve-id "-ob-fluid-vlayer"))
		       (set! event-name (string-append "make-" zone-name "-deforming-at-closing"))
		     (set! event-list(list-add event-list (define-deforming-zone-event (thread-name->id  zone-name) smoothing-start event-name ))) 

		     (set! zone-name (string-append "intf-" valve-id "-ib-fluid-ob-quad"))
		       (set! event-name (string-append "make-" zone-name "-deforming-at-closing"))
		     (set! event-list(list-add event-list (define-deforming-zone-event (thread-name->id  zone-name) smoothing-start event-name ))) 
		     
		      (set! zone-name (string-append "intf-int-" valve-id "-ob-fluid-vlayer"))
		      (set! adj-cell-zone (string-append "fluid-" valve-id "-vlayer"))
		     (set! event-name (string-append "reduce-stationary-" zone-name "-cell-height-at-closing"))
		      (set! event-list (list-add event-list (define-stationary-zone-event (thread-name->id  zone-name) smoothing-start event-name (thread-name->id  adj-cell-zone) 0) ) )
		    
		     (set! zone-name (string-append valve-id "-seat"))
		      (set! adj-cell-zone (string-append "fluid-" valve-id "-vlayer"))
		     (set! event-name (string-append "reduce-stationary-" zone-name "-cell-height-at-closing"))
		     (set! event-list (list-add event-list (define-stationary-zone-event (thread-name->id  zone-name) smoothing-start event-name (thread-name->id  adj-cell-zone) 0) )) 

		       (set! zone-name (string-append "fluid-" valve-id "-vlayer"))
		       (set! event-name (string-append "delete-" zone-name "-at-closing" ))
		     (set! event-list (list-add event-list (delete-dynamic-zone-event (thread-name->id  zone-name) smoothing-start event-name )))

	     
		     (rpsetvar 'dynamesh/in-cyn/event-list (append pre-event-list event-list))

		     ))

		   )
		)
	      

	      )
	   


(set! define-deforming-zone-event
  (lambda (zone-id start-angle event-name)
    (let ((deform-event '()))

      ;(pp "in define deforming zone")
      (set! deform-event (list '(active . #t) (cons 'name event-name) (cons 'start start-angle) (list 'action 15 (format #f "define dynamic-mesh zones create ~s deforming , , , , no ~s ~s 0.6" zone-id  (* 0.4 (caddr (%get-min-max-length-scale-on-zone zone-id #f)))  (* 1.4 (caddr (%get-min-max-length-scale-on-zone zone-id #f))) ))))
      deform-event
      )

    )
  )

  (set! delete-dynamic-zone-event
  (lambda (zone-id start-angle event-name)
    (let ((deform-event '()))
      ;(pp "in define deforming zone 2")

      (set! deform-event (list '(active . #t) (cons 'name event-name) (cons 'start start-angle) (list 'action 15 (format #f "define dynamic-mesh zones delete ~s" zone-id ))))
      deform-event
      )

    )
  )


(set! define-stationary-zone-event
  (lambda (zone-id start-angle event-name adjacent-zone cell-height)
    (let ((deform-event '()))

      ;(pp "in define deforming zone 3 3")
      (set! deform-event (list '(active . #t) (cons 'name event-name) (cons 'start start-angle) (list 'action 15 (format #f "define dynamic-mesh zones create ~s stationary ~s , ~s" zone-id  adjacent-zone cell-height ))))
     deform-event
      )

    )
  )

(set! define-rigidbody-zone-event
  (lambda (zone-id start-angle event-name motion-profile motion-axis)
    (let ((deform-event '()))

      ;(pp "in define deforming zone 4")
      (set! deform-event (list '(active . #t) (cons 'name event-name) (cons 'start start-angle) (list 'action 15 (format #f "define dynamic-mesh zones create ~s rigid-body  ~s  ~s ~s  ~s " zone-id  motion-profile (car motion-axis)  (cadr motion-axis)  (caddr motion-axis) ))))
      deform-event
      )

    )
  )

(set! read-valve-profile
  (lambda (valve-data)
    (let( (profile-name #f))
	      (if (assoc 'profile-path valve-data) (set! profile-name (cadr (assoc 'profile-path valve-data))))
	      (if  (and profile-name (file-exists? profile-name))

		(ti-menu-load-string (string-append  "/file/read-profile  " profile-name))
		)
     
      )))


(set! wb-ice-urf-setting

  (lambda ()
  (format  "\n\n\t\t ___________Setting Up URF_____________\t\n\n")  
    (ti-menu-load-string "/solve set p-v-coupling 22")
    (ti-menu-load-string "/solve set p-v-controls 1 1 no")
    (ti-menu-load-string "/solve set gradient-scheme yes")
    (ti-menu-load-string "/solve set discretization-scheme pressure 14")
    (ti-menu-load-string "/solve set discretization-scheme density 1")
    (ti-menu-load-string "/solve set discretization-scheme mom 1")
    (ti-menu-load-string "/solve set discretization-scheme k 1")
	(ti-menu-load-string "/solve set discretization-scheme epsilon 0")
	(ti-menu-load-string "/solve set discretization-scheme temp 1")
	
	(ti-menu-load-string "/solve set under-rel pressure 0.3")
	(ti-menu-load-string "/solve set under-rel mom 0.5")
	(ti-menu-load-string "/solve set under-rel k 0.4")
	(ti-menu-load-string "/solve set under-rel epsilon 0.4")
	(ti-menu-load-string "/solve set under-rel density 1")
	(ti-menu-load-string "/solve set under-rel b-f 1")
	(ti-menu-load-string "/solve set under-rel temp 1")
	
	(ti-menu-load-string "/sol/set/limits 1000 5e7 250 1000 ,,,")
	(ti-menu-load-string "(rpsetvar 'max-iterations-per-step 30)")
	;(ti-menu-load-string "/define d-m events e-e-f ice-event.event")

    )
)
(set! wb-ice-create-planes

(lambda ()
(format  "\n\n\t\t ___________Defining User planes_____________\t\n\n")  
    (define post-plane)
	(set! post-plane  (read-var solver-set-file 'post_planes))
	(do ((i 0 (+ i 1))) ((= i (length post-plane)))

  (let 
  	(
	 (point1_x (list-ref (list-ref (list-ref post-plane i) 1) 0))
	 (point1_y (list-ref (list-ref (list-ref post-plane i) 1) 1))
	 (point1_z (list-ref (list-ref (list-ref post-plane i) 1) 2))
	 (point2_x (list-ref (list-ref (list-ref post-plane i) 2) 0))
	 (point2_y (list-ref (list-ref (list-ref post-plane i) 2) 1))
	 (point2_z (list-ref (list-ref (list-ref post-plane i) 2) 2))
	 (point3_x (list-ref (list-ref (list-ref post-plane i) 3) 0))
	 (point3_y (list-ref (list-ref (list-ref post-plane i) 3) 1))
	 (point3_z (list-ref (list-ref (list-ref post-plane i) 3) 2))
	(name (list-ref (list-ref post-plane i) 0))
	
	) 
	(ti-menu-load-string (format #f "/sur plane ~s ~s ~s ~s ~s ~s ~s ~s ~s ~s" (string->symbol name) point1_x point1_y point1_z point2_x point2_y point2_z point3_x point3_y point3_z))
	
    )
    )
  )
)

(set! wb-ice-create-default-images-on-planes

(lambda ()
(format  "\n\n\t\t ___________Defining Images on User planes_____________\t\n\n")  
    (define post-plane)
	(set! post-plane  (read-var solver-set-file 'post_planes))
	(do ((i 0 (+ i 1))) ((= i (length post-plane)))

  (let* 
  	(
	 
	(name (list-ref (list-ref post-plane i) 0))
	(s1 "time-step")
	(s2 (format #f "/dis set cont surf ( ~s )" (string->symbol name)))
	(s3 "/dis cont v-m , ,")
	(s4 (format #f "/dis view rv my-~s-view" (string->symbol name)))
	(s5 (format #f "dis hc ../Report/ice-anim-vel-mag-on-~s-%t.jpg" (string->symbol name)))
	(s6 (format #f "/dis set ms (~s)" (string->symbol name)))
	(s7 "/dis s-m , ,")
	(s8 (format #f "dis hc ../Report/ice-anim-mesh-on-~s-%t.jpg" (string->symbol name)))
	) 
	(ti-menu-load-string "/dis set pic driver jpeg")
	(ti-menu-load-string (format #f "/dis set cont surf ~s ,"(string->symbol name)))
	(ti-menu-load-string "/dis set cont filled-cont yes")
	;;(ti-menu-load-string "view auto")
	
	(ti-menu-load-string (format #f  "/solv exe-comm a-e , 15 ~s ~s" s1 s2))
	(ti-menu-load-string (format #f  "/solv exe-comm a-e , 15 ~s ~s" s1 s3))
	(ti-menu-load-string (format #f  "/solv exe-comm a-e , 15 ~s ~s" s1 s4))
	(ti-menu-load-string (format #f  "/solv exe-comm a-e , 15 ~s ~s" s1 s5))
	(ti-menu-load-string (format #f  "/solv exe-comm a-e , 15 ~s ~s" s1 s6))
	(ti-menu-load-string (format #f  "/solv exe-comm a-e , 15 ~s ~s" s1 s7))
	(ti-menu-load-string (format #f  "/solv exe-comm a-e , 15 ~s ~s" s1 s8))


    )
    )
  )
)

(set! wb-ice-create-view-for-planes

(lambda ()
(format  "\n\n\t\t ___________ Defining Views for User Planes _____________\t\n\n") 
    (define post-plane)
	

	(set! post-plane  (read-var solver-set-file 'post_planes))
	(if post-plane
	  (begin 
	
	 (define cyl-axis  (read-var ICDyMeshSetFile 'cylinder_axis))
	(define up-vectorx (car cyl-axis))
	(define up-vectory (cadr cyl-axis))
	(define up-vectorz (caddr cyl-axis))

	 (define cyl-origin  (read-var ICDyMeshSetFile 'cylinder_axis_origin))
	(define ox (car cyl-origin))
	(define oy (cadr cyl-origin))
	(define oz (caddr cyl-origin))
	
	 (wb-ice-report-dir)
	 (ti-menu-load-string (format #f  "/dis view w-v ice-view1 front ()"))
	(define ip (open-input-file "ice-view1.vw"))
	(define sexp (read ip))
	(define l1 (cadr (caadr (caaadr sexp))))
	(define position (car l1))
	(define target (cadr l1))
	(define px (car position))
	(define py (cadr position))
	(define pz (caddr position))
	(define tx (car target))
	(define ty (cadr target))
	(define tz (cadr target))
	(do ((i 0 (+ i 1))) ((= i (length post-plane)))
	
  (let 
  	(
	 (x0 (list-ref (list-ref (list-ref post-plane i) 1) 0))

	 (y0 (list-ref (list-ref (list-ref post-plane i) 1) 1))
	 (z0 (list-ref (list-ref (list-ref post-plane i) 1) 2))
	 (x1 (list-ref (list-ref (list-ref post-plane i) 2) 0))
	 (y1 (list-ref (list-ref (list-ref post-plane i) 2) 1))
	 (z1 (list-ref (list-ref (list-ref post-plane i) 2) 2))
	 (x2 (list-ref (list-ref (list-ref post-plane i) 3) 0))
	 (y2 (list-ref (list-ref (list-ref post-plane i) 3) 1))
	 (z2 (list-ref (list-ref (list-ref post-plane i) 3) 2))
	 (name (list-ref (list-ref post-plane i) 0))
	)
	(define x-coeff)
	(define y-coeff)
	(define z-coeff)
	(set! x-coeff (- (* (- y2 y0) (- z1 z0)) (* (- y1 y0) (- z2 z0))))
	(set! y-coeff (- (* (- x1 x0) (- z2 z0)) (* (- x2 x0) (- z1 z0))))
	(set! z-coeff (- (* (- x2 x0) (- y1 y0)) (* (- x1 x0) (- y2 y0))))
	(define unit-x)
	(define unit-y)
	(define unit-z)
	(set! unit-x (/ x-coeff (sqrt (+ (* x-coeff x-coeff) (* y-coeff y-coeff) (* z-coeff z-coeff)))))
	(set! unit-y (/ y-coeff (sqrt (+ (* x-coeff x-coeff) (* y-coeff y-coeff) (* z-coeff z-coeff)))))
	(set! unit-z (/ z-coeff (sqrt (+ (* x-coeff x-coeff) (* y-coeff y-coeff) (* z-coeff z-coeff)))))

	(define pos-dis (sqrt (+ (* (- px tx) (- px tx)) (* (- py ty) (- py ty)) (* (- pz tz) (- pz tz)))) )
	(define px1 (+ x0 (* pos-dis unit-x)))
	(define py1 (+ y0 (* pos-dis unit-y)))
	(define pz1 (+ z0 (* pos-dis unit-z)))

	(ti-menu-load-string (format #f "/dis view camera po ~s ~s ~s" px1 py1 pz1))
	(ti-menu-load-string (format #f "/dis view camera target ~s ~s ~s" x0 y0 z0))
	(ti-menu-load-string (format #f "/dis view camera u-v ~s ~s ~s" up-vectorx up-vectory up-vectorz))
	(ti-menu-load-string (format #f "/dis set ms (~s)" (string->symbol name)))
	(ti-menu-load-string "/dis s-m , ,")
	(ti-menu-load-string "view auto")
	(ti-menu-load-string (format #f "/dis view camera target ~s ~s ~s" ox oy oz))
	(ti-menu-load-string (format #f "/dis view s-v my-~s-view" (string->symbol name)))
	)
  )))
  )
)
(set! wb-ice-postproc-monitors-default

	(lambda ()
	(pp "inside wb-ice-postproc-monitors-default")
	(define engine-type)
	(set! engine-type  (read-var ICDyMeshSetFile 'layering))
	(let 
		((s1 "Volume-Average")
		(s2 "Mass-Average"))
	  (if engine-type 
		;(pp "engine type is straight valve engine")
	   (begin
		(ti-menu-load-string (format #f "/sol/mon vol s-m pressure-mon1 ~s pressure fluid-ch-upper fluid-ch-lower , no no yes  pressure-mon1.out 1 no f-t" s1))
		(ti-menu-load-string (format #f "/sol/mon vol s-m temp-mon1 ~s temperature fluid-ch-upper fluid-ch-lower , no no yes  pressure-mon1.out 1 no f-t" s2))
		(ti-menu-load-string (string-append "/sol/mon vol s-m mesh-quality , cell-equivolume-skew fluid-ch-upper fluid-ch-lower , no no yes  mesh-mon.out 1 no f-t"))
	   )
		;(pp "Engine is canted valve engine")
	  (begin 
		(ti-menu-load-string (string-append "/sol/mon vol s-m mesh-quality , cell-equivolume-skew fluid-ch , no no yes  mesh-mon.out 1 no f-t"))
		(ti-menu-load-string (format #f "/sol/mon vol s-m temp-mon1 ~s temperature fluid-ch , no no yes  temp-mon.out 1 no f-t" s2))
		(ti-menu-load-string (format #f "/sol/mon vol s-m pressure-mon1 ~s pressure fluid-ch , no no yes  pressure-mon.out 1 no f-t" s1))
	  )
	
	 )
	)
	)
)
(set! wb-ice-postproc-monitors-user

	(lambda ()
	(let ((monitor-def  (read-var solver-set-file 'monitor_def))
	      )
	(if monitor-def
	  (begin
	    (format  "\n\n\t\t ___________Setting Up Monitors_____________\t\n\n")  
	(do ((i 0 (+ i 1))) ((= i (length monitor-def)))
	
	(let 
		(
		(mon-name (list-ref (list-ref monitor-def i) 0))
		(report-type (list-ref (list-ref monitor-def i) 2))
		(field-var (list-ref (list-ref monitor-def i) 3))
		(zones (list-ref (list-ref monitor-def i) 4))
		)
				
	(if (string=? (list-ref (list-ref monitor-def i) 1) "surface")
		(ti-menu-load-string (format #f "/sol/mon/surface s-m ~s ~s ~s ~s no no yes ~s 1 yes f-t" (string->symbol mon-name) report-type (string->symbol field-var) zones (string-append "ice-" mon-name ".out")))
		)
	(if (string=? (list-ref (list-ref monitor-def i) 1) "volume")
	(ti-menu-load-string (format #f "/sol/mon/volume s-m ~s ~s ~s ~s no no yes ~s 1 yes f-t" (string->symbol mon-name) report-type (string->symbol field-var) zones (string-append "ice-" mon-name ".out")))
	)
	(format "\n\n")
	)
	)
	(wb-ice-report-dir)
	(if (file-exists? "ice-monitors") 
	(ti-menu-load-string (format #f "/file start-trans ice-monitors yes"))
	(ti-menu-load-string (format #f "/file start-trans ice-monitors"))
	)
	(ti-menu-load-string (format #f "/sol mon surf li"))
	(ti-menu-load-string (format #f "/sol mon vol li"))
	(ti-menu-load-string (format #f "/file stop-trans"))
	)
	
	)
	)
	)
	)

(set! wb-ice-initialization

  (lambda ()
  (format  "\n\n\t\t ___________Initializing the FLow_____________\t\n\n")  

	(if (not (data-valid?))
		(ti-menu-load-string "/solve init init")
	)
	
	(begin
	(set! patching-values (read-var solver-set-file 'patching_values))
	(if patching-values
	(do ((i 0 (+ i 1))) ((= i (length patching-values)))
	(let
	(
	(patch-parts1 ())
	(patch-value (list-ref (list-ref patching-values i) 1))
	(patch-variable (list-ref (list-ref patching-values i) 0))
	(patch-parts (list-ref (list-ref patching-values i) 2))
	)

	(do ((j 0 (+ j 1))) ((= j (length patch-parts)))
	(let 
	((p (thread-name->id (list-ref patch-parts j))))
	(set! patch-parts1 (append patch-parts1 (list p)))
	)
	)

	(ti-menu-load-string (format #f "/solve patch ~s ~s  ~s " patch-parts1 (string->symbol patch-variable) patch-value))
	))
	))

	;(wb-ic-update-time-steps)

	
  )
)
	
(set! wb-ic-update-time-steps
  (lambda ()

    (rpsetvar 'number-of-time-steps 1)
    (ti-menu-load-string "/solve it 0")
    )
  )
(set! wb-ice-swirl-data

  (lambda ()
    (define swirl-data)
    (set! swirl-data (read-var solver-set-file 'swirl_data))
    (set! cylinder_axis_origin       (read-var ICDyMeshSetFile 'cylinder_axis_origin))
    (set! setup-valve-list           (read-var ICDyMeshSetFile 'setup_valve_list))
    (set! valve_origin (list-ref (list-ref setup-valve-list 0) 7) )
    (let* (
	   (unitx (- (car valve_origin) (car cylinder_axis_origin)))
	   (unity (- (cadr valve_origin) (cadr cylinder_axis_origin)))
	   (unitz (- (caddr valve_origin) (caddr cylinder_axis_origin)))
	   (unitm (sqrt (+ (* unitx unitx) (* unity unity) (* unitz unitz))))
	   ( a1 (/ unitx unitm))
	   ( a2 (/ unity unitm))
	   ( a3 (/ unitz unitm))
	   ( b1 (car cylinder_axis))
	   ( b2 (cadr cylinder_axis))
	   ( b3 (caddr cylinder_axis))
	   ( c1 (- (* a2 b3) (* a3 b2)))
	   ( c2 (- (* a3 b1) (* a1 b3)))
	   ( c3 (- (* a1 b2) (* a2 b1)))
	   (unitm1 (sqrt (+ (* c1 c1) (* c2 c2) (* c3 c3))))
	   )
      (define tumbx1 (/ c1 unitm1))
      (define tumbx2 (/ c2 unitm1))
      (define tumbx3 (/ c3 unitm1))
      (define tumby1 (- (* b2 tumbx3) (* b3 tumbx2)))
      (define tumby2 (- (* b3 tumbx1) (* b1 tumbx3)))
      (define tumby3 (- (* b1 tumbx2) (* b2 tumbx1)))
      (define unit1 (sqrt (+ (* tumby1 tumby1) (* tumby2 tumby2) (* tumby3 tumby3))))
      (set! tumby1 (/ tumby1 unit1))
      (set! tumby2 (/ tumby2 unit1))
      (set! tumby3 (/ tumby3 unit1))
      (ti-menu-load-string (format #f "/def d-m c ico yes 1 ~s ..\Report\ice-incylinder-output.txt cg ~s ~s ~s ~s ~s ~s ~s ~s ~s" swirl-data (car cylinder_axis) (cadr cylinder_axis) (caddr cylinder_axis) tumbx1 tumbx2 tumbx3 tumby1 tumby2 tumby3))
      ) 
    )
  )
(set! wb-ice-inlet-outlet-bc

  (lambda ()
    (set! inlet-outlet-bc-list           (read-var solver-set-file 'inlet_outlet_bc_list))
    (if inlet-outlet-bc-list
      (begin
	(format  "\n\n\t\t ___________Setting Up Inlet and Outlet BC_____________\t\n\n")  
	(do ((i 0 (+ i 1))) ((= i (length inlet-outlet-bc-list)))
	  (let
	    (
	     (value 0)
	     (temp 0)
	     (turb1 1)
	     (turb2 5)
	     (turbtype "int-hyd-dia")
	     (zonename #f)
	     (data-profile? "no")
	     (temp-profile? "no")
	     (data-profile-name "")
	     (data-profile-field "")
	     (temp-profile-name "")
	     (temp-profile-field "")
	     (temp-string #f)
	     (data-string #f)
	     (turb-string #f)
	     )

	    (if (assoc 'zone (list-ref inlet-outlet-bc-list i))
	      (set! zonename  (cadr (assoc 'zone (list-ref inlet-outlet-bc-list i))))
	      )

	    (if (assoc 'data-profile (list-ref inlet-outlet-bc-list i))
	      (begin
		(set! data-profile? (cadr(assoc 'data-profile (list-ref inlet-outlet-bc-list i))))
		(if (string-ci=? data-profile? "yes")
		  (begin

		    (set! data-profile-name (caddr (assoc 'data-profile (list-ref inlet-outlet-bc-list i))))
		    (set! data-profile-field (cadddr (assoc 'data-profile (list-ref inlet-outlet-bc-list i))))
		    (set! data-string (string-append " yes no \"" data-profile-name "\" \"" data-profile-field "\" " ))
		    )
		  (begin 
		    (set! value (caddr (assoc 'data-profile (list-ref inlet-outlet-bc-list i))))
		    (set! data-string (string-append " no " (number->string value) " " ))
		    )

		  )
		)
	      )


	    (if  (assoc 'turbulance (list-ref inlet-outlet-bc-list i))
	      (begin
		(set! turbtype (cadr (assoc 'turbulance (list-ref inlet-outlet-bc-list i))))
		(set! turb1 (caddr (assoc 'turbulance (list-ref inlet-outlet-bc-list i))))
		(set! turb2 (cadddr (assoc 'turbulance (list-ref inlet-outlet-bc-list i))))
		(if (string-ci=? turbtype "int-vis-ratio")
		  (set! turb-string  (string-append " no no yes " (number->string turb1) " " (number->string turb2) " " ))
		  (set! turb-string  (string-append " no no no yes " (number->string turb1) " " (number->string turb2) " " ))
		  )
		)
	      )

	    (if (assoc 'temperature-profile (list-ref inlet-outlet-bc-list i))
	      (begin 
		(set! temp-profile? (cadr (assoc 'temperature-profile (list-ref inlet-outlet-bc-list i))))
		(if (string-ci=? temp-profile? "yes")
		  (begin
		    (set! temp-profile-name (caddr (assoc 'temperature-profile (list-ref inlet-outlet-bc-list i))))
		    (set! temp-profile-field (cadddr (assoc 'temperature-profile (list-ref inlet-outlet-bc-list i))))
		    (set! temp-string (string-append " yes no \"" temp-profile-name "\" \"" temp-profile-field "\" " ))

		    )
		  (begin 
		    (set! temp (caddr (assoc 'temperature-profile (list-ref inlet-outlet-bc-list i))))
		    (set! temp-string (string-append " no " (number->string temp) " " ))
		    ))

		)
	      )


	    (format  "\n")
	    (if (and (zone-name->id zonename) data-string temp-string turb-string)
	      (begin
		(if (string=? (cadr (assoc 'type (list-ref inlet-outlet-bc-list i))) "pressure-inlet")
		  (begin 
		    (ti-menu-load-string (string-append "/def b-c m-z z-t " zonename " p-i"))
		    (format "\n")
		    (ti-menu-load-string (string-append "/def b-c p-i " zonename  " yes " data-string " no 0 " temp-string " no yes " turb-string ))
		    (format "\n")
		    )
		  )

		(if (string=? (cadr (assoc 'type (list-ref inlet-outlet-bc-list i))) "pressure-outlet")
		  (begin 
		   (ti-menu-load-string (string-append "/def b-c m-z z-t " zonename " p-o"))

		    (ti-menu-load-string (string-append "/def b-c p-o " zonename data-string  temp-string " no yes " turb-string " no no no " ))
		    )
		  )

		(if (string=? (cadr (assoc 'type (list-ref inlet-outlet-bc-list i))) "mass-flow-inlet")
		  (begin 
		    (ti-menu-load-string (string-append "/def b-c m-z z-t " zonename " mfi"))
		    (format  "\n")
		    (ti-menu-load-string (string-append "/def b-c mfi " zonename " yes yes " data-string temp-string " no 0 no yes " turb-string ))
		    (format  "\n")
		    )
		  )
		)
	      (format  "\n\n--Warning: Unable to setup BC for Zone ~s ---\n\n" zonename)
	      )

	    )
	  )
	)
      )

    ;;(ti-menu-load-string (format #f  "/file wbc ice-bc_rpm"))
    )
  )

(set! wb-ice-wall-bc

  (lambda ()
    ;(format  "\n\n\t\t ___________Setting Up Wall BC_____________\t\n\n")  
    (define wall-bc-list)
    (set! wall-bc-list           (read-var solver-set-file 'wall_bc_list))
    (if wall-bc-list
	(begin
	(format  "\n\n\t\t ___________Setting Up User Wall BC_____________\t\n\n") 
      (do ((i 0 (+ i 1))) ((= i (length wall-bc-list)))
	(let
	  (  
	     (zonename #f)
	     (value #f)
	     (data-profile? #f )
	     (data-profile-name #f)
	     (data-profile-field #f)
	     (data-string #f)
	     )
	   (if (assoc 'zone (list-ref wall-bc-list i))
	      (set! zonename  (cadr (assoc 'zone (list-ref wall-bc-list i))))
	      )

	    (if (assoc 'data-profile (list-ref wall-bc-list i))
	      (begin
		(set! data-profile? (cadr(assoc 'data-profile (list-ref wall-bc-list i))))
		(if (string-ci=? data-profile? "yes")
		  (begin

		    (set! data-profile-name (caddr (assoc 'data-profile (list-ref wall-bc-list i))))
		    (set! data-profile-field (cadddr (assoc 'data-profile (list-ref wall-bc-list i))))
		     (set! data-string (string-append " yes no \"" data-profile-name "\" \"" data-profile-field "\" " ))
		    )
		  (begin 
		    (set! value (caddr (assoc 'data-profile (list-ref wall-bc-list i))))
		    (set! data-string (string-append " no " (number->string value) " " ))
		    )

		  )
		)
	      )
	    (if (and (zone-name->id zonename) data-string)
	  (ti-menu-load-string (format #f "/def bc wall ~s 0 no 0 no yes temperature ~s no no no no , no , " (string->symbol zonename) (string->symbol data-string)))
	  (format  "\n\n--Warning: Unable to setup Wall BC for Zone ~s ---\n\n" zonename)
	  )
	  ))
      )
	  (begin
	  (format  "\n\n\t\t ___________Setting Up Default Wall BC_____________\t\n\n") 
	  (define initial-walls (map thread-name (sort-threads-by-name (get-threads-of-type 'wall))))
	  (do ((j 0 (+ j 1))) ((= j (length initial-walls)))
	  (ti-menu-load-string (format #f "/def bc wall ~s 0 no 0 no yes temperature no 300 no no no no , no , " (list-ref initial-walls j)))
	  )
	  )
      )
	  
    )
  )

(set! wb-ice-viscous-model-setup
  (lambda ()
    (format  "\n\n\t\t ___________Setting Up Models_____________\t\n\n")  
    (ti-menu-load-string "/def mode vis ke-standa yes")
    (ti-menu-load-string (format #f  "/def dm eve eef event"))
    (ti-menu-load-string "/dis set pic driver jpeg")
    (gui-events)
    (cx-gui-do cx-activate-item "Dynamic Mesh Events*PanelButtons*PushButton1(Cancel)")	
		(ti-menu-load-string (format #f  "/def dm eve eef ice-event_file"))
	)
)
(set! wb-ice-report-mesh

	(lambda ()
	  (format  "\n\n\t\t ___________Setting Up Mesh Report_____________\t\n\n")  
	(define cur-event-list)
	(define pre-event-list)
(let 
	(
	(my1 (list '(active . #t) (cons 'name "report-event-for-cell-count-at-0CA") (cons 'start 0) (list 'action 15 (format #f "/file start-tran ../Report/ice-cell-count-0")))) 
	(my2 (list '(active . #t) (cons 'name "report-event-for-cell-count-at-0CA") (cons 'start 0) (list 'action 15 (format #f "/mesh size-info")))) 
	(my3 (list '(active . #t) (cons 'name "report-event-for-cell-count-at-0CA") (cons 'start 0) (list 'action 15 (format #f "/file stop-tran"))))
	(my4 (list '(active . #t) (cons 'name "report-event-for-cell-count-at-180CA") (cons 'start 180) (list 'action 15 (format #f "/file start-tran ../Report/ice-cell-count-180")))) 
	(my5 (list '(active . #t) (cons 'name "report-event-for-cell-count-at-180CA") (cons 'start 180) (list 'action 15 (format #f "/mesh size-info")))) 
	(my6 (list '(active . #t) (cons 'name "report-event-for-cell-count-at-180CA") (cons 'start 180) (list 'action 15 (format #f "/file stop-tran"))))
	)
	(set! cur-event-list (list my1 my2 my3 my4 my5 my6)) 
	(set! pre-event-list (rpgetvar 'dynamesh/in-cyn/event-list))
	;)
	(rpsetvar 'dynamesh/in-cyn/event-list (append pre-event-list cur-event-list))
	;(pp "in here")
	
	)
	)
)
(set! wb-ice-report-urf
(lambda ()
(define valve-list)
(define cur-event-list)
(define pre-event-list)
(define list1)
(set! list1 ())
  (format  "\n\n\t\t ___________Reporting URF _____________\t\n\n")  

  ;set report dir
(wb-ice-report-dir)

(set! valve-list           (read-var ICDyMeshSetFile 'setup_valve_list))
(ti-menu-load-string "/dis set pic driver jpeg")
(do ((i 0 (+ i 1))) ((= i (length valve-list)))
(let 
	(
	;(list1 ())
	(my1 (list '(active . #t) (cons 'name "report-event-for-solver-setting") (cons 'start 0) (list 'action 15 (format #f "/report summ yes ../Report/ice-report-solver-settings-0")))) 
	(my2 (list '(active . #t) (cons 'name "report-event-for-solver-setting") (cons 'start (list-ref (list-ref setup-valve-list i) 3)) (list 'action 15 (format #f "/report summ yes ../Report/ice-report-solver-settings-~s" (list-ref (list-ref setup-valve-list i) 3))))) 
	(my3 (list '(active . #t) (cons 'name "report-event-for-solver-setting") (cons 'start (list-ref (list-ref setup-valve-list i) 4)) (list 'action 15 (format #f "/report summ yes  ../Report/ice-report-solver-settings-~s" (list-ref (list-ref setup-valve-list i) 4))))) 
	)
	(let
	(	
	(p (list-ref (list-ref setup-valve-list i) 2))
	)
	(set! list1 (append list1 (list (string->symbol p))))
	)
	
	
	(set! cur-event-list (list my1 my2 my3))
	(set! pre-event-list (rpgetvar 'dynamesh/in-cyn/event-list))
	(rpsetvar 'dynamesh/in-cyn/event-list (append pre-event-list cur-event-list))
	))

	)
) 
(define mesh-replace-initial 
  (lambda ()
    (set! initial-fluids    (map thread-name (sort-threads-by-name (get-threads-of-type 'fluid))))
    (set! initial-interiors (map thread-name (sort-threads-by-name (get-threads-of-type 'interior))))

    ;; saving original fluid-ch id to call it during copy cell information
    (set! fluid-ch-original (car (inquire-adjacent-threads (thread-name->id "cyl-tri"))))

    (set! initial-mesh-walls (map thread-name (sort-threads-by-name (get-threads-of-type 'wall))))
    ;; clear all monitors
    (ti-menu-load-string "\nsolve/monitors/volume/clear-monitors\n")
    (ti-menu-load-string "solve/monitors/surface/clear-monitors\n")
    ;; turn off in-cylinder output
    (ti-menu-load-string "define/dynamic-mesh/controls/ in-cylinder-out n\n")
    ;; adjust wall boundary conditions
    (if (or (string=? mduc-model "rif") (string=? mduc-model "saci"))
        (set-bc-wall-hflux initial-mesh-walls #f) ) ; turn off enthalpy flux profile (set-bc-wall-hflux define in custom-funcs.scm)
    ;; let the solver know that the models have changed
    (models-changed)
    ;;  ----- end custom ----- ;;;
    ) ; end lambda
  ) ; end mesh-replace-initial

(define mesh-replace-final
  (lambda ()
    (let ( (new-fluid-names '())
           (new-wall-names '())
           (new-interior-names '() )
           (fluid-ch-id  0)
           (fluid-piston-layer1-id 0)
           (fluid-piston-layer2-id 0)
           (zone1-name "")
           (zone2-name "")
           (no-uds 0)
           (strng-no-uds "")
           (res-new-fluid-names '()) ;used to copy cell information
           
           (cylinder-axis     (cadr (assoc 'geom-axis (assoc (thread-name->id "cyl-tri") (rpgetvar 'dynamesh/dynamic-zones)))) )
           (lpiston-full-prof (cadr (assoc 'udf       (assoc (thread-name->id "piston")  (rpgetvar 'dynamesh/dynamic-zones)))))
           ) ;; end let
      
      (define count 0)
      (set! final-fluids     (map thread-name (sort-threads-by-name (get-threads-of-type 'fluid))))
      (set! new-interiors    (map thread-name (sort-threads-by-name (get-threads-of-type 'interior))))
      (set! final-mesh-walls (map thread-name (sort-threads-by-name (get-threads-of-type 'wall))))
      
      ;; get new fluid names
      (do ((i 0 (+ i 1))) ((= i (length final-fluids)))
        (if (not (member (list-ref final-fluids i)  initial-fluids))
            (set! new-fluid-names (append  new-fluid-names (list (list-ref final-fluids i))))
            )
        )
      (display (format #f "new fluid zones found: ~d --> ~a\n" (length new-fluid-names) new-fluid-names))

      ;; store the new fluid names to be used to copy cell information
      (define new-cyl-names '())
      (define new-cyl-ids "")
      (do ((i 0 (+ i 1))) ((= i (length new-fluid-names)))
        (set! res-new-fluid-names (append res-new-fluid-names (list (thread-name->id (list-ref new-fluid-names i)))))
        )

      ;; get new interior zone names
      (do ((i 0 (+ i 1))) ((= i (length  new-interiors )))
        (if (not (member (list-ref new-interiors  i)   initial-interiors))
            (set! new-interior-names (append new-interior-names (list (list-ref  new-interiors i))))
            )
        )
      (display (format #f "new interior zones found: ~d --> ~a\n" (length new-interior-names) new-interior-names))

      ;; store new walls to adjust the boundary conditions (if necessary)
      (do ((i 0 (+ i 1))) ((= i (length final-mesh-walls)))
        (if (not (member (list-ref final-mesh-walls i)  initial-mesh-walls))
            (set! new-wall-names (append new-wall-names (list (list-ref final-mesh-walls i))))
            )
        )
      (display (format #f "new walls found: ~d --> ~a\n" (length new-wall-names) new-wall-names))
      
      ;; only do the this for piston insertion
      (if (> (length new-fluid-names) 0)
          (begin ;;
            (display "checking zone names after piston layer insertion:\n")
            
            ;; get the thread ids of zones adjacent to cylinder head and piston
            (set! fluid-ch-id            (car (inquire-adjacent-threads (thread-name->id "cyl-tri"))))
            (set! fluid-piston-layer2-id (car (inquire-adjacent-threads (thread-name->id "piston"))))
            (display (format #f "fluid-ch-id: ~d\n" fluid-ch-id))
            (display (format #f "fluid-piston-layer2-id: ~d\n" fluid-piston-layer2-id))

            ;; check to see if zone adjacent to cylinder head is newly created
            (display (format #f "new-fluid-names before cylinder head check: ~a\n" new-fluid-names))
            (if (member (thread-id->name fluid-ch-id) new-fluid-names)
                (begin
                  (set! new-fluid-names (list-remove new-fluid-names (thread-id->name fluid-ch-id)))
                  (set! new-fluid-names (append new-fluid-names '(fluid-ch)));;old name assigned to the new zone
                  )
                )
            (display (format #f "                after cylinder head check: ~a\n" new-fluid-names))

            (if ( not (= 2 (length new-fluid-names)))
                (display "some error in zone assignment; there should be only 2 new fluid zones\n") )
            
            ;; get the id of the layering zone
            (set! fluid-piston-layer1-id (thread-name->id (car (list-remove new-fluid-names (thread-id->name fluid-piston-layer2-id)))))
            (display (format #f "fluid-piston-layer1-id: ~d\n" fluid-piston-layer1-id))

            (do ((i 0 (+ i 1))) ((= i (length  new-interior-names )))
              (define adjthreads)
              (set! adjthreads (inquire-adjacent-threads (thread-name->id  (list-ref new-interior-names  i))))
              (if (and (member fluid-ch-id adjthreads) (member  fluid-piston-layer1-id adjthreads))
                  (set! zone2-name  (list-ref new-interior-names  i))
                  )
              (if (and (member fluid-piston-layer2-id adjthreads) (member  fluid-piston-layer1-id adjthreads))
                  (set! zone1-name  (list-ref new-interior-names  i))
                  )
              )

            (if (or (eq? zone1-name "") (eq? zone2-name ""))
                (display "some error in zone assignment; unable to get the new interior zones\n")
                )
            
            (do ((i 0 (+ i 1))) ((= i (length new-fluid-names)))
              (ti-menu-load-string (format #f "\n/def d-m zones create ~s rigid-body ~s ~s ~s ~s" (list-ref new-fluid-names i) (string->symbol lpiston-full-prof) (list-ref cylinder-axis 0) (list-ref cylinder-axis 1) (list-ref cylinder-axis 2)))
              ) ; end loop over new-fluid-names
            
            (ti-menu-load-string (format #f "\n/def d-m zones create ~s stationary , , 0 , , 0.001" zone2-name))
            (ti-menu-load-string (format #f "\n/def d-m zones delete ~s" zone1-name))

            ;; ----- Commands to copy cell zone information to new cell zones--Start
            ;; Commands to update UDS for all cell zones
            (set! no-uds (rpgetvar 'uds/n-uds))
            
            (do ((i 0 (+ i 1))) ((= i no-uds))
              (set! strng-no-uds (string-append strng-no-uds " , , ,")) )
            
            (if (> no-uds 0)
                (ti-menu-load-string(string-append "/def/ud/uds/ , , , " strng-no-uds)) )
          
            (ti-menu-load-string(format #f "\n/def b-c copy-bc ~s ~s ~s ," fluid-ch-original (list-ref res-new-fluid-names 0) (list-ref res-new-fluid-names 1)))
            ;; ----- Commands to copy cell zone information to new cell zones--Finished

            ;; ----- additional functions EMD ----- ;
            ;; rename the piston and cylinder layers
            (display (format #f "\nchanging zone name of ~a to fluid-ch-piston\n" (thread-id->name fluid-piston-layer2-id)))
            (ti-menu-load-string (format #f "/def bc zone-name ~s fluid-ch-piston\n" fluid-piston-layer2-id))
            (display (format #f "changing zone name of ~a to fluid-ch-layer\n" (thread-id->name fluid-piston-layer1-id)))
            (ti-menu-load-string (format #f "/def bc zone-name ~s fluid-ch-layer\n"  fluid-piston-layer1-id))
            (set! new-cyl-names '(fluid-ch-layer fluid-ch-piston))
            (set! new-cyl-ids    (format #f "~d ~d" fluid-piston-layer1-id fluid-piston-layer2-id))

            ;; rename the interior zones accoridingly
            (do ((i 0 (+ i 1))) ((= i (length  new-interiors )))
              (define adjthreads)
              (set! adjthreads (inquire-adjacent-threads (thread-name->id (list-ref new-interiors i))))
              (cond ( (and (= fluid-piston-layer1-id (list-ref adjthreads 0)) (= fluid-piston-layer1-id (list-ref adjthreads 1)))
                      (ti-menu-load-string (format #f "/def bc zone-name ~a interior-fluid-ch-layer\n"  (list-ref new-interiors i))) )
                    ( (and (= fluid-piston-layer2-id (list-ref adjthreads 0)) (= fluid-piston-layer2-id (list-ref adjthreads 1)))
                      (ti-menu-load-string (format #f "/def bc zone-name ~a interior-fluid-ch-piston\n" (list-ref new-interiors i))) )
                    ) ;; end cond
              )

            ;; rename the cylinder walls accordingly
            (do ((i 0 (+ i 1))) ((= i (length  new-wall-names )))
              (define adjthreads)
              (set! adjthreads (inquire-adjacent-threads (thread-name->id (list-ref new-wall-names i))))
              (display (format #f "threads adjacent to ~a: ~a\n" (list-ref new-wall-names i) adjthreads))
              (if (member fluid-piston-layer1-id adjthreads)
                  (ti-menu-load-string (format #f "/def bc zone-name ~a cyl-layer\n"  (list-ref new-wall-names i)))
                  )
              (if (member fluid-piston-layer2-id adjthreads)
                  (ti-menu-load-string (format #f "/def bc zone-name ~a cyl-piston\n" (list-ref new-wall-names i)))
                  )
              )
            ;; ----- end additional functions EMD ----- ;
            ) ;; end if true
          ) ;; end if for piston insertion
  
      ;; ensure that the top fluid cell is always names fluid-ch (do for both insertion and deletion)
      (define rename-fluid-ch-id 0)
      (set! rename-fluid-ch-id (car (inquire-adjacent-threads (thread-name->id "cyl-tri"))))
      (ti-menu-load-string (format #f "/def bc zone-name ~s fluid-ch\n" rename-fluid-ch-id))

      ;; now rename the interior-fluid-ch 
      (set! new-interiors (map thread-name (sort-threads-by-name (get-threads-of-type 'interior))))
      (do ((i 0 (+ i 1))) ((= i (length new-interiors)))
        (define adjthreads)
        (set! adjthreads (inquire-adjacent-threads (thread-name->id (list-ref new-interiors i))))
        (if (and (= rename-fluid-ch-id (list-ref adjthreads 0)) (= rename-fluid-ch-id (list-ref adjthreads 1)))
            (ti-menu-load-string (format #f "/def bc zone-name ~a interior-fluid-ch\n" (list-ref new-interiors i))) )
        )
      
      ;; ----- additional functions EMD ----- ;
      ;; update the mesh zones and walls again after name change
      (set! final-fluids     (map thread-name (sort-threads-by-name (get-threads-of-type 'fluid))))
      (set! final-mesh-walls (map thread-name (sort-threads-by-name (get-threads-of-type 'wall))))

      ;; get all the thread ids for renaming purposes
      (define final-fluid-id (map thread-id   (sort-threads-by-name (get-threads-of-type 'fluid))))

      ;; re-setup monitors
      (define cyl-fluid-base '(fluid-ch))        ; cell volumes that already exist and that needs to be included in the current monitors
      ;; append new zones/walls to base set
      ;; NOTE: ensure that eng-zone is defined somewhere in the case setup
      (set! eng-cyl        (append cyl-fluid-base new-cyl-names))
      (display (format #f "cylinder zones set to: ~a\n" eng-cyl))
      ;; pass the new threads to the MDUC library for renaming
      (rpsetvar 'mduc/thread-names (format #f "~a" final-fluids))
      (rpsetvar 'mduc/thread-ids   (format #f "~a" final-fluid-id))
      (ti-menu-load-string "(ti-udf-on-demand) \"updateSolverThreadNames::mduc\"")
      ;; reset any custom monitors
      (ti-read-scheme (format #f "~a/journal/monitor.scm" case-root ))
      (ti-menu-load-string
      (format #f "\ndefine/dynamic-mesh/controls/in-cylinder-output y 1 ~a
                  ~a/monitor/cylinder-swirl.~a cg 0 0 1 1 0 0 0 1 0"
               eng-cyl case-dir case-sub-id))
      ;; add boundary conditions back
      (if (or (string=? mduc-model "rif") (string=? mduc-model "saci"))
          (set-bc-wall-hflux final-mesh-walls #t) )
      ;; update models
      (models-changed)
      ;; ----- end additional functions EMD ----- ;

      ) ; end let
    ) ; end lambda
  ) ; end mesh-replace final

;;
;;Function to rename cylinder's top fluid zone to "fluid-ch"
;;
(define name-update-inter-del
  (lambda ()
    (define rename-fluid-ch-id 0)
    (set! rename-fluid-ch-id (car (inquire-adjacent-threads (thread-name->id "cyl-tri"))))
    (ti-menu-load-string (format #f "/def bc zone-name ~s fluid-ch" rename-fluid-ch-id))
    ) ;; end lambda
  ) ;; end name-update-inter-del
;;

(set! wb-ice-plot-profile
(lambda ()
(define valve-list)
(define valve-profile-list)
(define profile-list)
(set! valve-profile-list ())
(set! profile-list ())

  (format  "\n\n\t\t ___________Plotting the Profile_____________\t\n\n")  
;set report dir
(wb-ice-report-dir)
(set! valve-list           (read-var ICDyMeshSetFile 'setup_valve_list))
(ti-menu-load-string "/dis set pic driver jpeg")
(do ((i 0 (+ i 1))) ((= i (length valve-list)))
	(let
	(	
	(p (list-ref (list-ref setup-valve-list i) 2))
	)
	(set! valve-profile-list (append valve-profile-list (list (string->symbol p))))
	)
	)
	(set! profile-list (append valve-profile-list (list (string->symbol piston-full-prof))))
	(if (= piston_type 1)

	  	(set! profile-list (append profile-list (list (string->symbol piston-limit-prof))))
	  )

	(ti-menu-load-string (format #f "/def dm c icp ppl ~s ~s 720 10 yes" profile-list start-crank-angle))
	(ti-menu-load-string "/dis hc ice-profile-plot1.jpg")
	(ti-menu-load-string (format #f "/def dm c icp ppl ~s 0 720 10 yes" valve-profile-list))
	(ti-menu-load-string "/dis hc ice-profile-plot.jpg")
	)
)

(set! wb-ice-report-dir

  (lambda ()
  (define icedir)
  (define report-dir)
  (set! icedir (directory (directory (client-default-basename "r"))))
  (set! report-dir (string-append icedir "\Report"))
  (if (file-exists? report-dir)
    (chdir report-dir)
    (begin
      (chdir icedir)
      (system "mkdir Report")
      (chdir "Report")
      )
    )
  )
)

(set! wb-ice-report-dir-fluent

  (lambda ()
  (define icedir)
  (set! icedir (directory (directory (client-default-basename "r"))))
  (set! icedir (string-append icedir "\Fluent"))
  (chdir icedir)
  
  )
)


  
  
(set! wb-ice-solver-setup

  (lambda ()
    (let(
	 (initialize-flow? (read-var solver-set-file 'initialize-flow?))
	 (report-mesh?  (read-var solver-set-file 'report-mesh? ))
	 (set-default-urf? (read-var solver-set-file 'set-default-urf? ))
	 (set-default-models?  (read-var solver-set-file 'set-default-models? ))
	 (set-bc?  (read-var solver-set-file 'set-bc? ))
	 (report-urf?  (read-var solver-set-file 'report-urf? ))
	 (post-planes?  (read-var solver-set-file 'post_planes ))
	 (swirl-data?  (read-var solver-set-file 'swirl_data ))
	 (user-monitor?  (read-var solver-set-file 'monitor_def ))
	 (auto-save-freq (read-var solver-set-file 'auto-save-freq))
		      )
  (if (OR report-mesh? report-urf?  swirl-data? user-monitor?)
      (wb-ice-report-dir)
      )
  (if set-default-models?
	(wb-ice-viscous-model-setup)
	)
  (if set-default-urf?
	(wb-ice-urf-setting)
	)
  (if set-bc?
    (begin
	  (wb-ice-inlet-outlet-bc)
	  (wb-ice-wall-bc)
	  )
    )

  (if report-urf?
    (begin	
      (wb-ice-report-dir)
      (wb-ice-report-urf)
      (wb-ice-plot-profile)
      (ti-menu-load-string "/plot res-set ptf ice-residual")
      )

    )

  (if initialize-flow?
    (wb-ice-initialization)
    )

 (if report-mesh? 
	(wb-ice-report-mesh)
	)

 (if (and post-planes?)
   (begin
    (wb-ice-create-planes)
    (wb-ice-create-view-for-planes)
    (wb-ice-create-default-images-on-planes)
   
     )
     )

 (if swirl-data?
     (wb-ice-swirl-data)
     )
 (if user-monitor?
     (wb-ice-postproc-monitors-user)
     )
 (if auto-save-freq

   (ti-menu-load-string (format #f "/file/auto-save/data-frequency ~s " auto-save-freq))
   )
 (wb-ice-report-dir)

	)
	)

)


(set! wb-setup-ic-case
  (lambda (ic-prop-file)

    (format  "\n\n\t\t ___________Setting Up WB IC Cold Flow Case ..._____________\t\n\n")  
    (set! ICDyMeshSetFile ic-prop-file)
    (ice-read-data)
    (wb-ice-driver)

    )
)


(set! wb-ic-register-post-journal
  (lambda()
    (let()
     (cx-add-dependent 'calculation-complete #f wb-ice-post-iteration-method);register the post-iteration journal file

     )

    )
  )


(set! wb-ic-write-custom-files
  (lambda()
    (let((outputport #f))
     (wb-ice-report-dir)
     (set! outputport (open-output-file "ice-custom-file-list"))
     (if outputport
       (begin
	 (format outputport "user-bc-file		~a~%" wb-ice-user-bc-file )
	 (format outputport "pre-iter-file		~a~%" wb-ice-pre-iteration-journal )
	 (format outputport "post-iter-file		~a~%" wb-ice-post-iteration-journal )

	 (close-output-port outputport)
	 )
       )

     )

    )
  )

(set! wb-setup-ic-solver
  (lambda (ic-solver-set-file)

    (if (and ic-solver-set-file (file-exists? ic-solver-set-file))
      (begin
	(set! solver-set-file ic-solver-set-file)
	(wb-ice-solver-setup)
	(format  "\n\n\t\t ___________WB IC Cold Flow Solver Setup Complete ______________\t\n\n")
	)
      )

    )
)

(set! wb-ice-post-iteration-method
  (lambda ()

    (if (and wb-ice-post-iteration-journal (file-exists? wb-ice-post-iteration-journal))
      (begin
	(format "\n\n\t\t _______________ Executing Post Iteration Journal  ~s  _______________\n\n" wb-ice-post-iteration-journal)
	(ti-menu-load-string (string-append "/file/read-journal \"" wb-ice-post-iteration-journal "\""))


	(format "\n\n\t\t _______________ Executed Post Iteration Journal  ~s  _______________\n\n" wb-ice-post-iteration-journal)
	(wb-ice-report-dir-fluent)

	)
      )


    )

  )
