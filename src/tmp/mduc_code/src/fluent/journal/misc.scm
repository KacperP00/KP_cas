(define (retrieve-zones match-type)
  (for-each
   (lambda (name)
     (if (string=? (format #f "~a" (zone-type (get-zone name))) match-type)
         (display (format #f "~a " name)) ())
     )
   (inquire-zone-names)
   )
  )
