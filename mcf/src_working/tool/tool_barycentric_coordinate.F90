  SUBROUTINE tool_barycentric_coordinate(this, xp, pp, lambda, stat_info)
        !----------------------------------------------------
        ! Subroutine  : tool_barycentric_coordinate
        !----------------------------------------------------
        !
        ! Purpose     : Calculate barycentric coordinates
        !
        ! Routines    :
        !
        ! Remarks     : line, triangle
        !
        ! References  :
        !
        ! Revisions   : V0.1 05.02 2013, original version.
        !
        !----------------------------------------------------
        ! Author       : Babak Gholami
        ! Contact      : babak.gholami@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments     
        !----------------------------------------------------
        
        TYPE(Tool), INTENT(IN)                  :: this
        REAL(MK), DIMENSION(:,:),INTENT(IN)     :: xp
        REAL(MK), DIMENSION(:),INTENT(IN)       :: pp
        REAL(MK), DIMENSION(:),INTENT(OUT)      :: lambda
        INTEGER, INTENT(OUT)                    :: stat_info

        INTEGER                                 :: num_dim, shape
        REAL(MK), ALLOCATABLE, DIMENSION(:)     :: v0, v1, v2
        REAL(MK)                                :: dot00, dot01, dot02, dot11, dot12 
        REAL(MK)                                :: invDenom
        
        !----------------------------------------------------
        ! Initialization
        !
        ! This is supposed to be used, otherwise,
        ! compiler complains that it is not used.
        !----------------------------------------------------
        
        stat_info = this%flag
        stat_info = 0
        
        
        num_dim = SIZE(xp,1)
        shape = SIZE(xp,2)
        
        IF (shape == 2) THEN ! line


           allocate(v0(num_dim), v1(num_dim))
           
           ! Compute vectors (assumption: line sections (what else could it be?)   
           v0 = xp(1:num_dim, 2) - xp(1:num_dim, 1)
           v1 = pp(1:num_dim)    - xp(1:num_dim, 1)
           
           ! Compute dot products
           dot00 = DOT_PRODUCT(v0, v0)
           dot01 = DOT_PRODUCT(v0, v1) 

           ! Compute barycentric coordinates (!)
           invDenom = 1.0_MK / dot00
           lambda(2) = dot01 * invDenom
           lambda(1) = 1.0_MK - lambda(2)
           
           
        ELSE ! triangle
        

           allocate(v0(num_dim), v1(num_dim), v2(num_dim))
           

           ! Compute vectors (assumption: facets are triangles)        
           v0 = xp(1:num_dim, 2) - xp(1:num_dim, 1)
           v1 = xp(1:num_dim, 3) - xp(1:num_dim, 1)
           v2 = pp(1:num_dim)    - xp(1:num_dim, 1)

           ! Compute dot products
           dot00 = DOT_PRODUCT(v0, v0)
           dot01 = DOT_PRODUCT(v0, v1)
           dot02 = DOT_PRODUCT(v0, v2)
           dot11 = DOT_PRODUCT(v1, v1)
           dot12 = DOT_PRODUCT(v1, v2)

           ! Compute barycentric coordinates
           invDenom = 1.0_MK / (dot00 * dot11 - dot01 * dot01)
           lambda(3) = (dot11 * dot02 - dot01 * dot12) * invDenom ! u
           lambda(2) = (dot00 * dot12 - dot01 * dot02) * invDenom ! v
           lambda(1) = 1.0_MK - lambda(2) - lambda(3)

        END IF

        !----------------------------------------------------
        ! Return.
        !----------------------------------------------------
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE tool_barycentric_coordinate
      
      
