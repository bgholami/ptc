      SUBROUTINE colloid_check_boundary_particle(this,&
           p_x,l_sur,l_out,l_in,c_sid,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_check_boundary_particle
        !----------------------------------------------------
        !
        ! Purpose     : check if a particle is inside a
        !               colloid geometry or 
        !               surrounding a colloid surface;
        !               
        !
        ! Reference   : ellipse in Wikipedia.
        !
        ! Remark      : 1 :
        !               In case of periodic or Lees-Edwards
        !               boundaries, the images of colloid's
        !               center has to be taken into account
        !               to decide if a potential boundary
        !               particle is inside the geometry of
        !               a colloid.
        !               For 2D, maximum 3**2=9->8 images;
        !               For 3D, maximum 3**3=27->26 images.
        !              
        !               2 :
        !               In order to prevent the fluid 
        !               particle being too close to surface
        !               of a colloid, 'dout' restricts distance
        !               for the fluid particle from surface.    
        !
        !
        ! Revisions   : V0.2 29.03 2010, rectangular colloid
        !               (Gholami)   
        !
        !               V0.1 10.12 2009, original version.
        !
        !----------------------------------------------------
        ! Author      : Xin Bian
        ! Contact     : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments :
        !
        ! Input
        !
        ! this      : object of colloid class.
        ! p_x       : potential boundary particle's position.
        !
        ! Output  
        !
        ! l_sur     : particle is inside geometry surface.
        ! l_out     : particle is inside geometry surface+dout.
        ! l_in      : particle is inside inner ring, which
        !             is(probaly more than) cut off far 
        !             from surface.
        ! c_sid      : colloid's ID.
        !
        ! stat_info : status of this routine.
        !----------------------------------------------------
        
        TYPE(colloid), INTENT(IN)               :: this
        REAL(MK),DIMENSION(:),INTENT(IN)        :: p_x
        LOGICAL, INTENT(OUT)                    :: l_sur
        LOGICAL, INTENT(OUT)                    :: l_out
        LOGICAL, INTENT(OUT)                    :: l_in
        INTEGER, INTENT(OUT)                    :: c_sid
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here.
        !
        ! rp_x : relative position vector to center of
        !        the colloid.
        ! d_pc : distance of potential boundary particle to
        !        the center of the colloid.
        ! d_sc : distance of from surface to the center of
        !        the colloid.
        ! d_ps : shortest distance of particle to the surface.
        ! theta: the angel between vertor from particle to
        !        center of colloid and x+ direction.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub

        INTEGER                                 :: dim
        REAL(MK), DIMENSION(3)                  :: x_coll
        REAL(MK), DIMENSION(3)                  :: rp_x
        REAL(MK), DIMENSION(3)                  :: v_coll

        REAL(MK)                                :: d_pc
        REAL(MK)                                :: d_sc
        REAL(MK)                                :: d_ps
        REAL(MK)                                :: d_ar 
        REAL(MK), DIMENSION(:),ALLOCATABLE      :: n_ar
        INTEGER                                 :: f_ar, col_id
        REAL(MK)                                :: theta
        REAL(MK), DIMENSION(2)                  :: s_x
        
        INTEGER                                 :: i
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     =  0
        stat_info_sub =  0
        
        l_sur    = .FALSE.
        l_out    = .FALSE.
        l_in     = .FALSE.
        c_sid    = 0
        dim      = this%num_dim
        rp_x(:)  = 0.0_MK

        ALLOCATE(n_ar(dim))
        
        DO i = 1, this%num_colloid
           
           !-------------------------------------------------
           ! Consider the relative positions of nearest 
           ! image centers, including box(cell) itself.
           ! Get relative displacement of potential boundary 
           ! particle to the center of nearest image.
           !-------------------------------------------------
           
           CALL colloid_nearest_image(this,p_x(1:dim),i,&
                x_coll(1:dim),rp_x(1:dim),v_coll(1:dim),&
                stat_info_sub)
           
           !----------------------------------------------
           ! Calculate the distance of p from the center
           ! of a colloid, or its image.
           ! Then calculate the angle from x+ direction.
           !----------------------------------------------
           
           d_pc = SQRT(DOT_PRODUCT(rp_x(1:dim),rp_x(1:dim)))
           
           theta = polar_angle(rp_x(1),rp_x(2))
           
           !-------------------------------------------------
           ! According to different shapes.
           !-------------------------------------------------
           
           SELECT CASE ( this%shape(i) )
              
           CASE ( mcf_colloid_shape_sphere )
                 
              !----------------------------------------------
              ! 2D Disk/ 3D sphere.
              !----------------------------------------------
              
              IF ( d_pc <= this%radius(1,i) + this%dout ) THEN
                 
                 l_out = .TRUE.
                 
                 IF ( d_pc <= this%radius(1,i) ) THEN
                    
                    l_sur = .TRUE.
                    
                    IF ( d_pc <= this%radius(1,i) - this%din ) THEN
                       
                       l_in = .TRUE.
                       
                    END IF ! d_pc <= ra - din
                    
                 END IF ! d_pc <= ra
                 
              END IF ! d_pc <= ra + dout
              
              
           CASE ( mcf_colloid_shape_ellipsoid )
              
              !----------------------------------------------
              ! 2D Ellipse for now. 3D not implemented.
              !----------------------------------------------
              
              !----------------------------------------------
              ! Get distance of the point at angle theta
              ! on the surface to the center.
              !----------------------------------------------
              
              d_sc = polar_ellipse_r(this%radius(1,i),&
                      this%radius(2,i),theta,this%theta(3,i))
                 
              !----------------------------------------------
              ! particles is closer than the surface point.
              !----------------------------------------------
              
              IF ( d_pc <= d_sc ) THEN
                 
                 l_sur = .TRUE.
                 
                 !-------------------------------------------
                 ! Get distance of the point at angle theta
                 ! on the inner ring surface to the center.
                 !-------------------------------------------
                 
                 CALL cartesian_ellipse_shortestD(this%radius(1,i),&
                      this%radius(2,i),this%theta(3,i), &
                      rp_x(1),rp_x(2),s_x(1),s_x(2),d_ps)
                 
                 !-------------------------------------------
                 ! particles is closer than the inner ring
                 ! surface point.
                 !-------------------------------------------
                 
                 IF ( d_ps > this%din ) THEN
                    
                    l_in = .TRUE.
                    
                 END IF
                 
              END IF  ! l_sur
              
           CASE ( mcf_colloid_shape_star )
              
              !----------------------------------------------
              ! 2D Star with different frequency.
              !
              ! At angel theta, caculate the point on the
              ! surface to the center.
              !----------------------------------------------
              
              d_sc = polar_star_r(this%radius(1,i),&
                   this%radius(2,i), &
                   REAL(this%freq(i),MK),theta,this%theta(3,i))
              
              !----------------------------------------------
              ! IF the particle is closer than the point on 
              ! the surface at same angel theta, it is inside.
              !----------------------------------------------
              
              IF ( d_pc <= d_sc ) THEN
                 
                 l_sur = .TRUE.
                 
                 !-------------------------------------------
                 ! At angel theta, caculate the point on the
                 ! inner ring surface to the center.
                 !-------------------------------------------
                 
                 CALL  polar_star_shortestD(this%radius(1,i), &
                      this%radius(2,i),REAL(this%freq(i),MK),&
                      this%theta(3,i),rp_x(1),rp_x(2),&
                      s_x(1),s_x(2),d_ps,stat_info_sub )
                 
                 IF ( stat_info_sub /= 0 ) THEN
                    PRINT *, "colloid_check_boundary_particle : ", &
                         "Finding star shortest D wrong !"
                    stat_info = -1
                    GOTO 9999                    
                 END IF
                 
                 !-------------------------------------------
                 ! Particle is close than the point on the 
                 ! inner ring surface at same angel theta.
                 !-------------------------------------------
                 
                 IF ( d_ps >= this%din ) THEN
                    
                    l_in = .TRUE.
                    
                 END IF
                 
              END IF  ! l_sur

           
           CASE ( mcf_colloid_shape_rectangle )
                 
              !----------------------------------------------
              ! 2D rectangle/ 3D rectangular cube. (just 2D for now)
              !----------------------------------------------
              
              ! you can do a smarter impl.
              IF ( ABS(rp_x(1)) <= this%radius(1,i) / 2.0_MK + this%dout &
                  .AND. ABS(rp_x(2)) <= this%radius(2,i) / 2.0_MK + this%dout ) THEN
                 
                 l_out = .TRUE.
                 
                 IF ( ABS(rp_x(1)) <= this%radius(1,i) / 2.0_MK  &
                     .AND. ABS(rp_x(2)) <= this%radius(2,i) / 2.0_MK ) THEN
                    
                    l_sur = .TRUE.
                    
                    IF ( ABS(rp_x(1)) <= this%radius(1,i) / 2.0_MK - this%din &
                        .AND. ABS(rp_x(2)) <= this%radius(2,i) / 2.0_MK - this%din ) THEN
                       
                       l_in = .TRUE.
                       
                    END IF
                    
                 END IF
                 
              END IF


           CASE ( mcf_colloid_shape_arbitrary )
                 
              !----------------------------------------------
              ! arbitrary colloiod shape
              !---------------------------------------------- 

              !col_id = i
              col_id = 0
              
              CALL colloid_arbitrary_distance(this, &
                   .FALSE., col_id, p_x, d_ar, n_ar, f_ar, stat_info_sub)
              

              IF ( d_ar <= this%dout ) THEN
                 
                 l_out = .TRUE.
                 
                 IF ( d_ar <= 0.0_MK ) THEN
                    
                    l_sur = .TRUE.
                    
                    IF ( d_ar <= -this%din ) THEN
                       
                       l_in = .TRUE.
                       
                    END IF
                    
                 END IF
                 
              END IF

           END SELECT ! shape
           
           !-------------------------------------------------
           ! If p belong to one colloid image, 
           ! found and quit !
           !-------------------------------------------------
           
           IF ( l_sur ) THEN
              
              c_sid  = i
              
              EXIT
              
           END IF
           
        END DO ! i = 1, num_colloid
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  colloid_check_boundary_particle
      
      
