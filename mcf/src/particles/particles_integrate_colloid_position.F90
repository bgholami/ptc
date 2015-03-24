      SUBROUTINE particles_integrate_colloid_position(this,&
           dt,accuracy_order,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_integrate_colloid_position
        !----------------------------------------------------
        !
        ! Purpose     : Integrate the positions of colloidal
        !               boundary particles.
        !
        ! Reference   :
        !
        ! Remark      : Colloid are modelled as rigid body,
        !               therefore, position of colloidal 
        !               boundary particles are integrated
        !               according to colloid's translation
        !               and rotation velocity.
        !
        ! Revision    : V0.1  23.08.2010, original version.
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
        ! Arguments
        !
        ! this           : an object of Particles Class.
        ! dt             : time step.
        ! accuracy_order : accuracy required.
        ! stat_info      : return flag of status.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        REAL(MK), INTENT(IN)                    :: dt
        INTEGER, INTENT(IN)                     :: accuracy_order
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables
	!----------------------------------------------------

        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim, i
        INTEGER                                 :: ip, sid
        TYPE(Colloid),POINTER                   :: colloids
        LOGICAL                                 :: translate
        LOGICAL                                 :: rotate
        REAL(MK), DIMENSION(3)                  :: coll_x, coll_v
        REAL(MK), DIMENSION(:,:,:), POINTER     :: coll_rot_matrix
        REAL(MK), DIMENSION(3)                  :: tx, rx_a, rx_b
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        translate     = .FALSE.
        rotate        = .FALSE.
        
        dim = physics_get_num_dim(this%phys,stat_info_sub)
        
        NULLIFY(colloids)
        CALL physics_get_colloid(this%phys,colloids,stat_info_sub)
        translate = &
             colloid_get_translate(colloids,stat_info_sub)
        rotate    = &
             colloid_get_rotate(colloids,stat_info_sub)
        NULLIFY(coll_rot_matrix)
        
#ifdef __POSITION_FIXED
#else
        
        IF ( rotate ) THEN
           
           !-------------------------------------------------
           ! Get rotation matrix of colloids.
           !-------------------------------------------------
           
           CALL colloid_get_rotation_matrix(colloids,&
                coll_rot_matrix,stat_info_sub)
           
           IF ( stat_info_sub /= 0 ) THEN
              PRINT *, "particles_integrate_colloid_position : " ,&
                   "Getting colloid rot_matrix failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF ! rotate
        
        !----------------------------------------------------
        ! Select accuracy oder
        ! 1 velocity contribution;
        ! 2 velcoity + acceleration contribution. 
        !----------------------------------------------------
        
        SELECT CASE (accuracy_order)            
           
        CASE (1)
           
           IF ( translate .OR. rotate ) THEN
              
              !----------------------------------------------
              ! Loop over all colloidal boundary particles.
              !----------------------------------------------
              
              DO i = 1, this%num_part_colloid
                 
                 !-------------------------------------------
                 ! Get index of this boundary particle
                 ! and its species ID.
                 !-------------------------------------------
                 
                 ip  = this%part_colloid_list(1,i)
                 sid = this%part_colloid_list(2,i)
                 
                 !-------------------------------------------
                 ! Get the nearest image of colloid center.
                 !-------------------------------------------

                 CALL colloid_in_nearest_image(colloids, &
                      this%x(1:dim,ip),sid, coll_x,rx_a,&
                      coll_v,stat_info_sub)
                 
                 IF ( stat_info_sub /= 0 ) THEN
                    PRINT *, "particles_integrate_colloid_position : " ,&
                         "Getting colloid nearest image failed !"
                    stat_info = -1
                    GOTO 9999
                 END IF
                 
                 !-------------------------------------------
                 ! record displacement according to
                 ! translation velocity.
                 !-------------------------------------------
                 
                 tx(1:dim) = coll_v(1:dim) * dt
                 
                 IF ( rotate ) THEN
                    
                    !----------------------------------------
                    ! get new relative position according to 
                    ! rotation motion and add it up.e
                    !----------------------------------------
                    
                    rx_b(1:dim) = &
                         MATMUL(coll_rot_matrix(1:dim,1:dim,sid),rx_a(1:dim))
                 
                    this%x(1:dim,ip) = &
                         this%x(1:dim,ip) + rx_b(1:dim) - rx_a(1:dim)

                 END IF ! rotate
                 
                 !-------------------------------------------
                 ! Add up translation.
                 !-------------------------------------------
                 
                 this%x(1:dim,ip) = this%x(1:dim,ip) + tx(1:dim)
                 
              END DO ! i =1, num_part_colloid
              
           END IF ! translate OR rotate
           
           !-------------------------------------------------
           ! Other integrations.
           !-------------------------------------------------
           
        CASE DEFAULT
           PRINT *, "particles_integrate_position : ", &
                "Order of accuracy not available !"
           stat_info = -1
           GOTO 9999
           
        END SELECT ! accuracy_order

#endif        
        
9999    CONTINUE

        IF ( ASSOCIATED(coll_rot_matrix) ) THEN
           DEALLOCATE(coll_rot_matrix)
        END IF
        
        RETURN
        
      END SUBROUTINE particles_integrate_colloid_position
      
      
