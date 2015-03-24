      SUBROUTINE physics_adjust_parameters(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : physics_adjust_parameters
        !----------------------------------------------------
        !
        ! Purpose     : Adjust physics parameters after
        !               they have values, and resolve
        !               the small conflicts.
        !      
        ! Reference   :
        !
        ! Remark      :
        !              1)This is a rouinte of Class Physics,
        !               therefore, the variables of physics
        !               object can be accessed directly, 
        !               however, some pointer variables
        !               are easily allocated if we call 
        !               "_get_" routines.
     
        !
        ! Revisions   : V0.1 23.07.2009, original version.
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
        !----------------------------------------------------
        
        TYPE(Physics), INTENT(INOUT)    :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        LOGICAL                         :: relax_run
        LOGICAL                         :: read_external
        LOGICAL                         :: Brownian
        INTEGER                         :: cc_lub_type
        INTEGER                         :: cc_repul_type
        INTEGER                         :: cw_repul_type
        INTEGER                         :: cw_lub_type
        INTEGER                         :: adaptive_dt
        INTEGER                         :: num_species
        INTEGER                         :: num_dim
        REAL(MK), DIMENSION(:), POINTER :: min_phys
        REAL(MK), DIMENSION(:), POINTER :: max_phys
        INTEGER                         :: wall_rho_type
        REAL(MK)                        :: wall_layer 
        REAL(MK)                        :: IOghost_layer
        REAL(MK)                        :: ewidth
        REAL(MK), DIMENSION(:), POINTER :: min_phys_t
        REAL(MK), DIMENSION(:), POINTER :: max_phys_t      
        REAL(MK), DIMENSION(:), POINTER :: dx
        INTEGER                         :: lattice_type
        INTEGER, DIMENSION(:), POINTER  :: num_part_dim
        INTEGER, DIMENSION(:), POINTER  :: num_part_dim_t
        INTEGER                         :: num_part_tot
        INTEGER, DIMENSION(6)           :: num_part_wall
        REAL(MK)                        :: cut_off
        INTEGER, DIMENSION(:), POINTER  :: bcdef
        INTEGER                         :: num_colloid

        INTEGER                         :: i
        REAL(MK)                        :: coff  

        !----------------------------------------------------
        ! Colloid parameters (and related):
        !----------------------------------------------------

        TYPE(Colloid), POINTER                  :: colloids
        INTEGER                                 :: coll_arbitrary_num 
        INTEGER , DIMENSION(:)    , POINTER     :: coll_f_num
        REAL(MK), DIMENSION(:,:,:), POINTER     :: coll_v 
        INTEGER , DIMENSION(:,:,:), POINTER     :: coll_f_vlist
        REAL(MK), DIMENSION(:,:,:), POINTER     :: coll_n

        REAL(MK), DIMENSION(:), ALLOCATABLE     :: ext_pos, pos_min, pos_max
        REAL(MK)                                :: final_ext
        INTEGER                                 :: final_num
        INTEGER                                 :: vi, fi, ii, j
        INTEGER                                 :: num_vperf


        !----------------------------------------------------
        ! Boundary parameters (and related):
        !----------------------------------------------------
        TYPE(Boundary), POINTER                 :: d_boundary  
        INTEGER                                 :: num_inout
        INTEGER , DIMENSION(:)  , POINTER       :: num_iopoints
        REAL(MK), DIMENSION(:,:), POINTER       :: iopatch_n   
        REAL(MK), DIMENSION(:,:,:), POINTER     :: iopatch_x
        INTEGER, DIMENSION(6)                   :: num_part_IOghost
        


        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        NULLIFY(min_phys)
        NULLIFY(max_phys)
        NULLIFY(min_phys_t)
        NULLIFY(max_phys_t)      
        NULLIFY(dx)
        NULLIFY(num_part_dim)
        NULLIFY(num_part_dim_t)
        NULLIFY(bcdef)  

        NULLIFY(colloids)
        NULLIFY(coll_f_num) 
        NULLIFY(coll_v) 
        NULLIFY(coll_f_vlist) 
        NULLIFY(coll_n)  


        NULLIFY(num_iopoints)
        NULLIFY(iopatch_n)  
        NULLIFY(iopatch_x)     

        num_part_wall(1:6) = 0
        
        coff = 100.0_MK        
        
        !----------------------------------------------------
        ! Get min_phys_t, max_phys_t values from
        ! min_phys and max_phys as basics.        
        !----------------------------------------------------
        
        relax_run    = &
             control_get_relax_run(this%ctrl,stat_info_sub)
        read_external= &
             control_get_read_external(this%ctrl,stat_info_sub)
        Brownian     = &
             control_get_Brownian(this%ctrl,stat_info_sub)
        adaptive_dt  = &
             control_get_adaptive_dt(this%ctrl,stat_info_sub)
       
        num_species  = this%num_species
        num_dim      = this%num_dim
        CALL physics_get_min_phys(this,min_phys,stat_info_sub)
        CALL physics_get_max_phys(this,max_phys,stat_info_sub)
        CALL physics_get_min_phys(this,min_phys_t,stat_info_sub)
        CALL physics_get_max_phys(this,max_phys_t,stat_info_sub) 
        lattice_type = this%lattice_type
        CALL physics_get_num_part_dim(this,num_part_dim, stat_info_sub)
        CALL physics_get_num_part_dim(this,num_part_dim_t, stat_info_sub)
        CALL physics_get_dx(this,dx,stat_info_sub)
        cut_off      = this%cut_off
        CALL physics_get_bcdef(this,bcdef,stat_info_sub)
        
        wall_rho_type = &
             boundary_get_wall_rho_type(this%boundary,stat_info_sub)
        
        !----------------------------------------------------
        ! If one species, should be no colloid.
        !----------------------------------------------------
        
        IF ( num_species == 1 ) THEN
           
           this%num_colloid = 0
           
        END IF
        
        num_colloid  = this%num_colloid
        
        !----------------------------------------------------
        ! Adjust according to different lattice type.
        !----------------------------------------------------
        
        SELECT CASE( lattice_type )
           
           !-------------------------------------------------
           ! If we are using square 2D or simple cubic 3D,
           ! the numbers of particles in each dimension are 
           ! known, then we need to calculate dx.
           ! 
           ! For 2D staggered grid or 3D body center grid
           ! it is the same.
           !-------------------------------------------------
           
        CASE (1:2)
           
           dx(1:num_dim) = &
                (max_phys(1:num_dim) - min_phys(1:num_dim))/ &
                num_part_dim(1:num_dim)
           
           CALL physics_set_dx(this,dx,stat_info_sub)
           
        CASE (3)
           
           !-------------------------------------------------
           ! 2D Hexagonal lattice.
           !------------------------------------------------
           
           IF ( num_dim == 2 ) THEN
              
              dx(2) = &
                   (max_phys(2) - min_phys(2)) / num_part_dim(2)
              
              dx(1) = dx(2) * SQRT(3.0_MK)
              
              !num_part_dim(1) = &
              !     CEILING((max_phys(1)-min_phys(1)) / dx(1))
              
              !----------------------------------------------
              ! Trying to shrink a little bit in x direction,
              ! in order to make sure particles are evenly 
              ! distributed.
              !----------------------------------------------
              
              !dx(1) = &
              !     (max_phys(1)-min_phys(1))/num_part_dim(1)
              
              CALL physics_set_dx(this,dx,stat_info_sub)
              
              !----------------------------------------------
              ! Re-calculate the number of particles in each
              ! direction in case we have wall,
              ! take "Ceiling" of the integer, which
              ! might be bigger than the real particles
              ! number.
              !----------------------------------------------
              
              num_part_dim(1:2) = CEILING ( &
                   (max_phys_t(1:2) - min_phys_t(1:2)) / &
                   dx(1:2) )
              
              num_part_dim_t(1:2) = num_part_dim(1:2)

              !----------------------------------------------
              ! 3D face centered lattice.
              !----------------------------------------------
              
           ELSE IF(num_dim == 3) THEN
              
           END IF
           
        CASE DEFAULT
           
           PRINT *, "physics_adjust_parameters : ", &
                "lattice type not available !"
           stat_info = -1
           GOTO 9999
           
        END SELECT ! lattice_type
        
        
        !-------------------------------------------------
        ! 1 Calculate the nubmer of particles in solid
        !   wall, which is handeld by mcf, not by ppm.
        !
        ! 2 Extend the total physical domain to
        !   include the wall particles.
        !
        ! 3 Add up the wall particles to total
        !   particles in each dimension.
        !-------------------------------------------------
        
        DO i = 1, num_dim
           
           IF ( bcdef(2*i-1) == mcf_bcdef_wall_solid ) THEN
              
              wall_layer = mcf_wall_layer_coeff * cut_off
              
              IF ( wall_rho_type == mcf_wall_rho_type_dynamic ) THEN
                 
                 wall_layer = wall_layer * 2.0_MK
                 
              END IF
              
              num_part_wall(2*i-1) = &
                   CEILING(wall_layer / dx(i))
              
              min_phys_t(i) = min_phys_t(i) - &
                   num_part_wall(2*i-1) * dx(i)
              
              num_part_dim_t(i) = num_part_dim_t(i) + &
                   num_part_wall(2*i-1)
              
           END IF ! bcdef(2i-1)
           
           IF ( bcdef(2*i) == mcf_bcdef_wall_solid ) THEN
              
              wall_layer = mcf_wall_layer_coeff * cut_off
              
              IF ( wall_rho_type == mcf_wall_rho_type_dynamic ) THEN
                 
                 wall_layer = wall_layer * 2.0_MK
                 
              END IF
              
              num_part_wall(2*i) = &
                   CEILING(wall_layer / dx(i))
              
              max_phys_t(i) = max_phys_t(i) + &
                   num_part_wall(2*i) * dx(i)
              
              num_part_dim_t(i) = num_part_dim_t(i) + &
                   num_part_wall(2*i)
              
           END IF ! bcdef(2i)

           
        END DO ! i = 1, num_dim


        !-------------------------------------------------
        ! Repeat the with arbitrary geometry
        !-------------------------------------------------  
        
        ALLOCATE(ext_pos(num_dim))
        ALLOCATE(pos_min(num_dim), pos_max(num_dim))
        pos_min = min_phys_t
        pos_max = max_phys_t

        ! Get colloid
        num_colloid = &
             physics_get_num_colloid(this, stat_info_sub)
        
        coll_arbitrary_num = 0 

        IF ( num_colloid > 0 ) THEN
           
           CALL physics_get_colloid(this, colloids, stat_info_sub) 
           
           coll_arbitrary_num = colloid_get_arbitrary_num(colloids, stat_info_sub)
           
        END IF


        IF (coll_arbitrary_num > 0) THEN  

           !----------------------------------------------------
           ! initialize colloid variables
           !----------------------------------------------------  

           CALL colloid_get_f_num(colloids, coll_f_num, stat_info_sub)
           CALL colloid_get_coll_v(colloids, coll_v, stat_info_sub)  
           CALL colloid_get_coll_f_vlist(colloids, coll_f_vlist, stat_info_sub)  
           CALL colloid_get_coll_n(colloids, coll_n, stat_info_sub) 

           num_vperf = num_dim     

           wall_layer = mcf_wall_layer_coeff * cut_off


           DO j = 1, coll_arbitrary_num

              DO fi = 1, coll_f_num(j)

                 DO ii = 1, num_vperf
                    
                    ! vertex
                    vi = coll_f_vlist(j, ii, fi)

                    ! compute extended position (with -coll_n to point outward)
                    ext_pos(1:num_dim) = coll_v(j, 1:num_dim, vi) - &
                         wall_layer * coll_n(j, 1:num_dim, fi)           

                    ! evaluate min and max in each direction
                    DO i = 1, num_dim
                       pos_min(i) = MIN(pos_min(i), ext_pos(i)) 
                       pos_max(i) = MAX(pos_max(i), ext_pos(i))
                    END DO ! i = 1, num_dim
                           

                 END DO

              END DO

           END DO


        END IF ! coll_arbitrary_num   
     

        !---------------------------------------------------- 
        ! update values (if changed by arbitrary colloid)
        !---------------------------------------------------- 
        DO i = 1, num_dim

           ! min_phys_t
           final_ext = min_phys_t(i) - MIN(min_phys_t(i), pos_min(i))

           final_num = CEILING(final_ext / dx(i))  

           min_phys_t(i) = min_phys_t(i) - &
                final_num * dx(i)

           num_part_dim_t(i) = num_part_dim_t(i) + &
                final_num


           ! max_phys_t
           final_ext =  MAX(max_phys_t(i), pos_max(i)) - max_phys_t(i)

           final_num =  CEILING(final_ext / dx(i))

           max_phys_t(i) = max_phys_t(i) + &
                final_num * dx(i)

           num_part_dim_t(i) = num_part_dim_t(i) + &
                final_num

        END DO ! i = 1, num_dim


        !-------------------------------------------------
        ! And for inflow/outflow patches
        !-------------------------------------------------   
        pos_min = min_phys_t
        pos_max = max_phys_t

        ! get boundary
        CALL physics_get_boundary(this, d_boundary, stat_info_sub)

        num_inout = boundary_get_num_inout(d_boundary, stat_info_sub)  

        IF (num_inout > 0) THEN    

           !----------------------------------------------------
           ! initialize boundary variables
           !----------------------------------------------------

           CALL boundary_get_num_iopoints(d_boundary, num_iopoints, stat_info_sub)
           CALL boundary_get_iopatch_x(d_boundary, iopatch_x, stat_info_sub) 
           CALL boundary_get_iopatch_n(d_boundary, iopatch_n, stat_info_sub)

           IOghost_layer = mcf_IOghost_layer_coeff * cut_off
           ewidth = cut_off

  
           DO j = 1, num_inout

              DO vi =  1, num_iopoints(j) 

                 ! compute extended position (with -iopatch_n to point outward)
                 ext_pos(1:num_dim) = iopatch_x(j, 1:num_dim, vi) - &
                      (IOghost_layer + ewidth) * iopatch_n(j, 1:num_dim)           

                 ! evaluate min and max in each direction
                 DO i = 1, num_dim
                    pos_min(i) = MIN(pos_min(i), ext_pos(i)) 
                    pos_max(i) = MAX(pos_max(i), ext_pos(i))
                 END DO ! i = 1, num_dim

              END DO

           END DO 

        END IF           


        !---------------------------------------------------- 
        ! update values (if changed by arbitrary colloid or inflow/outflow BC)   
        !---------------------------------------------------- 
        DO i = 1, num_dim

           ! min_phys_t
           final_ext = min_phys_t(i) - MIN(min_phys_t(i), pos_min(i))

           final_num = CEILING(final_ext / dx(i))  

           min_phys_t(i) = min_phys_t(i) - &
                final_num * dx(i)

           num_part_dim_t(i) = num_part_dim_t(i) + &
                final_num


           ! max_phys_t
           final_ext =  MAX(max_phys_t(i), pos_max(i)) - max_phys_t(i)

           final_num =  CEILING(final_ext / dx(i))

           max_phys_t(i) = max_phys_t(i) + &
                final_num * dx(i)

           num_part_dim_t(i) = num_part_dim_t(i) + &
                final_num

        END DO ! i = 1, num_dim

        CALL boundary_set_bwidth(d_boundary, IOghost_layer, stat_info_sub)
        CALL boundary_set_ewidth(d_boundary, ewidth, stat_info_sub)
        !----------------------------------------------------
        ! Reset the min_phys_t, max_phys_t.
        ! Reset the number of particles in each dimension.
        ! Reset the number of particles in total.
        !----------------------------------------------------

        CALL physics_set_min_phys_t(this,min_phys_t,stat_info_sub)
        CALL physics_set_max_phys_t(this,max_phys_t,stat_info_sub)
        num_part_tot = PRODUCT(num_part_dim_t(:),1)

        CALL physics_set_num_part_dim(this,num_part_dim, stat_info_sub)
        CALL physics_set_num_part_dim_t(this,num_part_dim_t, stat_info_sub)
        CALL physics_set_num_part_tot(this,num_part_tot, stat_info_sub)
        
        CALL boundary_set_min_phys_t(this%boundary,min_phys_t,stat_info_sub)
        CALL boundary_set_max_phys_t(this%boundary,max_phys_t,stat_info_sub)
        
        !----------------------------------------------------
        ! Set the total physics boundary for colloid.
        ! Set the minimal distance of a fluid particle 
        ! from surface of colloids, if there is any.
        !----------------------------------------------------
           
        IF ( num_colloid > 0 ) THEN
           
           cc_lub_type = &
                control_get_cc_lub_type(this%ctrl,stat_info_sub) 
           cc_repul_type = &
                control_get_cc_repul_type(this%ctrl,stat_info_sub)
           cw_lub_type = &
                control_get_cw_lub_type(this%ctrl,stat_info_sub)
           cw_repul_type = &
                control_get_cw_repul_type(this%ctrl,stat_info_sub)
           CALL colloid_set_cc_lub_type(this%colloids,&
                cc_lub_type,stat_info_sub)
           CALL colloid_set_cc_repul_type(this%colloids,&
                cc_repul_type,stat_info_sub)          
           CALL colloid_set_cw_lub_type(this%colloids,&
                cw_lub_type,stat_info_sub)
           CALL colloid_set_cw_repul_type(this%colloids,&
                cw_repul_type,stat_info_sub)
            CALL colloid_set_min_phys_t(this%colloids, &
                min_phys_t,stat_info_sub)
           CALL colloid_set_max_phys_t(this%colloids, &
                max_phys_t,stat_info_sub)
           CALL colloid_set_dout(this%colloids, &
                dx(1)/2.0_MK,stat_info_sub)
           
        END IF
        
        !----------------------------------------------------
        ! Get absolute value of body force and record it.
        !----------------------------------------------------
        
        this%fa_max = SQRT(DOT_PRODUCT(&
             this%body_force(1:num_dim),this%body_force(1:num_dim))) 


        !----------------------------------------------------
        ! Compute sound speed
        !----------------------------------------------------
        
        CALL physics_compute_cs(this, stat_info_sub)
    
        
        !----------------------------------------------------
        ! Initialize time step.
        !----------------------------------------------------
        
        CALL physics_initialize_dt(this,stat_info_sub)
        
        !-------------------------------------------------
        ! IF step is given from input,
        ! set time end to be negative for adaptive dt.
        !-------------------------------------------------
        
        IF( this%step_end >= 0 ) THEN
           
           IF ( adaptive_dt  > 0 ) THEN
           
              this%time_start = 0.0_MK
              this%time_end   = -1.0_MK
              
           ELSE
              
              this%time_start = this%step_start * this%dt
              this%time_end   = this%step_end * this%dt
              
           END IF
           
           !-------------------------------------------------
           ! IF time is given from input,
           ! set step end to negative for adaptive dt.
           !-------------------------------------------------
           
        ELSE IF ( this%time_end >=0.0_MK ) THEN
              
           IF ( adaptive_dt  > 0 ) THEN
              
              this%step_start = 0
              this%step_end   = -1.0
              
           ELSE
              
              this%step_start = INT(this%time_start / this%dt)
              this%step_end   = CEILING(this%time_end / this%dt)
              
           END IF
         
        END IF ! step_end >= 0.0
           
        
        !----------------------------------------------------
        ! For relax run type 1, set step_relax and time_relax.
        !----------------------------------------------------
        
        IF ( relax_run ) THEN
           
           SELECT CASE ( this%relax_type )
              
           CASE (1)
              
              IF ( this%time_relax > 0 ) THEN
                 
                 this%step_relax = CEILING(this%time_relax / this%dt)
                 
              ELSE IF ( this%step_relax > 0 ) THEN
                 
                 this%time_relax = this%step_relax * this%dt
                 
              END IF
              
              this%disorder_level = 1.0_MK
              
           CASE (2)
              
              this%step_relax = -1
              this%time_relax = -1.0_MK
              
           END SELECT
        
        END IF
        
        !----------------------------------------------------
        ! For non-Brownian simulation, set kt = 0.
        !----------------------------------------------------
        
        IF ( .NOT. Brownian ) THEN
           
           this%kt = 0.0_MK
           
        END IF
        
        !----------------------------------------------------
        ! In case there is periodic or Lees-Edwards
        ! boundary condition, we have to generate the
        ! images of each colloid.
        !----------------------------------------------------
        
        IF ( num_colloid > 0 ) THEN
           
           CALL colloid_adjust_parameters(this%colloids,stat_info_sub)
           CALL colloid_initialize_image(this%colloids,stat_info_sub)
           CALL colloid_compute_image(this%colloids,stat_info_sub)

        END IF
        
9999    CONTINUE       
        
        IF(ASSOCIATED(min_phys))THEN
           DEALLOCATE(min_phys)
        END IF
        
        IF(ASSOCIATED(max_phys))THEN
           DEALLOCATE(max_phys)
        END IF
        
        IF(ASSOCIATED(min_phys_t))THEN
           DEALLOCATE(min_phys_t)
        END IF
        
        IF(ASSOCIATED(max_phys_t))THEN
           DEALLOCATE(max_phys_t)
        END IF
        
        IF(ASSOCIATED(dx))THEN
           DEALLOCATE(dx)
        END IF
        
        IF(ASSOCIATED(num_part_dim)) THEN
           DEALLOCATE(num_part_dim)
        END IF

        IF(ASSOCIATED(num_part_dim_t)) THEN
           DEALLOCATE(num_part_dim_t)
        END IF

        IF(ASSOCIATED(bcdef)) THEN
           DEALLOCATE(bcdef)
        END IF

        IF(ASSOCIATED(coll_f_num)) THEN
           DEALLOCATE(coll_f_num)
        END IF

        IF(ASSOCIATED(coll_v)) THEN
           DEALLOCATE(coll_v)
        END IF

        IF(ASSOCIATED(coll_f_vlist)) THEN
           DEALLOCATE(coll_f_vlist)
        END IF

        IF(ASSOCIATED(coll_n)) THEN
           DEALLOCATE(coll_n)
        END IF 

        IF(ASSOCIATED(num_iopoints)) THEN
           DEALLOCATE(num_iopoints)
        END IF   

        IF(ASSOCIATED(iopatch_x)) THEN
           DEALLOCATE(iopatch_x)
        END IF  
        
        IF(ASSOCIATED(iopatch_n)) THEN
           DEALLOCATE(iopatch_n)
        END IF  
        
        RETURN
        
      END SUBROUTINE physics_adjust_parameters
      
      
