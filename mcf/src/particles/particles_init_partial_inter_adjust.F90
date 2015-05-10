SUBROUTINE particles_init_partial_inter_adjust(this, stat_info)
  !----------------------------------------------------
  ! Subroutine  : particles_init_partial_inter_adjust
  !----------------------------------------------------
  !
  ! Purpose     : adjusting particles after relax run 
  !               in case of arbitrary shape colloids
  !
  ! 
  ! Routines    : 
  !
  !
  ! Remarks     : the code is very similar to 
  !               particles_init_global_inter
  !
  ! References  :
  !
  ! Revisions   : V0.1  29.01 2013, original version.
  !
  !----------------------------------------------------
  ! Author      : Babak Gholami
  ! Contact     : babak.gholami@aer.mw.tum.de
  !
  ! Dr. Marco Ellero's Emmy Noether Group,
  ! Prof. Dr. N. Adams' Chair of Aerodynamics,
  ! Faculty of Mechanical Engineering,
  ! Technische Universitaet Muenchen, Germany.
  !----------------------------------------------------

  !----------------------------------------------------
  ! Modules :
  !----------------------------------------------------

!!$        USE ppm_module_find_duplicates

  !----------------------------------------------------
  ! Arguments :
  !----------------------------------------------------

  TYPE(Particles), INTENT(INOUT)          :: this
  INTEGER, INTENT(OUT)	                :: stat_info

  !----------------------------------------------------
  ! Local variables start here :
  !
  ! num_dim      : number of dimension.
  !----------------------------------------------------

  INTEGER                                 :: stat_info_sub
  INTEGER                                 :: num_dim,i,j

  !----------------------------------------------------
  ! Colloids parameters.
  !----------------------------------------------------

  INTEGER                                 :: num_colloid
  TYPE(Colloid), POINTER                  :: colloids
  INTEGER, DIMENSION(:), POINTER          :: coll_num_numerical_part
  INTEGER                                 :: coll_arbitrary_num        

  !----------------------------------------------------
  ! parameters to handle inflow/outflow for arbitrary geometry
  !----------------------------------------------------
!  TYPE(Colloid)                           :: temp_colloids        
  TYPE(Boundary), POINTER                 :: d_boundary
  INTEGER                                 :: num_inout
  REAL(MK), DIMENSION(:), POINTER         :: dx
  REAL(MK)                                :: avg_dx, cut_off
  INTEGER                                 :: num_layers
  REAL(MK)                                :: distance, bwidth
  INTEGER                                 :: patch_id


  !----------------------------------------------------
  ! Counters & Flags :
  !
  ! num_total : number of total particles.
  ! num       : counter for particles.
  ! num_extra :
  !             number of colloid boundary particles
  !             which are paralle layers to surface.
  !
  ! l_w   : particle is a solid wall boundary particle.
  ! l_sur : particle is inside a colloid surface.
  ! l_out : particle is inside outter ring of a colloid.
  ! l_in  : particle is inside inner right of a colloid.
  ! b_v   : boundary particle's velocity.
  ! b_sid : boundary Particle's species ID.
  !----------------------------------------------------

  INTEGER                                 :: num_total
  INTEGER                                 :: num  
  LOGICAL                                 :: counted_wall_solid
  LOGICAL                                 :: counted_colloid
  LOGICAL                                 :: counted_ignore

  LOGICAL                                 :: l_sur
  LOGICAL                                 :: l_out
  LOGICAL                                 :: l_in
  INTEGER                                 :: p_sid

  !----------------------------------------------------
  ! Intermediate variables.
  !----------------------------------------------------

  REAL(MK), DIMENSION(:,:), POINTER       :: tx
  REAL(MK), DIMENSION(:,:), POINTER       :: tv
  INTEGER, DIMENSION(:,:), POINTER        :: tid 

  !----------------------------------------------------
  ! Check parameters.
  !----------------------------------------------------        

  INTEGER , DIMENSION(:,:), POINTER       :: ide
  INTEGER                                 :: nid


  !----------------------------------------------------
  ! Initialization of variables.
  !----------------------------------------------------

  stat_info     = 0
  stat_info_sub = 0

  num_dim = this%num_dim


  !----------------------------------------------------
  ! Initialization of variables of colloids.
  !----------------------------------------------------

  NULLIFY(colloids)
  NULLIFY(coll_num_numerical_part)    

  NULLIFY(d_boundary)   
  NULLIFY(dx)


  !----------------------------------------------------
  ! Initialization of intermediates variables.
  !----------------------------------------------------

  NULLIFY(tx)
  NULLIFY(tv)
  NULLIFY(tid)


  NULLIFY(ide)

  !----------------------------------------------------
  ! Get physics variables,
  ! species, dimension and which lattice.
  !----------------------------------------------------

  num_dim      = this%num_dim

  !----------------------------------------------------
  ! Colloid parameters :
  !----------------------------------------------------

  num_colloid = &
       physics_get_num_colloid(this%phys,stat_info_sub)

  coll_arbitrary_num = 0

  IF ( num_colloid > 0 ) THEN

     CALL physics_get_colloid(this%phys, &
          colloids,stat_info_sub)           

     ! number of arbitrary shape colloid
     coll_arbitrary_num = colloid_get_arbitrary_num(colloids, &
          stat_info_sub)

     ! make a copy of colloids
     ! We want to temporarily edit parts of colloids, so we make
     ! a copy of it and use it in this routine ONLY. The copy will
     ! be (automatically) deleted at the end of this routine.
!     temp_colloids = colloids

     ! this was supposed to temporarily extend colloid beyond 
     ! min/max_phys to make sure assigment of particle IDs (
     ! i.e. fluid/colloid and inflow/outflow/normal are done
     ! accuratly. However, it seems that with the current 
     ! method to calculate wall distance, there is not need
     ! for extension. Therefore, it is left out for now until
     ! either otherwise is concluded or it is 100% made sure
     ! that the routine is not going to be used.

!!$           ! extend temp_colloids to inflow/outflow BC         
     CALL physics_get_boundary(this%phys,d_boundary,stat_info_sub)    
     bwidth = boundary_get_bwidth(d_boundary,stat_info_sub)
!!$           num_inout = &
!!$                boundary_get_num_inout(d_boundary,stat_info_sub)   
!!$
!!$
!!$           IF (num_inout > 0) THEN 
!!$
!!$              CALL physics_get_dx(this%phys,dx,stat_info_sub)
!!$              cut_off = physics_get_cut_off(this%phys, stat_info_sub)
!!$              avg_dx = SUM(dx) / DBLE(num_dim)
!!$
!!$              num_layers = CEILING(mcf_IOghost_layer_coeff  * cut_off / avg_dx)
!!$
!!$              CALL colloid_adjust_arbitrary_colloid(temp_colloids, &
!!$                   num_layers, stat_info_sub)
!!$
!!$           END IF


  END IF

  ! only and only IF arbitrary shape colloid
  IF (coll_arbitrary_num > 0) THEN

     !----------------------------------------------------
     ! Allocate temporary memory (large enough) of
     ! fluid particles to keep particles.
     !
     ! tx   : position
     ! tv   : velocity
     ! tid  : particle ID,  species ID
     !----------------------------------------------------

     num_total = this%num_part_real

     ALLOCATE(tx(num_dim,num_total), STAT=stat_info_sub)
     ALLOCATE(tv(num_dim,num_total), STAT=stat_info_sub)
     ALLOCATE(tid(this%num_id,num_total), &
          STAT=stat_info_sub)

     IF( stat_info_sub /= 0 ) THEN
        PRINT *, &
             "particles_init_partial_inter_adjust : ", &
             "Allocating memory for variables has problem !"
        stat_info = -1
        GOTO 9999
     END IF

     !----------------------------------------------------
     ! Reset all counters to zero,
     ! prepare for counting number of fluid particles and
     ! number of wall solid boundary particle.
     !----------------------------------------------------

     this%num_part_fluid = 0

     ALLOCATE(coll_num_numerical_part(1:num_colloid))
     coll_num_numerical_part(1:num_colloid) = 0

     num = 0

     !----------------------------------------------------
     ! Loop over all generated particles to decide 
     ! species ID, i.e., which type of particle it is.
     !----------------------------------------------------


     DO j = 1, num_total

        counted_wall_solid = .FALSE.
        counted_colloid    = .FALSE.
        counted_ignore     = .FALSE.

        p_sid          = 0


        !-------------------------------------------------
        ! check if particle is inside any colloid geometry.
        !-------------------------------------------------


        l_sur = .FALSE.
        l_out = .FALSE.
        l_in  = .FALSE.

        IF (this%id(this%sid_idx, j) < 0) THEN ! i.e. if wall 

           counted_wall_solid = .TRUE.

        ELSE

           IF (this%id(this%sid_idx, j) > 0) THEN ! i.e. if colloid...
              ! these are colloid particles on surface 
              ! of the arbitrary colloid, we know this!

              !l_sur = .TRUE.
              !l_out = .FALSE.
              !l_in  = .FALSE.

              counted_colloid = .TRUE.    

              p_sid = this%id(this%sid_idx, j)


           ELSE


              CALL colloid_check_boundary_particle(colloids, &
                   this%x(1:num_dim,j),l_sur,l_out,l_in, p_sid,&
                   stat_info_sub)

              IF ( l_sur ) THEN

                 !-------------------------------------------
                 ! particle j is inside surface of a colloid
                 ! and the placement of it is on lattice.
                 !-------------------------------------------

                 IF( .NOT. l_in ) THEN

                    counted_colloid = .TRUE.

                 ELSE

                    !counted_colloid = .TRUE.
                    counted_ignore = .TRUE.
                    
                 END IF ! .NOT. l_in

              ELSE IF( l_out  ) THEN

                 !-------------------------------------------
                 ! If j is outside of a colloid, but close
                 ! enough to the surface, it should not be 
                 ! created as fluid particle when colloid
                 ! boundary particles are required to be 
                 ! placed on paralle surfaces, since sx would
                 ! be too close to the boundary particles
                 ! on the surface.
                 !-------------------------------------------

                 !counted_ignore = .TRUE.

              END IF ! l_sur

           END IF ! p_id

        END IF ! not wall

        !-------------------------------------------------
        ! Save the particle's information.
        !-------------------------------------------------

        IF ( .NOT. counted_ignore ) THEN


           !----------------------------------------------
           ! Increase counters of different species.
           !----------------------------------------------
           IF (counted_wall_solid) THEN ! write the wall particle  

              num = num + 1
              tx(1:num_dim,num)     = this%x(1:num_dim,j)
              tv(1:num_dim,num)     = 0.0_MK
              tid(this%pid_idx,num) = this%id(this%pid_idx, j)
              tid(this%sid_idx,num) = p_sid
              tid(this%bid_idx,num) = this%id(this%bid_idx, j)

           ELSE ! wall is counted before

              IF ( counted_colloid )  THEN   

                 ! returns distance: 0 for interior, positive distance for inflow, negative distance for outflow 
                 distance = 0.0_MK
                 patch_id = 0
                 CALL boundary_check_particle_stat(d_boundary, &
                      this%x(1:num_dim, j), -1, distance, patch_id, stat_info_sub)


                 IF (ABS(distance) <= bwidth) THEN
                
                    num = num + 1
                    tx(1:num_dim,num)     = this%x(1:num_dim,j)
                    tv(1:num_dim,num)     = 0.0_MK
                    tid(this%pid_idx,num) = this%id(this%pid_idx, j)
                    tid(this%sid_idx,num) = p_sid
                    tid(this%bid_idx,num) = 0

                    coll_num_numerical_part(p_sid) = &
                         coll_num_numerical_part(p_sid) + 1  

                 END IF

              ELSE


                 ! returns distance: 0 for interior, positive distance for inflow, negative distance for outflow
                 distance = 0.0_MK
                 patch_id = 0
                 CALL boundary_check_particle_stat(d_boundary, &
                      this%x(1:num_dim, j), -1, distance, patch_id, stat_info_sub)

                 ! set boundary_id 
                 ! (corresponding d_boundary%patch_id inflow and outflow, 0 for interior)
                 this%id(this%bid_idx, j) = patch_id


                 IF (this%id(this%bid_idx, j) == 0) THEN ! interior fluid particle


                    num = num + 1
                    tx(1:num_dim,num)     = this%x(1:num_dim,j)
                    tv(1:num_dim,num)     = 0.0_MK
                    tid(this%pid_idx,num) = this%id(this%pid_idx, j)
                    tid(this%sid_idx,num) = p_sid
                    tid(this%bid_idx,num) = this%id(this%bid_idx, j)

                    this%num_part_fluid = &
                         this%num_part_fluid + 1

                 ELSE ! inflow/outflow fluid particle
   

                    ! write only if within a distance from its inlte/outlet patch
                    IF (ABS(distance) <= bwidth) THEN              

                       num = num + 1
                       tx(1:num_dim,num)     = this%x(1:num_dim,j)
                       tv(1:num_dim,num)     = 0.0_MK
                       tid(this%pid_idx,num) = this%id(this%pid_idx, j)
                       tid(this%sid_idx,num) = p_sid
                       tid(this%bid_idx,num) = this%id(this%bid_idx, j)

                       this%num_part_fluid = &
                            this%num_part_fluid + 1


                    END IF


                 END IF

              END IF ! counted

           END IF

        END IF ! .NOT. counted_ignore

     END DO ! j = 1, num


     !----------------------------------------------------
     ! Copy the data from temporary varaiable 
     ! to particles members.
     !----------------------------------------------------

     num_total = num

     IF (ASSOCIATED(this%x)) THEN
        DEALLOCATE(this%x,STAT=stat_info_sub)
     END IF
     IF (ASSOCIATED(this%v)) THEN
        DEALLOCATE(this%v,STAT=stat_info_sub)
     END IF
     IF (ASSOCIATED(this%id)) THEN
        DEALLOCATE(this%id,STAT=stat_info_sub)
     END IF

     ALLOCATE(this%x(num_dim,num_total), STAT=stat_info_sub)
     ALLOCATE(this%v(num_dim,num_total), STAT=stat_info_sub)
     ALLOCATE(this%id(this%num_id,num_total),STAT=stat_info_sub)


     this%x(1:num_dim,1:num) = tx(1:num_dim,1:num) 
     this%v(1:num_dim,1:num) = tv(1:num_dim,1:num)
     this%id(1:this%num_id,1:num) = &
          tid(1:this%num_id,1:num)
     !! §§§ Copying in a loop is slower but can handle higher number of particles
     !DO i = 1, num
     !   this%x(1:num_dim,i) = tx(1:num_dim,i) 
     !   this%v(1:num_dim,i) = tv(1:num_dim,i)
     !   this%id(1:this%num_id,i) = &
     !        tid(1:this%num_id,i)
     !END DO



     !----------------------------------------------------
     ! Count the number of numerical particles which 
     ! consititute colloids.
     !----------------------------------------------------

     this%num_part_colloid = 0

     DO j = 1, num_colloid
        this%num_part_colloid = &
             this%num_part_colloid + coll_num_numerical_part(j)
     END DO

     CALL colloid_set_num_numerical_part(colloids,&
          coll_num_numerical_part(:),stat_info_sub)



     !----------------------------------------------------
     ! Set number of real, all and ghost particles.
     !----------------------------------------------------

     this%num_part_real  = num_total
     this%num_part_all   = this%num_part_real
     this%num_part_ghost = 0


!!$           !----------------------------------------------------
!!$           ! Check if each particles is unique.
!!$           !----------------------------------------------------
!!$           
!!$           nid = 0
!!$
!!$           CALL ppm_find_duplicates(this%x,num_dim,&
!!$                this%num_part_real,nid,ide,stat_info_sub)
!!$
!!$           IF ( stat_info_sub /= 0 ) THEN
!!$              PRINT *, &
!!$                   "particles_init_global_inter : ", &
!!$                   "Error by checking duiplicates !"
!!$              stat_info = -1
!!$              GOTO 9999
!!$           END IF
!!$
!!$           IF ( nid /= 0) THEN
!!$              PRINT *, &
!!$                   "particles_init_global_inter : ", &
!!$                   "Found collocating particles !" &
!!$                   ,this%x(:,ide(1,1)), this%x(:,ide(2,1))
!!$              stat_info = -1
!!$              GOTO 9999
!!$           END IF           



     !-------------------------------------------------
     ! Create the rest of quntities, besides,
     !  x,  y(,  z);
     ! vx, vy(, vz);
     ! p_ID, s_ID on each process.
     !-------------------------------------------------

     CALL particles_init_partial_inter(this, &
          stat_info_sub)

     IF( stat_info_sub /=0 ) THEN
        PRINT *, "mcf : ", &
             "Generating quantities locally (after adjustment) failed !"
        stat_info = -1
        GOTO 9999
     END IF

     !-------------------------------------------------
     ! Calculate particles' mass initially,
     ! from density.
     !-------------------------------------------------

     CALL particles_compute_mass(this, &
          stat_info_sub)

     IF( stat_info_sub /=0 ) THEN
        PRINT *, "mcf : ", &
             "Computing mass (after adjustment) failed !"
        stat_info = -1
        GOTO 9999
     END IF



  END IF

9999 CONTINUE


  IF(ASSOCIATED(coll_num_numerical_part) ) THEN
     DEALLOCATE(coll_num_numerical_part)
  END IF

  IF(ASSOCIATED(dx)) THEN
     DEALLOCATE(dx)
  END IF

  IF(ASSOCIATEd(tx)) THEN
     DEALLOCATE(tx)
  END IF

  IF(ASSOCIATEd(tv)) THEN
     DEALLOCATE(tv)
  END IF

  IF(ASSOCIATEd(tid)) THEN
     DEALLOCATE(tid)
  END IF

  IF ( ASSOCIATED(ide) ) THEN
     DEALLOCATE(ide)
  END IF


  RETURN

END SUBROUTINE particles_init_partial_inter_adjust


