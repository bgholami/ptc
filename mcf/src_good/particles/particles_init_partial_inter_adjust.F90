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

                    IF (this%id(this%bid_idx, j) == 0) THEN ! only if a normal fluid particle 
                       ! (inflow and outflow fluid particles should be left alone!!)

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

                    END IF

                 END IF ! p_id

              END IF ! not wall

              !-------------------------------------------------
              ! Save the particle's information.
              !-------------------------------------------------

              IF ( .NOT. counted_ignore ) THEN

                 num = num + 1
                 tx(1:num_dim,num)     = this%x(1:num_dim,j)
                 tv(1:num_dim,num)     = 0.0_MK
                 tid(this%pid_idx,num) = this%id(this%pid_idx, j)
                 tid(this%sid_idx,num) = p_sid
                 tid(this%bid_idx,num) = this%id(this%bid_idx, j)

                 !----------------------------------------------
                 ! Increase counters of different species.
                 !----------------------------------------------
                 IF ( .NOT. counted_wall_solid )  THEN ! wall is counted before

                    IF ( counted_colloid )  THEN
                       
                       coll_num_numerical_part(p_sid) = &
                            coll_num_numerical_part(p_sid) + 1
                       
                    ELSE       
                       
                       this%num_part_fluid = &
                            this%num_part_fluid + 1
                       
                    END IF ! counted

                 END IF ! .NOT. wall

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

9999    CONTINUE


        IF(ASSOCIATED(coll_num_numerical_part) ) THEN
           DEALLOCATE(coll_num_numerical_part)
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
      
      
