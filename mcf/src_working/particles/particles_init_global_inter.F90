      SUBROUTINE particles_init_global_inter(this, d_rank,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_init_global_inter
        !----------------------------------------------------
        !
        ! Purpose     : Create particles on certain lattice 
        !               globally and internally,
        !               x, y,(z), vx,vy,(vz), p_id,s_id.
        !
        ! 
        ! Routines    : particles_init_global_inter_square()
        !               particles_init_global_inter_hexagonal()
        !               particles_init_global_inter_cubic()
        !
        !
        ! Remarks     : Currently velocity are set zero initially.
        !
        !               Currently this routine is only used 
        !               by root process to creat essenstial 
        !               quntities of particles globally on
        !               root process and then distribute 
        !               to all processes.
        !
        !               We generate particles using different
        !               lattice.
        !               However, for 2D colloid, a little bit 
        !               more attention has to be paid:
        !
        !               colloid boundary particles can be 
        !               on lattice, also can be on 
        !               the layers parallel to the surface,
        !               which further has two situations:
        !               a) each layer's particles have
        !                  equal distance, 
        !               b) or each layer has a fixed 
        !                  number of particles
        !                  (Zhu et al. 1999)
        !
        ! References  :
        !
        ! Revisions   : V0.1  22.07 2009, original version.
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
        ! Modules :
        !----------------------------------------------------
        
        USE ppm_module_find_duplicates
        
        !----------------------------------------------------
        ! Arguments :
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(IN)                     :: d_rank
        INTEGER, INTENT(OUT)	                :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here :
        !
        ! num_dim          : number of dimension.
        ! lattice_type     : lattice type.
        ! dx               : initial distance between two particles;        
        ! d_boundary       : pointer to a boundary obeject.
        ! num_wall_solid   : number of solid wall boundaries.
        ! num_inflow       : number of inflow boundaries.
        ! num_outflow      : number of outflow boundaries.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: num_dim,i,j
        INTEGER                                 :: lattice_type
        REAL(MK), DIMENSION(:), POINTER         :: dx
        TYPE(Boundary), POINTER                 :: d_boundary
        INTEGER                                 :: num_wall_solid
        INTEGER                                 :: num_inout
        
        !----------------------------------------------------
        ! Colloids parameters.
        !----------------------------------------------------
        
        INTEGER                                 :: num_colloid
        TYPE(Colloid), POINTER                  :: colloids
        INTEGER, DIMENSION(:), POINTER          :: coll_num_numerical_part
        INTEGER                                 :: coll_place
        REAL(MK), DIMENSION(:,:), POINTER       :: c_x
        REAL(MK), DIMENSION(:,:), POINTER       :: c_v
        INTEGER, DIMENSION(:), POINTER          :: c_sid

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
        INTEGER                                 :: num_extra
        LOGICAL                                 :: counted_wall_solid  
        LOGICAL                                 :: counted_colloid
        LOGICAL                                 :: counted_ignore
        
        LOGICAL                                 :: l_w  
        LOGICAL                                 :: l_sur
        LOGICAL                                 :: l_out
        LOGICAL                                 :: l_in
        REAL(MK), DIMENSION(3)                  :: p_v
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
        ! Hack variables: hck
        !----------------------------------------------------
        INTEGER                                 :: ie, ic, closest_i
        REAL(MK)                                :: closest_d, r_norm
        REAL(MK) , DIMENSION(:)                 :: r(3)

        
        
#ifdef __DEBUG
        REAL(MK), DIMENSION(:,:), POINTER       :: ppp
#endif 
        
        IF ( d_rank /= 0 ) THEN
           PRINT *,&
                "particles_init_global_inter : ", &
                "can only be called by root process !"
           stat_info = -1
           GOTO 9999
        END IF
        
        
        !----------------------------------------------------
    	! Initialization of variables.
    	!----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        num_dim = this%num_dim
        
        NULLIFY(dx)
        NULLIFY(d_boundary)
        
        !----------------------------------------------------
        ! Initialization of variables of colloids.
        !----------------------------------------------------
        
        NULLIFY(colloids)
        NULLIFY(coll_num_numerical_part)
       
        
        NULLIFY(c_x)
        NULLIFY(c_v)
        NULLIFY(c_sid) 

        
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
        lattice_type = &
             physics_get_lattice_type(this%phys,stat_info_sub)
        CALL physics_get_dx(this%phys,dx,stat_info_sub)
        CALL physics_get_boundary(this%phys,d_boundary,stat_info_sub)
        num_wall_solid = &
             boundary_get_num_wall_solid(d_boundary,stat_info_sub)   
        num_inout = &
             boundary_get_num_inout(d_boundary,stat_info_sub)   
        
        !----------------------------------------------------
        ! Check number of dimension and lattice type
        ! then choose particles generator accordingly.
        !----------------------------------------------------
        
        IF ( num_dim == 2 ) THEN
           
           !----------------------------------------------
           ! 2D lattice :
           !
           ! 1 square; 2 staggered; 3 hexagonal.
           !----------------------------------------------
           
           SELECT CASE( lattice_type )
              
           CASE (1)
              
              CALL particles_init_global_inter_square(this,&
                   stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "particles_init_global_inter : ", &
                      "Calling *_square failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
           
           CASE (2)
              
              CALL particles_init_global_inter_staggered(this,&
                   stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "particles_init_global_inter : ", &
                      "Calling *_staggerd failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           CASE (3)
              
              CALL particles_init_global_inter_hexagonal(this,&
                   stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "particles_init_global_inter : ", &
                      "Calling *_hexagonal failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           CASE DEFAULT
              
              PRINT * , "particles_init_global_inter : ",&
                   "lattcie type ", lattice_type, " not available !"
              stat_info = -1
              GOTO 9999
              
           END SELECT ! lattice_type
           
           !-------------------------------------------------
           ! 3D lattice, 
           ! 1 Simple Cubic; 2 Body centered;
           ! 3 Face centered lattice.
           !-------------------------------------------------
           
        ELSE IF( num_dim == 3 ) THEN
           
           SELECT CASE ( lattice_type )
              
           CASE (1)
              
              CALL particles_init_global_inter_cubic(this,&
                   stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "particles_init_global_inter : ", &
                      "Calling _cubic failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           CASE (2)
              
              PRINT * , "particles_init_global_inter : ",&
                   "3D body center lattice not available yet !"
              stat_info = -1
              GOTO 9999
              !CALL particles_init_global_inter_body_center(this,&
              !stat_info_sub)
              
           CASE (3)
              
              PRINT * , "particles_init_global_inter : ",&
                   "3D face center lattice not available yet !"
              stat_info = -1
              GOTO 9999                
              !CALL particles_init_global_inter_face_center(this,&
              !stat_info_sub)
              
           CASE DEFAULT
              
              PRINT * , "particles_init_global_inter : ",&
                   "lattcie type ", lattice_type, " not available yet !"
              stat_info = -1
              GOTO 9999
              
           END SELECT ! lattice_type
           
        END IF ! num_dim
        
#ifdef __DEBUG

#ifdef __DEBUG_INIT
        !----------------------------------------------------
        ! For debug purpose, could be written into files.
        !----------------------------------------------------
        
        NULLIFY(ppp)
        num = this%num_part_real
        
        ALLOCATE(ppp(num_dim+1,num))
        
        ppp(1:num_dim,1:num) = this%x(1:num_dim,1:num)
        ppp(num_dim+1,1:num) = this%id(2,1:num)
        
        CALL debug_write_output(global_debug,d_rank,&
             "particles_init_global_inter ", &
             "ppp_real",0,ppp,1,this%num_part_real,stat_info_sub)
        
#endif
#endif
        
        !----------------------------------------------------
        ! Colloid parameters :
        !----------------------------------------------------
        
        num_colloid = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        
        IF ( num_colloid > 0 ) THEN

#ifdef __COLLOID_ON_LATTICE
        
           CALL particles_set_colloid_on_lattice(this,stat_info_sub)
           
           IF ( stat_info_sub /=0 ) THEN
              
              PRINT * , "particles_init_global_inter : ", &
                   "calling particles_set_colloid_on_lattice failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
#endif           
           CALL physics_get_colloid(this%phys, &
                colloids,stat_info_sub)
           coll_place= &
                colloid_get_place(colloids,stat_info_sub)
           
           CALL colloid_compute_image(colloids,stat_info_sub)  


        END IF
        
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
                "particles_init_global_inter : ", &
                "Allocating memory for variables has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Reset all counters to zero,
        ! prepare for counting number of fluid particles,
        ! number of wall solid boundary particle, and
        ! number of inflow/outflow boundary particles.
        !----------------------------------------------------
        
        this%num_part_fluid = 0
        this%num_part_wall_solid = 0 
        
        IF ( num_colloid > 0 ) THEN
           
           ALLOCATE(coll_num_numerical_part(1:num_colloid))
           coll_num_numerical_part(1:num_colloid) = 0
           
        END IF
        
        num = 0
        
        !----------------------------------------------------
        ! Loop over all generated particles to decide 
        ! species ID, i.e., which type of particle it is.
        !----------------------------------------------------
        
        
        DO j = 1, num_total

           counted_wall_solid = .FALSE.
           counted_colloid    = .FALSE.
           counted_ignore     = .FALSE.
           
           p_v(1:num_dim) = 0
           p_sid          = 0
           
           !-------------------------------------------------
           ! Check if x is a wall(solid)/inflow/outflow boundary particle.
           !-------------------------------------------------

           IF ( (num_wall_solid > 0)  .OR. &
                (num_inout > 0) ) THEN
              
              l_w = .FALSE.
              
              CALL boundary_check_wallsolid_inoutflow_particle(d_boundary, &
                   this%x(1:num_dim,j), l_w, p_sid, stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "paricles_init_global_inter : ", &
                      "Checking wall solid / inoutflow particle failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              IF ( l_w ) THEN
                 
                 counted_wall_solid = .TRUE.

              END IF
              
           END IF ! num_wall_solid or num_inflow/outflow


           coll_arbitrary_num = 0
           IF ( num_colloid > 0 ) THEN
              coll_arbitrary_num = colloid_get_arbitrary_num(colloids, &
                   stat_info_sub)
           END IF
           
           ! only if NOT arbitrary shape colloids 
           IF (coll_arbitrary_num == 0) THEN
              
              
              !-------------------------------------------------
              ! If sx is neither a solid wall boundary particle
              ! nor an inflow/outflow boundary particle
              ! check if it is inside any colloid geometry.
              !-------------------------------------------------
              
              IF ( .NOT. counted_wall_solid .AND. &
                   num_colloid > 0 ) THEN
                 
                 !----------------------------------------------
                 ! loop over all colloids 
                 !----------------------------------------------
                 
                 l_sur = .FALSE.
                 l_out = .FALSE.
                 l_in  = .FALSE.
                 
                 CALL colloid_check_boundary_particle(colloids, &
                      this%x(1:num_dim,j),l_sur,l_out,l_in, p_sid,&
                      stat_info_sub)
                 
                 IF ( l_sur .AND. &
                      coll_place == mcf_colloid_place_lattice ) THEN
                    
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
                    
                 ELSE IF( l_out .AND. &
                      coll_place /= mcf_colloid_place_lattice ) THEN
                    
                    !-------------------------------------------
                    ! If j is outside of a colloid, but close
                    ! enough to the surface, it should not be 
                    ! created as fluid particle when colloid
                    ! boundary particles are required to be 
                    ! placed on paralle surfaces, since sx would
                    ! be too close to the boundary particles
                    ! on the surface.
                    !-------------------------------------------
                    
                    counted_ignore = .TRUE.
                    
                 END IF ! l_sur
                 
              END IF ! .NOT. counted_wall_solid and .NOT. counted_inflow/outflow
              
           END IF ! NO arbitrary colloids
           
           !-------------------------------------------------
           ! Save the particle's information.
           !-------------------------------------------------
           
           IF ( .NOT. counted_ignore ) THEN
              
              num = num + 1
              tx(1:num_dim,num)     = this%x(1:num_dim,j)
              tv(1:num_dim,num)     = p_v(1:num_dim)
              tid(this%pid_idx,num) = num
              tid(this%sid_idx,num) = p_sid
              tid(this%bid_idx,num) = 0

              !----------------------------------------------
              ! Increase counters of different species.
              !----------------------------------------------
              
              IF ( counted_wall_solid )  THEN
                 
                 this%num_part_wall_solid = &
                      this%num_part_wall_solid + 1 
                 
              ELSE IF ( counted_colloid )  THEN
                 
                 coll_num_numerical_part(p_sid) = &
                      coll_num_numerical_part(p_sid) + 1
                 
              ELSE       
                 
                 this%num_part_fluid = &
                      this%num_part_fluid + 1  


              END IF ! counted
              
           END IF ! .NOT. counted_ignore
           
        END DO ! j = 1, num
        
        
        !----------------------------------------------------
        ! Record number of solid wall boundary particles.
        !----------------------------------------------------
        
        CALL boundary_set_num_part_wall_solid(d_boundary, &
             this%num_part_wall_solid,stat_info_sub) 


        !----------------------------------------------------
        ! Check if we need to create colloid boundary 
        ! particles on parallel surfaces of colloid.
        ! If needed, we call colloid_* to create them.
        !----------------------------------------------------
        
        num_extra = 0
        
        IF ( num_colloid > 0 .AND. &
             coll_place /=  mcf_colloid_place_lattice ) THEN
           
           CALL colloid_create_boundary_particle(colloids,&
                dx(1:num_dim),c_x,c_v,c_sid,stat_info_sub)
           
           num_extra = SIZE(c_sid,1)
           
        END IF
   

        !----------------------------------------------------
        ! Copy the data from temporary varaiable 
        ! to particles members.
	!----------------------------------------------------
        
        num_total = num + num_extra
        
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
        
        !----------------------------------------------------
        ! If there is extra particles, i.e.,
        ! created on layers parallel to colloid surface.
        !----------------------------------------------------
        
        IF ( num_extra > 0 ) THEN
           
           this%x(1:num_dim,num+1:num_total) = &
                c_x(1:num_dim,1:num_extra)
           this%v(1:num_dim,num+1:num_total) = &
                c_v(1:num_dim,1:num_extra)
           
           DO i = 1, num_extra
              
              this%id(this%pid_idx,num+i) = (num + i)
              this%id(this%sid_idx,num+i) = c_sid(i)
              
              coll_num_numerical_part(c_sid(i)) = &
                   coll_num_numerical_part(c_sid(i)) + 1 
              
           END DO ! i = 1, num_extra
           
        END IF ! num_extra > 0


        !----------------------------------------------------
        ! Count the number of numerical particles which 
        ! consititute colloids.
        !----------------------------------------------------
        
        this%num_part_colloid = 0
        
        IF ( num_colloid > 0 ) THEN
           
           DO j = 1, num_colloid
              this%num_part_colloid = &
                   this%num_part_colloid + coll_num_numerical_part(j)
           END DO
           
           CALL colloid_set_num_numerical_part(colloids,&
                coll_num_numerical_part(:),stat_info_sub)
        
        END IF
        
        !----------------------------------------------------
        ! Set number of real, all and ghost particles.
        !----------------------------------------------------
        
        this%num_part_real  = num_total
        this%num_part_all   = this%num_part_real
        this%num_part_ghost = 0
        
#ifdef __DEBUG

#ifdef __DEBUG_INIT
        !----------------------------------------------------
        ! For debug purpose, could be written into files.
        !----------------------------------------------------
        
        NULLIFY(ppp)
        num = this%num_part_real
        
        ALLOCATE(ppp(num_dim+1,num))
        
        ppp(1:num_dim,1:num) = this%x(1:num_dim,1:num)
        ppp(num_dim+1,1:num) = this%id(2,1:num)
        
        CALL debug_write_output(global_debug,d_rank,&
             "particles_init_global_inter ", &
             "ppp_real",1,ppp,1,this%num_part_real,stat_info_sub)

#endif
        
#endif
        
        !----------------------------------------------------
    	! Check if each particles is unique.
    	!----------------------------------------------------
        
        nid = 0
        
        CALL ppm_find_duplicates(this%x,num_dim,&
             this%num_part_real,nid,ide,stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, &
                "particles_init_global_inter : ", &
                "Error by checking duiplicates !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF ( nid /= 0) THEN
           PRINT *, &
                "particles_init_global_inter : ", &
                "Found collocating particles !" &
                ,this%x(:,ide(1,1)), this%x(:,ide(2,1))
           stat_info = -1
           GOTO 9999
        END IF
        
        
9999    CONTINUE
        
        IF(ASSOCIATED(dx)) THEN
           DEALLOCATE(dx)
        END IF
        
        IF(ASSOCIATED(coll_num_numerical_part) ) THEN
           DEALLOCATE(coll_num_numerical_part)
        END IF
        
        IF(ASSOCIATEd(c_x)) THEN
           DEALLOCATE(c_x)
        END IF
        
        IF(ASSOCIATEd(c_v)) THEN
           DEALLOCATE(c_v)
        END IF
        
        IF(ASSOCIATEd(c_sid)) THEN
           DEALLOCATE(c_sid)
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
        
#if __DEBUG
        IF (ASSOCIATED(ppp)) THEN
           DEALLOCATE(ppp)
        END IF
#endif
        
        RETURN
        
      END SUBROUTINE particles_init_global_inter
      
#include "particles_init_global_inter_square.F90"
#include "particles_init_global_inter_staggered.F90"
#include "particles_init_global_inter_cubic.F90"
#include "particles_init_global_inter_hexagonal.F90"

      
      
