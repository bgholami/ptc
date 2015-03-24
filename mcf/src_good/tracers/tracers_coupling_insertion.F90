SUBROUTINE tracers_coupling_insertion(this,d_particles,stat_info)
  !----------------------------------------------------
  ! Subroutine  : tracers_coupling_insertion
  !----------------------------------------------------
  !
  ! Purpose     : Coupling of tracers and particles
  !               information are exchanged at 
  !               the interface (insertion)
  !                  
  !
  ! Routines    :  
  !
  ! References  :  
  !
  !
  ! Remarks     :  
  !
  !        
  ! Revisions   : V0.1 12.03.2014, change in code structure
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

  USE ppm_module_typedef, ONlY : ppm_t_clist

  !----------------------------------------------------
  ! Arguments :
  !----------------------------------------------------

  TYPE(Particles), INTENT(INOUT)          :: this   
  TYPE(Particles), INTENT(INOUT)          :: d_particles  
  INTEGER, INTENT(OUT)	                  :: stat_info	


  !----------------------------------------------------
  ! Physics parameters :
  !----------------------------------------------------

  INTEGER                                 :: num_dim 

  !----------------------------------------------------
  ! MPI parameters.
  !----------------------------------------------------

  INTEGER                                 :: rank
  INTEGER                                 :: COMM
  INTEGER                                 :: MPI_PREC   

  
  !----------------------------------------------------
  ! Technique parameters :
  !
  ! cell list        : cell list
  ! min_sub, max_sub : min and max extent of subdomain
  ! num_sub          : number of sumdomains on this process
  ! (the rest we don't need here)
  !---------------------------------------------------- 
 
  TYPE(ppm_t_clist),DIMENSION(:),POINTER       :: cell_list
  REAL(MK), DIMENSION(:), POINTER              :: min_sub
  REAL(MK), DIMENSION(:), POINTER              :: max_sub 
  INTEGER                                      :: num_sub
  INTEGER, DIMENSION(:,:), POINTER             :: sub_bcdef
  INTEGER, DIMENSION(:,:), POINTER             :: inp
  INTEGER, DIMENSION(:,:), POINTER             :: jnp
  INTEGER                                      :: nnp 


  !----------------------------------------------------
  ! cell counters and indices
  !---------------------------------------------------- 

  INTEGER                                      :: n1
  INTEGER                                      :: n2  
  INTEGER                                      :: icstart, icend, ci
  INTEGER                                      :: jcstart, jcend, cj
  INTEGER                                      :: kcstart, kcend, ck
  INTEGER                                      :: ccell   
  INTEGER                                      :: istart, iend, ipart
  INTEGER                                      :: jstart, jend, jpart
  INTEGER                                      :: il, ih, jl, jh, kl, kh


  !----------------------------------------------------
  ! Colloid parameters :
  !----------------------------------------------------

  TYPE(Colloid),POINTER                   :: colloids
  INTEGER                                 :: num_colloid 
  INTEGER                                 :: coll_arbitrary_num 
  INTEGER , DIMENSION(:)    , POINTER     :: coll_v_num 
  REAL(MK), DIMENSION(:,:,:), POINTER     :: coll_v  
  INTEGER , DIMENSION(:,:,:), POINTER     :: coll_f_vlist


  !----------------------------------------------------
  ! Number of tracers
  !----------------------------------------------------

  INTEGER                                 :: num_part_real  


  !----------------------------------------------------
  ! local variables:
  !----------------------------------------------------   

  REAL(MK), DIMENSION(:)   , ALLOCATABLE  :: r
  REAL(MK), DIMENSION(:)   , ALLOCATABLE  :: temp_prob 
  REAL(MK), DIMENSION(:)   , ALLOCATABLE  :: cell_dx  
  INTEGER , DIMENSION(:)   , ALLOCATABLE  :: cell_index, ngl
  INTEGER , DIMENSION(:)   , ALLOCATABLE  :: temp_neighbors
  INTEGER , DIMENSION(:,:) , ALLOCATABLE  :: local_insertion
  INTEGER                                 :: to_be_deleted, pos, pid
  REAL(MK)                                :: rand  
  REAL(MK)                                :: W, r_norm, uni_rnd
  REAL(MK)                                :: cut_off, cut_off2 
  INTEGER                                 :: rand_int
  INTEGER                                 :: vperf, nmh
  INTEGER                                 :: i, ip, vi, ic, ie    
  INTEGER                                 :: stat_info_sub



  !----------------------------------------------------
  ! Initialization of variables.
  !----------------------------------------------------

  stat_info     = 0
  stat_info_sub = 0  


  NULLIFY(colloids)
  NULLIFY(coll_v_num) 
  NULLIFY(coll_v) 
  NULLIFY(coll_f_vlist)  
  NULLIFY(min_sub)
  NULLIFY(max_sub)  
  NULLIFY(cell_list)
  NULLIFY(sub_bcdef)        
  NULLIFY(inp)
  NULLIFY(jnp) 


  !----------------------------------------------------
  ! Physics parameters :
  !----------------------------------------------------

  num_dim     = &
       physics_get_num_dim(this%phys,stat_info_sub) 
  cut_off     = &
       physics_get_cut_off(this%phys,stat_info_sub)
  cut_off2 = cut_off * cut_off


  !----------------------------------------------------
  ! MPI parameters.
  !----------------------------------------------------

  rank     = technique_get_rank(this%tech,stat_info_sub)
  MPI_PREC = technique_get_MPI_PREC(this%tech,stat_info_sub)
  COMM     = technique_get_comm(this%tech,stat_info_sub)


 !----------------------------------------------------
  ! Get the boundary of this sub-domain. 
  !----------------------------------------------------

  CALL technique_get_min_sub(this%tech,min_sub,stat_info_sub)
  CALL technique_get_max_sub(this%tech,max_sub,stat_info_sub)


  !----------------------------------------------------
  ! Get particles' cell list and prepare data to rank
  !----------------------------------------------------        

  CALL technique_get_cell_list(d_particles%tech, &
       num_sub, cell_list, &
       inp, jnp, nnp, sub_bcdef, stat_info_sub)


  IF (num_sub .GT. 1) THEN
     print*, 'Tracers_interpolate_velocity : ', &
          'technique_get_min_sub returns one set of dimension only.', &
          'More than one subdomain on each process is not supported!'
     stat_info = -1
     GOTO 9999
  END IF

  ! allocate arrays
  allocate(cell_index(num_dim), ngl(num_dim))
  allocate(cell_dx(num_dim))

  DO i = 1, num_dim 
     nmh = INT((max_sub(i) - min_sub(i)) / d_particles%tech%ghost_size)

     ngl(i)     = cell_list(num_sub)%nm(i) - nmh - 1
     cell_dx(i) = (max_sub(i) - min_sub(i)) / REAL(nmh, MK)
  END DO


  ! set cell index limits
  ! symmetry, so cell indices starts from 0
  icstart = 0              
  jcstart = 0
  kcstart = 0

  icend = cell_list(num_sub)%nm(1)-2
  jcend = cell_list(num_sub)%nm(2)-2 

  n1 = cell_list(num_sub)%nm(1)

  IF ( num_dim == 2 ) THEN

     kcend = kcstart
     n2 = 0

  ELSE

     kcend = cell_list(num_sub)%nm(3)-2      
     n2 = cell_list(num_sub)%nm(1) * &
          cell_list(num_sub)%nm(2)

  END IF


  !----------------------------------------------------
  ! Get colloid
  !----------------------------------------------------
  num_colloid = &
       physics_get_num_colloid(this%phys,stat_info_sub)

  IF ( num_colloid > 0 ) THEN

     CALL physics_get_colloid(this%phys, &
          colloids,stat_info_sub)

  END IF


  !----------------------------------------------------
  ! initialize colloid variables
  !----------------------------------------------------  
  coll_arbitrary_num = colloid_get_arbitrary_num(colloids, stat_info_sub)
  CALL colloid_get_v_num(colloids, coll_v_num, stat_info_sub)
  CALL colloid_get_coll_v(colloids, coll_v, stat_info_sub)    
  CALL colloid_get_coll_f_vlist(colloids, coll_f_vlist, stat_info_sub)  


  num_part_real =  this%num_part_real
  vperf = num_dim 
  ALLOCATE(local_insertion(coll_arbitrary_num, MAXVAL(coll_v_num)))



  ! is it neccesary to allocate again here?

  IF (ASSOCIATED(this%h_p)) THEN
     DEALLOCATE(this%h_p)
  END IF

  IF (ASSOCIATED(this%n_vec)) THEN
     DEALLOCATE(this%n_vec)
  END IF

  IF (ASSOCIATED(this%facet)) THEN
     DEALLOCATE(this%facet)
  END IF

  IF (ASSOCIATED(this%col_i)) THEN
     DEALLOCATE(this%col_i)
  END IF

  ALLOCATE(this%h_p(num_part_real))
  ALLOCATE(this%n_vec(num_dim, num_part_real))
  ALLOCATE(this%facet(num_part_real))  
  ALLOCATE(this%col_i(num_part_real)) 


  CALL tracers_compute_wall_distance(this, &
       2, stat_info_sub)


  ! reset counters
  to_be_deleted = 0
  local_insertion = 0
  this%insertion_list = 0
  DO ip = 1, num_part_real


     IF (this%id(1, ip) == 0) THEN ! if a to-be-deleted tracer 

        ! add one to the number of to-be-deleted ones
        to_be_deleted = to_be_deleted + 1

        ! insert into a SPH particle (only if inside the domain)
        IF ( this%id(2, ip) == 0) THEN 


           ! The tracer has left the NW region:
           ! add it to the insertion array at the interface.

           ! randomly choose one of the vertices of tracers current facet
           CALL RANDOM_NUMBER(rand)
           rand_int = INT(rand*DBLE(vperf)) + 1
           vi = coll_f_vlist(this%col_i(ip), rand_int, this%facet(ip))

           ! assign one tracer
           local_insertion(this%col_i(ip), vi) = local_insertion(this%col_i(ip), vi) + 1

        END IF


     ELSE ! if a normal or deposited tracer 


        ! find which wall the tracers belongs to

        IF (this%id(1, ip) > 0) THEN ! in existing tracers if normal trace 

           ! be happy
           CONTINUE

        ELSE  ! Add to the number of deposited tracers (i.e. this%id(1, ip) < 0)

           this%num_tracer_deposited(this%col_i(ip), this%facet(ip)) = &
                this%num_tracer_deposited(this%col_i(ip), this%facet(ip)) + 1

           ! and set the ID to zero to delete after counting all of the tracers
           this%id(1, ip) = 0
           ! ...and add one to the number of to-be-deleted ones
           to_be_deleted = to_be_deleted + 1

        END IF

     END IF

  END DO



  ! Allreduce local_insertion
  CALL MPI_ALLREDUCE(local_insertion(1:coll_arbitrary_num, 1:MAXVAL(coll_v_num)), &
       this%insertion_list(1:coll_arbitrary_num, 1:MAXVAL(coll_v_num)), &
       coll_arbitrary_num * MAXVAL(coll_v_num), &
       MPI_INTEGER, &
       MPI_SUM, &
       COMM, &
       stat_info_sub)



  ! Now actual insertion 

  ALLOCATE(temp_prob(d_particles%num_part_real))  
  ALLOCATE(temp_neighbors(d_particles%num_part_real))  
  ALLOCATE(r(num_dim))


  DO ic = 1, coll_arbitrary_num
     DO vi = 1, coll_v_num(ic)

        ! if this is a lcoal interface vertex
        IF ( this%local_interface(ic, vi) .AND. &
             (this%insertion_list(ic, vi) > 0)) THEN


           ! find SPH neighbors
           pos = 0
           temp_neighbors = 0
           temp_prob = 0.0_MK


           ! find the cell it belongs to
           cell_index(1:num_dim) = FLOOR((this%x_interface(ic, 1:num_dim, vi) - min_sub(1:num_dim)) / &
                cell_dx(1:num_dim)) + ngl(1:num_dim)

           ! enforce right limits for cell index
           il = MAX(cell_index(1)-1, icstart)  
           ih = MIN(cell_index(1)+1, icend)   
           jl = MAX(cell_index(2)-1, jcstart)  
           jh = MIN(cell_index(2)+1, jcend)
           IF ( num_dim == 2 ) THEN
              kl = kcstart
              kh = kcend
           ELSE
              kl = MAX(cell_index(3)-1, kcstart)  
              kh = MIN(cell_index(3)+1, kcend)
           END IF

           ! loop over neighboring cells
           DO ck = kl, kh
              DO cj = jl, jh
                 DO ci = il, ih

                    ! cell number
                    ccell = ci + 1 + n1 * cj + n2 * ck 

                    ! Get pointers of first and last particles in ccell
                    istart = cell_list(num_sub)%lhbx(ccell)
                    iend   = cell_list(num_sub)%lhbx(ccell+1)-1

                    ! loop over particles in ccell
                    DO ipart = istart, iend

                       ! index of the particles (in data array)
                       ie = cell_list(num_sub)%lpdx(ipart)


                       r = this%x_interface(ic, 1:num_dim, vi) - d_particles%x(1:num_dim, ie)
                       r_norm = dot_product(r, r)

                       IF ( (d_particles%id(2, ie) == 0) .AND. &
                            (r_norm <= cut_off2) ) THEN 

                          r_norm = SQRT(r_norm)

                          CALL kernel_kernel(d_particles%kern, &
                               r_norm, W, stat_info_sub)

                          pos = pos + 1
                          temp_neighbors(pos) = ie
                          temp_prob(pos) = W

                       END IF

                    END DO

                 END DO

              END DO

           END DO


           ! assign probability
           temp_prob = temp_prob / SUM(temp_prob)


           DO i = 1, this%insertion_list(ic, vi)

              ! insert this tracer to one particle
              uni_rnd = random_random_uniform(this%random, stat_info_sub)

              pos = 0
              DO WHILE (uni_rnd >= 0.0_MK)
                 pos = pos + 1
                 uni_rnd = uni_rnd - temp_prob(pos)
              END DO

              pid = temp_neighbors(pos)
              d_particles%c_tracers(pid) = d_particles%c_tracers(pid) + 1.0_MK

           END DO

        END IF

     END DO

  END DO

  !----------------------------------------------------
  ! delete those tracers that have left the near-wall region 
  !----------------------------------------------------
  ! do we need to copy hp as well?
  CALL tracers_delete_tracers(this, to_be_deleted, stat_info_sub)    

  !----------------------------------------------------
  ! Return.
  !----------------------------------------------------

9999 CONTINUE  


  IF(ASSOCIATED(coll_v))THEN
     DEALLOCATE(coll_v)
  END IF

  IF(ASSOCIATED(coll_v_num))THEN
     DEALLOCATE(coll_v_num)
  END IF

  IF(ASSOCIATED(coll_f_vlist))THEN
     DEALLOCATE(coll_f_vlist)
  END IF  

  IF(ASSOCIATED(min_sub)) THEN
     DEALLOCATE(min_sub)
  END IF

  IF(ASSOCIATED(max_sub)) THEN
     DEALLOCATE(max_sub)
  END IF

  IF(ASSOCIATED(cell_list)) THEN
     DEALLOCATE(cell_list)
  END IF

  IF(ASSOCIATED(sub_bcdef)) THEN
     DEALLOCATE(sub_bcdef)
  END IF

  IF(ASSOCIATED(inp)) THEN
     DEALLOCATE(inp)
  END IF

  IF(ASSOCIATED(jnp)) THEN
     DEALLOCATE(jnp)
  END IF


  RETURN

END SUBROUTINE tracers_coupling_insertion
