SUBROUTINE colloid_adjust_arbitrary_colloid(this, num_layers, stat_info)
  !----------------------------------------------------
  ! Subroutine  : colloid_adjust_arbitrary_colloid
  !----------------------------------------------------
  !
  ! Purpose     : Adjust arbitrary colloid geometery
  !               for inflow/outflow BC
  !
  ! Reference   :
  !
  ! Remark      :
  !
  ! Revisions   : 
  !
  !----------------------------------------------------
  ! Author      : 
  ! Contact     : 
  !
  ! Dr. Marco Ellero's Emmy Noether Group,
  ! Prof. Dr. N. Adams' Chair of Aerodynamics,
  ! Faculty of Mechanical Engineering,
  ! Technische Universitaet Muenchen, Germany.
  !----------------------------------------------------

  TYPE(colloid), INTENT(INOUT)    :: this
  INTEGER,  INTENT(IN)            :: num_layers
  INTEGER , INTENT(OUT)           :: stat_info    

  !----------------------------------------------------
  ! Boundary parameters (and related):
  !----------------------------------------------------
  INTEGER                                 :: num_inout
  INTEGER , DIMENSION(:)  , POINTER       :: num_iopoints
  REAL(MK), DIMENSION(:,:), POINTER       :: iopatch_n   
  REAL(MK), DIMENSION(:,:,:), POINTER     :: iopatch_x 
  INTEGER , DIMENSION(:)  , POINTER       :: patch_id


  !----------------------------------------------------
  ! Local variables start here :
  !---------------------------------------------------- 

  INTEGER                                 :: stat_info_sub 
  REAL(MK), DIMENSION(:,:), POINTER       :: temp_coll_v 
  INTEGER , DIMENSION(:,:), POINTER       :: temp_coll_v_flist   
  INTEGER , DIMENSION(:,:), POINTER       :: temp_coll_f_vlist  
  REAL(MK), DIMENSION(:,:), POINTER       :: temp_coll_n 
  REAL(MK), DIMENSION(:,:), POINTER       :: px, vv 
  INTEGER , DIMENSION(:)  , POINTER       :: ind 
  REAL(MK), DIMENSION(:,:), POINTER       :: dum
  REAL(MK), DIMENSION(3)                  :: x, d1, d2, dn
  REAL(MK)                                :: len_layer, dx
  INTEGER                                 :: num_dim, num_vperf, max_degree
  INTEGER                                 :: facet_id, index, old_size1, old_size2, old_size3
  INTEGER                                 :: i, j, k, ip, ti


  !----------------------------------------------------
  ! Initialization of variables.
  !----------------------------------------------------
  
  stat_info = 0       

  NULLIFY(num_iopoints)
  NULLIFY(iopatch_n)  
  NULLIFY(iopatch_x)  
  NULLIFY(patch_id)

  NULLIFY(temp_coll_v)
  NULLIFY(temp_coll_v_flist)  
  NULLIFY(temp_coll_f_vlist)
  NULLIFY(temp_coll_n)
  NULLIFY(px)  
  NULLIFY(vv)
  NULLIFY(ind)
  NULLIFY(dum)

  num_dim = this%num_dim
  num_vperf = num_dim


  !-------------------------------------------------
  ! Get inflow/outflow patch data
  !-------------------------------------------------       
  num_inout = boundary_get_num_inout(this%boundary, stat_info_sub)  

  CALL boundary_get_num_iopoints(this%boundary, num_iopoints, stat_info_sub)
  CALL boundary_get_iopatch_x(this%boundary, iopatch_x, stat_info_sub) 
  CALL boundary_get_iopatch_n(this%boundary, iopatch_n, stat_info_sub) 
  CALL boundary_get_patch_id(this%boundary, patch_id, stat_info_sub)  
  CALL boundary_get_patch_id(this%boundary, patch_id, stat_info_sub) 

  len_layer = boundary_get_bwidth(this%boundary, stat_info_sub)

  dx = len_layer / DBLE(num_layers)


  DO ip = 1, num_inout


     ! allocate patch
     IF ( ASSOCIATED(px) ) THEN
        DEALLOCATE(px)
     END IF

     IF ( ASSOCIATED(vv) ) THEN
        DEALLOCATE(vv)
     END IF

     IF ( ASSOCIATED(ind) ) THEN
        DEALLOCATE(ind)
     END IF  


     ALLOCATE(px(num_dim, num_iopoints(ip)))   
     ALLOCATE(vv(num_dim, num_iopoints(ip))) 
     ALLOCATE(ind(num_iopoints(ip)))  

     px = iopatch_x(ip, 1:num_dim, 1:num_iopoints(ip)) 
     vv(1:num_dim, 1:num_iopoints(ip)) = 0.0_MK

     DO k = 1, num_layers  


        IF ( ASSOCIATED(dum) ) THEN
           DEALLOCATE(dum)
        END IF
        ALLOCATE(dum(num_dim, this%coll_v_num(1)))

        ! match indices of patch and geometry (identify inlet patch)
        ind(1:num_iopoints(ip)) = 0
        DO i = 1, num_iopoints(ip)
           x = px(1:num_dim, i)
           dum(1:num_dim, 1:this%coll_v_num(1)) = 0.0_MK
           DO j = 1, this%coll_v_num(1)
              dum(1:num_dim, j) = this%coll_v(1, 1:num_dim, j) - x  
           END DO
           ind(i) = MINLOC(SUM(ABS(dum), 1), 1)
        END DO

        ! extend each patch

        ! copy vertices
        DO i = 1, num_iopoints(ip)
           vv(1:num_dim, i) = this%coll_v(1, 1:num_dim, ind(i)) + iopatch_n(ip, 1:num_dim) * dx
        END DO


        facet_id = this%coll_f_num(1)
        ! reallocate (resize)
        old_size1 = SIZE(this%coll_v, 3)
        old_size2 = SIZE(this%coll_v_flist, 3)
        old_size3 = SIZE(this%coll_f_vlist, 3)
        max_degree = SIZE(this%coll_v_flist, 2)


        IF ( ASSOCIATED(temp_coll_v) ) THEN
           DEALLOCATE(temp_coll_v)
        END IF
        IF ( ASSOCIATED(temp_coll_v_flist) ) THEN
           DEALLOCATE(temp_coll_v_flist)
        END IF 
        IF ( ASSOCIATED(temp_coll_f_vlist) ) THEN
           DEALLOCATE(temp_coll_f_vlist)
        END IF
        IF ( ASSOCIATED(temp_coll_n) ) THEN
           DEALLOCATE(temp_coll_n)
        END IF

        ALLOCATE(temp_coll_v(1:num_dim, old_size1))
        ALLOCATE(temp_coll_v_flist(max_degree, old_size2))   
        ALLOCATE(temp_coll_f_vlist(num_vperf, old_size3))   
        ALLOCATE(temp_coll_n(num_dim, old_size3))

        temp_coll_v(1:num_dim, 1:old_size1) = &
             this%coll_v(1, 1:num_dim, 1:old_size1)   
        temp_coll_v_flist(1:max_degree, 1:old_size2) = &
             this%coll_v_flist(1, 1:max_degree, 1:old_size2)
        temp_coll_f_vlist(1:num_vperf, 1:old_size3) = & 
             this%coll_f_vlist(1, 1:num_vperf, 1:old_size3) 
        temp_coll_n(1:num_dim, 1:old_size3) = & 
             this%coll_n(1, 1:num_dim, 1:old_size3)



        IF ( ASSOCIATED(this%coll_v) ) THEN
           DEALLOCATE(this%coll_v)
        END IF
        IF ( ASSOCIATED(this%coll_v_flist) ) THEN
           DEALLOCATE(this%coll_v_flist)
        END IF  
        IF ( ASSOCIATED(this%coll_f_vlist) ) THEN
           DEALLOCATE(this%coll_f_vlist)
        END IF
        IF ( ASSOCIATED(this%coll_n) ) THEN
           DEALLOCATE(this%coll_n)
        END IF

        ALLOCATE(this%coll_v(1, num_dim, old_size1 + num_iopoints(ip)))
        ALLOCATE(this%coll_v_flist(1, max_degree, old_size2 + 2 * num_iopoints(ip)))    
        ALLOCATE(this%coll_f_vlist(1, num_vperf, old_size3 + 2 * num_iopoints(ip)))   
        ALLOCATE(this%coll_n(1, num_dim, old_size3 + 2 * num_iopoints(ip)))

        this%coll_v(1, 1:num_dim, 1:old_size1) = &
             temp_coll_v(1:num_dim, 1:old_size1)
        this%coll_v(1, 1:num_dim, old_size1+1:old_size1 + num_iopoints(ip)) = &
             vv(1:num_dim, 1:num_iopoints(ip))

        this%coll_v_flist(1, 1:max_degree, 1:old_size2) = &
             temp_coll_v_flist(1:max_degree, 1:old_size2)
        this%coll_v_flist(1, 1:max_degree, old_size2+1:old_size2 + 2 * num_iopoints(ip)) = 0  

        this%coll_f_vlist(1, 1:num_vperf, 1:old_size3) = &
             temp_coll_f_vlist(1:num_vperf, 1:old_size3)
        this%coll_f_vlist(1, 1:num_vperf, old_size3+1:old_size3 + 2 * num_iopoints(ip)) = 0   

        this%coll_n(1, 1:num_dim, 1:old_size3) = &
             temp_coll_n(1:num_dim, 1:old_size3)
        this%coll_n(1, 1:num_dim, old_size3+1:old_size3 + 2 * num_iopoints(ip)) = 0


        ! update data structure
        DO i = 1, num_iopoints(ip)

           ti = MOD(i-1, num_iopoints(ip)) + 1
           facet_id = facet_id + 1

           ! add facet1
           this%coll_f_vlist(1, 1:num_vperf, facet_id) = &
                (/ ind(i), ind(MOD(i, num_iopoints(ip))+1), this%coll_v_num(1)+i /)
           ! facet normal
           d1 = this%coll_v(1, 1:num_dim, this%coll_f_vlist(1, 2, facet_id)) - &
                this%coll_v(1, 1:num_dim, this%coll_f_vlist(1, 1, facet_id))
           d2 = this%coll_v(1, 1:num_dim, this%coll_f_vlist(1, 3, facet_id)) - &
                this%coll_v(1, 1:num_dim, this%coll_f_vlist(1, 1, facet_id))
           CALL tool_cross_product(this%tool, d1, d2, dn, stat_info_sub)

           this%coll_n(1, 1:num_dim, facet_id) = SIGN(1, patch_id(ip)) * dn / SQRT(DOT_PRODUCT(dn, dn))

           ! update facet list
           max_degree = SIZE(this%coll_v_flist, 2) 
           DO j = 1, num_vperf
              index = this%coll_f_vlist(1, j, facet_id)
              this%coll_v_flist(1, 2:max_degree, index) = this%coll_v_flist(1, 1:max_degree-1, index)
              this%coll_v_flist(1, 1, index) = facet_id
           END DO


           facet_id = facet_id + 1
           ! add facet2
           this%coll_f_vlist(1, 1:num_vperf, facet_id) = &
                (/ ind(MOD(i, num_iopoints(ip))+1), &
                this%coll_v_num(1)+MOD(i, num_iopoints(ip))+1, &
                this%coll_v_num(1)+i /)

           ! facet normal
           d1 = this%coll_v(1, 1:num_dim, this%coll_f_vlist(1, 2, facet_id)) - &
                this%coll_v(1, 1:num_dim, this%coll_f_vlist(1, 1, facet_id))
           d2 = this%coll_v(1, 1:num_dim, this%coll_f_vlist(1, 3, facet_id)) - &
                this%coll_v(1, 1:num_dim, this%coll_f_vlist(1, 1, facet_id))
           CALL tool_cross_product(this%tool, d1, d2, dn, stat_info_sub)

           this%coll_n(1, 1:num_dim, facet_id) = SIGN(1, patch_id(ip)) * dn / SQRT(DOT_PRODUCT(dn, dn))

           ! update facet list
           max_degree = SIZE(this%coll_v_flist, 2) 
           DO j = 1, num_vperf
              index = this%coll_f_vlist(1, j, facet_id)
              this%coll_v_flist(1, 2:max_degree, index) = this%coll_v_flist(1, 1:max_degree-1, index)
              this%coll_v_flist(1, 1, index) = facet_id
           END DO

        END DO


        this%coll_v_num(1) = this%coll_v_num(1) + num_iopoints(ip)
        this%coll_f_num(1) = facet_id
        
        px(1:num_dim, 1:num_iopoints(ip)) = vv(1:num_dim, 1:num_iopoints(ip))

     END DO

  END DO



  IF(ASSOCIATED(num_iopoints)) THEN
     DEALLOCATE(num_iopoints)
  END IF

  IF(ASSOCIATED(iopatch_x)) THEN
     DEALLOCATE(iopatch_x)
  END IF

  IF(ASSOCIATED(iopatch_n)) THEN
     DEALLOCATE(iopatch_n)
  END IF  

  IF(ASSOCIATED(patch_id)) THEN
     DEALLOCATE(patch_id)
  END IF

  IF(ASSOCIATED(temp_coll_v)) THEN
     DEALLOCATE(temp_coll_v)
  END IF

  IF(ASSOCIATED(temp_coll_v_flist)) THEN
     DEALLOCATE(temp_coll_v_flist)
  END IF

  IF(ASSOCIATED(temp_coll_f_vlist)) THEN
     DEALLOCATE(temp_coll_f_vlist)
  END IF

  IF(ASSOCIATED(temp_coll_n)) THEN
     DEALLOCATE(temp_coll_n)
  END IF

  IF(ASSOCIATED(px)) THEN
     DEALLOCATE(px)
  END IF
  
  IF(ASSOCIATED(vv)) THEN
     DEALLOCATE(vv)
  END IF 

  IF(ASSOCIATED(ind)) THEN
     DEALLOCATE(ind)
  END IF
  
  IF(ASSOCIATED(dum)) THEN
     DEALLOCATE(dum)
  END IF


  RETURN


END SUBROUTINE colloid_adjust_arbitrary_colloid
