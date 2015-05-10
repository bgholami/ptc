 SUBROUTINE tracers_delete_tracers(this, to_be_deleted, stat_info)
        !----------------------------------------------------
        ! Subroutine  : tracers_delete_tracers
        !----------------------------------------------------
        !
        ! Purpose     : deletes tracers that leave the 
        !               near-wall region.
        !                  
        !
        ! Routines    :  
        !
        ! References  :  
        !
        !
        ! Remarks     :  
        !
        ! Revisions   : V0.1 27.07 2011, original
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
        
        
        
	!----------------------------------------------------
      	! Arguments :
      	!----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)  :: this
        INTEGER                         :: to_be_deleted
        INTEGER, INTENT(OUT)		:: stat_info	
        
      	!----------------------------------------------------
      	! Local variables start here :
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: ip

        !----------------------------------------------------
        ! Physics parameters :
        !
        ! num_dim        : number of dimension.
        !----------------------------------------------------
        
        INTEGER                         :: num_dim

        !----------------------------------------------------
        ! Number of all tracers
     	!----------------------------------------------------
        
        INTEGER                         :: num_part_real

        !----------------------------------------------------
        ! local variables:
        !  
        !----------------------------------------------------

        REAL(MK), DIMENSION(:,:),POINTER:: x_t, v_t, id_t
        REAL(MK), DIMENSION(:)  ,POINTER:: hp_t
        INTEGER                         :: num_of_deletes, delete_id, num

  	!----------------------------------------------------
      	! Initialization of variables.
      	!----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0    
        
        NULLIFY(x_t)    
        NULLIFY(v_t)        
        NULLIFY(hp_t)    
        NULLIFY(id_t)


        
        !----------------------------------------------------
        ! Physics parameters :
        !
        ! from a object of Physics class.
        !
        !----------------------------------------------------
        
        num_dim     = &
             physics_get_num_dim(this%phys,stat_info_sub) 


        !----------------------------------------------------
        ! Number of real tracers
        !----------------------------------------------------
        num_part_real = &
             tracers_get_num_part_real(this,stat_info_sub) 

        !----------------------------------------------------
        ! Allocate memory for force,
        !
        ! for non-symmetry : allocate num_part_real
        !
        !----------------------------------------------------
        

        ! id: positive >> normal
        !     zero     >> delete
        !     negative >> deposited
        !num_of_deletes = SUM(this%del_id)  
        ! those to be deleted have 0 for id, so...
        !num_of_deletes = sum(0**abs(this%id(1,1:num_part_real))) ! this works too
        num_of_deletes = to_be_deleted


!!$        print*, num_part_real, num_of_deletes

        IF (num_of_deletes > 0) THEN
           
           ALLOCATE(x_t(num_dim,num_part_real-num_of_deletes), &
                STAT=stat_info_sub)         
           ALLOCATE(v_t(num_dim,num_part_real-num_of_deletes), &
                STAT=stat_info_sub)   
           ALLOCATE(hp_t(num_part_real-num_of_deletes), &
                STAT=stat_info_sub) 
           ALLOCATE(id_t(this%num_id,num_part_real-num_of_deletes), &
                STAT=stat_info_sub)
           
           
           num = 1
           DO ip= 1, num_part_real
              
              !IF (this%del_id(ip) /= 1) THEN
              IF (this%id(1,ip) /= 0) THEN
                 
                 x_t(1:num_dim,num) = this%x(1:num_dim,ip)
                 v_t(1:num_dim,num) = this%v(1:num_dim,ip)
                 hp_t(num)   = this%h_p(ip)
                 id_t(1:this%num_id,num) = this%id(1:this%num_id,ip)
                 num = num + 1
                 
              END IF
              
           END DO
           
           num_part_real = num_part_real - num_of_deletes
           this%num_part_real = num_part_real 
           this%num_part_all = num_part_real       
           
           IF (ASSOCIATED(this%x)) THEN
              DEALLOCATE(this%x,STAT=stat_info_sub)
           END IF
           IF (ASSOCIATED(this%v)) THEN
              DEALLOCATE(this%v,STAT=stat_info_sub)
           END IF
           IF (ASSOCIATED(this%h_p)) THEN
              DEALLOCATE(this%h_p,STAT=stat_info_sub)
           END IF
           IF (ASSOCIATED(this%id)) THEN
              DEALLOCATE(this%id,STAT=stat_info_sub)
           END IF
           
           ALLOCATE(this%x(num_dim,num_part_real), &
                STAT=stat_info_sub)         
           ALLOCATE(this%v(num_dim,num_part_real), &
                STAT=stat_info_sub)   
           ALLOCATE(this%h_p(num_part_real), &
                STAT=stat_info_sub) 
           ALLOCATE(this%id(this%num_id,num_part_real), &
                STAT=stat_info_sub)
           
           
           
           this%x = x_t        
           this%v = v_t   
           this%h_p = hp_t    
           this%id = id_t

        END IF
        
        
        !----------------------------------------------------
        ! Return.
        !----------------------------------------------------
        
9999    CONTINUE  
        
        
        
        IF (ASSOCIATED(x_t)) THEN
           DEALLOCATE(x_t,STAT=stat_info_sub)
        END IF
        
        IF (ASSOCIATED(v_t)) THEN
           DEALLOCATE(v_t,STAT=stat_info_sub)
        END IF
        
        IF (ASSOCIATED(hp_t)) THEN
           DEALLOCATE(hp_t,STAT=stat_info_sub)
        END IF
        
        IF (ASSOCIATED(id_t)) THEN
           DEALLOCATE(id_t,STAT=stat_info_sub)
        END IF
        
        
        
        
        RETURN
        
        
      END SUBROUTINE tracers_delete_tracers
      
