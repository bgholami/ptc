      SUBROUTINE io_read_deposited(this,d_tracers,stat_info)
        !----------------------------------------------------
        ! Subroutine  : io_read_deposited
        !----------------------------------------------------
        !
        ! Purpose     : Reading deposited from file.
        !
        ! Revision    : V0.1 10.09 2013, original version.
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
        ! Arguments
        !----------------------------------------------------
        
        TYPE(IO), INTENT(IN)                    :: this
        TYPE(Particles), INTENT(INOUT)          :: d_tracers
        INTEGER, INTENT(OUT)                    :: stat_info 


        !----------------------------------------------------
        ! Colloid parameters :
        !
        ! coll_x     : center of the colloid
        ! coll_rad   : dimensions of the colloid
        !
        !----------------------------------------------------
        
        TYPE(Colloid),POINTER                   :: colloids
        INTEGER                                 :: num_colloid  
        INTEGER                                 :: coll_arbitrary_num 
        INTEGER , DIMENSION(:)    , POINTER     :: coll_f_num
        REAL(MK), DIMENSION(:,:,:), POINTER     :: coll_v 
        INTEGER , DIMENSION(:,:,:), POINTER     :: coll_f_vlist
        

        !----------------------------------------------------
        ! Local variables
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        REAL(MK), DIMENSION(:,:), POINTER       :: tx
        REAL(MK), DIMENSION(:)  , POINTER       :: td

        INTEGER                                 :: ilenread
        CHARACTER(MAX_CHAR)                     :: cbuf
        LOGICAL                                 :: lExist
        INTEGER                                 :: num_dim
        INTEGER                                 :: num_points, vperf
        INTEGER                                 :: i, k, ic, vi
        REAL(MK), ALLOCATABLE                   :: mid_point(:)
        REAL(MK)                                :: val
        INTEGER                                 :: iline
        
   
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        
        NULLIFY(tx)
        NULLIFY(td)

        NULLIFY(colloids)
        NULLIFY(coll_f_num) 
        NULLIFY(coll_v) 
        NULLIFY(coll_f_vlist) 
           
        
	!----------------------------------------------------
        ! num_dim : number of dimension
	!----------------------------------------------------

        num_dim = tracers_get_num_dim(d_tracers,stat_info_sub) 


        !----------------------------------------------------
        ! Get colloid
        !----------------------------------------------------
        
        num_colloid = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        
        coll_arbitrary_num = 0 

        IF ( num_colloid > 0 ) THEN

           CALL physics_get_colloid(this%phys, &
                colloids,stat_info_sub) 

           coll_arbitrary_num = colloid_get_arbitrary_num(colloids, stat_info_sub)
           
        END IF
        

        IF (coll_arbitrary_num > 0) THEN
           !----------------------------------------------------
           ! initialize colloid variables
           !----------------------------------------------------  
           
           CALL colloid_get_f_num(colloids, coll_f_num, stat_info_sub)
           CALL colloid_get_coll_v(colloids, coll_v, stat_info_sub)  
           CALL colloid_get_coll_f_vlist(colloids, coll_f_vlist, stat_info_sub)

        END IF
           
	!----------------------------------------------------
      	! Check if name of particle file is empty.
      	!----------------------------------------------------
        
        ilenread = LEN_TRIM(this%read_deposited_file)
        IF ( ilenread < 1 ) THEN
           PRINT *,'io_read_deposited : ',&
                'No file name is given !'
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
      	! Check if the file exists. 
      	!----------------------------------------------------
        
        INQUIRE(FILE=this%read_deposited_file,EXIST=lExist)
	
        IF (.NOT.lExist) THEN
           WRITE(cbuf,'(2A)')'No such file: ', &
                this%read_deposited_file(1:ilenread)
           PRINT *, 'io_read_deposited : "', cbuf
           stat_info = -1
           GOTO 9999
        END IF
        
      	!----------------------------------------------------
      	! Open the file.
      	!----------------------------------------------------
        
        OPEN(this%read_deposited_unit,&
             FILE=this%read_deposited_file,&
             IOSTAT=stat_info_sub,ACTION='READ')
        
        IF (stat_info_sub /= 0) THEN
           WRITE(cbuf,'(2A)')'Failed to open file: ',&
                this%read_deposited_file(1:ilenread)
           PRINT *,'io_read_deposited : ', cbuf
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! assign number of facets
        !----------------------------------------------------
        
        num_points = SUM(coll_f_num)

        !----------------------------------------------------
        ! Allocate memory according to tracer number.
        !----------------------------------------------------
        
        IF( num_points > 0) THEN
           
           ALLOCATE(tx(num_dim,num_points))
           ALLOCATE(td(num_points))
           
        END IF
        
        !----------------------------------------------------
        ! Scan the file line by line and read tracers.
        !----------------------------------------------------
        
        iline = 0
        
        DO
           !-------------------------------------------------
           ! Increase line counter.
           !-------------------------------------------------
           
           iline = iline + 1
           
           
           IF( iline > num_points ) THEN
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Read information of current line.
           !-------------------------------------------------
           
           READ(this%read_deposited_unit,*,END=9999,ERR=200) &
                tx(1:num_dim,iline), td(iline)
                

        END DO
        
200     CONTINUE
        
        PRINT *, "io_read_deposited : ", &
             "reading deposited has problem at line ",&
             iline
        stat_info = -1
        GOTO 9999
        
9999    CONTINUE
        
        !----------------------------------------------------
        ! If actual read lines is not equal to the number
        ! of tracers given in the beginning of file,
        ! then there is something some.
        ! Otherwise, start matching.
        !----------------------------------------------------

        IF ( iline /= num_points+1 ) THEN
           
           PRINT *, "io_read_deposited : ", &
                "actual number of lines is not equal to number given !"
           stat_info = -1
           
        ELSE  

           ! calculate mid-points and match  
           vperf = num_dim
           allocate(mid_point(num_dim))   
           DO ic = 1, coll_arbitrary_num
              
              DO i = 1, coll_f_num(ic)
                 
                 ! find the mid point
                 mid_point = 0.0_MK
                 DO k = 1, vperf
                    vi = coll_f_vlist(ic, k, i)
                    mid_point = mid_point + coll_v(ic, 1:num_dim, vi)
                 END DO
                 mid_point = mid_point / DBLE(vperf)
                 
                 
                 ! match and assign 
                 val = 1.0_MK
                 DO k = 1, num_points
                    IF (sum(abs(tx(:, k) - mid_point)) < val) THEN
                       val = sum(abs(tx(:, k) - mid_point))
                       d_tracers%num_tracer_deposited(ic, i) = INT(td(k))
                    END IF
                 END DO
                 
                 
              END DO
              
           END DO
           
        END IF
        
        !----------------------------------------------------
        ! Release dynamic memory.
        !----------------------------------------------------
        
        IF(ASSOCIATED(tx)) THEN
           DEALLOCATE(tx)
        END IF
        
        IF(ASSOCIATED(td)) THEN
           DEALLOCATE(td)
        END IF
                       
        
        RETURN
        
      END SUBROUTINE io_read_deposited
      
