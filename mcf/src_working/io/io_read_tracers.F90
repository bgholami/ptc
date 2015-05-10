      SUBROUTINE io_read_tracers(this,d_rank,d_tracers,stat_info)
        !----------------------------------------------------
        ! Subroutine  : io_read_tracers
        !----------------------------------------------------
        !
        ! Purpose     : Reading tracers' configuration,
        !               i.e. position, velocity, and ID 
        !               from file.
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
        
        TYPE(IO), INTENT(INOUT)                 :: this
        INTEGER, INTENT(IN)                     :: d_rank
        TYPE(Particles), INTENT(INOUT)          :: d_tracers
        INTEGER, INTENT(OUT)                    :: stat_info
        

        !----------------------------------------------------
        ! Local variables
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        REAL(MK), DIMENSION(:,:), POINTER       :: x
        REAL(MK), DIMENSION(:,:), POINTER       :: v
        INTEGER, DIMENSION(:)   , POINTER       :: id 
        REAL(MK)                                :: inter_id
        REAL(MK)                                :: dum

        INTEGER                                 :: ilenread
        CHARACTER(MAX_CHAR)                     :: cbuf
        LOGICAL                                 :: lExist
        TYPE(Physics),POINTER                   :: phys
        INTEGER                                 :: num_dim
        INTEGER                                 :: num_id
        INTEGER                                 :: num_part
        INTEGER                                 :: iline, count
        
   
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        
        NULLIFY(x)
        NULLIFY(v)
        NULLIFY(id)
        NULLIFY(phys) 
        
        IF( d_rank /= 0 ) THEN
           PRINT *, "io_read_tracers : " , &
                "can only be called by root process ! "
           stat_info  = -1
           GOTO 9999
        END IF
        
        
	!----------------------------------------------------
        ! num_dim : number of dimension
        ! num_id  : number of different type of IDs for tracers.
	!----------------------------------------------------

        num_dim = tracers_get_num_dim(d_tracers,stat_info_sub)
        num_id  = tracers_get_num_id(d_tracers,stat_info_sub)
        
	!----------------------------------------------------
      	! Check if name of particle file is empty.
      	!----------------------------------------------------
        
        ilenread = LEN_TRIM(this%read_tracers_file)
        
        IF ( ilenread < 1 ) THEN
           PRINT *,'io_read_tracers : ',&
                'No file name is given !'
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
      	! Check if the particle file exists. 
      	!----------------------------------------------------
        
        INQUIRE(FILE=this%read_tracers_file,EXIST=lExist)
	
        IF (.NOT.lExist) THEN
           WRITE(cbuf,'(2A)')'No such file: ', &
                this%read_tracers_file(1:ilenread)
           PRINT *, 'io_read_tracers : "', cbuf
           stat_info = -1
           GOTO 9999
        END IF
        
      	!----------------------------------------------------
      	! Open the file.
      	!----------------------------------------------------
        
        OPEN(this%read_tracers_unit,&
             FILE=this%read_tracers_file,&
             IOSTAT=stat_info_sub,ACTION='READ')
        
        IF (stat_info_sub /= 0) THEN
           WRITE(cbuf,'(2A)')'Failed to open file: ',&
                this%read_tracers_file(1:ilenread)
           PRINT *,'io_read_tracers : ', cbuf
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! assign number of tracers from physics config file
        !----------------------------------------------------
        
        num_part = physics_get_tracers_restart_num(this%phys, stat_info_sub)
        
        !----------------------------------------------------
        ! Allocate memory according to tracer number.
        !----------------------------------------------------
        
        IF( num_part > 0) THEN
           
           ALLOCATE(x(num_dim,num_part))
           ALLOCATE(v(num_dim,num_part))
           ALLOCATE(id(num_part)) 
           
        END IF
        
        !----------------------------------------------------
        ! Scan the file line by line and read tracers.
        !----------------------------------------------------
        
        iline = 0

        DO count = 1, num_part

           !-------------------------------------------------
           ! Increase line counter.
           !-------------------------------------------------
           
           iline = iline + 1
           

           !-------------------------------------------------
           ! Read information of current line.
           !-------------------------------------------------
           
           READ(this%read_tracers_unit,*) &
                x(1:num_dim,iline), v(1:num_dim,iline),&
                dum, inter_id
                
           
           !-------------------------------------------------
           ! Convert floating point number of ID to integer.
           !-------------------------------------------------
           
           id(iline) = INT(inter_id)
           ! check if it is a dummy output (that should be disregarded)
           IF (id(iline) <= 0) THEN
           !IF ((id(iline) <= 0) .OR. (dum < 0.0_MK) .OR. (dum > 2.0_MK * 7.1000000E-04)) THEN
              iline = iline - 1 ! i.e. ignore this line
           END IF
           
        END DO
        
9999    CONTINUE
        
        
        ! allocate memory
        IF (ASSOCIATED(d_tracers%x)) THEN
           DEALLOCATE(d_tracers%x,STAT=stat_info_sub)
        END IF
           
        IF (ASSOCIATED(d_tracers%v)) THEN
           DEALLOCATE(d_tracers%v,STAT=stat_info_sub)
        END IF
        
        IF (ASSOCIATED(d_tracers%id)) THEN
           DEALLOCATE(d_tracers%id,STAT=stat_info_sub)
        END IF
        
        ALLOCATE(d_tracers%x(SIZE(x,1), iline), STAT=stat_info_sub) 
        ALLOCATE(d_tracers%v(SIZE(v,1), iline), STAT=stat_info_sub)
        ALLOCATE(d_tracers%id(num_id, iline), STAT=stat_info_sub)

        d_tracers%x(:, :)  = x(:, 1:iline)
        d_tracers%v(:, :)  = v(:, 1:iline)
        d_tracers%id(1, :) = id(1:iline)
        IF (num_id > 1) THEN
           d_tracers%id(2:num_id, :) = 0
        END IF
        
                   
        d_tracers%num_part_real = iline
        
        CALL physics_set_tracers_restart_num(this%phys, &
             iline, stat_info_sub)
                
        !----------------------------------------------------
        ! Release dynamic memory.
        !----------------------------------------------------
        
        IF(ASSOCIATED(x)) THEN
           DEALLOCATE(x)
        END IF
        
        IF(ASSOCIATED(v)) THEN
           DEALLOCATE(v)
        END IF
        
        IF(ASSOCIATED(id)) THEN
           DEALLOCATE(id)
        END IF
               
        
        RETURN
        
      END SUBROUTINE io_read_tracers
      
