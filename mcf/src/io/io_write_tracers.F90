     SUBROUTINE io_write_tracers(this,rank,step,&
           d_tracers,num_part_in,stat_info)
        !----------------------------------------------------
        ! Subroutine  : io_write_tracers
        !----------------------------------------------------
        !
        ! Purpose     : Writing tracers' quantities into
        !               files, which is centralized.
        !               x, v, species ID (if more than one species),
        !               particle ID.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.1 09.05.2011, original version.
        !               The routine is derived from
        !               io_write_particles.F90: V0.2.
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
        !
        ! this           : an object of Marching Class.
        ! rank           : rank of process.
        ! step           : index of current step.
        ! d_tracers      : an object of Tracers Class
        ! num_part       : first num_part needed to be written.
        ! stat_info      : return flag of status.
        !----------------------------------------------------
        
        !----------------------------------------------------
      	! Modules
      	!----------------------------------------------------
        
        USE ppm_module_io
        
        
      	!----------------------------------------------------
      	! Arguments
      	!----------------------------------------------------
        
        TYPE(IO),INTENT(IN)                     :: this
        INTEGER, INTENT(IN)                     :: rank
        INTEGER, INTENT(IN)                     :: step
        TYPE(Particles), INTENT(IN)             :: d_tracers
        INTEGER, INTENT(IN)                     :: num_part_in
        INTEGER,INTENT(OUT)                     :: stat_info  


        !----------------------------------------------------
        ! Colloid parameters :
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
        LOGICAL                                 :: read_external
        
        REAL(MK), DIMENSION(:,:), POINTER       :: x
        REAL(MK), DIMENSION(:,:), POINTER       :: v 
        REAL(MK), DIMENSION(:)  , POINTER       :: h_p
        INTEGER,  DIMENSION(:,:), POINTER       :: id
        
        CHARACTER(LEN=MAX_CHAR)                 :: file_name
        INTEGER                                 :: ifmt,output_unit
        
        INTEGER                                 :: num_x, num_dim, num_part
        INTEGER                                 :: num_v
        INTEGER                                 :: num_id
        INTEGER                                 :: data_dim
        INTEGER                                 :: current_dim
        INTEGER                                 :: total_num_interface, write_index
        REAL(MK), DIMENSION(:,:), POINTER       :: output
        
        CHARACTER(LEN=MAX_CHAR)                 :: cbuf
        INTEGER					:: clen

#ifdef __DEBUG
        INTEGER                                 :: debug_flag
        REAL(MK)                                :: time_routine_start
#endif
        INTEGER                                 :: i, ic, k, vi, vperf
        REAL(MK), DIMENSION(:), ALLOCATABLE     :: mid_point, dumv
        ! temporary deposited output XX
        CHARACTER(LEN=MAX_CHAR)                 :: deposited_file
!!$        integer            :: ip
!!$        real(MK)           :: ft, fn, vt, vn, S, FS, eta
!!$        REAL(MK), DIMENSION(:), ALLOCATABLE     :: v_point 
!!$        INTEGER , DIMENSION(:)  , ALLOCATABLE   :: i_index
     
	!----------------------------------------------------
      	! Initialization of variables.
      	!----------------------------------------------------

	stat_info     = 0
        stat_info_sub = 0

        NULLIFY(x)
        NULLIFY(v)
        NULLIFY(h_p)
        NULLIFY(id)      
        NULLIFY(colloids)
        NULLIFY(coll_f_num)
        NULLIFY(coll_v)  
        NULLIFY(coll_f_vlist) 
        NULLIFY(output) 


#ifdef __DEBUG
        debug_flag = debug_get_flag(global_debug,stat_info_sub)
        IF( debug_flag == 2 ) THEN
           CALL debug_substart(global_debug,rank,"io_write_particles",&
                time_routine_start,stat_info_sub)
        END IF
#endif

	!----------------------------------------------------
        ! Get parameters.
	!----------------------------------------------------

        num_dim = d_tracers%num_dim

        num_part = num_part_in

        read_external   = &
             control_get_read_external(this%ctrl,stat_info_sub)

        if (num_part > 0) then
           CALL tracers_get_x(d_tracers,x,num_part,stat_info)
           CALL tracers_get_v(d_tracers,v,num_part,stat_info)

           IF(ASSOCIATED(h_p)) THEN 
              DEALLOCATE(h_p)
           END IF
           
           ALLOCATE(h_p(num_part))  
           h_p(1:num_part) = d_tracers%h_p

           CALL tracers_get_id(d_tracers,id,num_part,stat_info)

        else

           num_part = 1

           IF(ASSOCIATED(x)) THEN 
              DEALLOCATE(x)
           END IF

           ALLOCATE(x(num_dim, num_part)) 
           x = 0.0_MK

           IF(ASSOCIATED(v)) THEN 
              DEALLOCATE(v)
           END IF

           ALLOCATE(v(num_dim, num_part))  
           v = 0.0_MK  

           IF(ASSOCIATED(h_p)) THEN 
              DEALLOCATE(h_p)
           END IF

           ALLOCATE(h_p(num_part))  
           h_p = 0.0_MK

           IF(ASSOCIATED(id)) THEN 
              DEALLOCATE(id)
           END IF

           ALLOCATE(id(3, num_part))
           id(1, :) = 0
           id(2, :) = -1
           id(3, :) = 0

        end if

!!$        ! temporary output of force by parallel and normal components  
!!$        vperf = num_dim
!!$        ALLOCATE(i_index(vperf)) 
!!$        ALLOCATE(v_point(num_dim)) 
!!$        CALL physics_get_colloid(this%phys, &
!!$             colloids,stat_info_sub)
!!$        CALL colloid_get_coll_f_vlist(colloids, coll_f_vlist, stat_info_sub)  
!!$        DO ip = 1, d_tracers%num_part_real
!!$
!!$           fn = DOT_PRODUCT(d_tracers%f(1:num_dim, ip), d_tracers%n_vec(1:num_dim, ip))
!!$           ft = DOT_PRODUCT(d_tracers%f(1:num_dim, ip) - fn * d_tracers%n_vec(1:num_dim, ip), &
!!$                d_tracers%f(1:num_dim, ip) - fn * d_tracers%n_vec(1:num_dim, ip)) ** 0.5_MK
!!$           v(1, ip) = fn
!!$           v(2, ip) = ft
!!$           v(3, ip) = 0.0_MK
!!$
!!$
!!$           FS = 1.66810_MK
!!$           eta = physics_get_eta(this%phys,stat_info_sub) 
!!$           ! find and copy corresponding interface points
!!$           i_index(1:vperf) = coll_f_vlist(d_tracers%col_i(ip), 1:vperf, d_tracers%facet(ip))
!!$
!!$           v_point = 0.0_MK
!!$           DO i = 1, vperf
!!$              v_point(1:num_dim) = v_point(1:num_dim) + &
!!$                   d_tracers%vc_interface(d_tracers%col_i(ip), 1:num_dim, i_index(i))
!!$           END DO
!!$           v_point = v_point / dble(vperf)
!!$
!!$           vn = DOT_PRODUCT(v_point, d_tracers%n_vec(1:num_dim, ip))
!!$           vt = DOT_PRODUCT(v_point - vn * d_tracers%n_vec(1:num_dim, ip), &
!!$                v_point - vn * d_tracers%n_vec(1:num_dim, ip)) ** 0.5_MK
!!$
!!$           S = vt / d_tracers%near_wall  
!!$           v(3, ip) = 6.0_MK * mcf_pi * d_tracers%radius * d_tracers%h_p(ip) * eta * S * FS
!!$
!!$        END DO


        !----------------------------------------------------
      	! Define the output file name for this time step.
        ! If we have read tracers from file,
        ! the output should include "_r" to distingusih.
      	!----------------------------------------------------

      	!----------------------------------------------------
      	! Define format of output file.
      	!----------------------------------------------------

        IF( read_external ) THEN

           WRITE(file_name,'(2A,I8.8,A)') &
                TRIM(this%output_tracers_file),'_r',step,'.out'

        ELSE

           WRITE(file_name,'(A,I8.8,A)') &
                TRIM(this%output_tracers_file),step,'.out'

        END IF

        IF (this%output_tracers_fmt(1:1) .EQ. 'f' .OR. &
             this%output_tracers_fmt(1:1) .EQ. 'F') THEN

           ifmt = ppm_param_io_ascii

        ELSE

           ifmt = ppm_param_io_binary

        END IF

        output_unit = this%output_tracers_unit

        !----------------------------------------------------
        ! Allocate memory for output.
        !----------------------------------------------------

        num_x  = SIZE(x,1)
        num_v  = SIZE(v,1)

        ! id's dimension is in fact (d_tracers%num_id, :)
        ! however, only its first column is interesting for output.
        ! Therefore, num_id equals 1 eventhough id's dimension is larger.
        num_id = 1 


        !----------------------------------------------------
        ! Allocate memory for output data.
        ! x,y(,z), vx,vy(,vz), pid,sid.
        !----------------------------------------------------

        data_dim = num_x + num_v + num_id + 1 

        ALLOCATE(output(data_dim,num_part),STAT=stat_info_sub)

        IF (stat_info_sub /=0) THEN          
           PRINT *, 'io_write_tarcers : ',&
                'Allocating buffer for output failed !'
           stat_info = -1
           GOTO 9999          
        END IF

        !----------------------------------------------------
        ! Copy the data into output.
        !----------------------------------------------------

        output(1:num_x,1:num_part) = x(1:num_x,1:num_part)        

        current_dim = num_x

        output(current_dim+1:current_dim+num_v,1:num_part) = &
             v(1:num_v,1:num_part)        

        current_dim = current_dim + num_v

        output(current_dim+1, 1:num_part) = &
             h_p(1:num_part)
        current_dim = current_dim + 1

        output(current_dim+1:current_dim+num_id,1:num_part) = &
             id(1:num_id,1:num_part)        

        current_dim = current_dim + num_id


        !----------------------------------------------------
        ! Open ppm I/O unit for centralized I/O.
        !----------------------------------------------------

        CALL ppm_io_open(output_unit,file_name,&
             ppm_param_io_write, ppm_param_io_replace, &
             ifmt,ppm_param_io_centralized,stat_info_sub)

        IF (stat_info_sub /= 0) THEN          
           PRINT *, 'io_write_tracers : ', &
                'Failed to open unit !'          
           stat_info = -1
           GOTO 9999          
        END IF


        WRITE(cbuf,'(A1,I2,A6)') '(', data_dim ,'E16.8)'
        clen = LEN_TRIM(cbuf)         

        CALL ppm_io(output_unit, output, &
             ACTN=ppm_param_io_write, &
             DIST=ppm_param_io_concat, &
             IOFMT=cbuf(1:clen), &
             STAT=stat_info_sub)

        IF ( stat_info_sub /= 0 ) THEN          
           PRINT *, 'io_write_tracers : ',&
                'Writing tracers failed'	   
           stat_info = -1
           GOTO 9999
        END IF

        !----------------------------------------------------
        ! Close file.
        !----------------------------------------------------

        CALL ppm_io_close(output_unit,stat_info_sub)  


        !----------------------------------------------------
        ! Write the deposition data
        !----------------------------------------------------


        ! output for deposited tracers
        WRITE(deposited_file,'(A,I8.8,A)') &
             TRIM('deposited_tracers'),step,'.out'


        ! allocate memory (reuse output)
        data_dim = num_dim + 1

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

           ! determine size of output (local part of interface)
           total_num_interface = 0
           Do ic = 1, coll_arbitrary_num
              total_num_interface = total_num_interface + coll_f_num(ic)
           END Do


           DEALLOCATE(output)
           ALLOCATE(output(data_dim, total_num_interface),STAT=stat_info_sub)

           IF (stat_info_sub /=0) THEN          
              PRINT *, 'io_write_tarcers (deposited) : ',&
                   'Allocating buffer for output failed !'
              stat_info = -1
              GOTO 9999          
           END IF



           vperf = num_dim
           ALLOCATE(mid_point(num_dim)) 
           ALLOCATE(dumv(num_dim))
           ! copy data into output
           write_index = 0
           DO ic = 1, coll_arbitrary_num

              DO i = 1, coll_f_num(ic)
                 write_index = write_index + 1

                 ! find the mid point
                 mid_point = 0.0_MK
                 dumv = 0.0_MK
                 DO k = 1, vperf
                    vi = coll_f_vlist(ic, k, i)
                    mid_point = mid_point + coll_v(ic, 1:num_dim, vi)

                    !dumv = dumv + d_tracers%vn_interface(ic, 1:num_dim, vi)
                 END DO
                 mid_point = mid_point / DBLE(vperf)
                 !dumv      = dumv      / DBLE(vperf)

                 output(1:num_dim, write_index) = &
                      mid_point(1:num_dim)

                 !output(num_dim+1:2*num_dim, write_index) = &
                 !     dumv

                 output(data_dim, write_index) =  &
                      DBLE(d_tracers%num_tracer_deposited(ic, i))  
                 

                 

              END DO

           END DO


           ! not a great hack, we sum up output over all subdomains,
           ! but the first num_dim enteris are location.
           IF (rank /= 0) THEN
              output(1:num_dim, :) = 0.0_MK
           END IF


           !----------------------------------------------------
           ! Open ppm I/O unit
           !----------------------------------------------------

           CALL ppm_io_open(output_unit, deposited_file, &
                ppm_param_io_write, ppm_param_io_replace, &
                ifmt, ppm_param_io_centralized, stat_info_sub)

           IF (stat_info_sub /= 0) THEN          
              PRINT *, 'io_write_tracers (deposited) : ', &
                   'Failed to open unit !'          
              stat_info = -1
              GOTO 9999          
           END IF


           WRITE(cbuf,'(A1,I2,A6)') '(', data_dim ,'E16.8)'
           clen = LEN_TRIM(cbuf)         

           CALL ppm_io(output_unit, output, &
                ACTN=ppm_param_io_write, &
                DIST=ppm_param_io_sum, &
                IOFMT=cbuf(1:clen), &
                STAT=stat_info_sub)

           IF ( stat_info_sub /= 0 ) THEN          
              PRINT *, 'io_write_tracers (deposited) : ',&
                   'Writing tracers failed'	   
              stat_info = -1
              GOTO 9999
           END IF

           !----------------------------------------------------
           ! Close file.
           !----------------------------------------------------

           CALL ppm_io_close(output_unit, stat_info_sub) 

        END IF
        
        !----------------------------------------------------
        ! Print out for user to confirm.
        !----------------------------------------------------
        
        IF (rank == 0) THEN
           
           WRITE(cbuf,'(2A)') 'Tracers written to ',&
                TRIM(file_name)
           PRINT *, TRIM(cbuf)
           
        END IF   
        
        !----------------------------------------------------
        ! Return.
        !----------------------------------------------------
        
9999    CONTINUE
        
        IF (ASSOCIATED(x)) THEN
           DEALLOCATE(x)
        END IF
        
        IF (ASSOCIATED(v)) THEN
           DEALLOCATE(v)
        END IF 

        IF (ASSOCIATED(h_p)) THEN
           DEALLOCATE(h_p)
        END IF
        
        IF (ASSOCIATED(id)) THEN
           DEALLOCATE(id)
        END IF  

        IF(ASSOCIATED(coll_f_num)) THEN
           DEALLOCATE(coll_f_num)
        END IF
  
        IF(ASSOCIATED(coll_v))THEN
           DEALLOCATE(coll_v)
        END IF
        
        IF(ASSOCIATED(coll_f_vlist))THEN
           DEALLOCATE(coll_f_vlist)
        END IF
        
        IF (ASSOCIATED(output)) THEN
           DEALLOCATE(output)
        END IF
        
        
#ifdef __DEBUG
        IF( debug_flag == 2 ) THEN
           CALL debug_substop(global_debug,rank,&
                "io_write_tracers",&
                time_routine_start,stat_info_sub)
        END IF
#endif
        RETURN	
        
        
      END SUBROUTINE io_write_tracers   
     
