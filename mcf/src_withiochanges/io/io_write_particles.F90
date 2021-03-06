     SUBROUTINE io_write_particles(this,rank,step,&
           d_particles,num_part_in,stat_info)
        !----------------------------------------------------
        ! Subroutine  : io_write_particle
        !----------------------------------------------------
        !
        ! Purpose     : Writing particles' quantities into
        !               files, which is centralized.
        !               x, v, mass/number density, mass,
        !               species ID(if more than one species),
        !               particle ID.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.2 09.07.2009, 
        !               check again the work flow is correct
        !               and supply with more comments for code.
        !  
        !               V0.1 01.02.2009, original version.
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
        !
        ! this           : an object of Marching Class.
        ! rank           : rank of process.
        ! step           : index of current step.
        ! d_particles    : an object of Particles Class
        ! num_part_in    : first num_part needed to be written.
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
        TYPE(Particles), INTENT(IN)             :: d_particles
        INTEGER, INTENT(IN)                     :: num_part_in
        INTEGER,INTENT(OUT)                     :: stat_info
        
      	!----------------------------------------------------
      	! Local variables 
      	!----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        LOGICAL                                 :: read_external
        LOGICAL                                 :: p_energy
        
        REAL(MK), DIMENSION(:,:), POINTER       :: x
        REAL(MK), DIMENSION(:,:), POINTER       :: v
        REAL(MK), DIMENSION(:),   POINTER       :: rho
        REAL(MK), DIMENSION(:),   POINTER       :: m
#if __TRACER
        REAL(MK),  DIMENSION(:),   POINTER       :: c_tracers
#endif
        INTEGER,  DIMENSION(:,:), POINTER       :: id
        REAL(MK), DIMENSION(:,:), POINTER       :: f
        REAL(MK), DIMENSION(:),   POINTER       :: u
        
        CHARACTER(LEN=MAX_CHAR)                 :: file_name
        INTEGER                                 :: ifmt,output_unit
        
        INTEGER                                 :: num_x, num_dim, num_part
        INTEGER                                 :: num_v
        INTEGER                                 :: num_id
#if __WRITE_FORCE        
        INTEGER                                 :: num_f
#endif
        INTEGER                                 :: data_dim
        INTEGER                                 :: current_dim
        REAL(MK), DIMENSION(:,:), POINTER       :: output
        
        CHARACTER(LEN=MAX_CHAR)                 :: cbuf
        INTEGER					:: clen

#ifdef __DEBUG
        INTEGER                                 :: debug_flag
        REAL(MK)                                :: time_routine_start
#endif

	!----------------------------------------------------
      	! Initialization of variables.
      	!----------------------------------------------------
        
	stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(x)
        NULLIFY(v)
        NULLIFY(rho)
        NULLIFY(m)
#if __TRACER
        NULLIFY(c_tracers)
#endif
        NULLIFY(id)
        NULLIFY(f)
        NULLIFY(u)
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

        num_dim = d_particles%num_dim

        num_part = num_part_in

        read_external   = &
             control_get_read_external(this%ctrl,stat_info_sub)
        p_energy        = &
             control_get_p_energy(this%ctrl,stat_info_sub)

        if (num_part > 0) then
           CALL particles_get_x(d_particles,x,num_part,stat_info)
           CALL particles_get_v(d_particles,v,num_part,stat_info)    
           CALL particles_get_rho(d_particles,rho,num_part,stat_info)        
           CALL particles_get_m(d_particles,m,num_part,stat_info)

#if __TRACER
           ALLOCATE(c_tracers(num_part))
           c_tracers(1:num_part) =  d_particles%c_tracers(1:num_part)  ! XX     
#endif
           
#ifdef __WRITE_FORCE
           CALL particles_get_f(d_particles,f,num_part,stat_info)
#endif
           
           CALL particles_get_id(d_particles,id,num_part,stat_info)
           
           IF( p_energy ) THEN
              CALL particles_get_u(d_particles,u,num_part,stat_info)
           END IF

        else   

           num_part = 1

           ALLOCATE(x(num_dim, num_part)) 
           x = 0.0_MK  

           ALLOCATE(v(num_dim, num_part))  
           v = 0.0_MK 

           ALLOCATE(rho(num_part))  
           rho = 0.0_MK 
           
           ALLOCATE(m(num_part))  
           m = 0.0_MK

#if __TRACER
           ALLOCATE(c_tracers(num_part))
           c_tracers = 0.0_MK
#endif

#ifdef __WRITE_FORCE
           ALLOCATE(f(num_dim, num_part))  
           f = 0.0_MK 
#endif
           
           ALLOCATE(id(3, num_part))
           id(1, :) = 0
           id(2, :) = -1
           id(3, :) = 0
           
           IF( p_energy ) THEN            
              ALLOCATE(u(num_part))  
              u = 0.0_MK 
           END IF

        end if

        !----------------------------------------------------
      	! Define the output file name for this time step.
        ! If we have read particles from file,
        ! the output should include "_r" to distingusih.
      	!----------------------------------------------------
        
      	!----------------------------------------------------
      	! Define format of output file.
      	!----------------------------------------------------
        
        IF( read_external ) THEN
           
           WRITE(file_name,'(2A,I8.8,A)') &
                TRIM(this%output_particles_file),'_r',step,'.out'
           
        ELSE
           
           WRITE(file_name,'(A,I8.8,A)') &
                TRIM(this%output_particles_file),step,'.out'
           
        END IF
        
        IF (this%output_particles_fmt(1:1) .EQ. 'f' .OR. &
             this%output_particles_fmt(1:1) .EQ. 'F') THEN
           
           ifmt = ppm_param_io_ascii
           
        ELSE
           
           ifmt = ppm_param_io_binary
           
        END IF
        
        output_unit = this%output_particles_unit
        
        !----------------------------------------------------
        ! Allocate memory for output.
        !----------------------------------------------------
        
        num_x  = SIZE(x,1)
        num_v  = SIZE(v,1)

        num_id = 3
        
        !----------------------------------------------------
        ! Allocate memory for output data.
        ! x,y(,z), vx,vy(,vz), rho/d, m, pid,sid; f, u.
        !----------------------------------------------------
        
        data_dim = num_x + num_v + 1 + 1 + num_id
        
#ifdef __WRITE_FORCE
        num_f = SIZE(f,1)
        data_dim = data_dim + num_f
#endif
        
        IF (p_energy )  THEN
           data_dim = data_dim + 1
        END IF
        
        ALLOCATE(output(data_dim,num_part),STAT=stat_info_sub)
        
        IF (stat_info_sub /=0) THEN          
           PRINT *, 'io_write_particles : ',&
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
        
        output(current_dim+1,1:num_part) = rho(1:num_part)
        output(current_dim+2,1:num_part) = m(1:num_part)

#if __TRACER
output(current_dim+2,1:num_part) =  c_tracers(1:num_part)  ! XX     
#endif
        
        current_dim = current_dim + 2
           
        output(current_dim+1:current_dim+num_id,1:num_part) = &
             id(1:num_id,1:num_part)

!#if __TRACER
!output(current_dim+1,1:num_part) =  rank  ! XX 
!#endif
        
        current_dim = current_dim + num_id
        
#ifdef __WRITE_FORCE
        output(current_dim+1:current_dim+num_f,1:num_part) = &
             f(1:num_f,1:num_part)
        current_dim = current_dim + num_f
#endif
        
        IF( p_energy ) THEN
           
           output(current_dim+1,1:num_part) = u(1:num_part)
           current_dim = current_dim + 1
           
        END IF


        !----------------------------------------------------
        ! Open ppm I/O unit for centralized I/O.
        !----------------------------------------------------
        
        CALL ppm_io_open(output_unit,file_name,&
             ppm_param_io_write, ppm_param_io_replace, &
             ifmt,ppm_param_io_centralized,stat_info_sub)
        
        IF (stat_info_sub /= 0) THEN          
           PRINT *, 'io_write_particles : ', &
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
           PRINT *, 'io_write_particles : ',&
                'Writing particles failed'	   
           stat_info = -1
           GOTO 9999
        END IF
        

        !----------------------------------------------------
        ! Close file.
        !----------------------------------------------------
        
        CALL ppm_io_close(output_unit,stat_info_sub)
        
        !----------------------------------------------------
        ! Print out for user to confirm.
        !----------------------------------------------------
        
        IF (rank == 0) THEN
           
           WRITE(cbuf,'(2A)') 'Particles written to ',&
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
        
        IF (ASSOCIATED(rho)) THEN
           DEALLOCATE(rho)
        END IF
        
        IF (ASSOCIATED(m)) THEN
           DEALLOCATE(m)
        END IF

#if __TRACER
        IF (ASSOCIATED(c_tracers)) THEN
           DEALLOCATE(c_tracers)
        END IF
#endif
       
        IF (ASSOCIATED(id)) THEN
           DEALLOCATE(id)
        END IF
        
        IF (ASSOCIATED(f)) THEN
           DEALLOCATE(f)
        END IF
        
        IF (ASSOCIATED(u) ) THEN
           DEALLOCATE(u)
        END IF
        
        IF (ASSOCIATED(output)) THEN
           DEALLOCATE(output)
        END IF
        
        
#ifdef __DEBUG
        IF( debug_flag == 2 ) THEN
           CALL debug_substop(global_debug,rank,&
                "io_write_particles",&
                time_routine_start,stat_info_sub)
        END IF
#endif
        RETURN	
        
        
      END SUBROUTINE io_write_particles   
     
