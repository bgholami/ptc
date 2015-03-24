!------------------------------------------------------------
! Here are set functions/subroutines of Class Boundary
!------------------------------------------------------------

      SUBROUTINE boundary_set_num_dim(this,d_num_dim, stat_info)
        !-----------------------------------------
        ! Set the number of dimension of boundary.
        !-----------------------------------------
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        INTEGER, INTENT(IN)             :: d_num_dim
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        !---------------------------------------
        ! Only 2D, 3D are supported
        !---------------------------------------
        
        IF(d_num_dim < 2 .OR. d_num_dim > 3 ) THEN
           PRINT *, "boundary_set_num_dim : ", &
                "Dimension is not supported !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !---------------------------------------
        ! If reset the dimension, 
        ! the memory has to be reallocated
        !---------------------------------------
        
        IF( d_num_dim /= this%num_dim ) THEN
           
           this%num_dim = d_num_dim
           
           IF (ASSOCIATED(this%bcdef)) THEN
              DEALLOCATE(this%bcdef)
           END IF
           ALLOCATE(this%bcdef(2*d_num_dim)) 
           this%bcdef(:) = 0

           IF (ASSOCIATED(this%shear_rate)) THEN
              DEALLOCATE(this%shear_rate)
           END IF
           ALLOCATE(this%shear_rate(d_num_dim,d_num_dim))
           this%shear_rate(:,:) = 0.0_MK

           IF (ASSOCIATED(this%shear_length)) THEN
              DEALLOCATE(this%shear_length)
           END IF
           ALLOCATE(this%shear_length(d_num_dim,2*d_num_dim))
           this%shear_length(:,:) = 0.0_MK
           
           IF (ASSOCIATED(this%shear_type)) THEN
              DEALLOCATE(this%shear_type)
           END IF
           ALLOCATE(this%shear_type(2*d_num_dim))
           this%shear_type(:) = 0
    
           IF (ASSOCIATED(this%shear_v0)) THEN
              DEALLOCATE(this%shear_v0)
           END IF
           ALLOCATE(this%shear_v0(d_num_dim,2*d_num_dim))
           this%shear_v0(:,:) = 0.0_MK
           
           IF (ASSOCIATED(this%shear_v)) THEN
              DEALLOCATE(this%shear_v)
           END IF
           ALLOCATE(this%shear_v(d_num_dim,2*d_num_dim))
           this%shear_v(:,:) = 0.0_MK
           
           IF (ASSOCIATED(this%shear_freq)) THEN
              DEALLOCATE(this%shear_freq)
           END IF
           ALLOCATE(this%shear_freq(2*d_num_dim))
           this%shear_freq(:) = 0.0_MK
                      
           IF (ASSOCIATED(this%drag)) THEN
              DEALLOCATE(this%drag)
           END IF
           ALLOCATE(this%drag(d_num_dim,2*d_num_dim))
           this%drag(:,:) = 0.0_MK
           
           this%min_phys(:) = 0.0_MK
           this%max_phys(:) = 0.0_MK
           this%min_phys_t(:) = 0.0_MK
           this%max_phys_t(:) = 0.0_MK

           
        END IF
        
9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE boundary_set_num_dim
      
      
      SUBROUTINE boundary_set_bcdef(this,d_bcdef,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        INTEGER, DIMENSION(:)           :: d_bcdef
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        INTEGER                         :: i
        
        stat_info = 0
        
        dim = SIZE(d_bcdef,1)
        
        IF( dim /= 2*this%num_dim) THEN
           PRINT *, "boundary_set_bcdef : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%bcdef(1:dim) = d_bcdef(1:dim)
        
        this%num_peri       = 0
        this%num_sym        = 0
        this%num_wall_sym   = 0
        this%num_wall_solid = 0
        this%num_inflow     = 0
        this%num_outflow    = 0
        this%num_le         = 0
        
        DO i = 1, dim
           
           SELECT CASE ( d_bcdef(i) )
              
           CASE ( mcf_bcdef_periodic ) 
              
              this%num_peri = this%num_peri + 1
              
           CASE ( mcf_bcdef_symmetry )
              
              this%num_sym = this%num_sym + 1
              
           CASE( mcf_bcdef_wall_sym )
              
              this%num_wall_sym = this%num_wall_sym + 1
              
           CASE( mcf_bcdef_wall_solid )
              
              this%num_wall_solid = this%num_wall_solid + 1

           CASE( mcf_bcdef_inflow )
              
              this%num_inflow     = this%num_inflow + 1  

           CASE( mcf_bcdef_outflow )
              
              this%num_outflow     = this%num_outflow + 1

           CASE( mcf_bcdef_LE )
              
              this%num_le = this%num_le + 1
              
           END SELECT ! d_bcdef(i)
           
        END DO ! i = 1,  dim
       
        this%num_shear = this%num_wall_sym  + &
             this%num_wall_solid + &
             this%num_le
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_bcdef
      

      SUBROUTINE boundary_set_shear_type(this,d_shear_type,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        INTEGER, DIMENSION(:)           :: d_shear_type
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: i,dim
        
        stat_info = 0
        
        dim = SIZE(d_shear_type,1)
        
        IF( dim /= 2*this%num_dim) THEN
           PRINT *, "boundary_set_shear_type : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%shear_type(1:dim) = d_shear_type(1:dim)
        
        this%num_osci = 0
        
        DO i =1, dim
           IF ( d_shear_type(i) == 2 ) THEN
              this%num_osci = this%num_osci + 1
           END IF
        END DO
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_shear_type
      
      
      SUBROUTINE boundary_set_shear_rate(this,d_shear_rate,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), DIMENSION(:,:)        :: d_shear_rate
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim1
        INTEGER                         :: dim2

        stat_info = 0
        
        dim1 = SIZE(d_shear_rate,1)
        dim2 = SIZE(d_shear_rate,2)
        
        IF( dim1 /= this%num_dim .OR. &
             dim2 /= this%num_dim ) THEN
           PRINT *, "boundary_set_shear_rate : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%shear_rate(1:dim1,1:dim2) = &
             d_shear_rate(1:dim1,1:dim2) 
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_shear_rate


      SUBROUTINE boundary_set_shear_length(this,d_shear_length,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), DIMENSION(:,:)        :: d_shear_length
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim1
        INTEGER                         :: dim2

        stat_info = 0
        
        dim1 = SIZE(d_shear_length,1)
        dim2 = SIZE(d_shear_length,2)
        
        IF( dim1 /= this%num_dim .OR. &
             dim2 /= 2*this%num_dim ) THEN
           PRINT *, "boundary_set_shear_length : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%shear_length(1:dim1,1:dim2) = &
             d_shear_length(1:dim1,1:dim2) 
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_shear_length
      
      
      SUBROUTINE boundary_set_shear_v0(this,d_shear_v0,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), DIMENSION(:,:)        :: d_shear_v0
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim1
        INTEGER                         :: dim2

        stat_info = 0
        
        dim1 = SIZE(d_shear_v0,1)
        dim2 = SIZE(d_shear_v0,2)
        
        IF( dim1 /= this%num_dim .OR. &
             dim2 /= 2*this%num_dim ) THEN
           PRINT *, "boundary_set_shear_v0 : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%shear_v0(1:dim1,1:dim2) = &
             d_shear_v0(1:dim1,1:dim2)
        
        CALL boundary_set_shear_v(this,d_shear_v0(1:dim1,1:dim2),stat_info)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_shear_v0
      
      
      SUBROUTINE boundary_set_shear_v(this,d_shear_v,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), DIMENSION(:,:)        :: d_shear_v
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim1
        INTEGER                         :: dim2

        stat_info = 0
        
        dim1 = SIZE(d_shear_v,1)
        dim2 = SIZE(d_shear_v,2)
        
        IF( dim1 /= this%num_dim .OR. &
             dim2 /= 2*this%num_dim ) THEN
           PRINT *, "boundary_set_shear_v : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%shear_v(1:dim1,1:dim2) = &
             d_shear_v(1:dim1,1:dim2) 
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_shear_v
      
     
      SUBROUTINE boundary_set_shear_freq(this,d_shear_freq,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), DIMENSION(:)          :: d_shear_freq
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim

        stat_info = 0
        
        dim = SIZE(d_shear_freq,1)
        
        IF( dim /= 2*this%num_dim ) THEN
           PRINT *, "boundary_set_shear_freq : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%shear_freq(1:dim) = &
             d_shear_freq(1:dim) 
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_shear_freq

      
      SUBROUTINE boundary_set_wall_rho_type(this,d_rho_type,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        INTEGER                         :: d_rho_type
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF( d_rho_type < 0  .OR. d_rho_type > 1) THEN
           PRINT *, "boundary_set_wall_rho_type : ", &
                "Wrong rho type !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%rho_type =d_rho_type
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_wall_rho_type
      

      SUBROUTINE boundary_set_wall_noslip_type(this,d_noslip_type,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        INTEGER                         :: d_noslip_type
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%noslip_type =d_noslip_type
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_wall_noslip_type
      
      
      SUBROUTINE boundary_set_dout(this,d_dout,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), INTENT(IN)            :: d_dout
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%dout = d_dout
        
        RETURN
        
      END SUBROUTINE boundary_set_dout
      
      
      SUBROUTINE boundary_set_min_phys(this,d_min_phys,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), DIMENSION(:)          :: d_min_phys
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        
        stat_info = 0
        
        dim = SIZE(d_min_phys)
        
        IF( dim /= this%num_dim ) THEN
           PRINT *, "boundary_set_min_phys : ",&
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%min_phys(1:dim) =d_min_phys(1:dim)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_min_phys

      
      SUBROUTINE boundary_set_max_phys(this,d_max_phys,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), DIMENSION(:)          :: d_max_phys
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        
        stat_info = 0
        
        dim = SIZE(d_max_phys)
        
        IF( dim /= this%num_dim) THEN
           PRINT *, "boundary_set_max_phys : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%max_phys(1:dim) =d_max_phys(1:dim)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_max_phys
      

      SUBROUTINE boundary_set_min_phys_t(this,d_min_phys_t,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), DIMENSION(:)          :: d_min_phys_t
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        
        stat_info = 0
        
        dim = SIZE(d_min_phys_t)
        
        IF( dim /= this%num_dim) THEN
           PRINT *, "boundary_set_min_phys_t : ",&
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%min_phys_t(1:dim) =d_min_phys_t(1:dim)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_min_phys_t

      
      SUBROUTINE boundary_set_max_phys_t(this,d_max_phys_t,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), DIMENSION(:)          :: d_max_phys_t
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        
        stat_info = 0

        dim = SIZE(d_max_phys_t)

        IF( dim /= this%num_dim) THEN
           PRINT *, "boundary_set_max_phys_t : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF

        this%max_phys_t(1:dim) =d_max_phys_t(1:dim)

9999    CONTINUE

        RETURN

      END SUBROUTINE  boundary_set_max_phys_t


      SUBROUTINE boundary_set_num_part_wall_solid(this,num_part,stat_info)

        !---------------------------------------
        ! Set the num of wall boundary particles
        ! created by MCF.
        !---------------------------------------

        TYPE(Boundary), INTENT(INOUT)   :: this
        INTEGER, INTENT(IN)             :: num_part
        INTEGER, INTENT(OUT)            :: stat_info

        stat_info = 0

        this%num_part_wall_solid = num_part

        RETURN

      END SUBROUTINE  boundary_set_num_part_wall_solid


      SUBROUTINE boundary_set_inflow_outflow(this, in_file, stat_info)
        !-----------------------------------------
        ! Set up the inflow/outflow boundary patches
        !-----------------------------------------

        TYPE(Boundary), INTENT(INOUT)   :: this
        CHARACTER(LEN=MAX_CHAR)         :: in_file
        INTEGER, INTENT(OUT)            :: stat_info  

        INTEGER                         :: ilenread
        LOGICAL                         :: lExist
        INTEGER                         :: num_dim
        REAL(MK), DIMENSION(3)          :: pp
        REAL(MK)                        :: dd
        INTEGER                         :: i, j, it, t1, t2


        stat_info = 0


        ! check the file
        ilenread = LEN_TRIM(in_file)

        !----------------------------------------------------
        ! Check if a file name is given.
        !----------------------------------------------------   

        IF ( ilenread < 1 ) THEN
           PRINT *,'bcdef_inout_file : ',&
                'No file name is given !'
           stat_info = -1
           GOTO 9999

        ELSE

           !----------------------------------------------------
           ! Check if the file exists. 
           !----------------------------------------------------

           INQUIRE(FILE=in_file,EXIST=lExist)
           
           IF (.NOT.lExist) THEN 
              PRINT *,'NO such file : ',&
                   in_file
              stat_info = -1
              GOTO 9999
           END IF
        END IF


        NULLIFY(this%num_iopoints)  
        NULLIFY(this%patch_id)
        NULLIFY(this%iopatch_n)  
        NULLIFY(this%iopatch_x)
        NULLIFY(this%ref_point)
        NULLIFY(this%clength)
        NULLIFY(this%total_num_grid)
        NULLIFY(this%BC_time)
        NULLIFY(this%num_grid)
        NULLIFY(this%tvec)
        NULLIFY(this%dtvec)
        NULLIFY(this%grido)
        NULLIFY(this%BC_vel)

        num_dim = this%num_dim 


        ! read input file
        OPEN(10,  FILE=in_file) 

        READ(10, *)  this%num_inflow, this%num_outflow
        this%num_inout = this%num_inflow + this%num_outflow

        ! allocate  
        IF (this%num_inout > 0) THEN
           ALLOCATE(this%num_iopoints(this%num_inout))   
           ALLOCATE(this%patch_id(this%num_inout))
           ALLOCATE(this%iopatch_n(this%num_inout, num_dim))
           ALLOCATE(this%ref_point(this%num_inout, num_dim))
           ALLOCATE(this%clength(this%num_inout))
        END IF

        ! read number of point in all patches
        READ(10, *) this%num_iopoints(1:this%num_inout) 
        ! read patch id of all patches (positive for inlet, negative for outlet)
        READ(10, *) this%patch_id(1:this%num_inout)

        ! allocate
        ALLOCATE(this%iopatch_x(this%num_inout, 1:num_dim, MAXVAL(this%num_iopoints)))

        DO j = 1, this%num_inout

           ! read patch normal
           READ(10, *) this%iopatch_n(j, 1:num_dim)

           ! read patch points
           DO i = 1, this%num_iopoints(j)
              READ(10, *) this%iopatch_x(j, 1:num_dim, i)
           END DO

        END DO

        ALLOCATE(this%total_num_grid(2))

        ! read number of time and space grid points
        READ(10, *) this%num_time, this%total_num_grid(1:2)

        ! allocate 
        ALLOCATE(this%BC_time(this%num_time))
        ALLOCATE(this%num_grid(this%num_inout, 2))
        ALLOCATE(this%tvec(this%num_inout, 2, num_dim))
        ALLOCATE(this%dtvec(this%num_inout, 2))
        ALLOCATE(this%grido(this%num_inout, num_dim))
        ALLOCATE(this%BC_vel(this%num_inout, num_dim, &
             this%num_time, PRODUCT(this%total_num_grid)))

        ! read time instances
        DO j = 1, this%num_time
           READ(10, *) this%BC_time(j)
        END DO

        this%BC_vel = 0.0_MK
        DO j = 1, this%num_inout
           ! read number of grid points of patch
           READ(10, *) this%num_grid(j, 1:2)
           
           ! read unit vectors of patch
           DO i = 1, 2
              READ(10, *) this%tvec(j, i, 1:num_dim)
           END DO

           ! read grid size alogn unit vectors
           READ(10, *) this%dtvec(j, 1:2)
           ! read origin of grid
           READ(10, *) this%grido(j, 1:num_dim)

           ! read velocity over time and space
           DO it = 1, this%num_time
              i = 0
              DO t1 = 1, this%num_grid(j, 1)
                 DO t2 = 1, this%num_grid(j, 2)

                    i = i + 1
                    READ(10, *) this%BC_vel(j, 1:num_dim, it, i)

                 END DO
              END DO
           END DO

        END DO

        CLOSE(10)

        ! calculate ref_point
        DO j = 1, this%num_inout
           
           pp(1:num_dim) = 0.0_MK
           
           DO i = 1, this%num_iopoints(j)
              pp(1:num_dim) = pp(1:num_dim) + this%iopatch_x(j, 1:num_dim, i)
           END DO

           this%ref_point(j, 1:num_dim) = pp(1:num_dim) / DBLE(this%num_iopoints(j))

        END DO  

        ! calculate clength
        DO j = 1, this%num_inout

           dd = 0.0_MK
           DO i = 1, this%num_iopoints(j)
              pp(1:num_dim) = this%iopatch_x(j, 1:num_dim, i) - this%ref_point(j, 1:num_dim)
              dd = MAX(dd, SQRT(DOT_PRODUCT(pp(1:num_dim), pp(1:num_dim))))
           END DO

           this%clength(j) = dd

        END DO

        ! default tstart
        this%tstart = 0.0_MK


9999    CONTINUE

        RETURN

      END SUBROUTINE boundary_set_inflow_outflow  

      
      SUBROUTINE boundary_set_bwidth(this,d_bwidth,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), INTENT(IN)            :: d_bwidth
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%bwidth = d_bwidth
        
        RETURN
        
      END SUBROUTINE boundary_set_bwidth


      SUBROUTINE boundary_set_tstart(this, dt, h, drho, stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), INTENT(IN)            :: dt   
        REAL(MK), INTENT(IN)            :: h
        REAL(MK), INTENT(IN)            :: drho
        INTEGER, INTENT(OUT)            :: stat_info

        INTEGER                         :: num_dim, i, j
        REAL(MK)                        :: vc, v0, tstart, itr
        
        stat_info = 0
        num_dim = this%num_dim

        ! calculate a characteristic velocity
        vc = (0.250_MK * h * drho / dt) ** 0.50_MK

        ! evaluate the maximum magnitude of velocity at the beginning of pulse
        v0 = 0.0_MK
        DO j = 1, this%num_inout

           DO i = 1, this%num_grid(j, 1) * this%num_grid(j, 2)
              
              v0 = MAX(v0, DOT_PRODUCT(this%BC_vel(j, 1:num_dim, 1, i), &
                   this%BC_vel(j, 1:num_dim, 1, i))**0.50_MK)

           END DO

        END DO
        
        ! nominal start time
        tstart = v0 / vc

        ! adapt (and increase)
        itr = tstart * 1.50_MK / dt ! even 1.5 times larger
        ! at lease starts from iteration 10000 (with increaments of 10000)
        itr = CEILING(itr/10000.0_MK) * 10000.0_MK
        itr = 30000.0_MK
        tstart = itr * dt

        this%tstart = tstart
        PRINT*, 'Inflow/Outflow starting iteration = ', itr

        
        RETURN
        
      END SUBROUTINE boundary_set_tstart

      


