      SUBROUTINE io_write_boundary(this,&
           rank,step,time,d_boundary,stat_info)
        !----------------------------------------------------
        !  Subroutine   :  io_write_boundary
        !----------------------------------------------------
        !
        !  Purpose      :  Write information into output files.
        !
        !  Reference    :
        !
        !  Remark       :
        !
        !  Revisions    : 
        !                 V0.1 13.01.2009, original version,
        !
        !----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        !--------------------------------
	! Arguments
	!--------------------------------
        TYPE(IO), INTENT(IN)            :: this
        INTEGER, INTENT(IN)             :: rank
        INTEGER,  INTENT(IN)	        :: step
        REAL(MK), INTENT(IN)	        :: time
        TYPE(Boundary), INTENT(IN)      :: d_boundary
        INTEGER,  INTENT(OUT)	        :: stat_info
        
        !--------------------------------
        !  Local Variables
        !--------------------------------
        INTEGER                         :: stat_info_sub
        INTEGER                         :: num_dim
        INTEGER                         :: num_shear
        INTEGER,DIMENSION(:), POINTER   :: bcdef
        REAL(MK),DIMENSION(:,:),POINTER :: drag
        
        REAL(MK),DIMENSION(12)          :: output
        INTEGER                         :: j_start,j_end
        CHARACTER(len=MAX_CHAR)	        :: cbuf,fbuf
        INTEGER			        :: i

        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(bcdef)
        NULLIFY(drag)
        
        IF ( rank /= 0 ) THEN
           PRINT *, &
                "io_write_boundary can only be used by root process !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !--------------------------------
        ! Get 
        ! number of dimension;
        ! number of shear;
        ! boundary conditions;
        ! drags on the walls.
        !--------------------------------
        
        num_dim   = boundary_get_num_dim(d_boundary,stat_info_sub)
        num_shear = boundary_get_num_shear(d_boundary,stat_info_sub)

        !--------------------------------
        ! If there is no wall, nothing
        ! to write yet.
        !--------------------------------
        
        IF ( num_shear > 0) THEN
           
           CALL boundary_get_bcdef(d_boundary,bcdef,stat_info_sub)
           CALL boundary_get_drag(d_boundary,drag,stat_info_sub)
           
           !-------------------
           ! Record current time.
           !-------------------
           
           j_start = 1
           j_end   = 1
           
           output(j_start) = time
           
           DO i = 1, 2*num_dim
              
              IF ( bcdef(i) == mcf_bcdef_wall_sym .OR. &
                   bcdef(i) == mcf_bcdef_wall_solid .OR. &
                   bcdef(i) == mcf_bcdef_LE ) THEN 
                 
                 !-------------------
                 ! Record force/drag.
                 !-------------------
                 
                 j_start = j_end + 1
                 j_end   = j_start + num_dim -1
                 
                 output(j_start:j_end) = &
                      drag(1:num_dim,i)
                 
              END IF
              
           END DO
           
           
           WRITE(fbuf, '(A,I2,A)' ) , '(I9,', j_end, 'E16.8)'
           WRITE(cbuf,fbuf )  step, output(1:j_end)
           
           WRITE(UNIT=this%boundary_unit,FMT='(A)',&
                IOSTAT=stat_info_sub)  TRIM(cbuf)
           
           IF( stat_info_sub /= 0 ) THEN
              PRINT *,"io_write_boundary : ",&
                   "Writting into boundary file failed!"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        
9999    CONTINUE
        
        IF(ASSOCIATED(bcdef)) THEN
           DEALLOCATE(bcdef)
        END IF
        
        IF(ASSOCIATED(drag)) THEN
           DEALLOCATE(drag)
        END IF
        
        RETURN 
        
      END SUBROUTINE io_write_boundary
      
