!----------------------------------------------------------------
! This contains the routines which calculate the background
! pressure force and returns two accelerations between two 
! particles.
!----------------------------------------------------------------
      SUBROUTINE rhs_force_ff_background_pressure(this,&
           xi, xj, dij, numi, numj, mi, mj, gradw, &
           fi, fj, stat_info)
        !----------------------------------------------------
        ! Subroutine : rhs_force_ff_background_pressure
        !----------------------------------------------------
        !
        ! Purpose    : 
        !
        ! Reference  : 
        !
        ! Remark     :
        !
        ! Revisions  : V0.1 09.10.2013
        !
        !----------------------------------------------------
        ! Author       : Babak Gholami
        ! Contact      : babak.gholami@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments
        ! 
        ! numi : number density of i
        ! numj : number density of j
        !----------------------------------------------------
        
        TYPE(Rhs), INTENT(INOUT)                :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xi
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xj
        REAL(MK), INTENT(IN)                    :: dij
        REAL(MK), INTENT(IN)                    :: numi
        REAL(MK), INTENT(IN)                    :: numj
        REAL(MK), INTENT(IN)                    :: mi
        REAL(MK), INTENT(IN)                    :: mj
        REAL(MK), INTENT(IN)                    :: gradw
        REAL(MK), DIMENSION(:), INTENT(INOUT)   :: fi
        REAL(MK), DIMENSION(:), INTENT(INOUT)   :: fj
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: num_dim
        REAL(MK), DIMENSION(3)          :: eij
        REAL(MK)                        :: f_bp  
        REAL(MK)                        :: bp
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0

        fi(:) = 0.0_MK
        fj(:) = 0.0_MK
        
        num_dim = this%num_dim

        bp = physics_get_bp(this%phys, stat_info_sub)
        
        !----------------------------------------------------
        ! Check if any denominator is non-positive.
        !----------------------------------------------------
        
        IF ( mi < mcf_machine_zero .OR. &
             mj < mcf_machine_zero .OR. &
             numi < mcf_machine_zero .OR. &
             numj < mcf_machine_zero ) THEN
           
           PRINT *,"xi,xj,dij : ",  xi,xj,dij
           PRINT *,"numi,numj : ", numi,numj
           PRINT *,"mi,mj : ", mi,mj
           PRINT *,"gradw : ", gradw          
           
           stat_info = -1
           GOTO 9999
           
        END IF
        
        !----------------------------------------------------
        ! Calculate normalized vector pointing from j to i.
        !----------------------------------------------------
        
        eij(1:num_dim) = (xi(1:num_dim)  - xj(1:num_dim)) / dij
        
        
        !----------------------------------------------------
	! Calculate the background pressure force
	! per unit mass.
  	!----------------------------------------------------
        
        f_bp = ( bp/(numi**2.0_MK)+ bp/(numj**2.0_MK)) * gradw
        
        fi(1:num_dim) = fi(1:num_dim) - f_bp * eij(1:num_dim) / mi
        fj(1:num_dim) = fj(1:num_dim) + f_bp * eij(1:num_dim) / mj
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE rhs_force_ff_background_pressure
      
      
