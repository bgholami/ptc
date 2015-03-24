      SUBROUTINE colloid_compute_repulsion_cw(this,&
           x,sid,F,FB,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_repulsion_cw
        !----------------------------------------------------
        !
        ! Purpose     : The gap between near contacting
        !               wall and colloid can be too small
        !               and tend to overlap.
        !               Compute a extra repulsive force
        !               which prevents them to overlap
        !               or become too close.
        !
        ! Routines    :
        !
        ! References  : Ball and Melrose, 
        !               Adv. Colloid Interface Sci.
        !               59 (1995) 19 -30.
        !              
        !
        ! Remarks     : V0.2 19.11.2010, one pair version.
        !
        !               V0.1 5.11 2010, original version,
        !               for all pairs.
        !
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
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:),INTENT(IN)       :: x
        INTEGER, INTENT(IN)                     :: sid
        REAL(MK), DIMENSION(:),INTENT(OUT)      :: F
        REAL(MK), DIMENSION(:,:),INTENT(OUT)    :: FB
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !
        ! s_drag   : drag on this process
        ! t_drag   : total drag on all processes
        ! s_torque : torque on this process
        ! t_torque : total torque on all processes        
        !----------------------------------------------------
        
        INTEGER                                  :: dim,num
        INTEGER                                  :: i,j
        REAL(MK)                                 :: hn,hm,h1,h2
        REAL(MK)                                 :: F0,Ft
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        
        dim = this%num_dim
        num = dim * 2
        
        IF ( SIZE(F,1) /= dim ) THEN
           PRINT *, "colloid_compute_repulsion_cw : ", &
                "input F dimension does not match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF ( SIZE(FB,1) /= dim ) THEN
           PRINT *, "colloid_compute_repulsion_cw : ", &
                "input FB dimension does not match !"
           stat_info = -1
           GOTO 9999
        END IF

        IF ( SIZE(FB,2) /= num ) THEN
           PRINT *, "colloid_compute_repulsion_cw : ", &
                "input FB number does not match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        
        F(:)    = 0.0_MK
        FB(:,:) = 0.0_MK
        
        hn = this%cw_repul_cut_off
        hm = this%cw_repul_cut_on
        F0 = this%cw_repul_F0
         
        IF ( dim == 2  ) THEN
           
           !-------------------------------------------------
           ! Loop each dimension of solid wall.
           !-------------------------------------------------
           
           DO j = 1, dim
              
              IF ( this%bcdef(2*j-1) == mcf_bcdef_wall_solid ) THEN
                 
                 !-------------------------------------------
                 ! Loop each colloid
                 !-------------------------------------------
                 
                 h1 = x(j)-this%min_phys(j)-this%radius(1,sid)
                 h2 = this%max_phys(j)-x(j)-this%radius(1,sid)
                    
                 IF ( h1 < hn ) THEN
                    
                    !----------------------------------------
                    ! If gap smaller than minimal allowed gap,
                    ! set it to minimum.
                    !----------------------------------------
                    
                    IF ( h1 < hm ) THEN
                       
                       h1 = hm
                       
                    END IF
                    
                    SELECT CASE (this%cw_repul_type)
                       
                    CASE (mcf_cw_repul_type_Hookean)
                       
                       Ft = F0 - F0*h1/hn
                       
                       
                    CASE (mcf_cw_repul_type_DLVO)
                       
                       Ft = F0 * EXP(-h1/hn) /(1.0_MK-EXP(-h1/hn))
                       
                    END SELECT
                    
                    F(j) = F(j) +  Ft
                    
                    !----------------------------------------
                    ! Collect force on boundary
                    !----------------------------------------
                    
                    FB(j,2*j-1) =  -Ft
                    
                 END IF ! h1 < hn
                 
                 IF ( h2 < hn ) THEN
                    
                    !----------------------------------------
                    ! If gap smaller than minimal allowed gap,
                    ! set it to minimum.
                    !----------------------------------------
                    
                    IF ( h2 < hm ) THEN
                       
                       h2 = hm
                       
                    END IF
                    
                    SELECT CASE (this%cw_repul_type)
                       
                    CASE (mcf_cw_repul_type_Hookean)
                       
                       Ft = F0 - F0*h2/hn
                       
                    CASE (mcf_cw_repul_type_DLVO)
                       
                       Ft = F0 * EXP(-h2/hn) /(1.0_MK-EXP(-h2/hn))
                       
                    END SELECT
                    
                    F(j) = F(j) - Ft 
                    
                    !----------------------------------------
                    ! Collect force on boundary
                    !----------------------------------------
                    
                    FB(j,2*j) = Ft
                    
                 END IF ! h2 < hn
                 
              END IF ! bcdef
              
           END DO ! i = 1 , dim
           
        END IF ! dim = 2
        
        
9999    CONTINUE
        
        
        RETURN          
        
      END SUBROUTINE colloid_compute_repulsion_cw
      
      
      
