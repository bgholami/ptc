SUBROUTINE physics_compute_cs(this,stat_info)
  !----------------------------------------------------
  ! Subroutine  : physics_compute_cs
  !----------------------------------------------------
  !
  ! Purpose     : Compute the sound speed based
  !               on guidlines of Morris et al. JCP 1997.
  !      
  ! Reference   : Morris et al. JCP 1997.
  !
  ! Remark      :      
  !
  ! Revisions   : V0.1 10.10.2013, original version.
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

  TYPE(Physics), INTENT(INOUT)    :: this
  INTEGER, INTENT(OUT)            :: stat_info

  !----------------------------------------------------
  ! Local variables
  !----------------------------------------------------

  REAL(MK)                        :: cs_velo, cs_visc, cs_pres
  INTEGER                         :: stat_info_sub


  !----------------------------------------------------
  ! Initialization of variables.
  !----------------------------------------------------

  stat_info     = 0
  stat_info_sub = 0



  !----------------------------------------------------
  ! Set basic initial dts
  !----------------------------------------------------

  cs_velo  = -1.0_MK
  cs_visc  = -1.0_MK
  cs_pres  = -1.0_MK
  IF ( this%c   <= 0.0_MK ) THEN


     !----------------------------------------------------
     ! find values
     !----------------------------------------------------

     cs_velo = sqrt(this%vel_ref**2 / this%rho_var)
     cs_visc = sqrt(this%vel_ref * this%eta / this%rho / &
          this%len_ref / this%rho_var)
     cs_pres = sqrt(this%fa_max * this%len_ref / this%rho_var)



     !----------------------------------------------------
     ! assign largest as sound speed
     !----------------------------------------------------     

     IF ( this%c < cs_velo )  THEN

        this%c = cs_velo

     END IF

     IF ( this%c < cs_visc )  THEN

        this%c = cs_visc

     END IF

     IF ( this%c < cs_pres )  THEN

        this%c = cs_pres

     END IF

  END IF


  !----------------------------------------------------
  ! same sound speed for the relax run
  !----------------------------------------------------
  IF ( this%c_relax <= 0.0_MK ) THEN
     this%c_relax = this%c
  END IF



9999 CONTINUE

  RETURN

END SUBROUTINE physics_compute_cs


