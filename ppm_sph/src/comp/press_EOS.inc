!------------------------------------------------
!< Subroutine that calculates pressure from density
!-----------------------------------------------
subroutine press_EOS

  use mod_data_global
  use mod_data_physics
  use mod_data_prtl
  use mod_data_sph

  implicit none
  !-----------------------------
  !< Arguments
  !-----------------------------

  !-----------------------------
  !< Local variables
  !-----------------------------
  integer                     :: i
  integer                     :: phasei
  real(mk)                    :: dd

  !-----------------------------
  !< Pressure "Equation of State"
  !-----------------------------
  DO i = 1,Npart
    phasei = ap(phaseidx,i)

    SELECT CASE (phase(abs(phasei))%type)

      ! Fluid
      CASE  (param_type_fluid)
        dd = 1._mk/(phase(phasei)%density)
        pdata(pressidx,i) = phase(phasei)%press_ref*((pdata(rhoidx,i)*dd)**phase(phasei)%gamma - 1._mk) + &
        phase(phasei)%press_b

      ! Boundary particles
      CASE default
        pdata(pressidx,i) = 0._MK + phase(phasei)%press_b

    END SELECT
  ENDDO


end subroutine press_EOS
