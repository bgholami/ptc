  !-------------------------------------------------
  ! Calculates the wall lift force on tracers
  !-------------------------------------------------
  
  !shear_rate = u_t / this%h_p(ip)
  ! force cap is enforced by controling kappa, so we write shear_rate like this.
  shear_rate = u_t * kappa  / a_p 

  u_s = v_t - u_t


  IF (ABS(shear_rate) <= mcf_machine_zero) THEN

     lift_function = tracers_eval_lift_function(kappa, &
          0.0_MK)

  ELSE 

     IF (ABS(u_s) <= mcf_machine_zero) THEN
        u_s = shear_rate * a_p
     END IF
     
     lift_function = tracers_eval_lift_function(kappa, &
          shear_rate*a_p/u_s)
     
  END IF

  
  ! notice that rho is fluid's density
  f_lift = &
       -rho * a_p**2.0_MK * this%n_vec(1:num_dim, ip) / m_p * &
       (u_s**2 * lift_function)
  
  this%f(1:num_dim, ip) = this%f(1:num_dim, ip) + f_lift
     

  
  
