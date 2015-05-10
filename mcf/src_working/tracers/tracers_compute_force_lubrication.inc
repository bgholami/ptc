  !-------------------------------------------------
  ! Calculates the wall lubrication force on tracers
  !-------------------------------------------------
  
  f_lub = &
       sgn_ni(1:num_dim) * ABS(u_ni(1:num_dim) - v_ni(1:num_dim))  * &
       (1.10_MK / (1.0_MK/kappa - 1.0_MK)) + &
       sgn_ti(1:num_dim) * ABS(u_ti(1:num_dim) - v_ti(1:num_dim))  * &
       (0.70_MK * kappa)
  
  f_lub = f_lub / tau_p
  
  
  this%f(1:num_dim, ip) = this%f(1:num_dim, ip) + f_lub
