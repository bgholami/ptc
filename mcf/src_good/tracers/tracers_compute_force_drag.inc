  !-------------------------------------------------
  ! Calculates the flow drag force on tracers
  !-------------------------------------------------
  
  this%f(1:num_dim,1:num_part_real) = &
       this%f(1:num_dim,1:num_part_real) + &
       (VI(1:num_dim,1:num_part_real) - &
       this%v(1:num_dim,1:num_part_real)) / tau_p  
