  f0(1:val_dim) = f0(1:val_dim) + W / d_value * f_value(1:val_dim)

  do i = 1, num_dim
     f1(1:val_dim, i) = f1(1:val_dim, i) + W / d_value * f_value(1:val_dim) * r(i)
  end do

  r0 = r0 + W / d_value

  r1(1:num_dim) = r1(1:num_dim) + W / d_value * r(1:num_dim)

  do i = 1, num_dim
     r2(i,1:num_dim) = r2(i,1:num_dim) + W / d_value * r(i) * r(1:num_dim)
  end do
