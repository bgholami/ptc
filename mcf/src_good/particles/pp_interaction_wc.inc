  !----------------------------------------------------------
  ! particle-particle_interaction_wall_colloid
  !
  ! ip is wall boundary particle, 
  ! jp is colloidal boundary particle. 
  !
  ! ip and jp's pressure are from equation of state.
  !
  ! If fluid is non-Newtonian,
  ! not implemented yet.
  !
  !
  ! w_c interaction takes place the same way
  ! as f_f interaction.
  !----------------------------------------------------------

  
  v_ip(1:num_dim) = this%v(1:num_dim, ip)
  v_jp(1:num_dim) = this%v(1:num_dim, jp)
  
  
  IF ( Newtonian ) THEN
     
     CALL rhs_force_ff(this%rhs,&
          this%x(1:num_dim,ip),&
          this%x(1:num_dim,jp),dij,&
          v_ip(1:num_dim),&
          v_jp(1:num_dim),&
          this%rho(ip),this%rho(jp),&
          this%p(ip),this%p(jp),&
          this%m(ip),this%m(jp),&
          w,gradw,fip,fjp,fbpip,fbpjp,stat_info=stat_info_sub)
     
  ELSE ! non-Newtonian
     
     rho_ip = this%rho(ip)
     rho_jp = this%rho(jp)
     
#include "pp_vgt_ip.inc"
     
     IF ( symmetry ) THEN
        
#include "pp_vgt_jp.inc"
        
     END IF
     
     this%pt(1:num_dim,1:num_dim,jp) = &
          this%pt(1:num_dim,1:num_dim,ip)
     
     CALL rhs_force_ff(this%rhs, &
          this%x(1:num_dim,ip), &
          this%x(1:num_dim,jp),dij, &
          v_ip(1:num_dim), &
          v_jp(1:num_dim), &
          rho_ip,rho_jp, &
          this%pt(1:num_dim,1:num_dim,ip), &
          this%pt(1:num_dim,1:num_dim,jp), &
          this%m(ip),this%m(jp), &
          w,gradw,fip,fjp,fbpip,fbpjp,stat_info=stat_info_sub)
     
  END IF ! Newtonian
  
  IF( stat_info_sub /=0 ) THEN
     PRINT *, "pp_inteaction_wc : ",& 
          "Pairwise force has problem !"
     stat_info = -1
     GOTO 9999
  END IF
