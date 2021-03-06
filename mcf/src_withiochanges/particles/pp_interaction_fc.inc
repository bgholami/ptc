  !----------------------------------------------------------
  ! particle_particle_interaction_fluid_colloid
  !
  ! ip is fluid particle, 
  ! jp is colloid boundary particle. 
  !
  ! Substitute jp's velocity with an artificial
  ! velocity for no slip condition,
  ! pressure of jp with Neumann condition.
  !
  ! If Non-Newtonian fluids, 
  ! calculate velocity gradient tensor.
  ! pressure tensor of jp with Neumann condition.
  !
  !
  ! f_c interaction takes place the same way
  ! as f_f interaction.
  !----------------------------------------------------------
  
  
  v_ip(1:num_dim) = this%v(1:num_dim,ip)
  v_jp(1:num_dim) = this%v(1:num_dim,jp)
  
  
  !----------------------------------------------------------
  ! Morris no slip conditin for colloid.
  !----------------------------------------------------------
  
  IF ( coll_noslip == 2 ) THEN
     
     CALL colloid_noslip(colloids, &
          this%x(1:num_dim,ip),this%x(1:num_dim,jp), &
          this%v(1:num_dim,ip),v_jp(1:num_dim), &
          this%id(this%sid_idx,jp),&
          stat_info_sub)
     
     IF( stat_info_sub /=0 ) THEN
        PRINT *, "pp_interaction_fc : ",&
             "Colloid no slip jp has problem !"
        stat_info = -1
        GOTO 9999
     END IF
     
  END IF  ! coll_noslip
  
  !----------------------------------------------------------
  ! For Newtonian fluid, pressure using Neumann conditon.
  !----------------------------------------------------------

  IF ( Newtonian ) THEN
     
     this%p(jp) = this%p(ip)
     
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
     
     !-------------------------------------------------------
     ! For non-Newtonian fluid, density and pressure 
     ! using Neumann conditon.
     !-------------------------------------------------------
     
     rho_ip = this%rho(ip)
     rho_jp = this%rho(ip)
     
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
     
  END IF
  
  IF( stat_info_sub /=0 ) THEN
     PRINT *, "pp_inteaction_fc : ",& 
          "Pairwise force has problem !"
     stat_info = -1
     GOTO 9999
  END IF
  
  
