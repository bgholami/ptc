      SUBROUTINE particles_finalize(this,stat_info)

        TYPE(Particles),INTENT(INOUT)      :: this
        INTEGER,INTENT(OUT)             :: stat_info
        
        !-----------------------
        ! Local variables.
        !------------------------
        stat_info = 0
        
        IF(ASSOCIATED(this%x)) THEN
           DEALLOCATE(this%x)
        END IF
        
        IF(ASSOCIATED(this%v)) THEN
           DEALLOCATE(this%v)
        END IF
        
        IF(ASSOCIATED(this%rho)) THEN
           DEALLOCATE(this%rho)
        END IF

        IF(ASSOCIATED(this%rho_norm)) THEN
           DEALLOCATE(this%rho_norm)
        END IF
        
        IF(ASSOCIATED(this%m)) THEN
           DEALLOCATE(this%m)
        END IF
        
        IF(ASSOCIATED(this%p)) THEN
           DEALLOCATE(this%p)
        END IF
        
        IF(ASSOCIATED(this%id)) THEN
           DEALLOCATE(this%id)
        END IF
        
        IF(ASSOCIATED(this%f)) THEN
           DEALLOCATE(this%f)
        END IF   

        IF(ASSOCIATED(this%f_bp)) THEN
           DEALLOCATE(this%f_bp)
        END IF
        
        IF(ASSOCIATED(this%u)) THEN
           DEALLOCATE(this%u)
        END IF
        
        IF(ASSOCIATED(this%au)) THEN
           DEALLOCATE(this%au)
        END IF
        
        IF(ASSOCIATED(this%part_sym_list)) THEN
           DEALLOCATE(this%part_sym_list)
        END IF
       
        IF(ASSOCIATED(this%part_wall_sym_list)) THEN
           DEALLOCATE(this%part_wall_sym_list)
        END IF
        
        IF(ASSOCIATED(this%part_wall_solid_real_list)) THEN
           DEALLOCATE(this%part_wall_solid_real_list)
        END IF

        IF(ASSOCIATED(this%part_wall_solid_ghost_list)) THEN
           DEALLOCATE(this%part_wall_solid_ghost_list)
        END IF
        
        IF(ASSOCIATED(this%part_le_list)) THEN
           DEALLOCATE(this%part_le_list)
        END IF
        
        IF(ASSOCIATED(this%part_colloid_list)) THEN
           DEALLOCATE(this%part_colloid_list)
        END IF
        
        IF(ASSOCIATED(this%vgt)) THEN
           DEALLOCATE(this%vgt)
        END IF
        
        IF(ASSOCIATED(this%evgt)) THEN
           DEALLOCATE(this%evgt)
        END IF
        
        IF(ASSOCIATED(this%eval)) THEN
           DEALLOCATE(this%eval)
        END IF
        
        IF(ASSOCIATED(this%aeval)) THEN
           DEALLOCATE(this%aeval)
        END IF
        
        IF(ASSOCIATED(this%evec)) THEN
           DEALLOCATE(this%evec)
        END IF
        
        IF(ASSOCIATED(this%aevec)) THEN
           DEALLOCATE(this%aevec)
        END IF
        
        IF(ASSOCIATED(this%ct)) THEN
           DEALLOCATE(this%ct)
        END IF
        
        IF(ASSOCIATED(this%act)) THEN
           DEALLOCATE(this%act)
        END IF
        
        IF(ASSOCIATED(this%pt)) THEN
           DEALLOCATE(this%pt)
        END IF

#if __TRACER  

        IF(ASSOCIATED(this%x_interface)) THEN
           DEALLOCATE(this%x_interface)
        END IF   

        IF(ASSOCIATED(this%vo_interface)) THEN
           DEALLOCATE(this%vo_interface)
        END IF 

        IF(ASSOCIATED(this%vc_interface)) THEN
           DEALLOCATE(this%vc_interface)
        END IF   

        IF(ASSOCIATED(this%vn_interface)) THEN
           DEALLOCATE(this%vn_interface)
        END IF   

        IF(ASSOCIATED(this%num_tracer_deposited)) THEN
           DEALLOCATE(this%num_tracer_deposited)
        END IF 

        IF(ASSOCIATED(this%c_tracers)) THEN
           DEALLOCATE(this%c_tracers)
        END IF 

        IF(ASSOCIATED(this%dn)) THEN
           DEALLOCATE(this%dn)
        END IF
        
#endif
        
        
        PRINT *, "particles_finalize : ", "Finished!"
        
        RETURN          
        
      END SUBROUTINE particles_finalize
      
