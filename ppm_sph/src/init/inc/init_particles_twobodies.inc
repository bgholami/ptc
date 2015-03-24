#if __DIM==2
DO j = 1,Np(2)
  y = min_compbox(2) + ((REAL(j,MK)-0.5_MK)*dpx(2))
  DO i=1,Np(1)
     x = min_compbox(1) + ((REAL(i,MK)-0.5_MK)*dpx(1))
     !-------------------------------------------------
     !  Initialize particles in cube
     !-------------------------------------------------
     if ( abs(x - body1_pos(1)) .le. 0.5_mk*body1_size(1) .and. &
          abs(y - body1_pos(2)) .le. 0.5_mk*body1_size(2)) then

       Npart = Npart + 1

       IF (MOD(Npart,Nproc) .EQ. rank) THEN
         Npart_local = Npart_local + 1

         ! Position
         xp(1,Npart_local) = x
         xp(2,Npart_local) = y

         ! Phase ID
         ap(phaseidx,Npart_local) = body1_phase

         ! Particle ID
         ap(partIDidx,Npart_local) = Npart

         ! Global number of particles
         Ntotal(ap(phaseidx,Npart_local)) = Ntotal(ap(phaseidx,Npart_local)) + 1
       ENDIF
     elseif ( abs(x - body2_pos(1)) .le. 0.5_mk*body2_size(1) .and. &
              abs(y - body2_pos(2)) .le. 0.5_mk*body2_size(2)) then

       Npart = Npart + 1

       IF (MOD(Npart,Nproc) .EQ. rank) THEN
         Npart_local = Npart_local + 1

         ! Position
         xp(1,Npart_local) = x
         xp(2,Npart_local) = y

         ! Phase ID
         ap(phaseidx,Npart_local) = body2_phase

         ! Particle ID
         ap(partIDidx,Npart_local) = Npart

         ! Global number of particles
         Ntotal(ap(phaseidx,Npart_local)) = Ntotal(ap(phaseidx,Npart_local)) + 1
       ENDIF
     endif
   ENDDO
ENDDO
#elif __DIM==3  
!The particle volumes are constant
DO k=1,Np(3)
 z = min_compbox(3) + ((REAL(k,MK)-0.5_MK)*dpx(3))
 DO j=1,Np(2)
    y = min_compbox(2) + ((REAL(j,MK)-0.5_MK)*dpx(2))
    DO i=1,Np(1)
       x = min_compbox(1) + ((REAL(i,MK)-0.5_MK)*dpx(1))
       !-------------------------------------------------
       !  Initialize particles in cube
       !-------------------------------------------------
       if ( abs(x - body1_pos(1)) .le. 0.5_mk*body1_size(1) .and. &
            abs(y - body1_pos(2)) .le. 0.5_mk*body1_size(2) .and. &
            abs(z - body1_pos(3)) .le. 0.5_mk*body1_size(3)) then

         Npart = Npart + 1

         IF (MOD(Npart,Nproc) .EQ. rank) THEN
           Npart_local = Npart_local + 1

           ! Position
           xp(1,Npart_local) = x
           xp(2,Npart_local) = y
           xp(3,Npart_local) = z

           ! Phase ID
           ap(phaseidx,Npart_local) = body1_phase

           ! Particle ID
           ap(partIDidx,Npart_local) = Npart

           ! Global number of particles
           Ntotal(ap(phaseidx,Npart_local)) = Ntotal(ap(phaseidx,Npart_local)) + 1
         ENDIF
       elseif ( abs(x - body2_pos(1)) .le. 0.5_mk*body2_size(1) .and. &
                abs(y - body2_pos(2)) .le. 0.5_mk*body2_size(2) .and. &
                abs(z - body2_pos(3)) .le. 0.5_mk*body2_size(3)) then

         Npart = Npart + 1

         IF (MOD(Npart,Nproc) .EQ. rank) THEN
           Npart_local = Npart_local + 1

           ! Position
           xp(1,Npart_local) = x
           xp(2,Npart_local) = y
           xp(3,Npart_local) = z

           ! Phase ID
           ap(phaseidx,Npart_local) = body2_phase

           ! Particle ID
           ap(partIDidx,Npart_local) = Npart

           ! Global number of particles
           Ntotal(ap(phaseidx,Npart_local)) = Ntotal(ap(phaseidx,Npart_local)) + 1
         ENDIF
       endif
    ENDDO
 ENDDO
ENDDO

#endif
