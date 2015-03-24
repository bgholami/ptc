     !-------------------------------------------------------------------------
     !  Module   :                  ppm_neighlist
     !-------------------------------------------------------------------------
     ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
     !                    Center for Fluid Dynamics (DTU)
     !
     !
     ! This file is part of the Parallel Particle Mesh Library (PPM).
     !
     ! PPM is free software: you can redistribute it and/or modify
     ! it under the terms of the GNU Lesser General Public License
     ! as published by the Free Software Foundation, either
     ! version 3 of the License, or (at your option) any later
     ! version.
     !
     ! PPM is distributed in the hope that it will be useful,
     ! but WITHOUT ANY WARRANTY; without even the implied warranty of
     ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
     ! GNU General Public License for more details.
     !
     ! You should have received a copy of the GNU General Public License
     ! and the GNU Lesser General Public License along with PPM. If not,
     ! see <http://www.gnu.org/licenses/>.
     !
     ! Parallel Particle Mesh Library (PPM)
     ! ETH Zurich
     ! CH-8092 Zurich, Switzerland
     !-------------------------------------------------------------------------
#if   __ACTION == __COUNT
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE count_neigh_s(p_idx, clist, domain, xp, cutoff, skin, nvlist)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE count_neigh_d(p_idx, clist, domain, xp, cutoff, skin, nvlist)
#endif
#elif __ACTION == __GET
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE get_neigh_s(p_idx, clist, domain, xp, cutoff, skin, vlist, nvlist)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE get_neigh_d(p_idx, clist, domain, xp, cutoff, skin, vlist, nvlist)
#endif
#endif
      !!! Given the particle index, this subroutine locates the cell that this
      !!! cell is located in, gathers all particles in these cells and updates
      !!! verlet lists of all these particles.
          IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
          INTEGER,  INTENT(IN)                 :: p_idx
          !!! Particle index
          TYPE(ppm_clist)  , INTENT(IN)        :: clist
          !!! cell list
          REAL(MK), DIMENSION(2*ppm_dim)       :: domain
          !!! Pysical extent of whole domain including ghost layers
          REAL(MK), INTENT(IN), DIMENSION(:,:) :: xp
          !!! Particle coordinates
          REAL(MK), INTENT(IN), DIMENSION(:)   :: cutoff
          !!! Particle cutoff radii
          REAL(MK), INTENT(IN)                 :: skin
          !!! Skin parameter
#if __ACTION == __GET
          INTEGER,  DIMENSION(:,:)             :: vlist
          !!! Verlet list, where vlist(j, i) contains the jth neighbor of
          !!! ith particle
#endif
          INTEGER,  DIMENSION(:)               :: nvlist
          !!! Number of neighbors of particles. nvlist(i) contains number of
          !!! neighbors particle i has.
          INTEGER                              :: info

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
          INTEGER                              :: p_ref
          ! Reference particle
          INTEGER                              :: p_neigh
          ! Particle in the neighbor cell
          REAL(MK),    DIMENSION(ppm_dim)      :: offset_coor
          ! Offset from midpoint coordinates of reference array to get to
          ! neighbor cells.
          REAL(MK),    DIMENSION(ppm_dim)      :: c_coor
          ! Coordinates
          REAL(MK),    DIMENSION(ppm_dim)      :: n_coor
          ! Coordinates of neighbor cell
          INTEGER(ppm_kind_int64)              :: c_idx
          ! Cell index
          INTEGER(ppm_kind_int64)              :: n_idx
          ! Cell index
          INTEGER(ppm_kind_int64)              :: parentIdx
          ! Index of parent cell
          REAL(MK),    DIMENSION(ppm_dim)      :: p_coor
          ! Particle coordinates
          INTEGER                              :: p_depth
          ! Depth of the particle
          INTEGER                              :: c_depth
          ! Depth of cell

      !-------------------------------------------------------------------------
      !  Counters
      !-------------------------------------------------------------------------
          INTEGER                              :: i
          INTEGER                              :: j
          INTEGER                              :: m
          INTEGER                              :: n

          !---------------------------------------------------------------------
          !  If this particle was visited before, skip it.
          !---------------------------------------------------------------------
          IF(used(p_idx))   RETURN

          !---------------------------------------------------------------------
          !  Initialize position on empty list to 0.
          !---------------------------------------------------------------------
          empty_pos = 0

          !---------------------------------------------------------------------
          !  Get coordinates and depth of the particle
          !---------------------------------------------------------------------
          CALL getParticleCoorDepth(p_idx, domain, p_coor, p_depth, xp, cutoff, skin)

          !---------------------------------------------------------------------
          !  Get index of the cell that this particle is located in.
          !---------------------------------------------------------------------
          c_idx   = getCellIdx(p_coor, p_depth, domain)

          !---------------------------------------------------------------------
          !  Get coordinates and depth of the particle
          !---------------------------------------------------------------------
          CALL getCellCoor_Depth(c_idx, domain, c_coor, c_depth, & 
 &                               clist%max_depth, info)

          !---------------------------------------------------------------------
          !  Compute offset coordinates, which will be used to find neighbor cells
          !---------------------------------------------------------------------
          offset_coor(1) = (domain(2) - domain(1))/(2**(c_depth-1))
          offset_coor(2) = (domain(4) - domain(3))/(2**(c_depth-1))
          DO i = 3,ppm_dim
              offset_coor(i) = (domain(2*i) - domain(2*i-1))/(2**(c_depth-1))
          END DO

          !---------------------------------------------------------------------
          !  From first cell until center cell, check whether they are empty
          !  or not and store those in empty list that contain no particles in
          !  their region.
          !---------------------------------------------------------------------
          ! From first cell to (center-1)th cell
          DO i = 1, ((3**ppm_dim - 1)/2)
              ! Compute midpoints of neighbor cell
              n_coor(1) = c_coor(1) + (ncells(i,1)*offset_coor(1))
              n_coor(2) = c_coor(2) + (ncells(i,2)*offset_coor(2))
              DO j = 3, ppm_dim
                  n_coor(j) = c_coor(j) + (ncells(i,j)*offset_coor(j))
              END DO

              ! If found coordinates are not inside domain, skip!
              IF(.NOT. inDomain(n_coor, domain))   cycle

              ! Get index of the neighbor cell
              n_idx = getCellIdx(n_coor, c_depth, domain)

              ! Store in empty list if empty
              IF(isEmpty(n_idx,clist%lookup))     CALL putInEmptyList(n_idx)
          END DO

          !---------------------------------------------------------------------
          !  Get particles that are in the same cell with input particle.
          !---------------------------------------------------------------------
          CALL getParticlesInCell(c_idx, xp, clist, own_plist, own_nplist)

          !---------------------------------------------------------------------
          !  At the center cell, every particle in the same cell is compared
          !  with another once. If they are neighbors, nvlist entry of those
          !  particles that are real are incremented.
          !---------------------------------------------------------------------
          ! For each particle in the same cell
          DO i = 1, own_nplist
              ! Pick a reference particle
              p_ref = own_plist(i)
              ! For each particle up to the reference particle
              DO j = 1, (i - 1)
                  ! Pick a candidate for neighbor particle
                  p_neigh = own_plist(j)
                  ! If they are neighbors and ...
                  IF(isNeighbor(p_ref, p_neigh, xp, cutoff, skin))  THEN
                      ! If the reference particle is a real particle ...
                      IF(p_ref .LE. clist%n_real_p)   THEN
                          ! Store neighbor particle in verlet list of reference
                          ! particle
                          nvlist(p_ref) = nvlist(p_ref) + 1
#if __ACTION == __GET
                          vlist(nvlist(p_ref), p_ref) = p_neigh
#endif
                      ENDIF

                      !If the neighbor particle is a real particle ...
                      IF(p_neigh .LE. clist%n_real_p)   THEN
                          ! Store reference particle in verlet list of neighbor
                          ! particle
                          nvlist(p_neigh) = nvlist(p_neigh) + 1
#if __ACTION == __GET
                          vlist(nvlist(p_neigh), p_neigh) = p_ref
#endif
                      ENDIF
                  END IF
              END DO
          END DO

          !---------------------------------------------------------------------
          !  Starting from the cell after the center cell upto last neighbor cell,
          !  get particles from neighbor particles and compare reference particles
          !  with those from neighboring cells. If two particles are found,
          !  increment nvlist entry of those particles that are real.
          !---------------------------------------------------------------------
          ! For each neighbor cell
          DO i = ((3**ppm_dim + 1)/2) + 1, 3**ppm_dim
              ! Compute midpoint coordinates of the neighbor cell
              n_coor(1) = c_coor(1) + (ncells(i,1)*offset_coor(1))
              n_coor(2) = c_coor(2) + (ncells(i,2)*offset_coor(2))
              DO j = 3, ppm_dim
                  n_coor(j) = c_coor(j) + (ncells(i,j)*offset_coor(j))
              END DO

              ! If computed coordinates are not inside the domain, skip!
              IF(.NOT. inDomain(n_coor, domain))   cycle

              ! Get index of the neighbor cell.
              n_idx = getCellIdx(n_coor, c_depth, domain)

              ! Get particles in the neighbor cell
              CALL getParticlesInCell(n_idx, xp, clist, neigh_plist, neigh_nplist)

              ! For each particle in the reference cell
              DO m = 1, own_nplist
                  ! Pick a reference particle
                  p_ref = own_plist(m)
                  ! For each particle in the neighbor cell
                  DO n = 1, neigh_nplist
                      ! Pick a candidate for neighbor particle
                      p_neigh = neigh_plist(n)
                      ! If particles are neighbors and ...
                      IF(isNeighbor(p_ref, p_neigh, xp, cutoff, skin)) THEN
                          ! If reference particle is a real particle ...
                          IF(p_ref .LE. clist%n_real_p)   THEN
                              ! Store neighbor particle in verlet list of
                              ! reference particle
                              nvlist(p_ref)   = nvlist(p_ref)   + 1
#if __ACTION == __GET
                              vlist(nvlist(p_ref), p_ref) = p_neigh
#endif
                          END IF

                          ! If neighbor particle is a real particle ...
                          IF(p_neigh .LE. clist%n_real_p) THEN
                              ! Store reference particle in verlet list of
                              ! neighbor particle
                              nvlist(p_neigh) = nvlist(p_neigh) + 1
#if __ACTION == __GET
                              vlist(nvlist(p_neigh), p_neigh) = p_ref
#endif
                          END IF
                      END IF
                  END DO
              END DO
          END DO

          !---------------------------------------------------------------------
          !  For each reference particle, get coordinates and depth and increment
          !  depth by 1 in each iteration such that each reference particle finds
          !  its neighbors that have smaller cutoff radius, hence they are located
          !  in deeper cells. If the cell that is taken as reference cell does
          !  not contain any cells in its region, then this cell is added to
          !  empty list such that children cells will be skipped since they are
          !  also empty.
          !---------------------------------------------------------------------
          ! For each reference particle
          DO m = 1, own_nplist
              ! Pick a reference particle
              p_ref = own_plist(m)
              ! Get coordinates and depth of this particle
              CALL getParticleCoorDepth(p_ref, domain, p_coor, p_depth, xp, cutoff, skin)

              ! Increment the depth of the particle by 1 to look down at deeper cells
              p_depth = p_depth + 1

              ! Until we reach the maximum depth ...
              DO WHILE(p_depth .LE. clist%max_depth)
                  ! Get the child cell this particle is located in
                  c_idx = getCellIdx(p_coor, p_depth, domain)
                  ! Get index of parent cell to check whether parent is empty or not
                  parentIdx = parent(c_idx)

                  ! If the parent is empty, also put the child in the empty list.
                  IF (inEmptyList(parentIdx)) THEN
                     CALL putInEmptyList(c_idx)
                  ENDIF
                  ! Get cell coordinates and depth
                  CALL getCellCoor_Depth(c_idx, domain, c_coor, c_depth, &
 &                                       clist%max_depth, info)

                  ! Compute offset coordinates to compute midpoint coordinates
                  ! of neighbor cells at this depth.
                  offset_coor(1) = (domain(2) - domain(1))/(2**(c_depth-1))
                  offset_coor(2) = (domain(4) - domain(3))/(2**(c_depth-1))
                  DO i = 3, ppm_dim
                      offset_coor(i) = (domain(2*i) - domain(2*i-1))/&
 &                                     (2**(c_depth-1))
                  END DO

                  ! For each neighbor
                  DO i = 1, 3**ppm_dim
                      ! Compute midpoint coordinates of the neighbor
                      n_coor(1) = c_coor(1) + (ncells(i,1)*offset_coor(1))
                      n_coor(2) = c_coor(2) + (ncells(i,2)*offset_coor(2))
                      DO j = 3,ppm_dim
                          n_coor(j) = c_coor(j) + (ncells(i,j)*offset_coor(j))
                      END DO

                      ! If coordinates are not inside the domain, skip!
                      IF(.NOT. inDomain(n_coor, domain))   cycle

                      ! Get index of the neighbor cell
                      n_idx = getCellIdx(n_coor, c_depth, domain)
                      
                      ! if parent of neighbor is empty,
                      ! add to emptylist and skip
                      IF (inEmptyList(parent(n_idx))) THEN
                          CALL putInEmptyList(n_idx)
                          CYCLE
                      ENDIF

                      ! Get particles in the neighbor cell
                      !-----------------------------------------------------
                      ! REMARK: If this cell has no particles, it will be stored
                      !         in empty cell inside getParticlesInCell
                      !         subroutine. If it has particles only at this
                      !         depth, but no deeper particles, it will again
                      !         be saved in the empty list just  before returning.
                      !-----------------------------------------------------
                      CALL getParticlesInCell(n_idx, xp, clist, neigh_plist, neigh_nplist)

                      ! For each neighbor element
                      DO n = 1, neigh_nplist
                          ! Pick a candidate for neighbor particle
                          p_neigh = neigh_plist(n)
                          ! If they are neighbors and ...
                          IF(isNeighbor(p_ref, p_neigh, xp, cutoff, skin))  THEN
                              ! If reference particle is real ...
                              IF(p_ref .LE. clist%n_real_p) THEN
                                  ! Store neighbor particle in verlet list
                                  ! of reference particle
                                  nvlist(p_ref) = nvlist(p_ref) + 1
#if __ACTION == __GET
                                  vlist(nvlist(p_ref), p_ref) = p_neigh
#endif
                              END IF

                              ! If neighbor particle is real ...
                              IF(p_neigh .LE. clist%n_real_p) THEN
                                  ! Store reference particle in verlet list
                                  ! of neighbor particle
                                  nvlist(p_neigh) = nvlist(p_neigh) + 1
#if __ACTION == __GET
                                  vlist(nvlist(p_neigh), p_neigh) = p_ref
#endif
                              END IF
                          END IF
                      END DO
                 END DO
                 ! To search for deeper particles, increment depth of particle by 1.
                 p_depth = p_depth + 1
              END DO
          END DO

          ! Label all reference particles as USED so that they wont be visited again.
          DO i = 1, own_nplist
              used(own_plist(i)) = .TRUE.
          END DO
#if   __ACTION == __COUNT
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE count_neigh_s
#elif   __KIND == __DOUBLE_PRECISION
      END SUBROUTINE count_neigh_d
#endif
#elif __ACTION == __GET
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE get_neigh_s
#elif   __KIND == __DOUBLE_PRECISION
      END SUBROUTINE get_neigh_d
#endif
#endif

#if  __ACTION == __COUNT
#if   __KIND == __SINGLE_PRECISION
          SUBROUTINE count_neigh_sym_s(p_idx, clist, domain, actual_domain,   &
          &xp, cutoff, skin, nvlist)
#elif   __KIND == __DOUBLE_PRECISION
          SUBROUTINE count_neigh_sym_d(p_idx, clist, domain, actual_domain,   &
          &xp, cutoff, skin, nvlist)
#endif
#elif __ACTION == __GET
#if   __KIND == __SINGLE_PRECISION
          SUBROUTINE get_neigh_sym_s(p_idx, clist, domain, actual_domain,   &
          &xp, cutoff, skin, vlist, nvlist)
#elif   __KIND == __DOUBLE_PRECISION
          SUBROUTINE get_neigh_sym_d(p_idx, clist, domain, actual_domain,   &
          &xp, cutoff, skin, vlist, nvlist)
#endif
#endif
          IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
          INTEGER,  INTENT(IN)                 :: p_idx
          !!! Particle index
          TYPE(ppm_clist)  , INTENT(IN)        :: clist
          !!! cell list
          REAL(MK), DIMENSION(2*ppm_dim)       :: domain
          !!! Pysical extent of whole domain including ghost layers
          REAL(MK)   , DIMENSION(2*ppm_dim)    :: actual_domain
          !!! Physical extent of actual domain excluding ghost layers
          REAL(MK), INTENT(IN), DIMENSION(:,:) :: xp
          !!! Particle coordinates
          REAL(MK), INTENT(IN), DIMENSION(:)   :: cutoff
          !!! Particle cutoff radii
          REAL(MK), INTENT(IN)                 :: skin
          !!! Skin parameter
#if __ACTION == __GET
          INTEGER,  DIMENSION(:,:)             :: vlist
          !!! Verlet list, where vlist(j, i) contains the jth neighbor of
          !!! ith particle
#endif
          INTEGER,  DIMENSION(:)               :: nvlist
          !!! Number of neighbors of particles. nvlist(i) contains number of
          !!! neighbors particle i has.
          INTEGER                              :: info
          !!! 0 on success

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
          INTEGER                              :: p_ref
          ! Reference particle
          INTEGER                              :: p_neigh
          ! Particle in the neighbor cell
          REAL(MK),    DIMENSION(ppm_dim)      :: offset_coor
          ! Offset from midpoint coordinates of reference array to get to
          ! neighbor cells.
          REAL(MK),    DIMENSION(ppm_dim)      :: c_coor
          ! Coordinates
          REAL(MK),    DIMENSION(ppm_dim)      :: n_coor
          ! Coordinates of neighbor cell
          INTEGER(ppm_kind_int64)              :: c_idx
          ! Cell index
          INTEGER(ppm_kind_int64)              :: n_idx
          ! Cell index
          INTEGER(ppm_kind_int64)              :: parentIdx
          ! Index of parent cell
          REAL(MK),    DIMENSION(ppm_dim)      :: p_coor
          ! Particle coordinates
          INTEGER                              :: p_depth
          ! Depth of the particle
          INTEGER                              :: c_depth
          ! Depth of cell

      !-------------------------------------------------------------------------
      !  Counters
      !-------------------------------------------------------------------------
          INTEGER                              :: i
          INTEGER                              :: j
          INTEGER                              :: m
          INTEGER                              :: n

          !---------------------------------------------------------------------
          !  If this particle was visited before, skip it.
          !---------------------------------------------------------------------
          IF(used(p_idx))   RETURN

          !---------------------------------------------------------------------
          !  Initialize position on empty list to 0.
          !---------------------------------------------------------------------
          empty_pos = 0

          !---------------------------------------------------------------------
          !  Get coordinates and depth of the particle
          !---------------------------------------------------------------------
          CALL getParticleCoorDepth(p_idx, domain, p_coor, p_depth, xp, cutoff, skin)

          !---------------------------------------------------------------------
          !  Get index of the cell that this particle is located in.
          !---------------------------------------------------------------------
          c_idx   = getCellIdx(p_coor, p_depth, domain)

          !---------------------------------------------------------------------
          !  Get coordinates and depth of the particle
          !---------------------------------------------------------------------
          CALL getCellCoor_Depth(c_idx, domain, c_coor, c_depth, &
 &                               clist%max_depth, info)

          !---------------------------------------------------------------------
          !  Compute offset coordinates, which will be used to find neighbor cells
          !---------------------------------------------------------------------
          offset_coor(1) = (domain(2) - domain(1))/(2**(c_depth-1))
          offset_coor(2) = (domain(4) - domain(3))/(2**(c_depth-1))
          DO i = 3,ppm_dim
              offset_coor(i) = (domain(2*i) - domain(2*i-1))/(2**(c_depth-1))
          END DO

          !---------------------------------------------------------------------
          !  From first cell until center cell, check whether they are empty
          !  or not and store those in empty list that contain no cells in
          !  their region.
          !---------------------------------------------------------------------
          ! From first cell to (center-1)th cell
          DO i = 1, ((3**ppm_dim - 1)/2)
              ! Compute midpoints of neighbor cell
              n_coor(1) = c_coor(1) + (ncells(i,1)*offset_coor(1))
              n_coor(2) = c_coor(2) + (ncells(i,2)*offset_coor(2))
              DO j = 3, ppm_dim
                  n_coor(j) = c_coor(j) + (ncells(i,j)*offset_coor(j))
              END DO

              ! If found coordinates are not inside domain, skip!
              IF(.NOT. inDomain(n_coor, domain))   cycle

              ! Get index of the neighbor cell
              n_idx = getCellIdx(n_coor, c_depth, domain)

              ! Store in empty list if empty
              IF(isEmpty(n_idx,clist%lookup))     CALL putInEmptyList(n_idx)
          END DO

          !---------------------------------------------------------------------
          !  Get particles that are in the same cell with input particle.
          !---------------------------------------------------------------------
          CALL getParticlesInCell(c_idx, xp, clist, own_plist, own_nplist)

          !---------------------------------------------------------------------
          !  At the center cell, every particle in the same cell is compared
          !  with another once. If they are neighbors and:
          !  1. If reference particle is real, nvlist entry of it is incremented
          !     by 1
          !  2. If reference particle is ghost, but neighbor particle is real,
          !     nvlist entry of neighbor particle is incremented by 1.
          !  3. If both particles are ghost and they are cross-neighbors,
          !     nvlist entry of reference particle is incremented by 1.
          !---------------------------------------------------------------------
          ! For each particle in the same cell
          DO i = 1, own_nplist
              ! Pick a reference particle
              p_ref = own_plist(i)
              ! For each neighbor particle up to the reference particle ...
              DO j = 1, i - 1
                  ! Pick a candidate for neighbor particle
                  p_neigh = own_plist(j)
                  ! If particles are neighbors ...
                  IF(isNeighbor(p_ref, p_neigh, xp, cutoff, skin)) THEN
                      ! If reference particle is real, no matter what neighbor
                      ! particle is ...
                      IF(p_ref .LE. clist%n_real_p)   THEN
                          ! add neighbor particle to verlet list of reference
                          ! particle.
                          nvlist(p_ref) = nvlist(p_ref)   + 1
#if __ACTION == __GET
                          vlist(nvlist(p_ref), p_ref) = p_neigh
#endif
                      ! If reference particle is ghost and ...
                      ELSE
                          ! neighbor particle is real,
                          IF(p_neigh .LE. clist%n_real_p) THEN
                               ! add reference particle to verlet list of
                               ! neighbor particle.
                               nvlist(p_neigh) = nvlist(p_neigh)   + 1
#if __ACTION == __GET
                               vlist(nvlist(p_neigh), p_neigh) = p_ref
#endif
                          ! If both particles are ghost and ...
                          ELSE
                              ! they are cross-neighbors,
                              IF(cross_neighbor(p_ref, p_neigh, xp, actual_domain)) THEN
                                  ! add neighbor particle to verlet list of reference
                                  ! particle.
                                  nvlist(p_ref) = nvlist(p_ref)   + 1
#if __ACTION == __GET
                                  vlist(nvlist(p_ref), p_ref) = p_neigh
#endif
                              END IF
                          END IF
                      END IF
                  END IF
              END DO
          END DO

          !---------------------------------------------------------------------
          !  Starting from the cell after the center cell upto last neighbor cell,
          !  get particles from neighbor particles and compare reference particles
          !  with those from neighboring cells. If two particles are found to be
          !  neighbors and:
          !  1. If reference particle is real, nvlist entry of it is incremented
          !     by 1
          !  2. If reference particle is ghost, but neighbor particle is real,
          !     nvlist entry of neighbor particle is incremented by 1.
          !  3. If both particles are ghost and they are cross-neighbors,
          !     nvlist entry of reference particle is incremented by 1.
          !---------------------------------------------------------------------
          ! For each neighbor cell
          DO i = ((3**ppm_dim + 1)/2) + 1, 3**ppm_dim
              ! Compute midpoint coordinates of the neighbor cell
              n_coor(1) = c_coor(1) + (ncells(i,1)*offset_coor(1))
              n_coor(2) = c_coor(2) + (ncells(i,2)*offset_coor(2))
              DO j = 3, ppm_dim
                  n_coor(j) = c_coor(j) + (ncells(i,j)*offset_coor(j))
              END DO

              ! If computed coordinates are not inside the domain, skip!
              IF(.NOT. inDomain(n_coor,domain))   cycle

              ! Get index of the neighbor cell.
              n_idx = getCellIdx(n_coor, c_depth, domain)

              ! Get particles in the neighbor cell
              CALL getParticlesInCell(n_idx, xp, clist, neigh_plist, neigh_nplist)

              ! For each particle in the reference cell
              DO m = 1, own_nplist
                  ! Pick a reference particle
                  p_ref = own_plist(m)
                  ! For each particle in the neighbor cell
                  DO n = 1, neigh_nplist
                      ! Pick a candidate for neighbor particle
                      p_neigh = neigh_plist(n)
                      ! If particles are neighbors and ...
                      IF(isNeighbor(p_ref, p_neigh, xp, cutoff, skin)) THEN
                          ! If reference particle is real, no matter what neighbor
                          ! particle is ...
                          IF(p_ref .LE. clist%n_real_p)   THEN
                              ! add neighbor particle to verlet list of reference
                              ! particle.
                              nvlist(p_ref) = nvlist(p_ref) + 1
#if __ACTION == __GET
                              vlist(nvlist(p_ref), p_ref) = p_neigh
#endif
                          ! If reference particle is ghost and ...
                          ELSE
                              ! neighbor particle is real,
                              IF(p_neigh .LE. clist%n_real_p)  THEN
                                  ! add reference particle to verlet list of
                                  ! neighbor particle.
                                  nvlist(p_neigh) = nvlist(p_neigh) + 1
#if __ACTION == __GET
                                  vlist(nvlist(p_neigh), p_neigh) = p_ref
#endif
                              ! If both particles are ghost and ...
                              ELSE
                                  ! they are cross-neighbors,
                                  IF(cross_neighbor(p_ref, p_neigh, xp, actual_domain)) THEN
                                      ! add neighbor particle to verlet list of
                                      ! reference particle.
                                      nvlist(p_ref) = nvlist(p_ref) + 1
#if __ACTION == __GET
                                      vlist(nvlist(p_ref), p_ref) = p_neigh
#endif
                                  END IF
                              END IF
                          END IF
                      END IF
                  END DO
              END DO
          END DO

          !---------------------------------------------------------------------
          !  For each reference particle, get coordinates and depth and increment
          !  depth by 1 in each iteration such that each reference particle finds
          !  its neighbors that have smaller cutoff radius, hence they are located
          !  in deeper cells. If the cell that is taken as reference cell does
          !  not contain any cells in its region, then this cell is added to
          !  empty list such that children cells will be skipped since they are
          !  also empty. If two particles are found to be neighbors and:
          !  1. If reference particle is real, nvlist entry of it is incremented
          !     by 1
          !  2. If reference particle is ghost, but neighbor particle is real,
          !     nvlist entry of neighbor particle is incremented by 1.
          !  3. If both particles are ghost and they are cross-neighbors,
          !     nvlist entry of reference particle is incremented by 1.
          !---------------------------------------------------------------------
          ! For each reference particle
          DO m = 1, own_nplist
              ! Pick a reference particle
              p_ref = own_plist(m)
              ! Get coordinates and depth of this particle
              CALL getParticleCoorDepth(p_ref, domain, p_coor, p_depth, xp, cutoff, skin)

              ! Increment the depth of the particle by 1 to look down at deeper cells
              p_depth = p_depth + 1

              ! Until we reach the maximum depth ...
              DO WHILE(p_depth .LE. clist%max_depth)
                  ! Get the child cell this particle is located in
                  c_idx = getCellIdx(p_coor, p_depth, domain)
                  ! Get index of parent cell to check whether parent is empty or not
                  parentIdx = parent(c_idx)

                  ! If the parent is empty, also put the child in the empty list.
                  IF (inEmptyList(parentIdx)) THEN
                     CALL putInEmptyList(c_idx)
                  ENDIF
                  ! Get cell coordinates and depth
                  CALL getCellCoor_Depth(c_idx, domain, c_coor, c_depth, &
 &                                       clist%max_depth, info)

                  ! Compute offset coordinates to compute midpoint coordinates
                  ! of neighbor cells at this depth.
                  offset_coor(1) = (domain(2) - domain(1))/(2**(c_depth-1))
                  offset_coor(2) = (domain(4) - domain(3))/(2**(c_depth-1))
                  DO i = 1, ppm_dim
                      offset_coor(i) = (domain(2*i) - domain(2*i-1))/&
    &                                  (2**(c_depth-1))
                  END DO

                  ! For each neighbor
                  DO i = 1, 3**ppm_dim
                      ! Compute midpoint coordinates of the neighbor
                      n_coor(1) = c_coor(1) + (ncells(i,1)*offset_coor(1))
                      n_coor(2) = c_coor(2) + (ncells(i,2)*offset_coor(2))
                      DO j = 3,ppm_dim
                          n_coor(j) = c_coor(j) + (ncells(i,j)*offset_coor(j))
                      END DO

                      ! If coordinates are not inside the domain, skip!
                      IF(.NOT. inDomain(n_coor, domain))   cycle

                      ! Get index of the neighbor cell
                      n_idx = getCellIdx(n_coor, c_depth, domain)
                      
                      ! if parent of neighbor empty,
                      ! put neighbor in emptylist and skip
                      IF (inEmptyList(parent(n_idx))) THEN
                          CALL putInEmptyList(n_idx)
                          CYCLE
                      ENDIF

                      ! Get particles in the neighbor cell
                      !-----------------------------------------------------
                      ! REMARK: If this cell has no particles, it will be stored
                      !         in empty cell inside getParticlesInCell
                      !         subroutine. If it has particles only at this
                      !         depth, but no deeper particles, it will again
                      !         be saved in the empty list just  before returning.
                      !-----------------------------------------------------
                      CALL getParticlesInCell(n_idx, xp, clist, neigh_plist, neigh_nplist)

                      ! For each neighbor element
                      DO n = 1, neigh_nplist
                          ! Pick a candidate for neighbor particle
                          p_neigh = neigh_plist(n)
                          ! If they are neighbors and ...
                          IF(isNeighbor(p_ref, p_neigh, xp, cutoff, skin))  THEN
                          ! If reference particle is real, no matter what neighbor
                          ! particle is ...
                              IF(p_ref .LE. clist%n_real_p)   THEN
                                  ! add neighbor particle to verlet list of
                                  ! reference particle.
                                  nvlist(p_ref) = nvlist(p_ref) + 1
#if __ACTION == __GET
                                  vlist(nvlist(p_ref), p_ref) = p_neigh
#endif
                              ! If reference particle is ghost and ...
                              ELSE
                                  ! neighbor particle is real,
                                  IF(p_neigh .LE. clist%n_real_p)  THEN
                                      ! add reference particle to verlet
                                      ! list of neighbor particle.
                                      nvlist(p_neigh) = nvlist(p_neigh) + 1
#if __ACTION == __GET
                                      vlist(nvlist(p_neigh), p_neigh) = p_ref
#endif
                                  ! If both particles are ghost and ...
                                  ELSE
                                      ! they are cross-neighbors,
                                      IF(cross_neighbor(p_ref, p_neigh, xp, &
                                      &  actual_domain)) THEN
                                      ! add neighbor particle to verlet
                                      ! list of reference particle.
                                          nvlist(p_ref) = nvlist(p_ref) + 1
#if __ACTION == __GET
                                          vlist(nvlist(p_ref), p_ref) = p_neigh
#endif
                                      END IF
                                  END IF
                              END IF
                          END IF
                      END DO
                 END DO
                 ! To search for deeper particles, increment depth of particle by 1.
                 p_depth = p_depth + 1
             END DO
          END DO

          ! Label all reference particles as USED so that they wont be visited again.
          DO i = 1, own_nplist
              used(own_plist(i)) = .TRUE.
          END DO
#if   __ACTION == __COUNT
#if   __KIND == __SINGLE_PRECISION
        END SUBROUTINE count_neigh_sym_s
#elif   __KIND == __DOUBLE_PRECISION
        END SUBROUTINE count_neigh_sym_d
#endif
#elif   __ACTION == __GET
#if   __KIND == __SINGLE_PRECISION
        END SUBROUTINE get_neigh_sym_s
#elif   __KIND == __DOUBLE_PRECISION
        END SUBROUTINE get_neigh_sym_d
#endif
#endif
