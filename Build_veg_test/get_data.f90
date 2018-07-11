      SUBROUTINE get_data (ng)
!
!svn $Id: get_data.F 791 2016-05-05 22:39:42Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in forcing, climatology and other data from      !
!  NetCDF files.  If there is more than one time-record,  data is      !
!  loaded into global  two-time  record arrays. The interpolation      !
!  is carried elsewhere.                                               !
!                                                                      !
!  Currently, this routine is only executed in serial mode by the      !
!  main thread.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
      USE mod_clima
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_sources
      USE mod_stepping
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical, dimension(3) :: update =                                 &
     &         (/ .FALSE., .FALSE., .FALSE. /)
      integer :: ILB, IUB, JLB, JUB
      integer :: LBi, UBi, LBj, UBj
      integer :: i, ic, my_tile
!
!  Lower and upper bounds for nontiled (global values) boundary arrays.
!
      my_tile=-1                           ! for global values
      ILB=BOUNDS(ng)%LBi(my_tile)
      IUB=BOUNDS(ng)%UBi(my_tile)
      JLB=BOUNDS(ng)%LBj(my_tile)
      JUB=BOUNDS(ng)%UBj(my_tile)
!
!  Lower and upper bounds for tiled arrays.
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 3)
!
!=======================================================================
!  Read in forcing data from FORCING NetCDF file.
!=======================================================================
!
!-----------------------------------------------------------------------
!  Point Sources/Sinks time dependent data.
!-----------------------------------------------------------------------
!
!  Point Source/Sink vertically integrated mass transport.
!
      IF (LuvSrc(ng).or.LwSrc(ng)) THEN
        CALL get_ngfld (ng, iNLM, idRtra, SSF(ng)%ncid,                 &
     &                  1, SSF(ng), update(1),                          &
     &                  1, Nsrc(ng), 1, 2, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % QbarG(:,1))
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Tracer Sources/Sinks.
!
      DO i=1,NT(ng)
        IF (LtracerSrc(i,ng)) THEN
          CALL get_ngfld (ng, iNLM, idRtrc(i), SSF(ng)%ncid,            &
     &                    1, SSF(ng), update(1),                        &
     &                    1, Nsrc(ng), N(ng), 2, 1, Nsrc(ng), N(ng),    &
     &                    SOURCES(ng) % TsrcG(:,:,:,i))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
!
!=======================================================================
!  Read in open boundary conditions from BOUNDARY NetCDF file.  In
!  grid refinement, only the coarser grid (RefineScale(ng)=0) open
!  boundary conditions data is processed and needed.
!=======================================================================
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        IF (LBC(iwest,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(iwest),                      &
     &                    ncBRYid(idU3bc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(iwest,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(iwest),                      &
     &                    ncBRYid(idV3bc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_west)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(ieast,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(ieast),                      &
     &                    ncBRYid(idU3bc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(ieast,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(ieast),                      &
     &                    ncBRYid(idV3bc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_east)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(isouth,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(isouth),                     &
     &                    ncBRYid(idU3bc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(isouth,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(isouth),                     &
     &                    ncBRYid(idV3bc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_south)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(inorth,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(inorth),                     &
     &                    ncBRYid(idU3bc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (LBC(inorth,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(inorth),                     &
     &                    ncBRYid(idV3bc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng), update(1),           &
     &                    ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_north)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
!
      IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
        DO i=1,NT(ng)
          IF (LBC(iwest,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(iwest,i),                  &
     &                      ncBRYid(idTbry(iwest,i),ng),                &
     &                      nBCfiles(ng), BRY(1,ng), update(1),         &
     &                      JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_west(:,:,:,i))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
        DO i=1,NT(ng)
          IF (LBC(ieast,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(ieast,i),                  &
     &                      ncBRYid(idTbry(ieast,i),ng),                &
     &                      nBCfiles(ng), BRY(1,ng), update(1),         &
     &                      JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_east(:,:,:,i))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
        DO i=1,NT(ng)
          IF (LBC(isouth,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(isouth,i),                 &
     &                      ncBRYid(idTbry(isouth,i),ng),               &
     &                      nBCfiles(ng), BRY(1,ng), update(1),         &
     &                      ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_south(:,:,:,i))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
        DO i=1,NT(ng)
          IF (LBC(inorth,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(inorth,i),                 &
     &                      ncBRYid(idTbry(inorth,i),ng),               &
     &                      nBCfiles(ng), BRY(1,ng), update(1),         &
     &                      ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_north(:,:,:,i))
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
      END IF
!
!=======================================================================
!  Read in data from Climatology NetCDF file.
!=======================================================================
!
!  Free-surface.
!
      IF (LsshCLM(ng)) THEN
        CALL get_2dfld (ng, iNLM, idSSHc, ncCLMid(idSSHc,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
     &                  GRID(ng) % rmask,                               &
     &                  CLIMA(ng) % sshG)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  2D momentum.
!
      IF (Lm2CLM(ng)) THEN
        CALL get_2dfld (ng, iNLM, idUbcl, ncCLMid(idUbcl,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
     &                  GRID(ng) % umask,                               &
     &                  CLIMA(ng) % ubarclmG)
        IF (exit_flag.ne.NoError) RETURN
!
        CALL get_2dfld (ng, iNLM, idVbcl, ncCLMid(idVbcl,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
     &                  GRID(ng) % vmask,                               &
     &                  CLIMA(ng) % vbarclmG)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  3D momentum.
!
      IF (Lm3CLM(ng)) THEN
        CALL get_3dfld (ng, iNLM, idUclm, ncCLMid(idUclm,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
     &                  GRID(ng) % umask,                               &
     &                  CLIMA(ng) % uclmG)
        IF (exit_flag.ne.NoError) RETURN
!
        CALL get_3dfld (ng, iNLM, idVclm, ncCLMid(idVclm,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
     &                  GRID(ng) % vmask,                               &
     &                  CLIMA(ng) % vclmG)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Tracers.
!
      ic=0
      DO i=1,NT(ng)
        IF (LtracerCLM(i,ng)) THEN
          ic=ic+1
          CALL get_3dfld (ng, iNLM, idTclm(i),                          &
     &                    ncCLMid(idTclm(i),ng),                        &
     &                    nCLMfiles(ng), CLM(1,ng), update(1),          &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,           &
     &                    GRID(ng) % rmask,                             &
     &                    CLIMA(ng) % tclmG(:,:,:,:,ic))
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 3)
      RETURN
      END SUBROUTINE get_data
