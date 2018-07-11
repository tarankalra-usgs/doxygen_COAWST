      SUBROUTINE wrt_rst (ng)
!
!svn $Id: wrt_rst.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes fields into restart NetCDF file.                !
!                                                                      !
!  Notice that only momentum is affected by the full time-averaged     !
!  masks.  If applicable, these mask contains information about        !
!  river runoff and time-dependent wetting and drying variations.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_sedbed
      USE mod_sediment
      USE mod_vegetation
      USE mod_vegarr
      USE mod_stepping
!
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: Fcount, gfactor, gtype, i, itrc, status, varid
      integer :: ntmp(1)
      real(r8) :: scale
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
      SourceFile='wrt_rst.F'
!
!-----------------------------------------------------------------------
!  Write out restart fields.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
      gfactor=1
!
!  Set time record index.
!
      RST(ng)%Rindex=RST(ng)%Rindex+1
      Fcount=RST(ng)%Fcount
      RST(ng)%Nrec(Fcount)=RST(ng)%Nrec(Fcount)+1
!
!  If requested, set time index to recycle time records in restart
!  file.
!
      IF (LcycleRST(ng)) THEN
        RST(ng)%Rindex=MOD(RST(ng)%Rindex-1,2)+1
      END IF
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, RST(ng)%name,                     &
     &                      TRIM(Vname(idtime,ng)), time(ng:),          &
     &                      (/RST(ng)%Rindex/), (/1/),                  &
     &                      ncid = RST(ng)%ncid,                        &
     &                      varid = RST(ng)%Vid(idtime))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out wet/dry mask at PSI-points.
!
      scale=1.0_r8
      gtype=gfactor*p2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idPwet),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % pmask,                              &
     &                   GRID(ng) % pmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idPwet)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at RHO-points.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idRwet),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % rmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRwet)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at U-points.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUwet),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % umask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUwet)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at V-points.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVwet),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % vmask,                              &
     &                   GRID(ng) % vmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVwet)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out free-surface (m).
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idFsur),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % rmask,                              &
     &                   OCEAN(ng) % zeta(:,:,kstp(ng)),                    &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUbar),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % umask_full,                         &
     &                   OCEAN(ng) % ubar(:,:,kstp(ng)))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVbar),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % vmask_full,                         &
     &                   OCEAN(ng) % vbar(:,:,kstp(ng)))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idUvel),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
     &                   GRID(ng) % umask_full,                         &
     &                   OCEAN(ng) % u(:,:,:,nrhs(ng)))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVvel),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
     &                   GRID(ng) % vmask_full,                         &
     &                   OCEAN(ng) % v(:,:,:,nrhs(ng)))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Tid(itrc),   &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
     &                     GRID(ng) % rmask,                            &
     &                     OCEAN(ng) % t(:,:,:,nrhs(ng),itrc))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
!
!  Write out density anomaly.
!
      scale=1.0_r8
      gtype=gfactor*r3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idDano),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
     &                   GRID(ng) % rmask,                              &
     &                   OCEAN(ng) % rho)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idDano)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out vertical viscosity coefficient.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idVvis),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
     &                   GRID(ng) % rmask,                              &
     &                   MIXING(ng) % Akv,                              &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvis)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTdif),   &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
     &                   GRID(ng) % rmask,                              &
     &                   MIXING(ng) % Akt(:,:,:,itemp),                 &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTdif)), RST(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out vegetation properties 
! 
      DO i=1,NVEGP 
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid,                    &
     &                       RST(ng)%Vid(idvprp(i)),                    &
     &                       RST(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, NVEG, scale,        &
     &                       GRID(ng) % rmask,                          &
     &                       VEG(ng) % plant(:,:,:,i))
          IF (status.ne.nf90_noerr) THEN 
            IF (Master) THEN 
              WRITE (stdout,10) TRIM(Vname(1,idvprp(i))), RST(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
      END DO 
!
!
!  Write out wave dissipation due to vegetation 
! 
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idWdvg), &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     VEG(ng)%Dissip_veg )
        IF (status.ne.nf90_noerr) THEN 
          IF (Master) THEN 
            WRITE (stdout,10) TRIM(Vname(1,idWdvg)), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
! 
!
!  Write out initial masking for marshes 
! 
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTims), &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     VEG(ng)%marsh_mask)
        IF (status.ne.nf90_noerr) THEN 
          IF (Master) THEN 
            WRITE (stdout,10) TRIM(Vname(1,idTims)), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out wave thrust on marsh output 
! 
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTmsk), &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     VEG(ng)%mask_thrust)
        IF (status.ne.nf90_noerr) THEN 
          IF (Master) THEN 
            WRITE (stdout,10) TRIM(Vname(1,idTmsk)), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Define maximum thrust due to waves. 
!
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTmax), &
     &                     RST(ng)%Rindex, gtype,                       &    
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     VEG(ng)%Thrust_max)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmax)), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Define maximum thrust due to waves. 
!
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTmax), &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     VEG(ng)%Thrust_max)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTmax)), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Define Tonelli masking based thrust due to waves. 
!
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idTton), &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     VEG(ng)%Thrust_tonelli)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTton)), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out exposed sediment layer properties. Notice that only the
!  first four properties (mean grain diameter, mean grain density,
!  mean settling velocity, mean critical erosion stress,
!  ripple length and ripple height) are written.
!
      DO i=1,6
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid,                      &
     &                     RST(ng)%Vid(idBott(i)),                      &
     &                     RST(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
     &                     GRID(ng) % rmask,                            &
     &                     SEDBED(ng) % bottom(:,:,i))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idBott(i))), RST(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
!
!  Write out 2D U-momentum stokes velocity.
!
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idU2Sd), &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % umask,                              &
     &                   OCEAN(ng) % ubar_stokes)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idU2Sd)),          &
     &                                  RST(ng)%Rindex
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out 2D V-momentum stokes velocity.
!
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idV2Sd), &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
     &                   GRID(ng) % vmask,                              &
     &                   OCEAN(ng) % vbar_stokes)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV2Sd)),          &
     &                                  RST(ng)%Rindex
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out 3D U-momentum stokes velocity.
!
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idU3Sd), &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
     &                   GRID(ng) % umask,                              &
     &                   OCEAN(ng) % u_stokes)
        IF (status.ne.nf90_noerr) THEN
          IF(Master) WRITE (stdout,10) TRIM(Vname(1,idU3Sd)),           &
     &                                 RST(ng)%Rindex
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Write out 3D V-momentum stokes velocity.
!
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, RST(ng)%ncid, RST(ng)%Vid(idV3Sd), &
     &                   RST(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
     &                   GRID(ng) % vmask,                              &
     &                   OCEAN(ng) % v_stokes)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idV3Sd)),          &
     &                                  RST(ng)%Rindex
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Synchronize restart NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, RST(ng)%name, RST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
      IF (Master) WRITE (stdout,20) kstp(ng), nrhs(ng), RST(ng)%Rindex
!
  10  FORMAT (/,' WRT_RST - error while writing variable: ',a,/,11x,    &
     &        'into restart NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_RST   - wrote re-start fields (Index=', i1,       &
     &        ',',i1,') into time record = ',i7.7)
      RETURN
      END SUBROUTINE wrt_rst
