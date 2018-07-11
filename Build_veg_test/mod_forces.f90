      MODULE mod_forces
!
!svn $Id: mod_forces.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Surface momentum stresses.                                          !
!                                                                      !
!  sustr        Kinematic surface momentum flux (wind stress) in       !
!                 the XI-direction (m2/s2) at horizontal U-points.     !
!  sustrG       Latest two-time snapshots of input "sustr" grided      !
!                 data used for interpolation.                         !
!  svstr        Kinematic surface momentum flux (wind stress) in       !
!                 the ETA-direction (m2/s2) at horizontal V-points.    !
!  svstrG       Latest two-time snapshots of input "svstr" grided      !
!                 data used for interpolation.                         !
!  Taux         Surface stress in the XI-direction at rho points       !
!                 from atm model.                                      !
!  Tauy         Surface stress in the ETA-direction at rho points      !
!                 from atm model.                                      !
!                                                                      !
!  Bottom momentum stresses.                                           !
!                                                                      !
!  bustr        Kinematic bottom momentum flux (bottom stress) in      !
!                 the XI-direction (m2/s2) at horizontal U-points.     !
!  bvstr        Kinematic bottom momentum flux (bottom stress) in      !
!                 ETA-direction (m2/s2) at horizontal V-points.        !
!                                                                      !
!  Surface wind induced waves.                                         !
!                                                                      !
!  Hwave        Surface wind induced wave height (m).                  !
!  HwaveG       Latest two-time snapshots of input "Hwave" grided      !
!                 data used for interpolation.                         !
!  Dwave        Surface wind induced wave direction (radians).         !
!  DwaveG       Latest two-time snapshots of input "Dwave" grided      !
!                 data used for interpolation.                         !
!  Lwave        Mean surface wavelength read in from swan output       !
!  LwaveG       Latest two-time snapshots of input "Lwave" grided      !
!                 data used for interpolation.                         !
!  Lwavep       Peak surface wavelength read in from swan output       !
!  LwavepG      Latest two-time snapshots of input "Lwavep" grided     !
!                 data used for interpolation.                         !
!  Pwave_top    Wind induced surface wave period (s).                  !
!  Pwave_topG   Latest two-time snapshots of input "Pwave_top" grided  !
!                 data used for interpolation.                         !
!  Pwave_bot    Wind induced bottom wave period (s).                   !
!  Pwave_botG   Latest two-time snapshots of input "Pwave_bot" grided  !
!                 data used for interpolation.                         !
!  Uwave_rms    Bottom orbital velocity read in from swan output       !
!  Uwave_rmsG   Latest two-time snapshots of input "Uwave_rms" grided  !
!                 data used for interpolation.                         !
!  wave_dissip  Wave dissipation                                       !
!  wave_dissipG Latest two-time snapshots of input "wave_dissip"       !
!                 gridded data used for interpolation.                 !
!  Wave_break   Percent of wave breaking for use with roller model.    !
!  Wave_breakG  Latest two-time snapshots of input "wave_break"        !
!                 gridded data used for interpolation.                 !
!  Wave_ds      Wave directional spreading.                            !
!  Wave_qp      Wave spectrum peakedness.                              !
!                                                                      !
!  Solar shortwave radiation flux.                                     !
!                                                                      !
!  srflx        Kinematic surface shortwave solar radiation flux       !
!                 (Celsius m/s) at horizontal RHO-points               !
!  srflxG       Latest two-time snapshots of input "srflx" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Cloud fraction.                                                     !
!                                                                      !
!  cloud        Cloud fraction (percentage/100).                       !
!  cloudG       Latest two-time snapshots of input "cloud" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface heat fluxes, Atmosphere-Ocean bulk parameterization.        !
!                                                                      !
!  lhflx        Kinematic net latent heat flux (degC m/s).             !
!  lrflx        Kinematic net longwave radiation (degC m/s).           !
!  shflx        Kinematic net sensible heat flux (degC m/s).           !
!                                                                      !
!  Surface air humidity.                                               !
!                                                                      !
!  Hair         Surface air specific (g/kg) or relative humidity       !
!                 (percentage).                                        !
!  HairG        Latest two-time snapshots of input "Hair" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface air pressure.                                               !
!                                                                      !
!  Pair         Surface air pressure (mb).                             !
!  PairG        Latest two-time snapshots of input "Pair" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface air temperature.                                            !
!                                                                      !
!  Tair         Surface air temperature (Celsius)                      !
!  TairG        Latest two-time snapshots of input "Tair" grided       !
!                 data used for interpolation.                         !
!  PotT         Surface air potential temperature (Kelvin)             !
!  Surface Winds.                                                      !
!                                                                      !
!  Uwind        Surface wind in the XI-direction (m/s) at              !
!                 horizontal RHO-points.                               !
!  UwindG       Latest two-time snapshots of input "Uwind" grided      !
!                 data used for interpolation.                         !
!  Vwind        Surface wind in the ETA-direction (m/s) at             !
!                 horizontal RHO-points.                               !
!  VwindG       Latest two-time snapshots of input "Vwind" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Rain fall rate.                                                     !
!                                                                      !
!  evap         Evaporation rate (kg/m2/s).                            !
!  rain         Rain fall rate (kg/m2/s).                              !
!  rainG        Latest two-time snapshots of input "rain" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  Snow fall rate.                                                     !
!                                                                      !
!  snow         Snow fall rate (kg/m2/s).                              !
!  snowG        Latest two-time snapshots of input "snow" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface tracer fluxes.                                              !
!                                                                      !
!  stflx        Kinematic surface flux of tracer type variables        !
!                 (temperature: degC m/s; salinity: PSU m/s) at        !
!                 horizontal RHO-points.                               !
!  stflxG       Latest two-time snapshots of input "stflx" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Bottom tracer fluxes.                                               !
!                                                                      !
!  btflx        Kinematic bottom flux of tracer type variables         !
!                 (temperature: degC m/s; salinity: PSU m/s) at        !
!                horizontal RHO-points.                                !
!  btflxG       Latest two-time snapshots of input "btflx" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface heat flux correction.                                       !
!                                                                      !
!  dqdt         Kinematic surface net heat flux sensitivity to SST,    !
!                 d(Q)/d(SST), (m/s).                                  !
!  dqdtG        Latest two-time snapshots of input "dqdt" grided       !
!                 data used for interpolation.                         !
!  sst          Sea surface temperature (Celsius).                     !
!  sstG         Latest two-time snapshots of input "sst" grided        !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface freshwater flux correction.                                 !
!                                                                      !
!  sss          Sea surface salinity (PSU).                            !
!  sssG         Latest two-time snapshots of input "sss" grided        !
!                 data used for interpolation.                         !
!  sssflx       Sea surface salinity flux correction.                  !
!  sssflxG      Latest two-time snapshots of input "sssflx" grided     !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface spectral downwelling irradiance.                            !
!                                                                      !
!  SpecIr       Spectral irradiance (NBands) from 400-700 nm at        !
!                 5 nm bandwidth.                                      !
!  avcos        Cosine of average zenith angle of downwelling          !
!                 spectral photons.                                    !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
        implicit none
        TYPE T_FORCES
!
!  Nonlinear model state.
!
          real(r8), pointer :: sustr(:,:)
          real(r8), pointer :: svstr(:,:)
          real(r8), pointer :: bustr(:,:)
          real(r8), pointer :: bvstr(:,:)
          real(r8), pointer :: Dwave(:,:)
          real(r8), pointer :: DwaveG(:,:,:)
          real(r8), pointer :: Hwave(:,:)
          real(r8), pointer :: HwaveG(:,:,:)
          real(r8), pointer :: Lwave(:,:)
          real(r8), pointer :: LwaveG(:,:,:)
          real(r8), pointer :: Pwave_top(:,:)
          real(r8), pointer :: Pwave_topG(:,:,:)
          real(r8), pointer :: Pwave_bot(:,:)
          real(r8), pointer :: Pwave_botG(:,:,:)
          real(r8), pointer :: Uwave_rms(:,:)
          real(r8), pointer :: Uwave_rmsG(:,:,:)
          real(r8), pointer :: Dissip_break(:,:)
          real(r8), pointer :: Dissip_wcap(:,:)
          real(r8), pointer :: Dissip_breakG(:,:,:)
          real(r8), pointer :: Dissip_wcapG(:,:,:)
          real(r8), pointer :: Dissip_fric(:,:)
          real(r8), pointer :: Dissip_fricG(:,:,:)
          real(r8), pointer :: stflx(:,:,:)
          real(r8), pointer :: btflx(:,:,:)
        END TYPE T_FORCES
        TYPE (T_FORCES), allocatable :: FORCES(:)
      CONTAINS
      SUBROUTINE allocate_forces (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( FORCES(Ngrids) )
!
!  Nonlinear model state
!
      allocate ( FORCES(ng) % sustr(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % svstr(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % bustr(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % bvstr(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % Dwave(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % DwaveG(LBi:UBi,LBj:UBj,2) )
      allocate ( FORCES(ng) % Hwave(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % HwaveG(LBi:UBi,LBj:UBj,2) )
      allocate ( FORCES(ng) % Lwave(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % LwaveG(LBi:UBi,LBj:UBj,2) )
      allocate ( FORCES(ng) % Pwave_top(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % Pwave_topG(LBi:UBi,LBj:UBj,2) )
      allocate ( FORCES(ng) % Pwave_bot(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % Pwave_botG(LBi:UBi,LBj:UBj,2) )
      allocate ( FORCES(ng) % Uwave_rms(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % Uwave_rmsG(LBi:UBi,LBj:UBj,2) )
      allocate ( FORCES(ng) % Dissip_break(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % Dissip_wcap(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % Dissip_breakG(LBi:UBi,LBj:UBj,2) )
      allocate ( FORCES(ng) % Dissip_wcapG(LBi:UBi,LBj:UBj,2) )
      allocate ( FORCES(ng) % Dissip_fric(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % Dissip_fricG(LBi:UBi,LBj:UBj,2) )
      allocate ( FORCES(ng) % stflx(LBi:UBi,LBj:UBj,NT(ng)) )
      allocate ( FORCES(ng) % btflx(LBi:UBi,LBj:UBj,NT(ng)) )
      RETURN
      END SUBROUTINE allocate_forces
      SUBROUTINE initialize_forces (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k
      integer :: itrc
      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!  Set array initialization range.
!
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
       IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            FORCES(ng) % sustr(i,j) = IniVal
            FORCES(ng) % svstr(i,j) = IniVal
            FORCES(ng) % bustr(i,j) = IniVal
            FORCES(ng) % bvstr(i,j) = IniVal
            FORCES(ng) % Dwave(i,j) = IniVal
            FORCES(ng) % DwaveG(i,j,1) = IniVal
            FORCES(ng) % DwaveG(i,j,2) = IniVal
            FORCES(ng) % Hwave(i,j) = IniVal
            FORCES(ng) % HwaveG(i,j,1) = IniVal
            FORCES(ng) % HwaveG(i,j,2) = IniVal
            FORCES(ng) % Lwave(i,j) = IniVal
            FORCES(ng) % LwaveG(i,j,1) = IniVal
            FORCES(ng) % LwaveG(i,j,2) = IniVal
            FORCES(ng) % Pwave_top(i,j) = IniVal
            FORCES(ng) % Pwave_topG(i,j,1) = IniVal
            FORCES(ng) % Pwave_topG(i,j,2) = IniVal
            FORCES(ng) % Pwave_bot(i,j) = IniVal
            FORCES(ng) % Pwave_botG(i,j,1) = IniVal
            FORCES(ng) % Pwave_botG(i,j,2) = IniVal
            FORCES(ng) % Uwave_rms(i,j) = IniVal
            FORCES(ng) % Uwave_rmsG(i,j,1) = IniVal
            FORCES(ng) % Uwave_rmsG(i,j,2) = IniVal
            FORCES(ng) % Dissip_break(i,j) = IniVal
            FORCES(ng) % Dissip_wcap(i,j) = IniVal
            FORCES(ng) % Dissip_breakG(i,j,1) = IniVal
            FORCES(ng) % Dissip_breakG(i,j,2) = IniVal
            FORCES(ng) % Dissip_wcapG(i,j,1) = IniVal
            FORCES(ng) % Dissip_wcapG(i,j,2) = IniVal
            FORCES(ng) % Dissip_fric(i,j) = IniVal
            FORCES(ng) % Dissip_fricG(i,j,1) = IniVal
            FORCES(ng) % Dissip_fricG(i,j,2) = IniVal
            DO itrc=1,NT(ng)
              FORCES(ng) % stflx(i,j,itrc) = IniVal
              FORCES(ng) % btflx(i,j,itrc) = IniVal
            END DO
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE initialize_forces
      END MODULE mod_forces
