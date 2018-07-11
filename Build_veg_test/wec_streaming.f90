      MODULE wec_streaming_mod
!
!svn $Id: wec_streaming.F 1428 2008-03-12 13:07:21Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!                                                   Nirnimesh Kumar    !
!================================================== John C. Warner ====!
!                                                                      !
!  This routine computes the terms corresponding to vortex forces in   !
!  momentum equations.                                                 !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Uchiyama, Y., McWilliams, J.C., and Shchepetkin, A.F. (2010).       !
!  Wave current interacation in an oceanic circulation model with a    !
!  vortex-force formalism: Applications to surf zone, Ocean Modeling,  !
!  34, 16-35.
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: wec_streaming
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wec_streaming (ng, tile)
!***********************************************************************
!
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_stepping
      USE mod_diags
      USE mod_vegarr
      USE vegetation_stream_mod, ONLY : vegetation_stream_cal
!
      integer, intent(in) :: ng, tile
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
      CALL vegetation_stream_cal (ng, tile)
      CALL wclock_on (ng, iNLM, 21)
      CALL wec_streaming_tile (ng, tile, LBi, UBi, LBj, UBj, N(ng),     &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs(ng),                                &
     &                         GRID(ng) % angler,                       &
     &                         GRID(ng) % om_u,                         &
     &                         GRID(ng) % om_v,                         &
     &                         GRID(ng) % on_u,                         &
     &                         GRID(ng) % on_v,                         &
     &                         GRID(ng) % Hz,                           &
     &                         GRID(ng) % z_r,                          &
     &                         GRID(ng) % z_w,                          &
     &                         FORCES(ng) % Hwave,                      &
     &                         FORCES(ng) % Dwave,                      &
     &                         FORCES(ng) % Lwave,                      &
     &                         FORCES(ng) % Pwave_top,                  &
     &                         FORCES(ng) % Dissip_fric,                &
     &                         DIAGS(ng) % DiaRU,                       &
     &                         DIAGS(ng) % DiaRV,                       &
     &                         MIXING(ng) % rubst2d,                    &
     &                         MIXING(ng) % rvbst2d,                    &
     &                         VEG(ng) % BWDXL_veg,                     &
     &                         VEG(ng) % BWDYL_veg,                     &
     &                         MIXING(ng) % rustr3d,                    &
     &                         MIXING(ng) % rvstr3d)
      CALL wclock_off (ng, iNLM, 21)
      RETURN
      END SUBROUTINE wec_streaming
!
!***********************************************************************
      SUBROUTINE wec_streaming_tile (ng, tile, LBi, UBi, LBj, UBj, UBk, &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               nrhs,                              &
     &                               angler,                            &
     &                               om_u, om_v, on_u, on_v,            &
     &                               Hz, z_r, z_w,                      &
     &                               Hwave, Dwave, Lwave,               &
     &                               Pwave_top,                         &
     &                               Dissip_fric,                       &
     &                               DiaRU, DiaRV,                      &
     &                               rubst2d, rvbst2d,                  &
     &                               BWDXL_veg, BWDYL_veg,              &
     &                               rustr3d,  rvstr3d) 
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
      USE bc_2d_mod
      USE bc_3d_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
      real(r8), intent(in) :: Pwave_top(LBi:,LBj:)
      real(r8), intent(in) :: Dissip_fric(LBi:,LBj:)
      real(r8), intent(inout) :: DiaRU(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: DiaRV(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: rubst2d(LBi:,LBj:)
      real(r8), intent(inout) :: rvbst2d(LBi:,LBj:)
      real(r8), intent(inout) :: rustr3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: rvstr3d(LBi:,LBj:,:)
      real(r8), intent(inout) :: BWDXL_veg(LBi:,LBj:,:)
      real(r8), intent(inout) :: BWDYL_veg(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8) :: cff, cff1, cff2, cff3, cff4
      real(r8) :: fac2, sqrt2
      real(r8), parameter :: ks=0.03_r8
      real(r8), parameter :: awd=1.0_r8
      real(r8), parameter :: KWDmax=200.0_r8
      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8), parameter :: kDmax = 5.0_r8
      real(r8), parameter :: Lwave_min = 1.0_r8
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: oDstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: kD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waven
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavenx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveny
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: KWD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: owd
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: EWD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: sigma
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: osigma
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BWDXL
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: BWDYL
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
      sqrt2=SQRT(2.0_r8)
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
!
!  Compute total depth
!
          Dstp(i,j)=z_w(i,j,N(ng))-z_w(i,j,0)
          oDstp(i,j)=1.0_r8/Dstp(i,j)
!
!  Compute wave amplitude (0.5*Hrms), wave number, intrinsic frequency.
!
          waven(i,j)=2.0_r8*pi/MAX(Lwave(i,j),Lwave_min)
          cff=1.5_r8*pi-Dwave(i,j)-angler(i,j)
          wavenx(i,j)=waven(i,j)*COS(cff)
          waveny(i,j)=waven(i,j)*SIN(cff)
          sigma(i,j)=MIN(SQRT(g*waven(i,j)*TANH(waven(i,j)*Dstp(i,j))), &
     &                   2.0_r8)
          osigma(i,j)=1.0_r8/sigma(i,j)
!
!  Compute wave celerity and nonlinear water depth
!
          kD(i,j)=MIN(waven(i,j)*Dstp(i,j)+eps,kDmax)
!
!  Compute metrics for vertical bottom streaming distribution
!          
          owd(i,j)=0.0_r8
!  Bottom Orbital Velocity
          cff=0.25_r8*sqrt2*sigma(i,j)*Hwave(i,j)/                      &
     &        (SINH(kD(i,j))+eps)
          cff1=awd*0.09_r8*ks*(cff/                                     &
     &         (ks*sigma(i,j)))**0.82_r8
          KWD(i,j)=MIN(Dstp(i,j)/(cff1+eps),KWDmax)
          DO k=1,N(ng)
            cff2=(z_r(i,j,N(ng))-z_r(i,j,k))*oDstp(i,j)
            owd(i,j)=owd(i,j)+                                          &
     &               Hz(i,j,k)*COSH(KWD(i,j)*cff2)
          END DO
          owd(i,j)=1.0_r8/owd(i,j)
!
!  Wave friction factor (based on Soulsby, 1997).
!  Hold this constant for now.  Need to add logic for 
!  zo if no bbl on, or zo_apparent if use a bbl.
         cff3=MIN(1.39_r8*(sigma(i,j)*(ks/30.0_r8)/                     &
     &               cff)**0.52_r8,0.2_r8)
!
!  Wave dissipation rate due to wave bottom drag Reniers et al. (2004b)
!
          EWD(i,j)=Dissip_fric(i,j)
        END DO
      END DO
!
!  Initialize depth-avg forcing terms for subsequent computation.
!  For now, these are used for the diagnostics.  This logic will change
!  when 2D method is implemented.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          rubst2d(i,j)=0.0_r8
          rvbst2d(i,j)=0.0_r8
        END DO
      END DO
!
! Compute bottom streaming based acceleration terms
!
      K_LOOP : DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            fac2=(z_r(i,j,N(ng))-z_r(i,j,k))*oDstp(i,j)
            cff2=COSH(fac2*KWD(i,j))
            cff3=EWD(i,j)*osigma(i,j)
            BWDXL(i,j)=cff2*cff3*wavenx(i,j)*owd(i,j)
            BWDYL(i,j)=cff2*cff3*waveny(i,j)*owd(i,j)
          END DO
        END DO
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            BWDXL(i,j)=BWDXL_veg(i,j,k)
            BWDYL(i,j)=BWDYL_veg(i,j,k)
          END DO
        END DO
!
! Compute contribution to U-momentum
! 
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff1=0.5_r8*(Hz(i-1,j,k)+                                   &
     &                   Hz(i  ,j,k))
            cff=0.5_r8*(BWDXL(i  ,j)*Hz(i  ,j,k)+                       &
     &                  BWDXL(i-1,j)*Hz(i-1,j,k))
            rustr3d(i,j,k)=rustr3d(i,j,k)-cff
            rubst2d(i,j)=rubst2d(i,j)+cff*cff1*om_u(i,j)*on_u(i,j)
            DiaRU(i,j,k,nrhs,M3bstm)=cff*om_u(i,j)*on_u(i,j)
          END DO
        END DO
!
! Compute contribution to V-momentum
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff1=0.5_r8*(Hz(i,j  ,k)+                                   &
     &                   Hz(i,j-1,k))
            cff=0.5_r8*(BWDYL(i,j  )*Hz(i,j  ,k)+                       &
     &                  BWDYL(i,j-1)*Hz(i,j-1,k))
            rvstr3d(i,j,k)=rvstr3d(i,j,k)-cff
            rvbst2d(i,j)=rvbst2d(i,j)+cff*cff1*om_v(i,j)*on_v(i,j)
            DiaRV(i,j,k,nrhs,M3bstm)=cff*om_v(i,j)*on_v(i,j)
          END DO
        END DO
!
      END DO K_LOOP
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff=0.5_r8*(Dstp(i,j)+Dstp(i-1,j))
          rubst2d(i,j)=rubst2d(i,j)/cff
!
          IF (j.ge.JstrV) THEN
            cff=0.5_r8*(Dstp(i,j)+Dstp(i,j-1))
            rvbst2d(i,j)=rvbst2d(i,j)/cff
          END IF
        END DO
      END DO
!
!  Apply boundary conditions.
      CALL bc_u3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  rustr3d)
      CALL bc_v3d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  rvstr3d)
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  rubst2d)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  rvbst2d)
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rustr3d,  rvstr3d)
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rubst2d, rvbst2d)
      RETURN
      END SUBROUTINE wec_streaming_tile
      END MODULE wec_streaming_mod
