      MODULE mod_vegarr
!                                                                      !
!svn $Id: vegarr_mod.h 429 2009-12-20 17:30:26Z arango $               !
!================================================== Hernan G. Arango ==!
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!================================================== John C. Warner ====!
!==================================================== Neil K. Ganju  ==! 
!==================================================== Alexis Beudin  ==! 
!==================================================Tarandeep S. Kalra==!
!                                                                      !
!  Vegetation Model Kernel Variables:                                  !
!  plant         Vegetation variable properties:                       !
!                   plant(:,:,:,phght) => height                       !
!                   plant(:,:,:,pdens) => density                      !
!                   plant(:,:,:,pthck) => thickness                    !
!                   plant(:,:,:,pdiam) => diameter                     !
!                   plant(:,:,:,pabbm) => above ground biomass         !
!                   plant(:,:,:,pbgbm) => below ground biomass         !
!  ru_veg         Momentum term for x direction(takes account for all  !
!                 vegetation types)                                    !
!  rv_veg         Momentum term for x direction(takes account for all  !
!                 vegetation types)                                    !
!  ru_veg_loc     Momentum term for x direction(takes account for only !
!                 local vegetation type)                               !
!  rv_veg_loc     Momentum term for x direction(takes account for all  !
!                 local vegetation types)                              !
!  step2d_uveg    Momentum term for 2d x direction                     !
!  step2d_vveg    Momentum term for 2d y direction                     !
!  bend           Bending for each vegetation                          !
!  Lveg           Effective blade length                               ! 
!  tke_veg        Turbulent kinetic energy from vegetation             !
!  gls_veg        Length scale change from vegetation                  !
!  dissip_veg     Dissipation from the SWAN model due to vegetation    !
!  BWDXL_veg      Wave streaming effect due to vegetation              !
!  BWDYL_veg      Wave streaming effect due to vegetation              !
!  visc2d_r_veg   Effect of viscosity change at vegetation interface   ! 
!  visc3d_r_veg   Effect of viscosity change at vegetation interface   ! 
!  marsh_mask     User input of marsh masking at MSL                   ! 
!  mask_thrust    Tonellis masking for wave thrust on marshes          !
!  Thrust_max     Maximum thrust from wave to marshes                  !
!  Thrust_tonelli Reduced thrust from tonelli's masking                !
!                                                                      !
!======================================================================!
!
      USE mod_kinds
!
      implicit none
      TYPE T_VEG
!
!  Nonlinear model state.
!
        real(r8), pointer :: plant(:,:,:,:)
!  Momentum terms go back to act as sink in rhs
        real(r8), pointer :: ru_veg(:,:,:)
        real(r8), pointer :: rv_veg(:,:,:)
!  Momentum terms feed to the turbulence model 
        real(r8), pointer :: ru_loc_veg(:,:,:,:)
        real(r8), pointer :: rv_loc_veg(:,:,:,:)
        real(r8), pointer :: step2d_uveg(:,:)
        real(r8), pointer :: step2d_vveg(:,:)
        real(r8), pointer :: Lveg(:,:,:)
        real(r8), pointer :: bend(:,:,:)
        real(r8), pointer :: tke_veg(:,:,:)
        real(r8), pointer :: gls_veg(:,:,:)
        real(r8), pointer :: dissip_veg(:,:)
        real(r8), pointer :: BWDXL_veg(:,:,:)
        real(r8), pointer :: BWDYL_veg(:,:,:)
        real(r8), pointer :: marsh_mask(:,:)
        real(r8), pointer :: mask_thrust(:,:)
        real(r8), pointer :: Thrust_max(:,:)
        real(r8), pointer :: Thrust_tonelli(:,:) 
      END TYPE T_VEG
      TYPE (T_VEG), allocatable :: VEG(:)
      CONTAINS
      SUBROUTINE allocate_vegarr (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_vegetation 
      implicit none 
!                       
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate structure variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( VEG(Ngrids) )
!
!  Nonlinear model state.
!
      allocate ( VEG(ng) % plant(LBi:UBi,LBj:UBj,NVEG,NVEGP) )
      allocate ( VEG(ng) % ru_veg(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( VEG(ng) % rv_veg(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( VEG(ng) % ru_loc_veg(LBi:UBi,LBj:UBj,N(ng),NVEG) )
      allocate ( VEG(ng) % rv_loc_veg(LBi:UBi,LBj:UBj,N(ng),NVEG) )
      allocate ( VEG(ng) % step2d_uveg(LBi:UBi,LBj:UBj) )
      allocate ( VEG(ng) % step2d_vveg(LBi:UBi,LBj:UBj) ) 
      allocate ( VEG(ng) % Lveg(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( VEG(ng) % bend(LBi:UBi,LBj:UBj,NVEG) )
      allocate ( VEG(ng) % tke_veg(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( VEG(ng) % gls_veg(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( VEG(ng) % dissip_veg(LBi:UBi,LBj:UBj) )
      allocate ( VEG(ng) % BWDXL_veg(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( VEG(ng) % BWDYL_veg(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( VEG(ng) % marsh_mask(LBi:UBi,LBj:UBj) )
      allocate ( VEG(ng) % mask_thrust(LBi:UBi,LBj:UBj) )
      allocate ( VEG(ng) % Thrust_max(LBi:UBi,LBj:UBj) )
      allocate ( VEG(ng) % Thrust_tonelli(LBi:UBi,LBj:UBj) )
!
!-----------------------------------------------------------------------
!  Allocate various input variables for vegetation module.
!-----------------------------------------------------------------------
!
      RETURN
      END SUBROUTINE allocate_vegarr
      SUBROUTINE initialize_vegarr (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize structure variables in the module using     !
!  first touch distribution policy. In shared-memory configuration,    !
!  this operation actually performs the propagation of the "shared     !
!  arrays" across the cluster,  unless another policy is specified     !
!  to  override the default.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_vegetation 
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k, iveg, ivpr
!
      real(r8), parameter :: IniVal = 0.0_r8
!
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
!  Initialize vegetation structure variables.
!-----------------------------------------------------------------------
!
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO ivpr=1,NVEGP
          DO iveg=1,NVEG
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                VEG(ng) % plant(i,j,iveg,ivpr) = IniVal
              END DO
            END DO
          END DO 
        END DO
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              VEG(ng) % ru_veg(i,j,k) = IniVal
              VEG(ng) % rv_veg(i,j,k) = IniVal
            END DO 
          END DO 
        END DO 
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              VEG(ng) % Lveg(i,j,k) = IniVal
            END DO 
          END DO 
        END DO 
        DO iveg=1,NVEG
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                VEG(ng) % ru_loc_veg(i,j,k,iveg) = IniVal
                VEG(ng) % rv_loc_veg(i,j,k,iveg) = IniVal
              END DO 
            END DO 
          END DO 
	END DO 
	DO j=Jmin,Jmax
	  DO i=Imin,Imax
            VEG(ng) % step2d_uveg(i,j) = IniVal
            VEG(ng) % step2d_vveg(i,j) = IniVal
	  END DO 
	END DO 
        DO iveg=1,NVEG
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              VEG(ng) % bend(i,j,iveg) = IniVal
            END DO 
          END DO 
        END DO 
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              VEG(ng) % tke_veg(i,j,k) = IniVal
              VEG(ng) % gls_veg(i,j,k) = IniVal
            END DO 
          END DO
        END DO 
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            VEG(ng) % dissip_veg(i,j) = IniVal
          END DO 
        END DO 
        DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              VEG(ng) % BWDXL_veg(i,j,k) = IniVal
              VEG(ng) % BWDYL_veg(i,j,k) = IniVal
            END DO 
          END DO
        END DO 
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            VEG(ng) % marsh_mask(i,j) = IniVal
            VEG(ng) % mask_thrust(i,j) = IniVal
            VEG(ng) % Thrust_max(i,j) = IniVal
            VEG(ng) % Thrust_tonelli(i,j) = IniVal
          END DO 
        END DO
!
      END IF
! 
      RETURN   
      END SUBROUTINE initialize_vegarr
      END MODULE mod_vegarr
