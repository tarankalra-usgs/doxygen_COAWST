       MODULE vegetation_stream_mod
!
!svn $Id: vegetation_stream.F 429 2015-04-20 17:30:26Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!================================================John C. Warner=========
!================================================Neil K. Ganju  ========
!================================================Alexis Beudin  ========
!==============================================Tarandeep S. Kalra=======
!                                                                      ! 
!  References:                                                         !   
!                                                                      !
!=======================================================================
!                                                                      !
!=======================================================================
      implicit none
      PRIVATE
      PUBLIC  :: vegetation_stream_cal
      CONTAINS
!
!***********************************************************************
      SUBROUTINE vegetation_stream_cal (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_vegarr
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
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
!
      CALL wclock_on (ng, iNLM, 16)
      CALL vegetation_stream_tile  (ng, tile,                           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                         GRID(ng) % angler,                       &
     &                         GRID(ng) % z_w,                          &
     &                         FORCES(ng) % Dwave,                      &
     &                         FORCES(ng) % Lwave,                      &
     &                         VEG(ng) % dissip_veg,                    &
     &                         VEG(ng) % Lveg,                          &
     &                         VEG(ng) % BWDXL_veg,                     &
     &                         VEG(ng) % BWDYL_veg)
      CALL wclock_off (ng, iNLM, 16)
      RETURN
      END SUBROUTINE vegetation_stream_cal 
!***********************************************************************
      SUBROUTINE vegetation_stream_tile  (ng, tile,                     &
     &                              LBi, UBi, LBj, UBj,                 &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                              angler,                             &
     &                              z_w,                                &
     &                              Dwave,                              &
     &                              Lwave,                              &
     &                              dissip_veg, Lveg,                   &
     &                              BWDXL_veg, BWDYL_veg)  
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars
      USE mod_vegetation
      USE mod_vegarr
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
      real(r8), intent(in) :: dissip_veg(LBi:,LBj:)
      real(r8), intent(in) :: Lveg(LBi:,LBj:,:)
      real(r8), intent(inout) :: BWDXL_veg(LBi:,LBj:,:)
      real(r8), intent(inout) :: BWDYL_veg(LBi:,LBj:,:)
!  Local variable declarations.
!
      integer :: i, j, k, iveg
      real(r8) :: cff1, cff2
      real(r8) :: EWD_veg
      real(r8), parameter :: Lwave_min = 1.0_r8
      real(r8) :: Dstp
      real(r8) :: waven, wavenx, waveny
      real(r8) :: sigma, osigma
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
!----------------------------------------------------------------------
!----------Executing the code------------------------------------------
!----------------------------------------------------------------------
!     
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
            Dstp=z_w(i,j,N(ng))-z_w(i,j,0)
!----------------------------------------------------------------------
!  Compute wave amplitude (0.5*Hrms), wave number, intrinsic frequency.
!----------------------------------------------------------------------
            waven=2.0_r8*pi/MAX(Lwave(i,j),Lwave_min)
            cff1=1.5_r8*pi-Dwave(i,j)-angler(i,j)
            wavenx=waven*COS(cff1)
            waveny=waven*SIN(cff1)
            sigma=MIN(SQRT(g*waven*TANH(waven*Dstp)),2.0_r8)
            osigma=1.0_r8/sigma
! 
!----------------------------------------------------------------------
!   Note: Alexis - check if we need a local dissip_veg here 
!   Also Lveg is for 1 veg type only 
!----------------------------------------------------------------------
!
            EWD_veg=dissip_veg(i,j)
            cff2=EWD_veg*osigma*Lveg(i,j,k)
            BWDXL_veg(i,j,k)=cff2*wavenx
            BWDYL_veg(i,j,k)=cff2*waveny
!           
           END DO
        END DO
      END DO    
!
      END SUBROUTINE vegetation_stream_tile
      END MODULE vegetation_stream_mod
