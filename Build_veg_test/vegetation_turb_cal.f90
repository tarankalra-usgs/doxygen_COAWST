      MODULE vegetation_turb_mod
!
!svn $Id: vegetation_turb_cal.F 429 2015-06-10 12:30:26Z arango $
!=======================================================================
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                                              !
!==================================================== John C. Warner ===
!==================================================== Neil K. Ganju  ===
!==================================================== Alexis Beudin  ===
!==================================================Tarandeep S. Kalra===
!                                                                      !
!  This routine computes the turbulent kinetic energy and length scale !
!  modifications due to vegetation for gls_corstep.F                   !
!                                                                      !
!  References:                                                         !
!                                                                      !
!   Uittenbogaard R. (2003): Modelling turbulence in vegetated aquatic !
!   flows. International workshop on RIParian FORest vegetated         !
!   channels: hydraulic, morphological and ecological aspects,         !
!   20-22 February 2003, Trento, Italy.                                !
!                                                                      !
!   Warner J.C., C.R. Sherwood, H.G. Arango, and R.P. Signell (2005):  !
!   Performance of four turbulence closure models implemented using a  !
!   generic length scale method, Ocean Modelling 8: 81-113.            !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: vegetation_turb_cal
      CONTAINS
!
!***********************************************************************
      SUBROUTINE vegetation_turb_cal (ng, tile)
!***********************************************************************
!
      USE mod_stepping 
      USE mod_grid
      USE mod_ocean
      USE mod_mixing 
      USE mod_vegarr
      USE vegetation_drag_mod
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
      CALL vegetation_turb_tile  ( ng, tile,                            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nstp(ng), nnew(ng),                       &
     &                        OCEAN(ng) % u,                            &
     &                        OCEAN(ng) % v,                            &
     &                        VEG(ng) % ru_loc_veg,                     &
     &                        VEG(ng) % rv_loc_veg,                     &
     &                        VEG(ng) % plant,                          &
     &                        VEG(ng) % bend,                           &
     &                        MIXING(ng) % gls,                         &
     &                        MIXING(ng) % tke,                         &
     &                        VEG(ng) % gls_veg,                        &
     &                        VEG(ng) % tke_veg )
      CALL wclock_off (ng, iNLM, 16)
      RETURN
      END SUBROUTINE vegetation_turb_cal 
!
!***********************************************************************
      SUBROUTINE vegetation_turb_tile ( ng, tile,                       &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nstp, nnew,                         &
     &                              u, v,                               &
     &                              ru_loc_veg, rv_loc_veg,             &
     &                              plant,                              &
     &                              bend,                               &
     &                              gls, tke,                           &
     &                              gls_veg, tke_veg )
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_vegetation 
      USE mod_vegarr
      USE vegetation_drag_mod
!
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew 
!
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ru_loc_veg(LBi:,LBj:,:,:)
      real(r8), intent(in) :: rv_loc_veg(LBi:,LBj:,:,:)
      real(r8), intent(in) :: plant(LBi:,LBj:,:,:)
      real(r8), intent(in) :: bend(LBi:,LBj:,:)
      real(r8), intent(in) :: gls(LBi:,LBj:,0:,:)
      real(r8), intent(in) :: tke(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: gls_veg(LBi:,LBj:,0:)
      real(r8), intent(inout) :: tke_veg(LBi:,LBj:,0:)
!
!  Local variable declarations.
!
      integer :: i, j, k, iveg
!
      real(r8), parameter :: one_half=1.0_r8/2.0_r8 
      real(r8), parameter :: one_third=1.0_r8/3.0_r8 
      real(r8), parameter :: Inival=0.0_r8
      real(r8), parameter :: cl_veg=1.0_r8, ck=0.09_r8
      real(r8), parameter :: max_L=10.0e10_r8 
      real(r8), parameter :: eps=1.0e-12_r8 
      real(r8) :: wrku1, wrku2, wrku3, wrku4, wrku
      real(r8) :: wrkv1, wrkv2, wrkv3, wrkv4, wrkv
      real(r8) :: wrk, cff1, cff2, cff3, dissip, inverse_dissip
      real(r8) :: solid, L, eqvegT
      real(r8) :: taufree, tauveg, taueff
      real(r8) :: tke_loc_veg, gls_loc_veg
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vegu
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vegv 
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
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
            tke_veg(i,j,k)=Inival
            gls_veg(i,j,k)=Inival
          END DO 
        END DO 
      END DO 
      cff1=3.0_r8+gls_p(ng)/gls_n(ng)
      cff2=1.5_r8+gls_m(ng)/gls_n(ng)
      cff3=-1.0_r8/gls_n(ng)
      VEG_LOOP: DO iveg=1,NVEG
        DO k=1,N(ng)-1
          DO j=Jstr,Jend
            DO i=Istr,Iend
!
!-----------------------------------------------------------------------
! Additional turbulence generated by the vegetation = 
! work spent by the fluid against the plants (in m3/s3)
!-----------------------------------------------------------------------
!
              wrku1=ru_loc_veg(i,j,k,iveg)*u(i,j,k,nstp)
              wrku2=ru_loc_veg(i,j,k+1,iveg)*u(i,j,k+1,nstp)
              wrku3=ru_loc_veg(i+1,j,k,iveg)*u(i+1,j,k,nstp)
              wrku4=ru_loc_veg(i+1,j,k+1,iveg)*u(i+1,j,k+1,nstp)
              wrku=0.25_r8*(wrku1+wrku2+wrku3+wrku4)
              wrkv1=rv_loc_veg(i,j,k,iveg)*v(i,j,k,nstp)
              wrkv2=rv_loc_veg(i,j,k+1,iveg)*v(i,j,k+1,nstp)
              wrkv3=rv_loc_veg(i,j+1,k,iveg)*v(i,j+1,k,nstp)
              wrkv4=rv_loc_veg(i,j+1,k+1,iveg)*v(i,j+1,k+1,nstp)
              wrkv=0.25_r8*(wrkv1+wrkv2+wrkv3+wrkv4)
              tke_loc_veg=sqrt(wrku*wrku+wrkv*wrkv)
!
!-----------------------------------------------------------------------
! Dissipation due to vegetation
!-----------------------------------------------------------------------
! Dissipation in GLS (Eq. 12 in Warner et al., 2005)
!
              wrk=MAX(tke(i,j,k,nstp),gls_Kmin(ng))
              dissip=(gls_cmu0(ng)**cff1)*(wrk**cff2)*                  &
     &                 (gls(i,j,k,nstp)**cff3)
              inverse_dissip=1.0_r8/MAX(dissip,eps)
!
! Dissipation time-scale for free turbulence
!
              taufree=wrk*inverse_dissip
!
!
! Equivalent thickness: horizontal projection of the bending plant
! 
!              eqvegT=plant(i,j,iveg,pthck)+sin(bend(i,j,iveg))*         &
!     &                                       plant(i,j,iveg,phght)
              eqvegT=plant(i,j,iveg,pthck)
!
!
! Solidity:cross-sectional area of a plant the number of plants per m2
!
!
              solid=plant(i,j,iveg,pdiam)*eqvegT*plant(i,j,iveg,pdens)
!
! Eddies typical size constrained by distance in between the plants
!
              L=cl_veg*((1.0_r8-MIN(solid,1.0_r8))/                     &
     &                 plant(i,j,iveg,pdens))**one_half
              L=MIN(L,max_L)
!
! Dissipation time-scale of eddies in between the plants
!
              tauveg=(L**2.0_r8/(ck**2.0_r8*tke_loc_veg))**one_third
!
! Effective dissipation time-scale
! 
              taueff=MIN(taufree,tauveg)
              gls_loc_veg=gls_c2(ng)*tke_loc_veg/taueff
!
!-----------------------------------------------------------------------
! Add the tke and gls changes from all vegetation types
!-----------------------------------------------------------------------
! 
              tke_veg(i,j,k)=tke_loc_veg + tke_veg(i,j,k)
              gls_veg(i,j,k)=gls_loc_veg + gls_veg(i,j,k)
            END DO 
          END DO
        END DO 
      END DO VEG_LOOP
!
      RETURN
      END SUBROUTINE vegetation_turb_tile
      END MODULE vegetation_turb_mod
