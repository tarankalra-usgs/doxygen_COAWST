      MODULE mct_coupler_params
      USE mod_coupler_kinds
      implicit none
!
!  Number of coupling models.
!
      integer :: N_mctmodels
!
!  Number of parallel nodes assigned to each model in the coupled
!  system.
!
      integer :: NnodesATM
      integer :: NnodesWAV
      integer :: NnodesOCN
!
!  Parallel nodes assined to the wave model.
!
        integer :: peWAV_frst          ! first atmosphere parallel node
        integer :: peWAV_last          ! last  atmosphere parallel node
!
!  Parallel nodes assined to the ocean model.
!
        integer :: peOCN_frst          ! first ocean parallel node
        integer :: peOCN_last          ! last  ocean parallel node
        integer, dimension(:), pointer :: roms_fwcoup
        integer, dimension(:), pointer :: roms_2wcoup
        integer, dimension(:), pointer :: roms_facoup
        integer, dimension(:), pointer :: roms_2acoup
!
!  Time interval (seconds) between coupling of models.
!
      real(m8) :: TI_ATM2WAV           ! atmosphere to wave coupling interval
      real(m8) :: TI_ATM2OCN           ! atmosphere to ocean coupling interval
      real(m8) :: TI_WAV2ATM           ! wave to atmosphere coupling interval
      real(m8) :: TI_WAV2OCN           ! wave to ocean coupling interval
      real(m8) :: TI_OCN2WAV           ! ocean to wave coupling interval
      real(m8) :: TI_OCN2ATM           ! ocean to atmosphere coupling interval
!
!  Number of atmosphere model time-steps and atmosphere model ID.
!
      integer :: Natm_grids
      integer :: Nocn_grids
      integer :: Nwav_grids
      real(m8), dimension(:), pointer :: dtocn
      real(m8), dimension(:), pointer :: dtwav
      real(m8), dimension(:), pointer :: dtatm
!
      integer, dimension(:,:), pointer :: nOCN2WAV
      integer, dimension(:,:), pointer :: nWAV2OCN
      integer, dimension(:,:), pointer :: nOCNFWAV
      integer, dimension(:,:), pointer :: nWAVFOCN
!
!  Coupled model components IDs.
!
      integer, dimension(:), pointer :: ocnids
      integer, dimension(:), pointer :: wavids
      integer, dimension(:), pointer :: atmids
      integer :: OCNid
      integer :: WAVid
      integer :: ATMid
      CONTAINS
      SUBROUTINE allocate_coupler_params
!=======================================================================
!                                                                      !
!  This routine initialize all the coupler vars.                       !
!                                                                      !
!=======================================================================
      allocate (nOCN2WAV(Nocn_grids,Nwav_grids))
      allocate (nWAV2OCN(Nwav_grids,Nocn_grids))
      allocate (nOCNFWAV(Nocn_grids,Nwav_grids))
      allocate (nWAVFOCN(Nwav_grids,Nocn_grids))
      RETURN
      END SUBROUTINE allocate_coupler_params
      END MODULE mct_coupler_params
