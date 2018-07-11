      PROGRAM mct_coupler
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!==================================================== John C. Warner ===
!                                                                      !
!  Master program to couple ROMS/TOMS to other models using the Model  !
!  Coupling Toolkit (MCT) library.                                     !
!                                                                      !
!  The following models are coupled to ROMS/TOMS:                      !
!                                                                      !
!  SWAN, Simulating WAves Nearshore model:                             !
!        http://vlm089.citg.tudelft.nl/swan/index.htm                  !
!                                                                      !
!=======================================================================
!
      USE mod_iounits
      USE mod_scalars
      USE swan_iounits
      USE mct_coupler_params
      USE mod_coupler_iounits
!
      USE m_MCTWorld, ONLY : MCTWorld_clean => clean
      USE ocean_control_mod, ONLY : ROMS_initialize
      USE ocean_control_mod, ONLY : ROMS_run
      USE ocean_control_mod, ONLY : ROMS_finalize
      USE waves_control_mod, ONLY : SWAN_driver_init
      USE waves_control_mod, ONLY : SWAN_driver_run
      USE waves_control_mod, ONLY : SWAN_driver_finalize
      USE ocean_coupler_mod, ONLY : finalize_ocn2wav_coupling
!
      implicit none
      include 'mpif.h'
!
!  Local variable declarations.
!
      logical, save :: first
      integer :: MyColor, MyCOMM, MyError, MyKey, Nnodes
      integer :: MyRank, pelast
      integer :: Ocncolor, Wavcolor, Atmcolor
      integer :: ng, iw, io, ia, icc, lcm, gcdlcm
      real(m4) :: CouplingTime             ! single precision
!
!-----------------------------------------------------------------------
!  Initialize distributed-memory (1) configuration
!-----------------------------------------------------------------------
!
!  Initialize 1 execution environment.
! 
      CALL mpi_init (MyError)
!
!  Get rank of the local process in the group associated with the
!  comminicator.
!
      CALL mpi_comm_size (MPI_COMM_WORLD, Nnodes, MyError)
      CALL mpi_comm_rank (MPI_COMM_WORLD, MyRank, MyError)
!
!  Read in coupled model parameters from standard input.
!
      CALL read_coawst_par(1)
!
!  Now that we know the input file names and locations for each model,
!  for each model read in the number of grids and the grid time steps.
!
      CALL read_model_inputs
!
      CALL allocate_coupler_params
!
!
!  Compute the mct send and recv instances.
!
!  For each model grid, determine the number of steps it should
!  compute before it sends data out.
!  For example, nWAV2OCN(1,2) is the number of steps the wave model 
!  grid 1 should take before it sends data to the ocn grid 2.
!
      DO iw=1,Nwav_grids
        DO io=1,Nocn_grids
          lcm=gcdlcm(dtwav(iw),dtocn(io))
          IF (MOD(TI_WAV2OCN,REAL(lcm,m8)).eq.0) THEN
            nWAV2OCN(iw,io)=INT(TI_WAV2OCN/dtwav(iw))
          ELSE
            lcm=gcdlcm(TI_WAV2OCN,REAL(lcm,m8))
            nWAV2OCN(iw,io)=INT(REAL(lcm,m8)/dtwav(iw))
          END IF
        END DO
      END DO
!
      DO io=1,Nocn_grids
        DO iw=1,Nwav_grids
          lcm=gcdlcm(dtwav(iw),dtocn(io))
          IF (MOD(TI_OCN2WAV,REAL(lcm,m8)).eq.0) THEN
            nOCN2WAV(io,iw)=INT(TI_OCN2WAV/dtocn(io))
          ELSE
            lcm=gcdlcm(TI_OCN2WAV,REAL(lcm,m8))
            nOCN2WAV(io,iw)=INT(REAL(lcm,m8)/dtocn(io))
          END IF
        END DO
      END DO
!
!  Similarly, for each model grid, determine the number of steps 
!  it should compute before it recvs data from somewhere.
!  For example, nWAVFOCN(1,2) is the number of steps the wave model 
!  grid 1 should take before it gets data from ocn grid 2.
!
      DO iw=1,Nwav_grids
        DO io=1,Nocn_grids
          lcm=gcdlcm(dtwav(iw),dtocn(io))
          IF (MOD(TI_OCN2WAV,REAL(lcm,m8)).eq.0) THEN
            nWAVFOCN(iw,io)=INT(TI_OCN2WAV/dtwav(iw))
          ELSE
            lcm=gcdlcm(TI_OCN2WAV,REAL(lcm,m8))
            nWAVFOCN(iw,io)=INT(REAL(lcm,m8)/dtwav(iw))
          END IF
        END DO
      END DO
!
      DO io=1,Nocn_grids
        DO iw=1,Nwav_grids
          lcm=gcdlcm(dtwav(iw),dtocn(io))
          IF (MOD(TI_WAV2OCN,REAL(lcm,m8)).eq.0) THEN
            nOCNFWAV(io,iw)=INT(TI_WAV2OCN/dtocn(io))
          ELSE
            lcm=gcdlcm(TI_WAV2OCN,REAL(lcm,m8))
            nOCNFWAV(io,iw)=INT(REAL(lcm,m8)/dtocn(io))
          END IF
        END DO
      END DO
!
!  Allocate several coupling variables.
!
      allocate(ocnids(Nocn_grids))
      allocate(wavids(Nwav_grids))
!
      N_mctmodels=0
      DO ng=1,Nocn_grids
        N_mctmodels=N_mctmodels+1
        ocnids(ng)=N_mctmodels
      END DO
      DO ng=1,Nwav_grids
        N_mctmodels=N_mctmodels+1
        wavids(ng)=N_mctmodels
      END DO
!
!  Assign processors to the models.
!
      pelast=-1
      peOCN_frst=pelast+1
      peOCN_last=peOCN_frst+NnodesOCN-1
      pelast=peOCN_last
      peWAV_frst=pelast+1
      peWAV_last=peWAV_frst+NnodesWAV-1
      pelast=peWAV_last
      IF (pelast.ne.Nnodes-1) THEN
        IF (MyRank.eq.0) THEN
          WRITE (stdout,10) pelast+1, Nnodes
 10       FORMAT (/,' mct_coupler - Number assigned processors: '       &
     &            ,i3.3,/,15x,'not equal to spawned MPI nodes: ',i3.3)
        END IF
        STOP
      ELSE
        IF (MyRank.eq.0) THEN
          WRITE (stdout,19)
 19       FORMAT (/,' Model Coupling: ',/)
          WRITE (stdout,20) peOCN_frst, peOCN_last
 20       FORMAT (/,7x,'Ocean Model MPI nodes: ',i3.3,' - ', i3.3)
          WRITE (stdout,21) peWAV_frst, peWAV_last
 21       FORMAT (/,7x,'Waves Model MPI nodes: ',i3.3,' - ', i3.3)
        END IF
      END IF
!
!  Split the communicator into SWAN, WRF, and ROMS subgroups based 
!  on color and key.
!
      Atmcolor=1
      Ocncolor=2
      Wavcolor=3
      MyKey=0
      IF ((peOCN_frst.le.MyRank).and.(MyRank.le.peOCN_last)) THEN
        MyColor=OCNcolor
      END IF
      IF ((peWAV_frst.le.MyRank).and.(MyRank.le.peWAV_last)) THEN
        MyColor=WAVcolor
      END IF
      CALL mpi_comm_split (MPI_COMM_WORLD, MyColor, MyKey, MyCOMM,      &
     &                     MyError)
!
!-----------------------------------------------------------------------
!  Run coupled models according to the processor rank.
!-----------------------------------------------------------------------
!
      IF (MyColor.eq.WAVcolor) THEN
        CALL SWAN_driver_init (MyCOMM)
        CALL SWAN_driver_run
        CALL SWAN_driver_finalize
      END IF
      IF (MyColor.eq.OCNcolor) THEN
        first=.TRUE.
        Nrun=1
        IF (exit_flag.eq.NoError) THEN
          CALL ROMS_initialize (first, MyCOMM)
        END IF
        IF (exit_flag.eq.NoError) THEN
          run_time=0.0_m8
          DO ng=1,Ngrids
            run_time=MAX(run_time, dt(ng)*ntimes(ng))
          END DO
          CALL ROMS_run (run_time)
        END IF
        CALL ROMS_finalize
        CALL finalize_ocn2wav_coupling
      END IF
!
!-----------------------------------------------------------------------
!  Terminates all the mpi-processing and coupling.
!-----------------------------------------------------------------------
!
      CALL mpi_barrier (MPI_COMM_WORLD, MyError)
      CALL MCTWorld_clean ()
      CALL mpi_finalize (MyError)
      STOP
      END PROGRAM mct_coupler
      FUNCTION gcdlcm (dtAin, dtBin)
!
!=======================================================================
!                                                                      !
!  This function computes the greatest common denominator              !
!  and lowest common multiple.                                         !
!                                                                      !
!  On Input:                                                           !
!     dtA        time step of model A                                  !
!     dtB        time step of model B                                  !
!                                                                      !
!  On Output:                                                          !
!     lcm        least common multiple                                 !
!                                                                      !
!=======================================================================
!
      USE mod_coupler_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      real(m8), intent(in) :: dtAin, dtBin
      integer :: gcdlcm
!
!  Local variable declarations.
!
      logical :: stayin
      real(m8) :: r, m, n, p, gcd, dtA, dtB
!
!-----------------------------------------------------------------------
!  Compute greatest common denominator and least common multiplier.
!-----------------------------------------------------------------------
      dtA=dtAin
      dtB=dtBin
      m=dtA
      n=dtB
      IF (dtA.gt.dtB) THEN
        p=dtA
        dtA=dtB
        dtB=p
      END IF
      stayin=.true.
      DO WHILE (stayin)
        r=mod(dtB,dtA)
        IF (r.eq.0) THEN
          gcd=dtA
          stayin=.false.
        ELSE
          dtB=dtA
          dtA=r
        END IF
      END DO
      gcdlcm=m*n/dtA
      RETURN
      END FUNCTION gcdlcm
