subroutine saconly(n_hrus, dt, sim_length, year, month, day, hour, &
    sac_pars, &
    peadj, pxadj, &
    map, etd, &
    return_states, &
    tci, aet, uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc, &
    roimp, sdro, ssur, sif, bfs, bfp )

! !     Subroutine Description
! !     -----------------------------------
! !     The sacsnow subroutine is a wrapper to run SNOW17 and
! !     SAC-SMA models, returning total channel inflow and
! !     optionally model states
! !
! !     Arguments
! !     -----------------------------------
! !     INPUTS
! !     n_hrus:  Number of zones(integer)
! !     dt:  model timestep in seconds (integer)
! !     sim_length: length of simulation in days (integer)
! !     year:  The year associated with each time step (integer array)
! !     month:  The month associated with each time step (integer array)
! !     day:  The day associated with each time step (integer array)
! !     hour:  The hour associated with each time step (integer array)
! !     sac_pars: (double array): SAC-SMA parameters for each zone.  See line below for order of parameters (double array)
! !     uztwm, uzfwm, lztwm, lzfpm, lzfsm, adimp, uzk, lzpk, lzsk, zperc, rexp, pctim, pfree, riva, side, rserv, efc
! !     peadj, pxadj: (double array):  zone specific etd (peadj) and map (pxadj) multiplication factor
! !     init_swe: (double array):  zone specific initial SWE in mm (double array)
! !     scf, mfmax, mfmin, uadj, si, nmf, tipm, mbase, plwhc, daygm, adc_a, adc_b, adc_c,
! !     map, ptps, mat,etd: (double array):  precipitation as mm, precent precipitation as snow as decimal,
! !                                          temperature as DegC, and evaporation demand as mm (double array)
! !     return_states: option to return SNOW17 and SAC-SMA states array or return TCI only(logical)
! !     OUTPUTS
! !     tci:  total channel inflow for each zone in mm (double array)
! !     aet:  actual evapotranspiration for each zone as mm (double array)
! !     uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc:  SAC-SMA content or state for each zone as mm (double array)
! !     roimp, sdro, ssur, sif, bfs, bfp:  tci contribution from each SAC-SMA runoff source form each zone
! !     swe, aesc, neghs, liqw, raim, psfall, prain:  Snow17 state for each zone (double array)

    ! ! zone info 
    ! ! sac-sma params in a matrix, see the variable declaration
    ! sac_pars, &
    ! ! zone specific etd (peadj) and map (pxadj) adjustments 
    ! peadj, pxadj, &
    ! ! snow17 params in a matrix, see the variable declaration
    ! ! initial state value for swe 
    ! init_swe, & 
    ! ! forcings 
    ! map, ptps, mat, etd, &
    ! ! outputs
    ! tci, aet, uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc, &
    ! swe, aesc, neghs, liqw, raim, psfall, prain)

  use utilities

  implicit none

  double precision, parameter:: pi=3.141592653589793238462643383279502884197d0
  double precision, parameter:: sec_day = 86400.     !seconds in a day
  double precision, parameter:: sec_hour = 3600.     !seconds in an hour
  integer, parameter:: sp = KIND(1.0)
  integer:: k

  logical:: return_states

  integer, intent(in):: n_hrus ! number of zones

  ! sac pars matrix 
  ! uztwm, uzfwm, lztwm, lzfpm, lzfsm, adimp, uzk, lzpk, lzsk, zperc, rexp, pctim, pfree, riva, side, rserv, efc
  double precision, dimension(17,n_hrus), intent(in):: sac_pars 

  ! snow pars marix 
  ! scf, mfmax, mfmin, uadj, si, nmf, tipm, mbase, plwhc, daygm, adc_a, adc_b, adc_c

  ! this code is currently not set up to do any timestep less than 1 hour, 
  ! nor could it do fractional hour timesteps.
  integer, intent(in):: dt    ! model timestep in seconds
  integer:: dt_hours          ! model timestep in hours
  integer:: ts_per_day, ts_per_year

  ! initial states for a cold start run
  ! used in all model HRUs
  ! model state variables not listed start at 0
  double precision, dimension(6):: spin_up_start_states, spin_up_end_states
  integer:: spin_up_counter, spin_up_max_iter
  double precision:: pdiff
  double precision, dimension(n_hrus):: init_uztwc, init_uzfwc, init_lztwc, init_lzfsc, &
          init_lzfpc, init_adimc

  ! SAC_model params & other key inputs in the sace param file
  !character(len = 20), dimension(n_hrus) :: hru_id   ! local hru id
  double precision, dimension(n_hrus):: uztwm, uzfwm, uzk, pctim, adimp, zperc, rexp, &
                                lztwm, lzfsm, lzfpm, lzsk, lzpk, pfree, &
                                riva, side, rserv, efc, peadj, pxadj

  ! local variables
  integer:: nh,i          ! AWW index for looping through areas
  integer:: sim_length   ! length of simulation (days)

  ! single precision sac-sma and snow variables
  ! these are single precision so as to be supplied to 
  ! NWS f77 models 
  real(sp):: uztwc_sp, uzfwc_sp, lztwc_sp, lzfsc_sp, lzfpc_sp, adimc_sp
  real(sp):: roimp_sp, sdro_sp, ssur_sp, sif_sp, bfs_sp, bfp_sp

  ! sac-sma state variables
  double precision, dimension(sim_length ,n_hrus), intent(out):: uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc
  double precision, dimension(sim_length ,n_hrus), intent(out):: roimp, sdro, ssur, sif, bfs, bfp

  ! sac-sma output variables and channel inflow
  real(sp):: aet_sp, tci_sp
  double precision, dimension(sim_length ,n_hrus), intent(out):: tci, aet

  ! date variables
  integer, dimension(sim_length), intent(in):: year, month, day, hour
  integer:: houri

  ! atmospheric forcing variables
  !f2py intent(in,out) map, etd
  double precision, dimension(sim_length, n_hrus), intent(inout):: map, etd
  double precision:: map_step, etd_step

  ! initilize outputs 
  tci = 0
  if(return_states)then
    aet = 0 
    uztwc = 0 
    uzfwc = 0 
    lztwc = 0 
    lzfsc = 0 
    lzfpc = 0 
    adimc = 0 
    roimp = 0
    sdro = 0
    ssur = 0
    sif = 0 
    bfs = 0 
    bfp =0
  end if

  ! pull out sac params to separate variables
  uztwm = sac_pars(1,:)
  uzfwm = sac_pars(2,:)
  lztwm = sac_pars(3,:)
  lzfpm = sac_pars(4,:)
  lzfsm = sac_pars(5,:)
  adimp = sac_pars(6,:)
    uzk = sac_pars(7,:)
   lzpk = sac_pars(8,:)
   lzsk = sac_pars(9,:)
  zperc = sac_pars(10,:)
   rexp = sac_pars(11,:)
  pctim = sac_pars(12,:)
  pfree = sac_pars(13,:)
   riva = sac_pars(14,:)
   side = sac_pars(15,:)
  rserv = sac_pars(16,:)
    efc = sac_pars(17,:)

  ts_per_day = 86400/dt
  dt_hours = dt/3600
  
  ! write(*,*)'Timesteps per day:',ts_per_day

  ! ========================= ZONE AREA LOOP ========================================================
  !   loop through the zones, running the lumped model code for each

  do nh=1,n_hrus
    print*, 'Running area',nh,'out of',n_hrus

    !print run dates
    write(*,*)'  start:',year(1), month(1), day(1), hour(1)
    write(*,*)'    end:',year(sim_length), month(sim_length), day(sim_length), hour(sim_length)

    ! set the areal depletion curve based on parameters ax^b+(1-a)x^c
    ! 0 < a < 1; b, c > 0 
    ! if b < 1 & c < 1 curve is concave up
    ! if b > 1 & c > 1 curve is concave up
    ! if b < 1 & c > 1 OR b > 1 & c < 1 curve is s-shaped
    ! "A value of As = 0.05 is used for a W/Ai = 0.0 ratio so that small amounts of snow
    ! don't continue to exist well past the time when all the snow is gone in nature."
    ! - snow 17 manual
    ! adc = adc_a(nh)*adc_y**adc_b(nh)+(1.-adc_a(nh))*adc_y**adc_c(nh)
    ! do i=1,11
    !   if(adc(i) < 0.05) adc(i) = 0.05
    ! end do
    
    ! =============== Spin up procedure =====================================

    ! starting values
    spin_up_start_states = 1d0 
    spin_up_end_states = 0d0
    pdiff = 1d0
    ts_per_year = ts_per_day * 365
    spin_up_counter = 0
    !spin_up_counter = 50
    spin_up_max_iter = 50

    do while (pdiff > 0.01 .and. spin_up_counter < spin_up_max_iter)

      spin_up_counter = spin_up_counter + 1

      ! put the ending states from the previous iteration as the starting states 
      uztwc_sp = real(spin_up_end_states(1))
      uzfwc_sp = real(spin_up_end_states(2))
      lztwc_sp = real(spin_up_end_states(3))
      lzfsc_sp = real(spin_up_end_states(4))
      lzfpc_sp = real(spin_up_end_states(5))
      adimc_sp = real(spin_up_end_states(6))

      ! run for 1 year 
      do i = 1,ts_per_year

        ! apply pe and px adjustments (zone-wise) for the current timestep
        map_step = map(i,nh) * pxadj(nh)
        etd_step = etd(i,nh) * peadj(nh)

!        call exsnow19(int(dt/sec_hour,4),int(day(i),4),int(month(i),4),int(year(i),4),&
!            !SNOW17 INPUT AND OUTPUT VARIABLES
!            real(map_step), real(ptps(i,nh)), real(mat(i,nh)), &
!            raim_sp, sneqv_sp, snow_sp, snowh_sp, psfall_sp, prain_sp, aesc_sp,&
!            !SNOW17 PARAMETERS
!            !ALAT,SCF,MFMAX,MFMIN,UADJ,SI,NMF,TIPM,MBASE,PXTEMP,PLWHC,DAYGM,ELEV,PA,ADC
!            real(latitude(nh)), real(scf(nh)), real(mfmax(nh)), real(mfmin(nh)), &
!            real(uadj(nh)), real(si(nh)), real(nmf(nh)), &
!            real(tipm(nh)), real(mbase(nh)), real(pxtemp(nh)), real(plwhc(nh)), real(daygm(nh)),&
!            real(elev(nh)), real(pa), real(adc_x), &
!            !SNOW17 CARRYOVER VARIABLES
!            cs, taprev_sp) 

        ! modify ET demand using the effective forest cover 
        ! Anderson calb manual pdf page 232
        !etd_step = efc(nh)*etd_step+(1d0-efc(nh))*(1d0-dble(aesc_sp))*etd_step
    
        ! write(*,*) "call exsac:", "map_step=", map_step, "etd_step = ", etd_step
        call exsac(real(dt), map_step, real(etd_step), &
            !SAC PARAMETERS
            !UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC, &
            !REXP,LZTWM,LZFSM,LZFPM,LZSK,LZPK,PFREE, &
            !SIDE,RSERV, &
            real(uztwm(nh)), real(uzfwm(nh)), real(uzk(nh)), real(pctim(nh)), &
            real(adimp(nh)), real(riva(nh)), real(zperc(nh)), &
            real(rexp(nh)), real(lztwm(nh)), real(lzfsm(nh)), real(lzfpm(nh)), &
            real(lzsk(nh)), real(lzpk(nh)), real(pfree(nh)),&
            real(side(nh)), real(rserv(nh)), &
            !SAC State variables
            uztwc_sp, uzfwc_sp, lztwc_sp, lzfsc_sp, lzfpc_sp, adimc_sp, &
            !SAC Runoff variables
            roimp_sp,sdro_sp,ssur_sp,sif_sp,bfs_sp,bfp_sp, &
            !SAC OUTPUTS
            tci_sp, aet_sp)

      end do  ! spin up 1 year loop 

      spin_up_end_states(1) = dble(uztwc_sp)
      spin_up_end_states(2) = dble(uzfwc_sp)
      spin_up_end_states(3) = dble(lztwc_sp)
      spin_up_end_states(4) = dble(lzfsc_sp)
      spin_up_end_states(5) = dble(lzfpc_sp)
      spin_up_end_states(6) = dble(adimc_sp)

      pdiff = 0.0
      do k=1,6
        ! avoid divide by zero 
        if(spin_up_start_states(k) < 0.000001)then
          cycle
        end if
        pdiff = pdiff + abs(spin_up_start_states(k)-spin_up_end_states(k))/spin_up_start_states(k)
      end do
      ! on the first iteration all the states are at zero so 
      ! artificially set pdiff and keep going
      if(spin_up_counter .eq. 1) pdiff = 1.0

      spin_up_start_states = spin_up_end_states

      !write(*,'(7f10.3)')pdiff, spin_up_start_states

    end do 
    ! write(*,*)

    ! Save the spun up states to use for init in the full run
    init_uztwc(nh) = spin_up_end_states(1)
    init_uzfwc(nh) = spin_up_end_states(2)
    init_lztwc(nh) = spin_up_end_states(3)
    init_lzfsc(nh) = spin_up_end_states(4)
    init_lzfpc(nh) = spin_up_end_states(5)
    init_adimc(nh) = spin_up_end_states(6)

    ! =============== End spin up procedure =====================================

    ! set single precision sac state variables to initial values
    uztwc_sp = real(init_uztwc(nh))
    uzfwc_sp = real(init_uzfwc(nh))
    lztwc_sp = real(init_lztwc(nh))
    lzfsc_sp = real(init_lzfsc(nh))
    lzfpc_sp = real(init_lzfpc(nh))
    adimc_sp = real(init_adimc(nh))

!ADIMC=46.04635792983159
!FGIX=0.0
!LZFPC=1.743847376956119
!LZFSC=1.2986819249148157E-4
!LZTWC=30.5294362947644
!UZFWC=0.0
!UZTWC=13.566566157631149

! DSBT2 states
!    uztwc_sp = 13.566566157631149
!    uzfwc_sp = 0.0
!    lztwc_sp = 30.5294362947644
!    lzfsc_sp = 1.2986819249148157E-4
!    lzfpc_sp = 1.743847376956119
!    adimc_sp = 46.04635792983159

! DNGT2 states
!    uztwc_sp = 12.840677368638566
!    uzfwc_sp = 0.005736753812190701
!    lztwc_sp = 40.622529644160586
!    lzfsc_sp = 0.3676078556792409
!    lzfpc_sp = 7.633354536427002
!    adimc_sp = 60.83319949305382
! DNGT2 cold states
!     uztwc_sp = 2.689043543766021
!     uzfwc_sp = 0.008559990638373783
!     lztwc_sp = 17.124192588147444
!     lzfsc_sp = 0.036227178579246656
!     lzfpc_sp = 0.8110268652692039
!     adimc_sp = 17.730247869515264
! DNGT2 debug output states
!     uztwc_sp = 121.21354579420947
!     uzfwc_sp = 0.0
!     lztwc_sp = 0.6240390544095995
!     lzfsc_sp = 0.0
!     lzfpc_sp = 0.0
!     adimc_sp = 122.75349678884632

    ! =============== START SIMULATION TIME LOOP =====================================
    do i = 1,sim_length,1

      ! dummy use of the hour variable to shut the compiler up, 
      ! we may want to use the hour as an input in the future
      if(i .eq. 1) houri = hour(i)

      ! apply adjustments (zone-wise) for the current timestep
      map_step = map(i,nh) * pxadj(nh)
      etd_step = etd(i,nh) * peadj(nh) 

!      call exsnow19(int(dt/sec_hour,4),int(day(i),4),int(month(i),4),int(year(i),4),&
!          !SNOW17 INPUT AND OUTPUT VARIABLES
!          real(map_step), real(ptps(i,nh)), real(mat(i,nh)), &
!          raim_sp, sneqv_sp, snow_sp, snowh_sp, psfall_sp, prain_sp, aesc_sp,&
!          !SNOW17 PARAMETERS
!          !ALAT,SCF,MFMAX,MFMIN,UADJ,SI,NMF,TIPM,MBASE,PXTEMP,PLWHC,DAYGM,ELEV,PA,ADC
!          real(latitude(nh)), real(scf(nh)), real(mfmax(nh)), real(mfmin(nh)), &
!          real(uadj(nh)), real(si(nh)), real(nmf(nh)), &
!          real(tipm(nh)), real(mbase(nh)), real(pxtemp(nh)), real(plwhc(nh)), real(daygm(nh)),&
!          real(elev(nh)), real(pa), real(adc_x), &
!          !SNOW17 CARRYOVER VARIABLES
!          cs, taprev_sp) 


      ! modify ET demand using the effective forest cover 
      ! Anderson calb manual pdf page 232
      !etd_step = efc(nh)*etd_step+(1d0-efc(nh))*(1d0-dble(aesc_sp))*etd_step
  
!      write(*,*) "call exsac ***********************"
!      write(*,*) "uztwm = ", real(uztwm(nh))
!      write(*,*) "uzfwm = ", real(uzfwm(nh))
!      write(*,*) "uzk = ", real(uzk(nh))
!      write(*,*) "pctim = ",  real(pctim(nh))
!      write(*,*) "adimp = ",  real(adimp(nh))
!      write(*,*) "riva = ",  real(riva(nh))
!      write(*,*) "zperc = ", real(zperc(nh))
!      write(*,*) "rexp = ",  real(rexp(nh))
!      write(*,*) "lztwm = ", real(lztwm(nh))
!      write(*,*) "lzfsm = ", real(lzfsm(nh))
!      write(*,*) "lzfpm = ", real(lzfpm(nh))
!      write(*,*) "lzsk = ",  real(lzsk(nh))
!      write(*,*) "lzpk = ",  real(lzpk(nh))
!      write(*,*) "pfree = ",  real(pfree(nh))
!      write(*,*) "side = ",  real(side(nh))
!      write(*,*) "rserv = ",  real(rserv(nh))
!      write(*,*) "uztwc_sp = ",  uztwc_sp
!      write(*,*) "uzfwc_sp = ",  uzfwc_sp
!      write(*,*) "lztwc_sp = ",  lztwc_sp
!      write(*,*) "lzfsc_sp = ",  lzfsc_sp
!      write(*,*) "lzfpc_sp = ",  lzfpc_sp
!      write(*,*) "adimc_sp = ",  adimc_sp
!      write(*,*) "dt = ",  real(dt)
!      write(*,*) "map_step = ",  map_step

      call exsac(real(dt), real(map_step), real(etd_step), &
          !SAC PARAMETERS
          !UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC, &
          !REXP,LZTWM,LZFSM,LZFPM,LZSK,LZPK,PFREE, &
          !SIDE,RSERV, &
          real(uztwm(nh)), real(uzfwm(nh)), real(uzk(nh)), real(pctim(nh)), &
          real(adimp(nh)), real(riva(nh)), real(zperc(nh)), &
          real(rexp(nh)), real(lztwm(nh)), real(lzfsm(nh)), real(lzfpm(nh)), &
          real(lzsk(nh)), real(lzpk(nh)), real(pfree(nh)),&
          real(side(nh)), real(rserv(nh)), &
          !SAC State variables
          uztwc_sp, uzfwc_sp, lztwc_sp, lzfsc_sp, lzfpc_sp, adimc_sp, &
          !SAC Runoff variables
          roimp_sp,sdro_sp,ssur_sp,sif_sp,bfs_sp,bfp_sp, &
          !SAC OUTPUTS
          tci_sp, aet_sp)
    
      !write(*,*) "i =", i, " map_step = ", map_step, " etd_step = ", etd_step, " tci_sp = ", tci_sp 
      ! place state variables in output arrays
      tci(i,nh) = dble(tci_sp)

      if(return_states)then
        uztwc(i,nh) = dble(uztwc_sp)
        uzfwc(i,nh) = dble(uzfwc_sp)
        lztwc(i,nh) = dble(lztwc_sp)
        lzfsc(i,nh) = dble(lzfsc_sp)
        lzfpc(i,nh) = dble(lzfpc_sp)
        adimc(i,nh) = dble(adimc_sp)
        roimp(i,nh)=dble(roimp_sp)
        sdro(i,nh)=dble(sdro_sp)
        ssur(i,nh)=dble(ssur_sp)
        sif(i,nh)=dble(sif_sp)
        bfs(i,nh)=dble(bfs_sp)
        bfp(i,nh)=dble(bfp_sp)
        aet(i,nh) = dble(aet_sp)

        ! inout forcings to capture the pe/pxadj and efc
        map(i,nh) = map_step
        etd(i,nh) = etd_step

      end if

      
      ! PQNET
      ! PRAIN
      ! PROBG
      ! PSNWRO
      ! SNSG
      ! TINDEX
      ! SWE

      ! SNOW=SXFALL
      ! RAIM=RM(1)
      ! SNEQV=TWE/1000.
      ! SNOWH=SNDPT/100.

    end do  ! ============ end simulation time loop ====================

  end do   ! ========== end of simulation areas loop   ====================

end subroutine
