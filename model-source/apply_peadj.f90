subroutine  apply_peadj(n_hrus, dt, sim_length, year, month, day, hour, &
    peadj_m, mape, etd)

  use utilities
  implicit none

  ! this code is currently not set up to do any timestep less than 1 hour,
  ! nor could it do fractional hour timesteps.

  ! ! INPUTS ARGUMENTS ! !
  integer, intent(in):: n_hrus ! number of zones
  integer, intent(in):: dt    ! model timestep in seconds
  integer, intent(in):: sim_length   ! length of simulation for dt
  integer, dimension(sim_length), intent(in):: year, month, day, hour !TS of date component
  double precision, dimension(sim_length, n_hrus), intent(in):: mape  ! atmospheric forcing variables
  ! static pet adjustment
  double precision, dimension(12, n_hrus), intent(in):: peadj_m

  double precision :: dayn, dayi,interp_day, decimal_day !used to interpolate fa per ts
  integer:: mo  !current month for ts 

  ! ! LOCAL VARIBLES ! !
  integer:: i,k, nh           ! AWW index for looping through areas
  integer:: dt_hours, ts_per_day    ! model timestep in hours and number of timesteps per day
  integer, dimension(12) :: mdays, mdays_prev !lookup tables with number of days in month
  double precision, dimension(12, n_hrus):: peadj_m_prev, peadj_m_next
  double precision:: peadj_step, pe_step !adj for current ts in loop

  ! ! OUTPUT ARGUMENTS ! !
  ! forcing adjusted timeseries
  double precision, dimension(sim_length, n_hrus), intent(out):: etd

  mdays =      (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /) 
  mdays_prev = (/ 31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30 /) 

  etd = 0

  ts_per_day = 86400/dt
  dt_hours = dt/3600

  do nh=1,n_hrus
     ! there are monthly crop coefficients per zone 
     peadj_m_prev(1, nh) = peadj_m(12, nh)
     peadj_m_prev(2:12, nh) = peadj_m(1:11, nh)
     peadj_m_next(12, nh) = peadj_m(1, nh)
     peadj_m_next(1:11, nh) = peadj_m(2:12, nh)


     ! disintegrate the pe timeseries

     do i = 1, sim_length, ts_per_day
      !write(*,*) year(i), month(i), day(i), hour(i), mape(i, nh)
      pe_step = mape(i, nh) / ts_per_day
      do k = 0,  ts_per_day -1
        if ( i + k .le.  sim_length ) then
           etd(i + k, nh) = pe_step 
        else
           exit
        end if
      end do    
     end do

     do i = 1, sim_length
      ! adjust days in february if the year is a leap year
      if(mod(year(i),100) .ne. 0 .and. mod(year(i),4) .eq. 0) then
        mdays(2) = 29 ! leap year
      else if(mod(year(i),400).eq.0) then
        mdays(2) = 29 ! leap year
      else
        mdays(2) = 28 ! not leap year
      endif

      ! interpolate between (x0,y0) and (x1,y1)
      ! y = y0 + (x-x0)*(y1-y0)/(x1-x0)
      ! interpolate between (day0,limit0)=(0,limit0) and (day1,limit1)=(dayn,limit1)
      ! y = limit0 + dayi/dayn*(limit1-limit0)

      ! decimal day to start interpolation
      interp_day = 16. + dble(dt_hours)/24
      ! current decimal day
      decimal_day = dble(day(i)) + dble(dt_hours)/24.

      mo = month(i)

      if(decimal_day >= interp_day)then
        dayn = dble(mdays(mo))
        dayi = decimal_day - interp_day 
        peadj_step = peadj_m(mo, nh) + dayi/dayn*(peadj_m_next(mo, nh)-peadj_m(mo, nh))
      else 
        dayn = dble(mdays_prev(mo))
        dayi = decimal_day - interp_day + mdays_prev(mo) 
        peadj_step = peadj_m_prev(mo, nh) + dayi/dayn*(peadj_m(mo, nh)-peadj_m_prev(mo, nh))
      end if 
      etd(i, nh) = etd(i, nh) * peadj_step
      !write(*,*) year(i), month(i), day(i), hour(i), mape(i,nh), peadj_step, etd(i,nh)

    end do  ! ============ end simulation time loop ====================
  end do 


end subroutine
