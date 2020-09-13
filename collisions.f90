! this module performes collisions between representative bodies
! in this module all the quantities should be computed at the center of current cell
! in the matrices: first index -> representative particle, 2nd -> physical
module collisions

  use constants
  use parameters
  use types

  implicit none

  public :: mc_collisions
  private

  contains

  ! the routine performes collisional evolution on swarms located in the cell nr with the aid of the MC algorithm
  subroutine mc_collisions(swarms, volume, rzone, dtime, dtzone)
    implicit none
    type(swarm), pointer, dimension(:)                      :: swarms      ! local rps array
    real(kind=8), intent(in)                                :: rzone       ! radial distance of zone center
    real(kind=8), intent(in)                                :: dtime       ! time step that is used now
    real(kind=8), intent(out)                               :: dtzone      ! a typical collision time step the zone communicates back
    real(kind=8), intent(in)                                :: volume      ! volume of the cell
    integer                                                 :: nsws        ! number of representative particles in given cell
    integer                                                 :: nri, nrk    ! indices of physical and representative particles choosen to the next collision
    real(kind=8), dimension(:,:), allocatable               :: colrates    ! collision rates matrix
    real(kind=8), dimension(:,:), allocatable               :: relvels     ! relative velocities matrix
    real(kind=8), dimension(:,:), allocatable               :: accelncol   ! coagulation acceleration matrix (in the case of high mass ratio,
                                                                           ! instead of performing every collision separately, we group the collisions)
    real(kind=8), dimension(:), allocatable                 :: colri       ! collision rates for rps
    real(kind=8), dimension(:), allocatable                 :: vterm       ! terminal velocity of each particle calculated in zone center
    real(kind=8)                                            :: vff         ! free-fall velocity of each particle
    real(kind=8)                                            :: menczone
    real(kind=8), dimension(:,:), allocatable               :: vsph         ! spherical coordinates velocities
    real(kind=8)                                            :: local_time, dt ! physical time spent in the cell, time step between two collisions
    real(kind=8)                                            :: rand        ! random number
    real(kind=8)                                            :: totr        ! total collision rate
    integer                                                 :: i, k
    real(kind=8)                                            :: dustdens

    !-------------------------------------------------------------------------------------

    nsws = size(swarms)          ! number of swarms in the current cell

    dustdens =  real(nsws) * mswarm / volume

    allocate (colrates(nsws,nsws), accelncol(nsws,nsws), relvels(nsws,nsws), colri(nsws),vsph(nsws,3))
    allocate (vterm(nsws))

    menczone = sum(swarms(:)%menc) / real(nsws)
    vff = dsqrt(2.*Ggrav*menczone/R0)*dsqrt((R0-rzone)/rzone)
    vterm(:) = vff
    if (gas) vterm(:) = min(swarms(:)%stnr / omegaK * Ggrav * menczone / rzone**2.,vff)

    ! velocities cartesian -> spherical
    !vr
    vsph(:,1) = swarms(:)%loc(1)*swarms(:)%vel(1)/swarms(:)%rdis+swarms(:)%loc(2)*swarms(:)%vel(2)/swarms(:)%rdis + &
              swarms(:)%loc(3)*swarms(:)%vel(3)/swarms(:)%rdis
    !vtheta*r*cos(phi)
    vsph(:,2) = swarms(:)%vel(1)*swarms(:)%loc(2)/dsqrt(swarms(:)%loc(1)**2.+swarms(:)%loc(2)**2.)&
              -swarms(:)%loc(1)*swarms(:)%vel(2)/dsqrt(swarms(:)%loc(1)**2.+swarms(:)%loc(2)**2.)
    !vphi*r
    vsph(:,3) = swarms(:)%loc(3)*(swarms(:)%loc(1)*swarms(:)%vel(1)+swarms(:)%loc(2)*swarms(:)%vel(2))&
              /(swarms(:)%rdis*dsqrt(swarms(:)%loc(1)**2.+swarms(:)%loc(2)**2.)) - &
              (swarms(:)%loc(1)**2.+swarms(:)%loc(2)**2.)*swarms(:)%vel(3)&
              /(swarms(:)%rdis*dsqrt(swarms(:)%loc(1)**2.+swarms(:)%loc(2)**2.))

    ! calculates relative velocities and collision rates matrix
    do i = 1, nsws
      call rel_vels(i, vterm, vsph, relvels)
      call col_accel(i ,swarms, accelncol)
      call col_rates(i, swarms, relvels, colrates, accelncol, volume)
    enddo

    ! collision rate for representative particles
    colri(:) = sum(colrates(:,:),dim=2)

    !------------ MAIN COLLISIONS LOOP ----------------------------------------------------
    local_time = 0.0
    k = 0
    dtzone = 0
    do while (local_time < dtime)
      totr = sum(colri)                 ! total collision rate
      call random_number(rand)
      dt = -1. * log(rand) / totr       ! time step between this and the next collision
      if (dt > dtime) then ! 0 or 1 collisions, decided by a random number
        call random_number(rand)
        if (rand> dtime/dt) then
          local_time = dtime
          dtzone = dt
          cycle
        endif
      endif
      local_time = local_time + dt      ! update of the local time
      dtzone = dtzone + dt
      call choose_swarms(nri,nrk,colrates,colri,totr)
      swarms(nri)%nrcl = swarms(nri)%nrcl + 1
      call collision(nri,nrk,swarms,relvels, accelncol, vsph)
      if (gas) vterm(nri) = min(swarms(nri)%stnr / omegaK * Ggrav * menczone / rzone**2.,vff)
      call rel_vels(nri, vterm, vsph, relvels)
      relvels(:,nri) = relvels(nri,:)
      call col_accel(nri ,swarms, accelncol)
      call col_accel_r(nri ,swarms, accelncol)
      colri(:) = colri(:) - colrates(:,nri)
      call col_rates(nri, swarms, relvels, colrates, accelncol, volume)
      call col_rates_r(nri, swarms, relvels, colrates, accelncol, volume)
      colri(:) = colri(:) + colrates(:,nri)
      colri(nri) = sum(colrates(nri,:))
      k = k + 1
    enddo
    !-------------------------------------------------------------------------------------
    dtzone = dtzone / (real(k)+1.e-10)

    ! velocities spherical -> cartesian
    swarms(:)%vel(1) = swarms(:)%loc(1)*vsph(:,1)/swarms(:)%rdis + &
                      swarms(:)%loc(2)*vsph(:,2)/dsqrt(swarms(:)%loc(1)**2+swarms(:)%loc(2)**2) + &
                      swarms(:)%loc(3)*swarms(:)%loc(1)*vsph(:,3) / &
                      (swarms(:)%rdis*dsqrt(swarms(:)%loc(1)**2+swarms(:)%loc(2)**2))

    swarms(:)%vel(2) = swarms(:)%loc(2)*vsph(:,1)/swarms(:)%rdis - &
                      swarms(:)%loc(1)*vsph(:,2)/dsqrt(swarms(:)%loc(1)**2+swarms(:)%loc(2)**2) + &
                      swarms(:)%loc(3)*swarms(:)%loc(2)*vsph(:,3) / &
                      (swarms(:)%rdis*dsqrt(swarms(:)%loc(1)**2+swarms(:)%loc(2)**2))

    swarms(:)%vel(3) = swarms(:)%loc(3)*vsph(:,1)/swarms(:)%rdis - &
                      dsqrt(swarms(:)%loc(1)**2.+swarms(:)%loc(2)**2.) * vsph(:,3) / swarms(:)%rdis

    deallocate (colrates, accelncol, relvels, colri, vsph, vterm)

    return
  end subroutine mc_collisions

  subroutine rel_vels(nri, vterm, vsph, relvels)
    implicit none
    real(kind=8), dimension(:,:), allocatable      :: relvels
    integer                                        :: nri
    real(kind=8), dimension(:,:), allocatable      :: vsph
    real(kind=8), dimension(:), allocatable        :: vterm

    relvels(nri,:)  = dsqrt((vterm(nri) - vterm(:))**2. + (vsph(nri,2) - vsph(:,2))**2. + &
                      (vsph(nri,3) - vsph(:,3))**2.)

    return
  end subroutine rel_vels

  ! calculation of the collision rates between representative particle nl and all physical particles
  subroutine col_rates(nl, swarms, relvels, colrates, accelncol, vol)
    implicit none
    integer, intent(in)                                     :: nl       ! number of line of colrates to update
    type(swarm), pointer, dimension(:)                      :: swarms
    real(kind=8), dimension(:,:), allocatable               :: colrates
    real(kind=8), dimension(:,:), allocatable               :: accelncol
    real(kind=8), dimension(:,:), allocatable               :: relvels
    real(kind=8), intent(in)                                :: vol      ! volume of the cell

    ! basic eq. (see e.g. Drazkowska et al 2013, Eqs 19-20)
    colrates(nl,:) = swarms(:)%npar * relvels(nl,:) * con1 * (swarms(nl)%mass**(1./3.) + swarms(:)%mass**(1./3.))**2./vol
    ! instead of perform 1000 identical collisions, we divide the collision rate by 1000 but if the collision happens,
    ! we stick 1000 small particles at once
    where(accelncol(nl,:) > 0.0) colrates(nl,:) = colrates(nl,:) / accelncol(nl,:)
    ! if the representative particle represents less than 1 particle, the method is not valid anymore,
    ! so the collision rate is supressed
    if (swarms(nl)%npar <= 1.0) colrates(nl,:) = 0.0
    if (swarms(nl)%set) colrates(nl,:) = 0.0

    return
  end subroutine col_rates

  ! for updating the collision rates after collision: calculating the collsion rates between physical particle nl
  ! and all the representative particles
  subroutine col_rates_r(nl, swarms, relvels, colrates, accelncol, vol)
    implicit none
    integer, intent(in)                                     :: nl       ! number of column of colrates to update
    type(swarm), pointer, dimension(:)                      :: swarms
    real(kind=8), dimension(:,:), allocatable               :: colrates
    real(kind=8), dimension(:,:), allocatable               :: accelncol
    real(kind=8), dimension(:,:), allocatable               :: relvels
    real(kind=8), intent(in)                                :: vol      ! volume of the cell

    colrates(:,nl) = swarms(nl)%npar * relvels(nl,:) * con1 * (swarms(nl)%mass**(1./3.) + swarms(:)%mass**(1./3.))**2./vol
    where(accelncol(:,nl) > 0.0) colrates(:,nl) = colrates(:,nl) / accelncol(:,nl)
    where(swarms(:)%npar <= 1.0) colrates(:,nl) = 0.0
    where(swarms(:)%set) colrates(:,nl) = 0.0

    return
  end subroutine col_rates_r

  ! calculating the number of collisions in the case of collision acceleration
  ! only in the case when the representative particle is much larger than the physical particle
  subroutine col_accel(nri, swarms, accelncol)
    implicit none
    integer, intent(in)                                     :: nri
    real(kind=8), dimension(:), pointer                     :: accelncol_c
    type(swarm), pointer, dimension(:)                      :: swarms
    real(kind=8), allocatable, dimension(:,:), target       :: accelncol

    accelncol_c => accelncol(nri,:)

    accelncol_c(:) = 0.0
    where ((swarms(:)%mass / swarms(nri)%mass) < dmmax)
      accelncol_c(:) = swarms(nri)%mass * dmmax / swarms(:)%mass
    endwhere

    return
  end subroutine col_accel

  subroutine col_accel_r(nri, swarms, accelncol)
    implicit none
    integer, intent(in)                                     :: nri
    real(kind=8), dimension(:), pointer                     :: accelncol_r
    type(swarm), pointer, dimension(:)                      :: swarms
    real(kind=8), allocatable, dimension(:,:), target       :: accelncol

    accelncol_r => accelncol(:,nri)

    accelncol_r(:) = 0.0
    where ((swarms(nri)%mass / swarms(:)%mass) < dmmax)
        accelncol_r(:) = swarms(:)%mass * dmmax / swarms(nri)%mass
    endwhere

    return
  end subroutine col_accel_r

  ! choosing particles to the next collision
  ! nri -> representative
  ! nrk -> physical
  subroutine choose_swarms(nri, nrk, colrates, ri, totrate)
    implicit none
    integer, intent(out)                                    :: nri, nrk
    real(kind=8), dimension(:,:), allocatable, intent(in)   :: colrates
    integer                                                 :: i
    real(kind=8), dimension(:), allocatable, intent(in)     :: ri
    real(kind=8), intent(in)                                :: totrate
    real(kind=8), dimension(2)                              :: rand
    real(kind=8)                                            :: fin

    call random_number(rand) ! drawing 2 random numbers

    ! choosing the representative particle: the higher ri, the higher chance that it's particle "i"
    rand(1) = rand(1) * totrate
    i = 1
    fin = ri(1)
    do while(rand(1) > fin)
      fin = fin + ri(i+1)
      i = i + 1
    enddo
    nri = i

    ! choosing the physical particle
    i = 1
    rand(2) = rand(2) * ri(nri)
    fin = colrates(nri,1)
    do while (rand(2) > fin)
        fin = fin + colrates(nri,i+1)
        i = i + 1
    enddo
    nrk = i

    return
  end subroutine choose_swarms

  ! performing the collision: deciding the collision outcome - put your collision model here
  ! only the representative particle is updated
  subroutine collision(nri,nrk,swarms,relvels, accelncol, vsph)
    implicit none
    integer, intent(in)                                     :: nri, nrk
    type(swarm), pointer, dimension(:)                      :: swarms
    real(kind=8), dimension(:,:), allocatable               :: accelncol
    real(kind=8), dimension(:,:), allocatable               :: relvels
    real(kind=8), dimension(:,:), allocatable               :: vsph
    real(kind=8)                                            :: rvel
    real(kind=8)                                            :: fmass, redmass

    rvel = relvels(nri,nrk)
    redmass = swarms(nri)%mass*swarms(nrk)%mass/(swarms(nri)%mass+swarms(nrk)%mass)
    fmass = swarms(nrk)%mass/swarms(nri)%mass

    if (rvel > 1000. .and. fragm) then
      call fragmentation(nri,swarms)
    else
      if (rvel > sqrt(5.*pi*a0*Froll/redmass) .and. bounc) then
        call bouncing(nri,nrk,swarms, accelncol, vsph)
      else
        call hit_and_stick(nri,nrk,swarms, accelncol, vsph)
      endif
    endif
    swarms(nri)%npar = mswarm / swarms(nri)%mass ! update the number of particles in case the mass changed
    swarms(nri)%stnr = pi*0.5*con2*swarms(nri)%mass**(1./3.)*rhos/sigmag

    return
  end subroutine collision

  ! sticking collision
  subroutine hit_and_stick(nri,nrk,swarms, accelncol, vsph)
    implicit none
    integer, intent(in)                                     :: nri, nrk
    type(swarm), pointer, dimension(:)                      :: swarms
    real(kind=8), dimension(:,:), allocatable               :: accelncol
    real(kind=8), dimension(:,:), allocatable               :: vsph

    if (accelncol(nri, nrk) > 0.0) then
      vsph(nri,1:3) = (swarms(nri)%mass*vsph(nri,1:3) + accelncol(nri, nrk)*swarms(nrk)%mass*vsph(nrk,1:3)) /  &
                      (swarms(nri)%mass + accelncol(nri, nrk) * swarms(nrk)%mass)
      swarms(nri)%mass = swarms(nri)%mass + accelncol(nri, nrk) * swarms(nrk)%mass
    else
      vsph(nri,1:3) = (swarms(nri)%mass*vsph(nri,1:3) + swarms(nrk)%mass*vsph(nrk,1:3)) /  &
                      (swarms(nri)%mass + swarms(nrk)%mass)
      swarms(nri)%mass = swarms(nri)%mass + swarms(nrk)%mass
    endif

    return
  end subroutine hit_and_stick

  ! bouncing collision
  subroutine bouncing(nri,nrk,swarms, accelncol, vsph)
    implicit none
    integer, intent(in)                                     :: nri, nrk
    type(swarm), pointer, dimension(:)                      :: swarms
    real(kind=8), dimension(:,:), allocatable               :: accelncol
    real(kind=8), dimension(:,:), allocatable               :: vsph

    swarms(nri)%nbou = swarms(nri)%nbou + 1

    if (accelncol(nri, nrk) > 0.0) then
      vsph(nri,1:3) = (swarms(nri)%mass*vsph(nri,1:3) + accelncol(nri, nrk)*swarms(nrk)%mass*vsph(nrk,1:3) &
                      + accelncol(nri, nrk)*swarms(nrk)%mass*COR*(vsph(nrk,1:3) - vsph(nri,1:3))) /  &
                      (swarms(nri)%mass + accelncol(nri, nrk) * swarms(nrk)%mass)
    else
      vsph(nri,1:3) = (swarms(nri)%mass*vsph(nri,1:3) + swarms(nrk)%mass*vsph(nrk,1:3) &
                      + swarms(nrk)%mass * COR * (vsph(nrk,1:3) - vsph(nri,1:3))) /  &
                      (swarms(nri)%mass + swarms(nrk)%mass)
    endif

    return
  end subroutine bouncing

  ! fragmentation collision
  subroutine fragmentation(nri,swarms)
    implicit none
    integer, intent(in)                                     :: nri
    real(kind=8)                                            :: ran
    type(swarm), pointer, dimension(:)                      :: swarms
    real, parameter                                         :: kappa = 1./6.  ! n(m) ~ m^(kappa - 2)

    !write(*,*) 'F'
    swarms(nri)%nfra = swarms(nri)%nfra + 1

    call random_number(ran)
    swarms(nri)%mass = (ran * (swarms(nri)%mass**kappa - monomer**kappa) +  monomer**kappa )**(1./kappa)
    swarms(nri)%mass = max(swarms(nri)%mass, monomer)

    return
  end subroutine fragmentation

end
