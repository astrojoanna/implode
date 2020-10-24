program collapse
  use constants
  use types
  use maxwell
  use parameters
  use collisions
  use sort
  use eom
  use constantsrkf

  implicit none

  character(len=100)  :: ctrl_file        ! parameter file
  character(len=10) :: numberp
  character(len=50) :: file_name_appendix
  character(len=90) :: energy_fn, output_fn, sizedistr_fn

  ! ----- superparticles (swarms) ---
  type(swarm), dimension(:), allocatable, target :: allswarms   ! list of all superparticles, sorted by distance
  type(swarm), dimension(:), allocatable, target :: idallswarms ! list of all superparticles, sorted by id (used for outputs only)

  ! -- cloud-global variables -------
  real(kind=8), dimension(:), allocatable :: xvels, yvels, zvels ! initial random velocities
  real(kind=8), dimension(:), allocatable :: menc                ! enclosed mass for selfgravity calculation
  real(kind=8), dimension(:,:), allocatable :: vsph
  real(kind=8), dimension(:), allocatable :: vff
  real(kind=8), allocatable :: listy(:,:),rd(:),abstol(:,:),vd(:),ers(:,:),ynew(:,:),ynewfifth(:,:),acc(:,:)

  real(kind=8) :: ss
  real(kind=8) :: U0,T0,Uexp,Texp ! initial and current potential and kinetic energies
  real(kind=8) :: vvir0  ! virial speed calculated from the initial kinetic energy
  real(kind=8) :: sigma  ! velocity dispersiton
  real(kind=8) :: tff    ! free-fall time
  real(kind=8) :: mcore  ! mass of the core
  real(kind=8) :: dt, dteom  ! time step
  real(kind=8) :: time   ! current time
  real(kind=8) :: timenextout ! time of the next output (maximum, if no output before)
  real(kind=8), dimension(3) :: ran ! random number
  real(kind=8) :: phi, theta
  integer :: switch = 0
  integer, parameter :: s = 2 ! determines initial distribution of particles (here: homogenous density)
  integer :: iter        ! number of iterations
  integer :: i, j, sbc
  logical :: step

  ! -- radial grid ------------------
  integer :: NZ, NZ0 ! number of radial zones
  real(kind=8), dimension(:), allocatable :: rwall ! zone walls
  real(kind=8), dimension(:), allocatable :: rcent   ! zone center
  real(kind=8), dimension(:), allocatable :: dr      ! zone width
  real(kind=8), dimension(:), allocatable :: dens    ! zone density
  real(kind=8), dimension(:), allocatable :: volume  ! zone volume
  real(kind=8), dimension(:), allocatable :: zonedt  ! zone collisional timestep
  integer :: nperzone  ! number of superparticles per zone
  integer :: k1, k2    ! indices of first and last particle in zone
  type(listofswarms), dimension(:), allocatable :: zonelist ! list of swarms in zone

  ! -- settling ----------------------
  real(kind=8), dimension(:), allocatable ::  rdismin ! minimum distance of each particle to keep core density at solid density
  real(kind=8) :: rinner ! settling threshold
  real(kind=8) :: tau1, tau   ! optical depth of the innermost zone
  integer :: nonepercent
  integer :: nrset, nrsetold  ! number of particles which are settled

  ! - data for histogram of particle masses
  integer, parameter :: nbins = 80
  real(kind=8), dimension(nbins+1) :: mgrid
  real(kind=8), dimension(nbins) :: m2fm
  real(kind=8) :: ll, lll, nord
  !-----------------------------------

  ! read parameter file
  call get_command_argument(1, ctrl_file)
  write(*,*) 'Reading parameters...'
  call read_parameters(ctrl_file)

  write(*,*) 'Initializing the simulation.......'

  ! allocate memory
  allocate(allswarms(NN),xvels(NN),yvels(NN),zvels(NN),idallswarms(NN),rdismin(NN),menc(NN),vsph(NN,3),vff(NN))

  ! set output filenames
  write(numberp, '(i0)') NN
  file_name_appendix = '-N'//trim(adjustl(numberp))
  if (fdisp>0.0) then
    file_name_appendix = trim(file_name_appendix) // '-disp'
  else
    file_name_appendix = trim(file_name_appendix) // '-nodisp'
  endif
  if (gas) then
    file_name_appendix = trim(file_name_appendix) // '-gas'
  else
    file_name_appendix = trim(file_name_appendix) // '-nogas'
  endif
  if (colls) then
    file_name_appendix = trim(file_name_appendix) // '-colls'
  else
    file_name_appendix = trim(file_name_appendix) // '-nocolls'
  endif
  file_name_appendix = trim(file_name_appendix) // '.dat'

  energy_fn = 'energy' // trim(file_name_appendix)
  output_fn = 'output' // trim(file_name_appendix)
  sizedistr_fn = 'sizedistr' // trim(file_name_appendix)

  ! initialization of representative particles: all swarms represent the same physical mass and particles of identical mass
  do i = 1, NN
    allswarms(i)%idnr = i
    allswarms(i)%mass = mmax
    call random_number(ran)
    if (initsizedistr) then
      allswarms(i)%mass = (ran(1) * (mmax**kappa - mini**kappa) + mini**kappa )**(1./kappa) ! size distribution between mmax and mini
    endif
    allswarms(i)%mass = max(allswarms(i)%mass, monomer)
    allswarms(i)%npar = mswarm/allswarms(i)%mass
    allswarms(i)%rdis = ((R0**(s+1.))*real(i)/real(NN-1))**(1./(s+1.))+1.d-16 ! this distributes particles with an uniform density if in a sphere
    theta = 2.d0*pi*ran(2)
    phi = asin(2.d0*ran(3)-1.0)
    allswarms(i)%loc(1) = allswarms(i)%rdis * cos(phi) * cos(theta)
    allswarms(i)%loc(2) = allswarms(i)%rdis * cos(phi) * sin(theta)
    allswarms(i)%loc(3) = allswarms(i)%rdis * sin(phi)
    allswarms(i)%stnr = 1.d+6
    if (gas) allswarms(i)%stnr = pi*0.5*con2*allswarms(i)%mass**(1./3.)*rhos/sigmag
    allswarms(i)%set = .false.
    allswarms(i)%nrcl = 0
    allswarms(i)%nbou = 0
    allswarms(i)%nfra = 0
  enddo

  ! make histogram of initial pebbles size distribution
  mgrid(1) = 0.5*monomer
  nord = 20
  ll = real(nord) / real(nbins)
  lll = ll
  do i = 2, nbins+1
    mgrid(i) = monomer * 10.0 ** lll
    lll = ll * i
  enddo

  m2fm(:) = 0.0
  do i = 1, NN
    j = 1
    do while (allswarms(i)%mass > mgrid(j))
      j = j + 1
    enddo
    m2fm(j) = m2fm(j) + ((allswarms(i)%npar) * allswarms(i)%mass**2) / &
              ((mgrid(j+1) - mgrid(j)) * mswarm * NN)
  enddo

  open(22,file='sizedistr-init' // trim(file_name_appendix),status='new')
  do j = 1, nbins
    write(22,*) sqrt(mgrid(j+1)*mgrid(j))/monomer, con2*sqrt(mgrid(j+1)*mgrid(j))**(1./3.), m2fm(j)
  enddo
  close(22)

  ! calculate minimum distance of each consecutive paricle so that density=rhos
  do i = 1, NN
    rdismin(i) = 0.999*con2*(mswarm*real(i))**(1./3.)+1.e-9
    menc(i) = real(i)*mswarm
    allswarms(i)%menc = menc(i)
  enddo

  ! build the initial grid: allswarms is already sorted by distance
  NZ = NN/nparts
  NZ0 = NZ
  allocate(rwall(NZ+1), rcent(NZ), dr(NZ), dens(NZ), volume(NZ), zonedt(NZ))
  allocate(zonelist(NZ))
  nperzone = nparts

  do i = 1, NZ
    k1 = (i-1)*nperzone + 1
    k2 = i * nperzone
    zonelist(i)%s(1:) => allswarms(k1:k2) ! binned swarms lists are pointers so allswarms is modified automatically
  enddo

  rwall(1) = 0.0
  do i = 2, NZ
    rwall(i) = 0.5*(zonelist(i)%s(1)%rdis + zonelist(i-1)%s(nperzone)%rdis)
  enddo
  rwall(NZ+1) = zonelist(NZ)%s(nperzone)%rdis + 0.5*(zonelist(NZ)%s(nperzone)%rdis-zonelist(NZ)%s(nperzone-1)%rdis)
  dr(:) = rwall(2:NZ+1) - rwall(1:NZ)
  do i = 1, NZ
    rcent(i) = sum(zonelist(i)%s(:)%rdis)/real(size(zonelist(i)%s))
  enddo

  ! calculate volumes and densities in each radial zone
  volume(:) = 4./3.*pi*(rwall(2:NZ+1)**3-rwall(1:NZ)**3)
  dens(:) = nperzone*mswarm/volume(:)

  !write initial density of each zone
  if (screenout) then
    write(*,*) 'Densities in each zone:'
    do j = 1, NZ
      write(*,*) j, rcent(j)/R0, dens(j)
    enddo
  endif

  ! calculate initial random velocities assuming virial equilibrium
  U0 = -0.6 * Ggrav * M0**2. / R0
  T0 = -0.5*U0
  tff = 0.5*pi*R0*sqrt(R0/(2.*Ggrav*M0)) ! free-fall time
  vvir0 = sqrt(2.*T0 / M0)
  sigma = vvir0/sqrt(3.)

  write(*,*) 'Free fall timescale of the cloud is:', tff/year,' years'
  if (gas) write(*,*) 'Average Stopping time of particles is', sum(allswarms(:)%stnr)/real(NN)/omegaK/year,' years'
  if (gas) write(*,*) 'Average Stokes number is', sum(allswarms(:)%stnr)/real(NN)
  if (gas) write(*,*)

  ! draw the initial speeds of particles
  call drawvels(NN,sigma,xvels,yvels,zvels)
  allswarms(:)%vel(1) = fdisp*xvels(:)
  allswarms(:)%vel(2) = fdisp*yvels(:)
  allswarms(:)%vel(3) = fdisp*zvels(:)

  if (initvff) then
    ! velocities cartesian -> spherical
    vsph(:,1) = allswarms(:)%loc(1)*allswarms(:)%vel(1)/allswarms(:)%rdis+&
          allswarms(:)%loc(2)*allswarms(:)%vel(2)/allswarms(:)%rdis + &
          allswarms(:)%loc(3)*allswarms(:)%vel(3)/allswarms(:)%rdis

    vsph(:,2) = allswarms(:)%vel(1)*allswarms(:)%loc(2)/dsqrt(allswarms(:)%loc(1)**2.+allswarms(:)%loc(2)**2.)&
          -allswarms(:)%loc(1)*allswarms(:)%vel(2)/dsqrt(allswarms(:)%loc(1)**2.+allswarms(:)%loc(2)**2.)

    vsph(:,3) = allswarms(:)%loc(3)*(allswarms(:)%loc(1)*allswarms(:)%vel(1)+allswarms(:)%loc(2)*allswarms(:)%vel(2))&
          /(allswarms(:)%rdis*dsqrt(allswarms(:)%loc(1)**2.+allswarms(:)%loc(2)**2.)) - &
          (allswarms(:)%loc(1)**2.+allswarms(:)%loc(2)**2.)*allswarms(:)%vel(3)&
          /(allswarms(:)%rdis*dsqrt(allswarms(:)%loc(1)**2.+allswarms(:)%loc(2)**2.))

    ! replace the radial component with analytical prediction
    if (gas) then
      vff(:) = allswarms(:)%stnr / omegaK * Ggrav * menc(:) / allswarms(:)%rdis**2.
      where (allswarms(:)%stnr > 0.25) vff(:) = 0.d0
    else
      vff(:) = 0.d0
    endif
    vsph(:,1) = -1.d0*vff(:)

    ! velocities spherical -> cartesian
    allswarms(:)%vel(1) = allswarms(:)%loc(1)*vsph(:,1)/allswarms(:)%rdis + &
                     allswarms(:)%loc(2)*vsph(:,2)/dsqrt(allswarms(:)%loc(1)**2+allswarms(:)%loc(2)**2) + &
                     allswarms(:)%loc(3)*allswarms(:)%loc(1)*vsph(:,3) / &
                     (allswarms(:)%rdis*dsqrt(allswarms(:)%loc(1)**2+allswarms(:)%loc(2)**2))

    allswarms(:)%vel(2) = allswarms(:)%loc(2)*vsph(:,1)/allswarms(:)%rdis - &
                     allswarms(:)%loc(1)*vsph(:,2)/dsqrt(allswarms(:)%loc(1)**2+allswarms(:)%loc(2)**2) + &
                     allswarms(:)%loc(3)*allswarms(:)%loc(2)*vsph(:,3) / &
                     (allswarms(:)%rdis*dsqrt(allswarms(:)%loc(1)**2+allswarms(:)%loc(2)**2))

    allswarms(:)%vel(3) = allswarms(:)%loc(3)*vsph(:,1)/allswarms(:)%rdis - &
                     dsqrt(allswarms(:)%loc(1)**2.+allswarms(:)%loc(2)**2.) * vsph(:,3) / allswarms(:)%rdis
  endif

  ! >>>>>>>>>>>>>>>>>>>>>>>
  ! calculate the energies explicitly
  Uexp = 0.
  Texp = 0.
  do i = 1, NN
    Uexp = Uexp - Ggrav/allswarms(i)%rdis * mswarm * menc(i)
    Texp = Texp + 0.5*mswarm*(allswarms(i)%vel(1)**2.+allswarms(i)%vel(2)**2.+allswarms(i)%vel(3)**2.)
  enddo
  write(*,*) 'Initial energy budget:'
  write(*,*) 'U0=',U0, ' Uexp=',Uexp
  write(*,*) 'T0=',T0, ' Texp=',Texp
  write(*,*)
  ! <<<<<<<<<<<<<<<<<<<<<<<

  ! initialization
  time = 0.d0
  iter = 0
  sbc = 0
  mcore = 0.0
  zonedt(:) = 100.d0/omegaK
  timenextout = outputdt
  dteom = delt
  nrset = 0
  if (gas) switch = 1
  rinner = rdismin(1) ! location of the settling threshold initially, then updated to tau=1 location

  ! check if the inner 1% of mass is already optically thick
  nonepercent = nint(0.01*NN) ! particles with 1% of total mass
  tau1 = con1*sum(allswarms(1:nonepercent)%npar*allswarms(1:nonepercent)%mass**(2./3.) / &
        (4.*pi*allswarms(1:nonepercent)%rdis**2.))
  ! if yes, check where the tau=1 line is and set accretion threshold there
  if (tau1 > 1.) then
    tau = 0.
    i = 1
    do while (tau < 1.)
      tau = tau + allswarms(i)%npar*con1*allswarms(i)%mass**(2./3.)/(4.*pi*allswarms(i)%rdis**2.)
      i = i + 1
    enddo
    rinner = allswarms(i-1)%rdis
  endif

  open(99,file=energy_fn,status='new')
  write(99,*) time/year, Texp, abs(Uexp), mcore/M0, sum(allswarms(:)%nrcl), &
              sum(allswarms(:)%nbou), sum(allswarms(:)%nfra), tau1, rinner/R0

  open(11,file=output_fn,status='new')
  call write_output(time,allswarms)
  if (.not.fullout) close(11)


  ! >> MAIN LOOP >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  do while (.not.all(allswarms(:)%set) .and. time<tmax) ! the code runs until either all particles are settled or tmax is exceeded
    ! calculate the timestep
    dt = min(dteom, outputdt)
    if (colls) dt = min(dt,0.5*real(nperzone)*minval(zonedt(:)))
    if (gas) dt = min(dt,minval(allswarms(:)%stnr,.not.allswarms(:)%set)/omegaK)
    if (dt<=0.d0) then
      write(*,*) 'SOMETHING WENT WRONG, HERE IS THE LIST OF SWARMS: '
      do i = 1, NN
        write(*,*) allswarms(i)%idnr, allswarms(i)%rdis/R0, allswarms(i)%vel(:), allswarms(i)%stnr, allswarms(i)%set
      enddo
      stop
    endif

    ! advecton part starts here ---->
    step = .false.

    allocate(listy(NN-nrset,6),rd(NN-nrset),abstol(NN-nrset,6),vd(NN-nrset), &
            ers(NN-nrset,6),ynew(NN-nrset,6),ynewfifth(NN-nrset,6),acc(NN-nrset,6))

    listy(:,1) =allswarms(nrset+1:NN)%vel(1)
    listy(:,2) =allswarms(nrset+1:NN)%vel(2)
    listy(:,3) =allswarms(nrset+1:NN)%vel(3)
    listy(:,4)=allswarms(nrset+1:NN)%loc(1)
    listy(:,5)=allswarms(nrset+1:NN)%loc(2)
    listy(:,6)=allswarms(nrset+1:NN)%loc(3)
    rd = dsqrt(allswarms(nrset+1:NN)%loc(1)**2. + allswarms(nrset+1:NN)%loc(2)**2. + allswarms(nrset+1:NN)%loc(3)**2.)

    do while (.not.step)
      call rkf(listy,acc,NN-nrset,softl,allswarms(nrset+1:NN)%stnr/omegaK,dt,ynew,ynewfifth,rd,menc(nrset+1:NN),switch)

      rd = dsqrt(ynew(:,4)**2 + ynew(:,5)**2 + ynew(:,6)**2)
      vd = dsqrt(ynew(:,1)**2 + ynew(:,2)**2 + ynew(:,3)**2)

      abstol(:,1) = rd* tol
      abstol(:,2) = abstol(:,1)
      abstol(:,3) = abstol(:,1)
      abstol(:,4) = vd* tol
      abstol(:,5) = abstol(:,4)
      abstol(:,6) = abstol(:,4)
      ers = abs(ynew - ynewfifth) + 1e-40
      ss = 0.84 * (minval(abstol / ers))**0.25

      if (maxval(ers / abstol)<= 1.) then
        step = .true.
        time = time + dt
        listy = ynew

        allswarms(nrset+1:NN)%vel(1) = ynew(:,1)
        allswarms(nrset+1:NN)%vel(2) = ynew(:,2)
        allswarms(nrset+1:NN)%vel(3) = ynew(:,3)
        allswarms(nrset+1:NN)%loc(1) = ynew(:,4)
        allswarms(nrset+1:NN)%loc(2) = ynew(:,5)
        allswarms(nrset+1:NN)%loc(3) = ynew(:,6)

        dteom = dt * ss
      else
        if (screenout) write(*,*) "Making timestep smaller at time ",time/year
        dt = dt * ss
      endif

    enddo

    deallocate(listy,rd,abstol,vd,ers,ynew,ynewfifth,acc)
    !<---- advection part ends here

    !update the radial distance of swarms
    allswarms(:)%rdis = dsqrt(allswarms(:)%loc(1)**2. + allswarms(:)%loc(2)**2. + allswarms(:)%loc(3)**2.)

    ! at this place sorting is neccesary if settling allows for overtaking
    call shell_sort_r(allswarms)
    do i = 1, NN
      allswarms(i)%menc = menc(i)
    enddo

    ! check the settling threshold: where the cloud center reaches the optical depth of 1
    ! (if this didn't happen previosuly = the inner 1% of mass just got to tau=1)
    nonepercent = min(nint(0.01*NN)+nrset,NN)
    tau1 = con1*sum(allswarms(nrset+1:nonepercent)%npar*allswarms(nrset+1:nonepercent)%mass**(2./3.)/ &
          (4.*pi*(allswarms(nrset+1:nonepercent)%rdis**2.)))
    if (tau1 >= 1. .and. rinner < rdismin(nrset+2)) rinner = allswarms(nonepercent)%rdis

    ! check for settling
    do i = 1,NN
      if (allswarms(i)%rdis < rinner) then
        allswarms(i)%rdis = rdismin(i)
        call random_number(ran)
        allswarms(i)%loc(1) = allswarms(i)%rdis * cos(phi) * cos(theta)
        allswarms(i)%loc(2) = allswarms(i)%rdis * cos(phi) * sin(theta)
        allswarms(i)%loc(3) = allswarms(i)%rdis * sin(phi)
        allswarms(i)%vel(:) = 1.e-20
        allswarms(i)%set = .true.
      endif
    enddo
    nrsetold = nrset
    rinner = max(rdismin(nrset+1),rinner)
    nrset = count(allswarms(:)%set) ! number of settled particles
    mcore = real(nrset)*mswarm ! mass of the core

    ! shell breaking
    if (gas .and. maxval(dens(:))>rhos) then
      do i = 1, NZ
        if (dens(i) > rhos .and. .not.zonelist(i)%s(1)%set) then
          if (screenout) write(*,*) 'Shell breaking in zone ',i, ' density ',dens(i), ' time ',time/year
          sbc = sbc + 1
          zonelist(i)%s(1:nperzone)%mass = boulder
          zonelist(i)%s(1:nperzone)%npar = mswarm/boulder
        endif
      enddo
    endif

    ! update grid - if more particles settled, rebuild the zones
    if (nrset .ne. nrsetold) then
      ! remove old grid
      deallocate(rwall, rcent, dr, dens, volume, zonedt, zonelist)

      ! how many zones now?
      do while ( (real(NN-nrset)/real(NZ)) < real(nparts) )
        NZ = NZ - 1
      enddo
      if (NZ < 1) NZ = 1

      ! how many particles per zone?
      nperzone = nint(real(NN-nrset)/real(NZ))

      allocate(rwall(NZ+1), rcent(NZ), dr(NZ), dens(NZ), volume(NZ), zonedt(NZ))
      allocate(zonelist(NZ))

      if (NZ > 1) then
        do i = 1, NZ-1
          k1 = nrset + (i-1)*nperzone + 1
          k2 = nrset + i * nperzone
          zonelist(i)%s(1:) => allswarms(k1:k2) ! binned swarms lists are pointers so allswarms is modified automatically
        enddo
        zonelist(NZ)%s(1:) => allswarms(k2+1:NN)
      else
        zonelist(1)%s(1:) => allswarms(nrset+1:NN)
      endif
    endif

    ! update the walls / cell centers
    if (nrset > 0) then
      rwall(1) = rinner
    else
      rwall(1) = 0.0
    endif
    do i = 2, NZ
      rwall(i) = 0.5*(zonelist(i)%s(1)%rdis + zonelist(i-1)%s(nperzone)%rdis)
    enddo
    rwall(NZ+1) = allswarms(NN)%rdis + 0.5*(allswarms(NN)%rdis-allswarms(NN-1)%rdis)
    dr(:) = rwall(2:NZ+1) - rwall(1:NZ)
    do i = 1, NZ
      rcent(i) = sum(zonelist(i)%s(:)%rdis)/real(size(zonelist(i)%s))
    enddo

    ! calculate volumes and densities
    volume(:) = 4.*pi*(rwall(2:NZ+1)**3-rwall(1:NZ)**3)/3.
    dens(:) = nperzone*mswarm/volume(:)

    ! collisions
    zonedt(:) = 100./omegaK
    Texp = 0.
    do i = 1, NN
      Texp = Texp + 0.5*mswarm*(allswarms(i)%vel(1)**2.+allswarms(i)%vel(2)**2.+allswarms(i)%vel(3)**2.)
    enddo

    if (colls) then
    !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(DYNAMIC)
    do i = 1, NZ
      if (dens(i) < rhos) call mc_collisions(zonelist(i)%s, volume(i), rcent(i), dt, zonedt(i))
    enddo
    !$OMP END PARALLEL DO
    endif
    Texp = 0.
    do i = 1, NN
      Texp = Texp + 0.5*mswarm*(allswarms(i)%vel(1)**2.+allswarms(i)%vel(2)**2.+allswarms(i)%vel(3)**2.)
    enddo

    iter = iter + 1

    ! WRITE OUTPUT (every outfreq iterations)
    if (mod(iter,outfreq)==0 .or. time>timenextout) then

      idallswarms = allswarms
      call shell_sort_id(idallswarms)
      if (.not.fullout) open(11,file=output_fn,status='replace')
      call write_output(time,idallswarms)
      if (.not.fullout) close(11)
      timenextout = time + outputdt

      ! calculate the energies
      Uexp = 0.
      Texp = 0.
      do i = 1, NN
        Uexp = Uexp - Ggrav/allswarms(i)%rdis * mswarm * menc(i)
        Texp = Texp + 0.5*mswarm*(allswarms(i)%vel(1)**2.+allswarms(i)%vel(2)**2.+allswarms(i)%vel(3)**2.)
      enddo
      write(99,*) time/year, Texp, abs(Uexp), mcore/M0, sum(allswarms(:)%nrcl), &
                  sum(allswarms(:)%nbou), sum(allswarms(:)%nfra),tau1, rinner/R0

      write(*,*) 'Output written at t=',time/year,'years, there are',nrset,&
                ' settled particles, there were ',sum(allswarms(:)%nrcl),' collisions'

      if (screenout) write(*,*) '    --  number of zones is now: ', NZ
      !write density in each zone
      if (screenout) then
        write(*,*) 'Densities in each zone:'
        do j = 1, NZ
          write(*,*) j, rcent(j)/R0, dens(j)
        enddo
        write(*,*)
      endif
    endif

  enddo
  !-------------- END OF THE MAIN LOOP -----------------------------------------

  ! final output
  idallswarms = allswarms
  call shell_sort_id(idallswarms)
  if (.not.fullout) open(11,file=output_fn,status='replace')
  call write_output(time,idallswarms)
  if (.not.fullout) close(11)
  ! calculate the energies
  Uexp = 0.
  Texp = 0.
  do i = 1, NN
    Uexp = Uexp - Ggrav/allswarms(i)%rdis * mswarm * menc(i)
    Texp = Texp + 0.5*mswarm*(allswarms(i)%vel(1)**2.+allswarms(i)%vel(2)**2.+allswarms(i)%vel(3)**2.)
  enddo
  write(99,*) time/year, Texp, abs(Uexp), mcore/M0, sum(allswarms(:)%nrcl), &
              sum(allswarms(:)%nbou), sum(allswarms(:)%nfra),tau1, rinner/R0

  ! make histogram of pebbles size distribution
  m2fm(:) = 0.0
  do i = 1, NN
    j = 1
    do while (allswarms(i)%mass > mgrid(j))
      j = j + 1
    enddo
    m2fm(j) = m2fm(j) + ((allswarms(i)%npar) * allswarms(i)%mass**2) / &
            ((mgrid(j+1) - mgrid(j)) * mswarm * NN)
  enddo

  open(22,file=sizedistr_fn,status='new')
  do j = 1, nbins
    write(22,*) sqrt(mgrid(j+1)*mgrid(j))/monomer, con2*sqrt(mgrid(j+1)*mgrid(j))**(1./3.), m2fm(j)
  enddo
  close(22)

  write(*,*) 'Run completed! The output can be found in `output...dat` and `energy...dat`.  &
              Mass/size distribution is in `sizedistr...dat`'

  deallocate(allswarms)
  deallocate(idallswarms)
  deallocate(xvels,yvels,zvels)
  deallocate(zonelist)
  deallocate(rdismin)
  deallocate(rwall, rcent, dr, dens, volume, zonedt)
  deallocate(menc)
  deallocate(vsph)

  if (fullout) close(11)
  close(99)

end program
