module parameters

use constants
use types
implicit none
public

logical      :: gas           ! gas switch
logical      :: colls         ! collisions switch
integer      :: outfreq       ! output frequency (nr of timesteps, higher nr = smaller file)
real(kind=8) :: outputdt      ! maximum timestep between outputs
logical      :: screenout     ! screen output switch
logical      :: bounc         ! bouncing on/off
logical      :: fragm         ! fragmentation on/off
logical      :: initsizedistr ! initial size distribution on/off
logical      :: initvff       ! initialize radial velocity with free-fall/terminal velocities on/off
logical      :: fullout       ! write full output on/off
integer      :: NN            ! total number of swarms in the simulation
integer      :: nparts        ! number of particles per zone
real(kind=8) :: delt          ! basic timestep for advection solver
real(kind=8) :: tmax          ! end time of the simulation in the case of no full collapse is reached
real(kind=8) :: Rsolid        ! final radius of the planetesimal [cm]
real(kind=8) :: ainimax       ! maximum initial size of a pebble [cm]
real(kind=8) :: ainimin       ! minimum initial size of a pebble [cm]
real(kind=8) :: rhos          ! solid density [g cm^-3]
real(kind=8) :: adis          ! semi-major axis of the cloud
real(kind=8) :: omegaK        ! Keplerian frequency
real(kind=8) :: sigmag        ! gas surface density
real(kind=8) :: M0            ! mass of the cloud [g]
real(kind=8) :: a0            ! monomer radius [cm]
real(kind=8) :: monomer       ! monomer mass
real(kind=8) :: mini, mmax    ! initial mass of a pebble (min and max)
real(kind=8) :: R0            ! initial radius of the cloud [cm]
real(kind=8) :: dmmax         ! collision acceleration parameter
real(kind=8) :: mswarm        ! physical mass of dust represented by each representative particle
real(kind=8) :: COR           ! coefficient of restitution for bouncing collisions
real(kind=8) :: con1          ! optimization
real(kind=8) :: con2          ! optimization
real(kind=8) :: boulder       ! mass of a boulder = crushed shell fragment
real(kind=8) :: softl         ! softening length for selfgravity calculation
real(kind=8) :: kappa         ! adjusting initial size distribution
real(kind=8) :: fdisp         ! initial velocity dispersion factor


contains

subroutine read_parameters(ctrl_file)
  implicit none
  character(len=100), intent(in)   :: ctrl_file
  character(len=100)               :: buffer, label
  integer                          :: pos
  integer                          :: ios = 0
  integer                          :: line = 0
  integer, parameter               :: fh = 1

  ! reading the parameter file
  open(fh,file=ctrl_file,action='read')

  do while (ios == 0)
    read(fh, '(A)', iostat=ios) buffer
    if (ios == 0) then
      line = line + 1

      pos = scan(buffer, ' 	')
      label = buffer(1:pos)
      buffer = buffer(pos+1:)

      select case (label)

      case ('gas')
        read(buffer, *, iostat=ios) gas
        print *, 'Gas switch: ', gas
      case ('collisions')
        read(buffer, *, iostat=ios) colls
        print *, 'Collsions switch: ', colls
      case ('output_freq')
        read(buffer, *, iostat=ios) outfreq
        print *, 'Output frequency: ', outfreq
      case ('output_dt_year')
        read(buffer, *, iostat=ios) outputdt
        print *, 'Maximum time between outputs: ', outputdt, 'years'
        outputdt = outputdt * year
      case('screen_out')
        read(buffer, *, iostat=ios) screenout
        print *, 'Screen output switch: ', screenout
      case('nr_parts')
        read(buffer, *, iostat=ios) NN
        print *, 'Total number of superparticles: ', NN
      case('nr_parts_per_zone')
        read(buffer, *, iostat=ios) nparts
        print *, 'Number of particles per zone: ', nparts
      case('delta_t_omegaK')
        read(buffer, *, iostat=ios) delt
        print *, 'Basic timestep for advection solver: ', delt, '/OmegaK'
      case('max_time_yr')
        read(buffer, *, iostat=ios) tmax
        print *, 'Maximum time of the simulation: ', tmax, 'years'
        tmax = tmax * year
      case('radial_distance_AU')
        read(buffer, *, iostat=ios) adis
        print *, 'Radial distance: ', adis,'AU'
        adis = adis*AU
      case('gas_surfde_cgs')
        read(buffer, *, iostat=ios) sigmag
        print *, 'Gas surface density: ', sigmag, 'g cm-2'
      case('dmmax')
        read(buffer, *, iostat=ios) dmmax
        print *, 'Collisions speedup parameter: ', dmmax
      case('final_radius_km')
        read(buffer, *, iostat=ios) Rsolid
        print *, 'Final planetesimal radius: ', Rsolid, 'km'
        Rsolid = Rsolid * 1.e+5
      case('init_max_size')
        read(buffer, *, iostat=ios) ainimax
        print *, 'Initial maximum size: ', ainimax, 'cm'
      case('init_min_size')
        read(buffer, *, iostat=ios) ainimin
        print *, 'Initial minimum size: ', ainimin, 'cm'
      case('material_density_cgs')
        read(buffer, *, iostat=ios) rhos
        print *, 'Material density: ', rhos, 'g cm-3'
      case('monomer_size_cm')
        read(buffer, *, iostat=ios) a0
        print *, 'Monomer size: ', a0, 'cm'
      case('COR')
        read(buffer, *, iostat=ios) COR
        print *, 'Coefficient of restitution:', COR
      case('bouncing')
        read(buffer, *, iostat=ios) bounc
        print *, 'Bouncing:', bounc
      case('fragmentation')
        read(buffer, *, iostat=ios) fragm
        print *, 'Fragmentation:', fragm
      case('initial_size_distr')
        read(buffer, *, iostat=ios) initsizedistr
        print *, 'Initial size distribution:', initsizedistr
      case('softening_length_Rsolid')
        read(buffer, *, iostat=ios) softl
        print *, 'Softening length for self-gravity calculation:', softl,' final radii'
        softl=softl*Rsolid
      case('initial_vel_disp_factor')
        read(buffer, *, iostat=ios) fdisp
        print *, 'Initial velocity dispersion factor:', fdisp
      case('init_free_fall')
        read(buffer, *, iostat=ios) initvff
        print *, 'Initialize radial velocity with analytical prediction:', initvff
      case('full_output')
        read(buffer, *, iostat=ios) fullout
        print *, 'Write full output:', fullout
      case('size_distr_kappa')
        read(buffer, *, iostat=ios) kappa
        print *, 'Size distribution parameter kappa:', kappa

      case default
        print *, 'Skipping invalid label at line', line
      end select
    end if
  end do

  close(fh)

  ! calculate other parameters from the input
  omegaK = sqrt(Ggrav*Msun/adis**3.)  ! Keplerian frequency
  M0 = 4.*pi*rhos*Rsolid**3./3.       ! mass of the cloud [g]
  monomer = 4.*pi*rhos*a0**3./3.      ! monomer mass
  mmax = 4.*pi*rhos*ainimax**3./3.    ! initial maximum mass of a pebble
  mini = 4.*pi*rhos*ainimin**3./3.    ! initial minimim mass of a pebble
  boulder = 4.*(2.*sigmag)**3./(3.*(pi*rhos)**2.)    ! mass of a St=1 particle
  R0 = adis*(M0/(3.*Msun))**(1./3.)   ! initial radius of the cloud [cm] (here: equal to the Hill radius)
  mswarm=M0/real(NN)                  ! physical mass of dust represented by each representative particle
  con1 = pi**(1./3.) * (0.75 / rhos)**(2./3.) ! optimization
  con2 = (0.75 / pi / rhos)**(1./3.)          ! optimization

  return
end subroutine read_parameters

subroutine write_output(time,list)
  implicit none
  type(swarm), dimension(:), allocatable :: list
  real(kind=8) :: time

  write(11,*) time/year, mswarm, R0
  write(11,*) list(:)%loc(1)
  write(11,*) list(:)%loc(2)
  write(11,*) list(:)%loc(3)
  write(11,*) list(:)%mass/mini
  write(11,*) list(:)%vel(1)
  write(11,*) list(:)%vel(2)
  write(11,*) list(:)%vel(3)
  write(11,*) list(:)%nrcl
  write(11,*) list(:)%set
  write(11,*) list(:)%nbou
  write(11,*) list(:)%nfra
  write(11,*)

end subroutine write_output

end module
