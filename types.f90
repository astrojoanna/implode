module types
  implicit none
  public

  type :: swarm
    integer                         :: idnr     ! id of the swarm
    real(kind=8)                    :: npar     ! number of bodies
    real(kind=8)                    :: mass     ! mass of one body
    real(kind=8)                    :: stnr     ! Stokes number
    real(kind=8),dimension(3)       :: loc      ! x,y,z coordinates
    real(kind=8)                    :: rdis     ! distance to the cloud center
    real(kind=8),dimension(3)       :: vel      ! z,y,z, components of velocity
    integer                         :: nrcl     ! number of collisions the partice underwent
    integer                         :: nbou     ! number of bouncing collisions the particle underwent
    integer                         :: nfra     ! number of fragmentation collisions the particle underwent
    logical                         :: set      ! is the particle already settled?
    real(kind=8)                    :: menc     ! enclosed mass inside of this particle
  end type swarm

  type :: listofswarms
    type(swarm), pointer, dimension(:) :: s
  end type listofswarms

end
