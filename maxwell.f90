! calculate the Maxwellian velocity distribution
module maxwell
  use constants
  implicit none
  public

  contains

  subroutine drawvels(NN,sigma,xvel,yvel,zvel)
    implicit none
    integer, intent(in) :: NN
    real(kind=8), intent(in) :: sigma  ! sigma for the final distribution
    real(kind=8), dimension(NN), intent(out) :: xvel, yvel, zvel
    real(kind=8), dimension(2)  :: temp
    real(kind=8), dimension(3)  :: gauss
    real(kind=8) :: r, theta
    integer :: i, j

    do i = 1, NN
      ! from uniform to gaussian distribution
      do j = 1, 3
        call random_number(temp)
        r=(-2.0*log(temp(1)))**0.5
        theta = 2.0*pi*temp(2)
        gauss(j) = sigma*r*sin(theta)
      enddo
      xvel(i) = gauss(1)
      yvel(i) = gauss(2)
      zvel(i) = gauss(3)
    enddo

  end subroutine

end module
