module sort
  use types
  implicit none
  public
  contains

! sorting particles by radial distance using shell sort algorithm
  subroutine shell_sort_r(swrm)
    implicit none
    type(swarm), dimension(:), allocatable :: swrm
    integer                                :: i, j, increment
    type(swarm)                            :: tempswarm  ! temp = a(i) = swrm(i)%rdis !

    increment = size(swrm) / 2
    do while (increment > 0)
      do i = increment+1, size(swrm)
        j = i
        tempswarm = swrm(j)
        do while (j >= increment+1 .and. swrm(j-increment)%rdis > tempswarm%rdis)
          swrm(j) = swrm(j-increment)
          j = j - increment
        enddo
        swrm(j) = tempswarm
      enddo
      if (increment == 2) then
        increment = 1
      else
        increment = increment * 5 / 11
      endif
    enddo

    return
  end subroutine shell_sort_r

  ! sorting particles by id
  subroutine shell_sort_id(swrm)
    implicit none
    type(swarm), dimension(:), allocatable :: swrm
    integer                                :: i, j, increment
    type(swarm)                            :: tempswarm  ! temp = a(i) = swrm(i)%idnr !

    increment = size(swrm) / 2
    do while (increment > 0)
      do i = increment+1, size(swrm)
        j = i
        tempswarm = swrm(j)
        do while (j >= increment+1 .and. swrm(j-increment)%idnr > tempswarm%idnr)
          swrm(j) = swrm(j-increment)
          j = j - increment
        enddo
        swrm(j) = tempswarm
      enddo
      if (increment == 2) then
        increment = 1
      else
        increment = increment * 5 / 11
      endif
    enddo

    return
  end subroutine shell_sort_id

end module
