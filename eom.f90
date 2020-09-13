module eom
  use constants
  use constantsrkf
  contains

  subroutine calcacc(listy,acc,ND,soft,stnr,rd,Mch,switch)

    implicit none

    integer, parameter :: dp1 = SELECTED_REAL_KIND(15)
    integer :: ND
    real(dp1):: soft

    real(dp1),intent(out):: acc(:,:)
    real(dp1),intent(in):: listy(:,:)
    real(dp1) :: Fdrx(ND),Fdry(ND),Fdrz(ND),Fplx(ND),Fply(ND),Fplz(ND),rd(ND),stnr(ND)
    real(dp1) Mch(ND)
    integer :: switch

    rd = dsqrt(listy(:,4)**2 + listy(:,5)**2 + listy(:,6)**2 + soft**2)

    if (switch == 1) then
      !drag force
      Fdrx = -listy(:,1) / stnr
      Fdry = -listy(:,2) / stnr
      Fdrz = -listy(:,3) / stnr
    else
      Fdrx = 0
      Fdry = 0
      Fdrz = 0
    endif
    !gravity
    Fplx = -Ggrav * Mch * listy(:,4)/ rd**3
    Fply =-Ggrav * Mch * listy(:,5) / rd**3
    Fplz = -Ggrav * Mch *listy(:,6)/ rd**3
    acc(:,1) = Fdrx + Fplx
    acc(:,2) = Fdry + Fply
    acc(:,3) = Fdrz + Fplz
    acc(:,4) = listy(:,1)
    acc(:,5) = listy(:,2)
    acc(:,6) = listy(:,3)
    return
  end subroutine calcacc

  subroutine rkf(listy,acc,ND,soft,stnr,h,ynew,ynewfifth,rd,Mch,switch)

    implicit none

    INTEGER, PARAMETER :: dp1 = SELECTED_REAL_KIND(15)
    real(dp1) :: k1(ND,6)
    real(dp1) :: k2(ND,6)
    real(dp1) :: k3(ND,6)
    real(dp1) :: k4(ND,6)
    real(dp1) :: k5(ND,6)
    real(dp1):: k6(ND,6)
    real(dp1), intent(inout) :: ynew(:,:),ynewfifth(:,:)
    real(dp1), intent(in) :: listy(:,:)
    integer :: ND
    real(dp1) :: h
    real(dp1) :: soft
    real(dp1) :: stnr(ND)
    real(dp1) :: acc(ND,6)
    real(dp1):: rd(ND)
    real(dp1):: Mch(ND)
    integer :: switch

    call calcacc(listy,acc,ND,soft,stnr,rd,Mch,switch)

    k1 = h * acc

    call calcacc(listy + a1 * k1,acc,ND,soft,stnr,rd,Mch,switch)

    k2 = h  * acc

    call calcacc(listy + a3 * k2 + a4 * k2,acc,ND,soft,stnr,rd,Mch,switch)
    k3 = h  * acc

    call calcacc(listy+ a6 * k1 - a7 * k2 + a8 * k3,acc,ND,soft,stnr,rd,Mch,switch)
    k4 = h * acc

    call calcacc(listy + a9 * k1 -  a10 * k2 + a11 * k3 - a12 * k4,acc,ND,soft,stnr,rd,Mch,switch)
    k5 = h * acc

    call calcacc(listy - a13 * k1 + a14 * k2 - a15 * k3 + a16 * k4 - a17 * k5,acc,ND,soft,stnr,rd,Mch,switch)
    k6 = h * acc

    ynew = listy+ (a18 * k1 + a19 * k3 + a20 * k4 - a21 * k5)

    ynewfifth = listy + (a22 * k1 + a23 * k3 + a24 * k4 - a25 * k5 + a26* k6)

    return
  end subroutine rkf

end module eom
