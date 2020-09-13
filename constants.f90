module constants
  implicit none
  public

  real(kind=8), parameter         :: pi = 3.14159265358979323846
  real(kind=8), parameter         :: AU = 1.49597870691e+13
  real(kind=8), parameter         :: Msun = 1.9891e+33
  real(kind=8), parameter         :: Ggrav = 6.67384e-8
  real(kind=8), parameter         :: Froll = 8.5e-5
  real(kind=8), parameter         :: year = 365.256363051*24.0*3600.0
end module
