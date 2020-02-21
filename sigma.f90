program main
  implicit none
  integer,parameter :: nt = 51677
  integer,parameter :: nw = 1024
  real(8),parameter :: ev = 1d0/27.2114d0
  real(8),parameter :: wi = 45d0*ev, wf = 65d0*ev
  real(8),parameter :: dw = (wf-wi)/nw
  real(8),parameter :: gamma = 0.2d0*ev
  complex(8),parameter :: zi = (0d0, 1d0)
  real(8) :: tt(0:nt),Apump(0:Nt),Eprobe(0:nt),jt(0:nt)
  real(8) :: ww, dt
  integer :: iw, it, it_ini
  complex(8) :: zsigma, zjw, zEw, zfact

  open(20,file="Ac_Et_jt.out")
  read(20,*); read(20,*)
  do it = 0, nt
    read(20,*)tt(it),Apump(it),Eprobe(it),jt(it)
  end do
  close(20)
  dt = tt(1)-tt(0)
  
  do it = 0, nt
    if(Eprobe(it) /=0)then
      it_ini = it
      exit
    end if
  end do

  it_ini = 0

  open(30,file="zsigma.out")
  do iw = 0, nw
    ww = wi + dw*iw
    zjw = 0d0
!    zEw = 0d0
    do it = it_ini, nt
      zfact = exp(zI*ww*tt(it)-gamma*(tt(it)-tt(it_ini)))
      zjw = zjw + zfact*jt(it)
!      zEw = zEw + zfact*Eprobe(it)
    end do
    zjw = zjw*dt
    zEw = 1d0*dt
    zsigma = zjw/zEw
    write(30,"(999e26.16e3)")ww,zsigma, abs(zEw)

  end do
  close(30)


end program main
