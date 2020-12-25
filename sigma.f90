! Copyright 2020 Shunsuke A. Sato
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
program main
  implicit none
  integer,parameter :: nt = 38758
  integer,parameter :: nw = 1024
  real(8),parameter :: ev = 1d0/27.2114d0
  real(8),parameter :: wi = 45d0*ev, wf = 80d0*ev
  real(8),parameter :: dw = (wf-wi)/nw
  real(8),parameter :: gamma = 0.1d0*ev
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

!  it_ini = 0 ! impulsive

  open(30,file="zsigma.out")
  do iw = 0, nw
    ww = wi + dw*iw
    zjw = 0d0
    zEw = 0d0
    do it = it_ini, nt
      zfact = exp(zI*ww*tt(it)-gamma*(tt(it)-tt(it_ini)))
      zjw = zjw + zfact*jt(it)
      zEw = zEw + zfact*Eprobe(it)
    end do
    zjw = zjw*dt
    zEw = zEw*dt
!    zEw = 1d0 !impulsive
    zsigma = zjw/zEw
    write(30,"(999e26.16e3)")ww,zsigma, abs(zEw)

  end do
  close(30)


end program main
