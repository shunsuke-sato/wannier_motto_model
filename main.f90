! Wannier-Motto model for 1D
module global_variables
  implicit none
  public
! math paramters
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! physica constants

! unit
  real(8),parameter :: ev = 1d0/27.2114d0
  real(8),parameter :: fs = 1d0/0.024189d0

! physical systems
  integer :: nkx
  real(8) :: dkx, kx_max
  real(8) :: eps_cutoff
  real(8) :: mu_mass
  real(8) :: v0_coulomb
  real(8) :: eps_gap, pvc
  real(8),allocatable :: kx0(:)
  real(8),allocatable :: ham0(:,:), ham_t(:,:)
  real(8),allocatable :: eigval_ham0(:),eigvec_ham0(:,:)
  complex(8),allocatable :: zpsi(:)
  real(8),allocatable :: epsk_t(:)

! time propagation
  real(8) :: Tprop, dt
  real(8) :: Tprop_fs
  integer :: nt

! laser
  real(8),allocatable :: Apump(:),Eprobe(:)
  real(8) :: Tpump_fs, Tpump
  real(8) :: omega_pump_ev, omega_pump
  real(8) :: Tprobe_fs, Tprobe
  real(8) :: omega_probe_ev, omega_probe
  real(8) :: Tdelay_fs, Tdelay
  real(8) :: A0_pump, E0_probe
  integer :: n_field_type
  integer,parameter :: N_FIELD_TYPE_LASER_PULSE = 0
  integer,parameter :: N_FIELD_TYPE_IMPULSIVE_KICK = 1

  integer,parameter :: N_INIT_WF_GS = 0
  integer,parameter :: N_INIT_WF_ECITON_1s = 1
  integer :: n_init_wf = N_INIT_WF_ECITON_1s
  real(8),allocatable :: rho_rs_gs(:)
  

end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call input
  call preparation

  call time_propagation

end program main
!-------------------------------------------------------------------------------
subroutine input
  use global_variables
  implicit none

  eps_cutoff = 15d0*ev
  mu_mass = 0.24d0 ! light electron
!  mu_mass = 0.58d0 ! heavy electron

  nkx = 512+1

  eps_gap = 55.8d0*ev
  pvc = 1d0
  v0_coulomb = 0.187d0

! lasers
!  A0_pump = 1d-2
!  E0_probe = 1d-4
!  Tpump_fs = 16d0
!  Tprobe_fs = 1d0
!  omega_pump_ev = 1.55d0
!  omega_probe_ev = 55d0
!  Tdelay_fs = 0d0

  Tprop_fs = 150d0
  dt = 0.16d0


  read(*,*)A0_pump,E0_probe
  read(*,*)Tpump_fs,Tprobe_fs
  read(*,*)omega_pump_ev, omega_probe_ev
  read(*,*)Tdelay_fs


  kx_max = sqrt(2d0*mu_mass*eps_cutoff)
  dkx = 2d0*kx_max/nkx
  
  
  Tprop = Tprop_fs*fs
  nt = aint(Tprop/dt)+1

  n_field_type = N_FIELD_TYPE_LASER_PULSE
!  n_field_type = N_FIELD_TYPE_IMPULSIVE_KICK


  write(*,*)"nt=",nt

end subroutine input
!-------------------------------------------------------------------------------
subroutine preparation
  use global_variables
  implicit none
  integer :: ikx, ikx2
  real(8) :: kdist
!==LAPACK==
  real(8),allocatable :: work(:)
  integer :: lwork, info

  lwork = nkx**2+8*nkx
  allocate(work(lwork))
!==LAPACK==

  allocate(kx0(nkx))
  allocate(ham0(nkx,nkx), ham_t(nkx,nkx))
  allocate(eigval_ham0(nkx),eigvec_ham0(nkx,nkx))
  allocate(zpsi(nkx))

  zpsi = 0d0

! k-grid
  do ikx = 1, nkx
    kx0(ikx) = dkx*(ikx-nkx/2)
  end do

  ham0 = 0d0
  do ikx = 1, nkx
    ham0(ikx,ikx) = eps_gap + 0.5d0*kx0(ikx)**2/mu_mass
    do ikx2 = ikx+1,nkx
      kdist = abs(kx0(ikx) -kx0(ikx2))
      ham0(ikx,ikx2) = -2d0*v0_coulomb*(dkx/(2d0*pi))*besk0(kdist)
      ham0(ikx2,ikx) = ham0(ikx,ikx2)
    end do
  end do


  eigvec_ham0 = ham0
  call dsyev('V','U',nkx,eigvec_ham0,nkx,eigval_ham0,work,lwork,info)
  if(info /= 0)then
    write(*,*)'Error in LAPACK, info=',info
    stop
  end if

  open(20,file='eigenvalue.out')
  do ikx = 1, nkx
    write(20,"(I7,2x,e26.16e3)")ikx,eigval_ham0(ikx)
  end do
  close(20)

  open(20,file='eigenvector.out')
  do ikx = 1, nkx
    write(20,"(I7,2x,99999e26.16e3)")ikx,eigvec_ham0(ikx,:)
  end do
  close(20)

  contains
    include "external_lib/specfun.f90"

end subroutine preparation
!-------------------------------------------------------------------------------
subroutine time_propagation
  use global_variables
  implicit none
  integer :: it
  real(8) :: x_ave, x2_ave, Et_tmp
  real(8),allocatable :: jt(:)

  allocate(jt(0:nt+1))

  call init_laser
  call init_wf


  jt(0) = 2d0*pvc*sum(zpsi)*dkx/(2d0*pi)
  if(n_init_wf == N_INIT_WF_ECITON_1s)then  
    open(202,file='exciton_td.out')
  end if

  do it = 0, nt
    if(n_init_wf == N_INIT_WF_ECITON_1s)then
      if(mod(it,100)==0)call compute_real_space_wf(it,x_ave,x2_ave)
      Etmp = -0.5d0*(Apump(it+1)-Apump(it-1))/dt
      write(202,"(999e26.16e3)")dt*it,Etmp,x_ave,x2_ave
    end if
    call dt_evolve(it)
    jt(it+1) = 2d0*pvc*sum(zpsi)*dkx/(2d0*pi)
  end do
  open(20,file="Ac_Et_jt.out")
  write(20,"(A,2x,I9)")"#nt = ",nt
  write(20,"(A,2x,e26.16e3)")"#dtt = ",dt
  do it = 0, nt
    write(20,"(999e26.16e3)")dt*it, Apump(it), Eprobe(it), jt(it)
  end do
  close(20)

  if(n_init_wf == N_INIT_WF_ECITON_1s)then  
    close(202)
  end if

end subroutine time_propagation
!-------------------------------------------------------------------------------
subroutine init_laser
  use global_variables
  implicit none
  integer :: it
  real(8) :: tt, xx
  real(8) :: t_offset

  allocate(Apump(-1:nt+1), Eprobe(-1:nt+1))

  Tpump = Tpump_fs*fs
  omega_pump = omega_pump_ev*ev
  Tprobe = Tprobe_fs*fs
  omega_probe = omega_probe_ev*ev
  Tdelay = Tdelay_fs*fs

  t_offset = min(-0.5d0*Tpump, -0.5d0*Tprobe + Tdelay)

  
  Apump = 0d0
  Eprobe = 0d0
  if(n_field_type /= N_FIELD_TYPE_LASER_PULSE)then
    return
  end if

  do it = 0, nt+1
    tt = dt*it
! A-pump
    xx = tt + t_offset
    if(abs(xx)<=0.5d0*Tpump)then
      Apump(it) = A0_pump*sin(omega_pump*xx) *cos(pi*xx/Tpump)**2
    end if

! Eprobe
    xx = tt-Tdelay + t_offset
    if(abs(xx)<=0.5d0*Tprobe)then
      Eprobe(it) = E0_probe*sin(omega_probe*xx)*cos(pi*xx/Tprobe)**4
    end if

  end do

  open(20,file='lasers.out')
  do it = 0,nt+1
    write(20,"(999e26.16e3)")dt*it,Apump(it),Eprobe(it)
  end do
  close(20)
  

end subroutine init_laser
!-------------------------------------------------------------------------------
subroutine init_wf
  use global_variables
  implicit none
  integer :: ikx
  real(8) :: eps_t, ss

  zpsi = 0d0
  if(n_field_type == N_FIELD_TYPE_IMPULSIVE_KICK)then
    do ikx = 1, nkx
      eps_t = eps_gap + 0.5d0*kx0(ikx)**2/mu_mass
      zpsi(ikx) = pvc/eps_t
    end do
  end if

  if(n_init_wf == N_INIT_WF_ECITON_1s)then
    zpsi = eigvec_ham0(:,1)
    ss = sum(abs(zpsi)**2)*dkx
    zpsi = zpsi/sqrt(ss)
  end if

end subroutine init_wf
!-------------------------------------------------------------------------------
subroutine dt_evolve(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  real(8) :: eps_t, kxt(nkx), dip_t(nkx)
  real(8) :: vpot_t
  integer :: ikx

! apply E_probe field for time, t
  if(Eprobe(it) /=0d0)then

    kxt = kx0 + Apump(it)
    do ikx = 1, nkx
      eps_t = eps_gap + 0.5d0*kxt(ikx)**2/mu_mass
      dip_t(ikx) = pvc*Eprobe(it)/eps_t
    end do

    zpsi = zpsi +0.5d0*dt*dip_t

  end if

! apply A_pump field for time, t
  if(Apump(it) /=0d0)then

    eps_t = 0.5d0*Apump(it)**2/mu_mass
    do ikx= 1, nkx
      vpot_t = kx0(ikx)*Apump(it)/mu_mass + eps_t
      zpsi(ikx) = exp(-zI*0.5d0*dt*vpot_t)*zpsi(ikx)
    end do

  end if


! apply bare Hamiltonian for propagation, t-> t+dt
  zpsi = matmul(transpose(eigvec_ham0),zpsi)
  zpsi = exp(-zI*dt*eigval_ham0)*zpsi
  zpsi = matmul(eigvec_ham0,zpsi)

! apply A_pump field for time, t+dt
  if(Apump(it+1) /=0d0)then

    eps_t = 0.5d0*Apump(it+1)**2/mu_mass
    do ikx= 1, nkx
      vpot_t = kx0(ikx)*Apump(it+1)/mu_mass + eps_t
      zpsi(ikx) = exp(-zI*0.5d0*dt*vpot_t)*zpsi(ikx)
    end do

  end if

! apply E_probe field for time, t
  if(Eprobe(it+1) /=0d0)then

    kxt = kx0 + Apump(it+1)
    do ikx = 1, nkx
      eps_t = eps_gap + 0.5d0*kxt(ikx)**2/mu_mass
      dip_t(ikx) = pvc*Eprobe(it+1)/eps_t
    end do

    zpsi = zpsi +0.5d0*dt*dip_t

  end if


end subroutine dt_evolve
!-------------------------------------------------------------------------------
subroutine compute_real_space_wf(int_in, x_ave, x2_ave)
  use global_variables
  implicit none
  integer,intent(in) :: int_in
  real(8),intent(out) :: x_ave, x2_ave
  integer :: ix, nx, ikx
  real(8) :: length_x, dx, xx, rr, ss
  character(256) :: cit, cfilename
  complex(8),allocatable :: zpsi_rs(:)
  real(8),allocatable :: rho_rs(:)
  complex(8) :: zs


  length_x = 50d0
  nx = 500
  dx = length_x/nx

  allocate(zpsi_rs(-nx:nx))
  allocate(rho_rs(-nx:nx))
    


  do ix = -nx, nx
    xx = dx*ix

    zs = 0d0
    do ikx = 1, nkx
      zs = zs + zpsi(ikx)*exp(zi*kx0(ikx)*xx)
    end do

    zs = zs*dkx/(2d0*pi)
    rr = abs(zs)**2
!    write(201,"(999e26.16e3)")xx,rr
    rho_rs(ix)

  end do


  if(.not. allocated(rho_rs_gs))then
    allocate(rho_rs_gs(-nx:nx))
    rho_rs_gs = rho_rs
  end if

  write(cit,"(I7.7)")int_in
  cfilename = trim(cit)//'_rho_exciton.out'
  
  open(201,file=cfilename)
  do ix = -nx, nx
    xx = dx*ix
    write(201,"(999e26.16e3)")xx,rho_rs(ix),rho_rs(ix)-rho_rs_gs(ix)
  end do
  close(201)


  x_ave = 0d0
  x2_ave = 0d0
  ss = 0d0
  do ix = -nx, nx
    xx = dx*ix
    ss = ss + rho_rs(ix)
    x_ave = x_ave + rho_rs(ix)*xx
    x2_ave = x2_ave + rho_rs(ix)*xx**2
  end do

  x_ave = x_ave/ss
  x2_ave = x2_ave/ss

end subroutine compute_real_space_wf
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
