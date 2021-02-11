module global_variables
  implicit none

! mathematical constants
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! physical constants
  real(8),parameter :: ev = 1d0/27.2114d0
  real(8),parameter :: fs = 1d0/0.024189d0

! model paramters
  real(8) :: Egap, Egap_23
  real(8) :: d_12, d_23

  complex(8) :: zrho_dm(3,3)
  real(8) :: Hmat(3,3)
!  complex(8) :: zHmat(3,3)

! relaxation paramters
  real(8) :: T2_12, T2_13, T2_23

! parameters for time-propagation
  integer :: nt
  real(8) :: dt, Tprop
  

! laser paraemter
  real(8) :: E0_1, omega0_1, tpulse_1
  real(8),allocatable :: Et_1(:),Et_1_dt2(:)
  real(8) :: E0_2, omega0_2, tpulse_2, tdelay
  real(8),allocatable :: Et_2(:),Et_2_dt2(:)

! Floquet decomposition
  integer,parameter :: ndim_F = 2
  real(8),allocatable :: Et_1_env(:), phi_1(:) ! Et_1 = Et_1_env*cos(omega0_1*t+phi_1)
  complex(8),allocatable :: zpsi_F(:,:),zpsi_F_old(:,:),zpsi_F_new(:,:)
  complex(8),allocatable :: zham_F(:,:)
  real(8) :: eps_F(2), eps_F_old(2), eps_F_new(2)

end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call input
  call initialize
!  call time_propagation
  call time_propagation_floquet_decomp

end program main
!-------------------------------------------------------------------------------
subroutine input
  use global_variables
  implicit none
  real(8) :: omega0_1_ev, tpulse_1_fs
  real(8) :: omega0_2_ev, tpulse_2_fs, tdelay_fs

  Egap = 0.1999162505324391E+001 
  Egap_23 = 0.2047462404435500E+001 -Egap

  T2_12 = 1d0/(0.25d0*ev)

  d_12 = 1d0 !sqrt(0.14d0/(4d0*pi*T2_12))
  d_23 = 0.3493516210248626d1



  Tprop = 180d0*fs
  dt = 0.1d0
  nt = aint(Tprop/dt) + 1
  write(*,*)'nt=',nt

  open(20,file='input')
  read(20,*)E0_1, omega0_1_ev, tpulse_1_fs
!  E0_1 = 0d-2
  omega0_1 = omega0_1_ev*ev
  tpulse_1 = tpulse_1_fs*fs

  read(20,*)E0_2, omega0_2_ev, tpulse_2_fs
  read(20,*)tdelay_fs
!  E0_2 = 1d-4
  omega0_2 = omega0_2_ev*ev
  tpulse_2 = tpulse_2_fs*fs

  tdelay = tdelay_fs*fs

  close(20)

end subroutine input
!-------------------------------------------------------------------------------
subroutine initialize
  use global_variables
  implicit none

  zrho_dm = 0d0
  zrho_dm(1,1) = 1d0

  call init_laser

end subroutine initialize
!-------------------------------------------------------------------------------
subroutine init_laser
  use global_variables
  implicit none
  integer :: it
  real(8) :: xx, tt
  real(8) :: t_offset

  t_offset = 20d0*fs

  allocate(Et_1(-1:nt+1),Et_1_dt2(-1:nt+1))
  allocate(Et_2(-1:nt+1),Et_2_dt2(-1:nt+1))
  Et_1 = 0d0; Et_1_dt2 = 0d0
  Et_2 = 0d0; Et_2_dt2 = 0d0

! pump pulse
  do it = 0, nt+1
    tt = dt*it
    xx = tt - 0.5d0*tpulse_1 -t_offset
    if(abs(xx)<0.5d0*tpulse_1)then
      Et_1(it) = E0_1/omega0_1*(&
        -2d0*pi/tpulse_1*cos(pi*xx/tpulse_1)*sin(pi*xx/tpulse_1)*sin(omega0_1*xx) &
        +omega0_1*cos(pi*xx/tpulse_1)**2*cos(omega0_1*xx) )
    end if

    tt = dt*it+0.5d0*dt
    xx = tt - 0.5d0*tpulse_1 -t_offset
    if(abs(xx)<0.5d0*tpulse_1)then
      Et_1_dt2(it) = E0_1/omega0_1*(&
        -2d0*pi/tpulse_1*cos(pi*xx/tpulse_1)*sin(pi*xx/tpulse_1)*sin(omega0_1*xx) &
        +omega0_1*cos(pi*xx/tpulse_1)**2*cos(omega0_1*xx) )
    end if

  end do

! probe pulse
  do it = 0, nt+1
    tt = dt*it
    xx = tt - 0.5d0*tpulse_1 - tdelay -t_offset
    if(abs(xx)<0.5d0*tpulse_2)then
      Et_2(it) = E0_2/omega0_2*( &
        -4d0*pi/tpulse_2*cos(pi*xx/tpulse_2)**3*sin(pi*xx/tpulse_2)*sin(omega0_2*xx) &
        +omega0_2*cos(pi*xx/tpulse_2)**4*cos(omega0_2*xx) )
    end if

    tt = dt*it+0.5d0*dt
    xx = tt - 0.5d0*tpulse_1 - tdelay -t_offset
    if(abs(xx)<0.5d0*tpulse_2)then
      Et_2_dt2(it) = E0_2/omega0_2*( &
        -4d0*pi/tpulse_2*cos(pi*xx/tpulse_2)**3*sin(pi*xx/tpulse_2)*sin(omega0_2*xx) &
        +omega0_2*cos(pi*xx/tpulse_2)**4*cos(omega0_2*xx) )
    end if

  end do


end subroutine init_laser
!-------------------------------------------------------------------------------
subroutine time_propagation
  use global_variables
  implicit none
  integer :: it
  real(8) :: dipole

  open(20,file='Et_dipole.out')
  write(20,"(A,2x,I7)")'#nt=',nt
  
  it = 0
  dipole = -2d0*d_12*real(zrho_dm(1,2))
  write(20,"(999e26.16e3)")dt*it,Et_1(it),Et_2(it),dipole

  do it = 0, nt

    call dt_evolve(it)
    dipole = -2d0*d_12*real(zrho_dm(1,2))
    write(20,"(999e26.16e3)")dt*(it+1),Et_1(it+1),Et_2(it+1),dipole


  end do


  close(20)

end subroutine time_propagation
!-------------------------------------------------------------------------------
subroutine time_propagation_floquet_decomp
  use global_variables
  implicit none
  integer :: it
  real(8) :: dipole

  call calc_envelope_and_phase
  call init_floquet

  open(20,file='Et_dipole.out')
  write(20,"(A,2x,I7)")'#nt=',nt
  
  it = 0
  dipole = -2d0*d_12*real(zrho_dm(1,2))
  write(20,"(999e26.16e3)")dt*it,Et_1(it),Et_2(it),dipole

  do it = 0, nt

    call dt_evolve_floquet(it)
    dipole = -2d0*d_12*real(zrho_dm(1,2))
    write(20,"(999e26.16e3)")dt*(it+1),Et_1(it+1),Et_2(it+1),dipole


  end do


  close(20)

end subroutine time_propagation_floquet_decomp
!-------------------------------------------------------------------------------
subroutine dt_evolve(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer :: irk
  real(8) :: H12, H23
  complex(8) :: zrho_rk(3,3,4)
  complex(8) :: zrho_t(3,3)

  Hmat = 0d0
  Hmat(1,1) = -0.5d0*Egap
  Hmat(2,2) =  0.5d0*Egap
  Hmat(3,3) =  0.5d0*Egap+Egap_23


! at time, t
  H23 = Et_1(it)*d_23
  H12 = Et_2(it)*d_12
  Hmat(1,2) = H12; Hmat(2,1) = H12
  Hmat(2,3) = H23; Hmat(3,2) = H23

!RK1
  irk = 1
  zrho_t = zrho_dm
  zrho_rk(:,:,irk) = -zi*(matmul(Hmat,zrho_t)-matmul(zrho_t,Hmat))
  zrho_rk(1,2,irk) = zrho_rk(1,2,irk) -zrho_dm(1,2)/T2_12
  zrho_rk(2,1,irk) = zrho_rk(2,1,irk) -zrho_dm(2,1)/T2_12

! at time, t+dt/2
  H23 = Et_1_dt2(it)*d_23
  H12 = Et_2_dt2(it)*d_12
  Hmat(1,2) = H12; Hmat(2,1) = H12
  Hmat(2,3) = H23; Hmat(3,2) = H23

!RK2
  irk = 2
  zrho_t = zrho_dm + 0.5d0*dt*zrho_rk(:,:,1)
  zrho_rk(:,:,irk) = -zi*(matmul(Hmat,zrho_t)-matmul(zrho_t,Hmat))
  zrho_rk(1,2,irk) = zrho_rk(1,2,irk) -zrho_dm(1,2)/T2_12
  zrho_rk(2,1,irk) = zrho_rk(2,1,irk) -zrho_dm(2,1)/T2_12

!RK3
  irk = 3
  zrho_t = zrho_dm + 0.5d0*dt*zrho_rk(:,:,2)
  zrho_rk(:,:,irk) = -zi*(matmul(Hmat,zrho_t)-matmul(zrho_t,Hmat))
  zrho_rk(1,2,irk) = zrho_rk(1,2,irk) -zrho_dm(1,2)/T2_12
  zrho_rk(2,1,irk) = zrho_rk(2,1,irk) -zrho_dm(2,1)/T2_12

! at time, t+dt
  H23 = Et_1(it+1)*d_23
  H12 = Et_2(it+1)*d_12
  Hmat(1,2) = H12; Hmat(2,1) = H12
  Hmat(2,3) = H23; Hmat(3,2) = H23

!RK4
  irk = 4
  zrho_t = zrho_dm + dt*zrho_rk(:,:,3)
  zrho_rk(:,:,irk) = -zi*(matmul(Hmat,zrho_t)-matmul(zrho_t,Hmat))
  zrho_rk(1,2,irk) = zrho_rk(1,2,irk) -zrho_dm(1,2)/T2_12
  zrho_rk(2,1,irk) = zrho_rk(2,1,irk) -zrho_dm(2,1)/T2_12


  zrho_dm = zrho_dm + dt/6d0*(zrho_rk(:,:,1) &
                         +2d0*zrho_rk(:,:,2) &
                         +2d0*zrho_rk(:,:,3) &
                         +    zrho_rk(:,:,4))

end subroutine dt_evolve
!-------------------------------------------------------------------------------
subroutine calc_envelope_and_phase
  use global_variables
  implicit none
  integer :: nt_pulse
  real(8),allocatable :: Et_tmp(:)
  complex(8),allocatable :: zEw_tmp(:),zEt_tmp(:)
  integer :: it, iw

  allocate(Et_1_env(-1:nt+1), phi_1(-1:nt+1))

  nt_pulse = 0
  do it = 0, nt+1
    if(Et_1(it)/= 0d0)nt_pulse = it
  end do

  allocate(Et_tmp(0:nt_pulse),zEw_tmp(0:nt_pulse),zEt_tmp(0:nt_pulse))
  Et_tmp(0:nt_pulse) = Et_1(0:nt_pulse)

  zEw_tmp = 0d0
  do iw = 0, nt_pulse
    do it = 0, nt_pulse
      zEw_tmp(iw) = zEw_tmp(iw) + exp(zi*2d0*pi*it*iw/dble(nt_pulse+1))*Et_tmp(it)
    end do
  end do

! Hilbert transform
  zEt_tmp = 0d0
  do it = 0, nt_pulse
    do iw = 0, nt_pulse/2
      zEt_tmp(it) = zEt_tmp(it) + exp(-zi*2d0*pi*it*iw/dble(nt_pulse+1))*zEw_tmp(iw)
    end do
  end do
  zEt_tmp = zEt_tmp/dble(nt_pulse+1)

  Et_1_env = 0d0
  phi_1 = 0d0
  do it = 0, nt_pulse
    Et_1_env(it) = abs(zEt_tmp(it))
!    if(Et_1_env(it)/=0d0) phi_1(it) = acos(real(zEt_tmp(it))/Et_1_env(it))-omega0_1*dt*it
    if(Et_1_env(it)/=0d0) phi_1(it) = -aimag(log(zEt_tmp(it)/Et_1_env(it)))-omega0_1*dt*it
  end do
  Et_1_env = Et_1_env*2d0 

  phi_1(-1) = phi_1(0)
  phi_1(nt_pulse+1:nt+1) = phi_1(nt_pulse)

  open(20,file='laser_chk.out')
  do it = 0, nt
    write(20,"(999e26.16e3)")dt*it,Et_1(it),Et_1_env(it)&
      ,Et_1_env(it)*cos(omega0_1*dt*it+phi_1(it)),phi_1(it)
  end do
  close(20)

end subroutine calc_envelope_and_phase
!-------------------------------------------------------------------------------
subroutine init_floquet
  use global_variables
  implicit none

  allocate(zpsi_F(2*(2*ndim_F+1),2) &
          ,zpsi_F_old(2*(2*ndim_F+1),2) &
          ,zpsi_F_new(2*(2*ndim_F+1),2) )

  allocate(zham_F(2*(2*ndim_F+1),2*(2*ndim_F+1)))

  eps_F_old(1) = 0.5d0*Egap
  eps_F_old(2) = 0.5d0*Egap+Egap_23

  zpsi_F_old = 0d0
  zpsi_F_old(2*ndim_F+1,1) = 1d0
  zpsi_F_old(2*ndim_F+2,2) = 1d0

  eps_F = eps_F_old
  zpsi_F = zpsi_F_old


end subroutine init_floquet
!-------------------------------------------------------------------------------
subroutine calc_floquet(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  real(8) :: H22, H33
  complex(8) :: zH23
  integer :: ifloquet, jfloquet
  complex(8) :: zvec(2,2), zvec_new(2,2)
  real(8) :: tt
  complex(8) :: zs
! lapack arrays
  integer :: lwork,info
  complex(8),allocatable :: work(:)
  real(8),allocatable    :: rwork(:),w(:)

  lwork = 8*(2*(2*ndim_f+1))**2
  allocate(w(2*(2*ndim_f+1)))
  allocate(work(lwork), rwork(3*2*(2*ndim_f+1)))
! lapack arrays


  H22  = 0.5d0*Egap
  H33  = 0.5d0*Egap+Egap_23
  zH23 = 0.5d0*d_23*Et_1_env(it)*exp(zi*phi_1(it))

  zham_F = 0d0
  do ifloquet = 1, 2*ndim_F+1
    do jfloquet = 1, 2*ndim_F+1
      if(ifloquet == jfloquet)then
        zham_F(2*(ifloquet-1)+1,2*(ifloquet-1)+1) = H22+omega0_1*(ifloquet-ndim_F-1)
        zham_F(2*(ifloquet-2)+1,2*(ifloquet-1)+2) = H33+omega0_1*(ifloquet-ndim_F-1)
      else if(ifloquet == jfloquet-1)then
        zham_F(2*(ifloquet-1)+1,2*(jfloquet-1)+2) = zH23
        zham_F(2*(ifloquet-1)+2,2*(jfloquet-1)+1) = zH23
      else if(ifloquet == jfloquet+1)then
        zham_F(2*(ifloquet-1)+1,2*(jfloquet-1)+2) = conjg(zH23)
        zham_F(2*(ifloquet-1)+2,2*(jfloquet-1)+1) = conjg(zH23)
      end if
    end do
  end do

  call zheev('V', 'U', 2*(2*ndim_F+1), zham_F, 2*(2*ndim_F+1), w, work, lwork, rwork, info)

  zpsi_F_new(:,1) = zham_F(:,2*ndim_F+1)
  zpsi_F_new(:,2) = zham_F(:,2*ndim_F+2)
  eps_F_new(1)    = w(2*ndim_F+1)
  eps_F_new(2)    = w(2*ndim_F+2)

  tt = it*dt
  zvec = 0d0
  zvec_new = 0d0
  do ifloquet = 1, 2*ndim_F+1
    zvec_new(1:2,1) = zvec_new(1:2,1) &
      + exp(-zi*(eps_F_new(1)+omega0_1*(ifloquet-1-ndim_F))*tt)&
      * zpsi_F_new(2*(ifloquet-1)+1:2*(ifloquet-1)+2,1)
    zvec_new(1:2,2) = zvec_new(1:2,2) &
      + exp(-zi*(eps_F_new(2)+omega0_1*(ifloquet-1-ndim_F))*tt)&
      * zpsi_F_new(2*(ifloquet-1)+1:2*(ifloquet-1)+2,2)

    zvec(1:2,1) = zvec(1:2,1) &
      + exp(-zi*(eps_F(1)+omega0_1*(ifloquet-1-ndim_F))*tt)&
      * zpsi_F(2*(ifloquet-1)+1:2*(ifloquet-1)+2,1)
    zvec(1:2,2) = zvec(1:2,2) &
      + exp(-zi*(eps_F(2)+omega0_1*(ifloquet-1-ndim_F))*tt)&
      * zpsi_F(2*(ifloquet-1)+1:2*(ifloquet-1)+2,2)
  end do

  zs = sum(conjg(zvec(:,1))*zvec_new(:,1))
  zpsi_F_new(:,1) = zpsi_F_new(:,1)*exp(-zi*aimag(log(zs)))
  zs = sum(conjg(zvec(:,1))*zvec_new(:,2))
  zpsi_F_new(:,2) = zpsi_F_new(:,2)*exp(-zi*aimag(log(zs)))


end subroutine calc_floquet
!-------------------------------------------------------------------------------
subroutine dt_evolve_floquet(it)
  use global_variables
  implicit none
  integer,intent(in):: it
  complex(8) :: zrho_t(3,3), zrho_F(3,3)
  complex(8) :: zvec_new(2,2), zvec(2,2), zvec_old(2,2), zdvec(2,2)
  complex(8) :: zUvec(3,3), zham_org(3,3), zham_eff(3,3)
  integer :: ifloquet
  real(8) :: H12, tt

! relaxation
  zrho_dm(1,2) = zrho_dm(1,2) -0.5d0*dt*zrho_dm(1,2)/T2_12
  zrho_dm(2,1) = zrho_dm(2,1) -0.5d0*dt*zrho_dm(2,1)/T2_12


! == START: propagation from t to t+dt/2 ==
  tt = dt*it

  do ifloquet = 1, 2*ndim_F+1
    zvec_new(1:2,1) = zvec_new(1:2,1) &
      + exp(-zi*(eps_F_new(1)+omega0_1*(ifloquet-1-ndim_F))*tt)&
      * zpsi_F_new(2*(ifloquet-1)+1:2*(ifloquet-1)+2,1)
    zvec_new(1:2,2) = zvec_new(1:2,2) &
      + exp(-zi*(eps_F_new(2)+omega0_1*(ifloquet-1-ndim_F))*tt)&
      * zpsi_F_new(2*(ifloquet-1)+1:2*(ifloquet-1)+2,2)

    zvec(1:2,1) = zvec(1:2,1) &
      + exp(-zi*(eps_F(1)+omega0_1*(ifloquet-1-ndim_F))*tt)&
      * zpsi_F(2*(ifloquet-1)+1:2*(ifloquet-1)+2,1)
    zvec(1:2,2) = zvec(1:2,2) &
      + exp(-zi*(eps_F(2)+omega0_1*(ifloquet-1-ndim_F))*tt)&
      * zpsi_F(2*(ifloquet-1)+1:2*(ifloquet-1)+2,2)

    zvec_old(1:2,1) = zvec_old(1:2,1) &
      + exp(-zi*(eps_F_old(1)+omega0_1*(ifloquet-1-ndim_F))*tt)&
      * zpsi_F_old(2*(ifloquet-1)+1:2*(ifloquet-1)+2,1)
    zvec_old(1:2,2) = zvec_old(1:2,2) &
      + exp(-zi*(eps_F_old(2)+omega0_1*(ifloquet-1-ndim_F))*tt)&
      * zpsi_F_old(2*(ifloquet-1)+1:2*(ifloquet-1)+2,2)
  end do
  zdvec = 0.5d0/dt*(zvec_new-zvec_old)


  zUvec = 0d0
  zUvec(2:3,2:3) = zvec(1:2,1:2)
  zUvec(1,1) = exp(-zi*(-0.5d0*Egap*tt))

  H12 = Et_2(it)*d_12
  zham_org(1,2) = H12
  zham_org(2,1) = H12

  zham_eff = matmul(transpose(conjg(zUvec)),matmul(zham_org,zUvec))
  zham_eff(2,3) = zham_eff(2,3) -zi*sum(conjg(zvec(:,1))*zdvec(:,2))
  zham_eff(3,2) = conjg(zham_eff(2,3))


  zrho_F = matmul(transpose(conjg(zUvec)),matmul(zrho_dm,zUvec))
  
  zrho_t = matmul(zham_eff,zrho_F)-matmul(zrho_F,zham_eff)
  zrho_F = zrho_F -zi*0.5d0*dt*zrho_t

! == END: propagation from t to t+dt/2 ==
! == START: propagation from t+dt/2 to t+dt ==
  zpsi_F_old = zpsi_F
  zpsi_F = zpsi_F_new
  call calc_floquet(it+1)

  tt = dt*(it+1)
  do ifloquet = 1, 2*ndim_F+1
    zvec_new(1:2,1) = zvec_new(1:2,1) &
      + exp(-zi*(eps_F_new(1)+omega0_1*(ifloquet-1-ndim_F))*tt)&
      * zpsi_F_new(2*(ifloquet-1)+1:2*(ifloquet-1)+2,1)
    zvec_new(1:2,2) = zvec_new(1:2,2) &
      + exp(-zi*(eps_F_new(2)+omega0_1*(ifloquet-1-ndim_F))*tt)&
      * zpsi_F_new(2*(ifloquet-1)+1:2*(ifloquet-1)+2,2)

    zvec(1:2,1) = zvec(1:2,1) &
      + exp(-zi*(eps_F(1)+omega0_1*(ifloquet-1-ndim_F))*tt)&
      * zpsi_F(2*(ifloquet-1)+1:2*(ifloquet-1)+2,1)
    zvec(1:2,2) = zvec(1:2,2) &
      + exp(-zi*(eps_F(2)+omega0_1*(ifloquet-1-ndim_F))*tt)&
      * zpsi_F(2*(ifloquet-1)+1:2*(ifloquet-1)+2,2)

    zvec_old(1:2,1) = zvec_old(1:2,1) &
      + exp(-zi*(eps_F_old(1)+omega0_1*(ifloquet-1-ndim_F))*tt)&
      * zpsi_F_old(2*(ifloquet-1)+1:2*(ifloquet-1)+2,1)
    zvec_old(1:2,2) = zvec_old(1:2,2) &
      + exp(-zi*(eps_F_old(2)+omega0_1*(ifloquet-1-ndim_F))*tt)&
      * zpsi_F_old(2*(ifloquet-1)+1:2*(ifloquet-1)+2,2)
  end do
  zdvec = 0.5d0/dt*(zvec_new-zvec_old)


  zUvec = 0d0
  zUvec(2:3,2:3) = zvec(1:2,1:2)
  zUvec(1,1) = exp(-zi*(-0.5d0*Egap*tt))

  H12 = Et_2(it)*d_12
  zham_org(1,2) = H12
  zham_org(2,1) = H12

  zham_eff = matmul(transpose(conjg(zUvec)),matmul(zham_org,zUvec))
  zham_eff(2,3) = zham_eff(2,3) -zi*sum(conjg(zvec(:,1))*zdvec(:,2))
  zham_eff(3,2) = conjg(zham_eff(2,3))


  zrho_F = matmul(transpose(conjg(zUvec)),matmul(zrho_dm,zUvec))
  
  zrho_t = matmul(zham_eff,zrho_F)-matmul(zrho_F,zham_eff)
  zrho_F = zrho_F -zi*0.5d0*dt*zrho_t
! == END: propagation from t+dt/2 to t+dt ==  

! relaxation
  zrho_dm(1,2) = zrho_dm(1,2) -0.5d0*dt*zrho_dm(1,2)/T2_12
  zrho_dm(2,1) = zrho_dm(2,1) -0.5d0*dt*zrho_dm(2,1)/T2_12

end subroutine dt_evolve_floquet
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
