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
  real(8),allocatable :: kx0(:)
  real(8),allocatable :: ham0(:,:), ham_t(:,:)
  real(8),allocatable :: eigval_ham0(:),eigvec_ham0(:,:)
  complex(8),allocatable :: zpsi(:)

! time propagation
  real(8) :: Tprop, dt
  integer :: nt

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
  mu_mass = 1d0
  nkx = 256+1
  v0_coulomb = 1d0

  kx_max = sqrt(2d0*mu_mass*eps_cutoff)
  dkx = 2d0*kx_max/nkx
  

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
    ham0(ikx,ikx) = 0.5d0*kx0(ikx)**2/mu_mass
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

  

end subroutine time_propagation
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
