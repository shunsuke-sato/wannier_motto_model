module global_variables
  implicit none
! math paramters
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! physica constants
  real(8),parameter :: ev = 1d0/27.2114d0

! physical systems
  integer :: nkx
  real(8) :: dkx, kx_max
  real(8) :: eps_cutoff
  real(8) :: mu_mass
  real(8) :: v0_coulomb
  real(8),allocatable :: kx0(:)
  real(8),allocatable :: ham0(:,:), ham_t(:,:)
  real(8),allocatable :: eigvec_ham0(:,:)
  complex(8),allocatable :: zpsi(:)

  contains
    include "external_lib/specfun.f90"

end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call input
  call preparation

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
  dkx = 2d0*kx_max/nk
  

end subroutine input
!-------------------------------------------------------------------------------
subroutine preparation
  use global_variables
  implicit none
  integer :: ikx, ikx2
  real(8) :: kdist

  allocate(ham0(nkx,nkx), ham_t(nkx,nkx))
  allocate(eigvec_ham0(nkx,nkx))
  allocate(zpsi(nkx))
  allocatable(kx0(nkx))
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



end subroutine preparation
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
