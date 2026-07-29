! mod_perturb_param.f90
! =====================================
!  Module for Parameter Perturbations
! ======================================

module mod_perturb_param
  use common_pom_var
  use MYNNF_lev25_2012, only: alp2
  implicit none
  private

  public :: init_perturb_alp2, update_perturb_alp2

  ! Flag for random seed initialization (for Update mode)
  logical, save :: is_seed_initialized = .false.

contains

  ! ===================================================================
  !  Initialization: Set initial parameter fields (Cold Start Only)
  ! ===================================================================
  subroutine init_perturb_alp2(im, jm, irestart, nens, iens)
    integer, intent(in) :: im, jm
    logical, intent(in) :: irestart
    integer, intent(in) :: nens, iens

    real(kind = r_size) :: assigned_val
    real(kind = r_size), parameter :: alp2_default = 0.53d0

    if (.not. allocated(alp2)) allocate(alp2(im, jm))
    
    ! Do nothing for restart runs to preserve LETKF-updated values
    if (irestart) return 

    ! Set initial values for cold start
    select case (i_init_alp2)
    case (1)
       ! Pattern 1: Linspace (Equally spaced spread)
       if (nens > 1) then
          assigned_val = alp2_range_min + &
               (alp2_range_max - alp2_range_min) * real(iens - 1, kind=r_size) / real(nens - 1, kind=r_size)
       else
          assigned_val = alp2_default 
       end if
       alp2(:,:) = assigned_val

    case default
       ! Pattern 0: Common default value for all members
       alp2(:,:) = alp2_default

    end select

  end subroutine init_perturb_alp2


  ! ===================================================================
  !  Update: Time evolution and Physical Constraints
  ! ===================================================================
  subroutine update_perturb_alp2(im, jm, iens)
    integer, intent(in) :: im, jm
    integer, intent(in) :: iens

    integer :: i, j, seed_size
    integer, allocatable :: seed_array(:)
    real(kind = r_size) :: rand_val
    real(kind = r_size), allocatable :: rand_field(:,:)
    
    real(kind = r_size), parameter :: pert_amp = 0.1d0
    real(kind = r_size), parameter :: alp2_min = 0.1d0
    real(kind = r_size), parameter :: alp2_max = 1.0d0

    ! --- 1. Add time-varying noise (Modes 1 and 2) ---
    if (i_update_alp2 == 1 .or. i_update_alp2 == 2) then
       
       ! Setup random seed (First time only)
       if (.not. is_seed_initialized) then
          call random_seed(size=seed_size)
          allocate(seed_array(seed_size))
          seed_array = iens * 1000000 
          call random_seed(put=seed_array)
          deallocate(seed_array)
          is_seed_initialized = .true.
       end if

       if (i_update_alp2 == 2) then
          ! Mode 2: 2D spatial random noise
          allocate(rand_field(im, jm))
          call random_number(rand_field)
          do j = 1, jm
             do i = 1, im
                alp2(i,j) = alp2(i,j) + (rand_field(i,j) - 0.5d0) * 2.0d0 * pert_amp
             end do
          end do
          deallocate(rand_field)
       else
          ! Mode 1: 1D uniform random noise
          call random_number(rand_val)
          alp2(:,:) = alp2(:,:) + (rand_val - 0.5d0) * 2.0d0 * pert_amp
       end if
    end if
    ! (Skip noise addition if i_update_alp2 == 0)

    ! --- 2. Physical constraints (Clipping) ---
    ! Ensure value safety at every step regardless of the mode
    do j = 1, jm
       do i = 1, im
          alp2(i,j) = max(alp2_min, min(alp2_max, alp2(i,j)))
       end do
    end do

  end subroutine update_perturb_alp2

end module mod_perturb_param
