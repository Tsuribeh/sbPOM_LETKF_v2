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

  ! Internal flag to remember if the perturbation has already been applied
  ! (Useful for time-constant perturbation mode)
  logical, save :: is_alp2_perturbed = .false.

contains

  ! ===================================================================
  !  Initialization: MYNN alp2 parameter
  ! ===================================================================
  subroutine init_perturb_alp2(im, jm, irestart)
    integer, intent(in) :: im, jm
    logical, intent(in) :: irestart ! .true. if this is a restart run

    ! Allocate alp2 array if not yet allocated
    if (.not. allocated(alp2)) allocate(alp2(im, jm))
    
    ! Apply default values ONLY for cold starts.
    ! During a restart run (e.g., in LETKF parameter estimation cycles), 
    ! alp2 should have already been read from the restart file, 
    ! so we must NOT overwrite it here.
    if (.not. irestart) then
       alp2(:,:) = 0.53d0
       is_alp2_perturbed = .false.
    else
       is_alp2_perturbed = .true.
    end if
    
  end subroutine init_perturb_alp2


  ! ===================================================================
  !  Update: Apply perturbations to MYNN alp2 parameter
  ! ===================================================================
  subroutine update_perturb_alp2(im, jm, nens, iens, step_count)
    integer, intent(in) :: im, jm
    integer, intent(in) :: nens, iens, step_count

    integer :: i, j, seed_size
    integer, allocatable :: seed_array(:)
    real(kind = r_size) :: rand_val, assigned_val
    real(kind = r_size), allocatable :: rand_field(:,:)
    
    ! Perturbation settings (Can be moved to namelist in the future)
    real(kind = r_size), parameter :: alp2_default = 0.53d0
    real(kind = r_size), parameter :: pert_amp     = 0.1d0
    real(kind = r_size), parameter :: alp2_min     = 0.1d0
    real(kind = r_size), parameter :: alp2_max     = 1.0d0

    ! ---------------------------------------------------------
    ! Condition 3: Is perturbation enabled at all?
    ! ---------------------------------------------------------
    if (.not. l_pert_alp2) return 

    ! ---------------------------------------------------------
    ! Linspace Mode (Equally spaced assignments)
    ! ---------------------------------------------------------
    if (l_pert_alp2_linspace) then
       if (.not. is_alp2_perturbed) then
          if (nens > 1) then
             ! Assign values linearly from min to max across ensemble members
             assigned_val = alp2_range_min + &
                  (alp2_range_max - alp2_range_min) * real(iens - 1, kind=r_size) / real(nens - 1, kind=r_size)
          else
             ! Failsafe for single runs
             assigned_val = alp2_default 
          end if
          
          alp2(:,:) = assigned_val
          is_alp2_perturbed = .true.
       end if
       return ! Exit here without adding random noise
    end if

    ! ---------------------------------------------------------
    ! Condition 2: Time-varying or Time-constant? (Random mode)
    ! ---------------------------------------------------------
    ! If it is set to time-constant and already perturbed, do nothing
    if (.not. l_pert_alp2_time .and. is_alp2_perturbed) return

    ! --- Random Seed Setup ---
    ! Ensure reproducibility across runs and members
    call random_seed(size=seed_size)
    allocate(seed_array(seed_size))
    
    if (l_pert_alp2_time) then
       ! Time-varying: Seed depends on member ID and current step
       seed_array = iens * 1000000 + step_count
    else
       ! Time-constant: Seed depends ONLY on member ID
       seed_array = iens * 1000000
    end if
    call random_seed(put=seed_array)

    ! ---------------------------------------------------------
    ! Condition 1: Spatially varying or Uniform? (Random mode)
    ! ---------------------------------------------------------
    if (l_pert_alp2_space) then
       ! (A) Spatially varying (Stochastic 2D)
       allocate(rand_field(im, jm))
       call random_number(rand_field)
       do j = 1, jm
          do i = 1, im
             ! Convert 0.0~1.0 random number to -1.0~1.0, then scale by pert_amp
             alp2(i,j) = alp2_default + (rand_field(i,j) - 0.5d0) * 2.0d0 * pert_amp
          end do
       end do
       deallocate(rand_field)
    else
       ! (B) Spatially uniform (Uniform 1D)
       call random_number(rand_val)
       alp2(:,:) = alp2_default + (rand_val - 0.5d0) * 2.0d0 * pert_amp
    end if

    ! --- Guardrails (Clipping) ---
    ! Prevent physical inconsistencies and model crashes
    do j = 1, jm
       do i = 1, im
          alp2(i,j) = max(alp2_min, min(alp2_max, alp2(i,j)))
       end do
    end do

    ! Mark as perturbed to prevent redundant updates in time-constant mode
    is_alp2_perturbed = .true.
    deallocate(seed_array)
    
  end subroutine update_perturb_alp2

end module mod_perturb_param
